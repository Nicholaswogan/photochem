
program test
  implicit none

  call test_Radtran()
  call test_custom_optical_properties()
  print *, "Radiative transfer and custom optical property checks passed."

contains

  subroutine test_Radtran()
    use iso_fortran_env, only: output_unit
    use clima, only: Radtran, TAU_MODE_DZ
    use clima_types, only: AtmosphereFile, unpack_atmospherefile
    use clima_eqns, only: vertical_grid
    use clima_const, only: k_boltz, dp
    use clima_test_paths, only: fixture_file, data_dir, output_file
    implicit none
  
    type(Radtran) :: rad
    type(AtmosphereFile) :: atm
    character(:), allocatable :: err
    
    integer :: nz, i
    integer :: num_zeniths
    real(dp) :: surface_albedo, T_shift
    real(dp), allocatable :: T(:), P(:), densities(:,:), mix(:,:), dz(:), z(:), density(:)
    real(dp), allocatable :: pdensities(:,:), radii(:,:)
    
    atm = AtmosphereFile(fixture_file('atmosphere.txt'), err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
  
    nz = size(atm%columns(1,:))

    num_zeniths = 8
    surface_albedo = 0.15_dp
    rad = Radtran(fixture_file('settings_radtran.yaml'), &
      fixture_file('sun.txt'), num_zeniths, surface_albedo, nz, data_dir, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    ! rad%diurnal_fac = 1.0_dp
    
    ! step up the atmosphere
    allocate(T(nz), P(nz), densities(nz,rad%ng), mix(nz,rad%ng), dz(nz), z(nz), density(nz))
    allocate(pdensities(nz,rad%np), radii(nz,rad%np))
    
    call vertical_grid(0.0_dp, 1.0e7_dp, nz, z, dz)
    call unpack_atmospherefile(atm, rad%species_names, z, mix, T, P, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    T_shift = 0.0_dp
    T = T + T_shift
    
    density = (P*1.0e6_dp)/(k_boltz*T)
    do i = 1,rad%ng
      densities(:,i) = mix(:,i)*density
    enddo

    pdensities(:,:) = 1.0_dp
    radii(:,:) = 1.0e-5_dp
    
    call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    print*,rad%wrk_sol%fdn_n(nz+1)*1.0e-3_dp

    open(unit=1,file=output_file('ModernEarth.dat'),form='unformatted',status='replace')
    write(1) rad%ir%freq
    write(1) rad%wrk_ir%fup_a(nz+1,:)
    write(1) rad%sol%freq
    write(1) rad%wrk_sol%fup_a(nz+1,:)
    write(1) rad%photons_sol
    close(1)

    write(output_unit,'(a)') rad%opacities2yaml()

    ! Test custom opacity
    block
      real(dp) :: dtau_dz(2,2), w0(2,2), g0(2,2)
      dtau_dz = 0.0_dp
      w0 = 0.0_dp
      g0 = 0.0_dp
      call rad%set_custom_optical_properties([1.0_dp, 2.0_dp],[1.0_dp, 0.1_dp], dtau_dz, w0, g0, &
                                            TAU_MODE_DZ, err)
      if (allocated(err)) then
        print*,err
        stop 1
      endif

      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, err=err)
      if (allocated(err)) then
        print*,err
        stop 1
      endif

      call rad%unset_custom_optical_properties()
    end block

  end subroutine
  
  subroutine check_error(err)
    character(:), allocatable, intent(in) :: err
    if (allocated(err)) then
      print *, err
      stop 1
    endif
  end subroutine

  subroutine test_custom_optical_properties()
    use clima, only: AdiabatClimate, Radtran, dp, TAU_MODE_DZ, TAU_MODE_DLOG10P
    use clima_const, only: Rgas, k_boltz
    use clima_test_paths, only: fixture_file, data_dir
    type(AdiabatClimate) :: c
    character(:), allocatable :: err
    real(dp), allocatable :: mix(:,:), expected(:), tau(:), w(:), g(:)
    real(dp) :: pressure(3), temperature(3), dtau(2,2), w0(2,2), g0(2,2), isr, olr
    integer :: mode, i, j

    ! Both prescriptions must produce the same optical depths and fluxes
    ! when related by the hydrostatic layer-center conversion.
    block
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      integer, parameter :: nz = 4
      type(Radtran) :: rad
      real(dp) :: T(nz), P(nz), dz(nz)
      real(dp), allocatable :: densities(:,:), pdensities(:,:), radii(:,:)
      real(dp) :: dtau(2,2), w0(2,2), g0(2,2)
      real(dp) :: height(nz), flux(nz+1), tau(nz), w(nz), g(nz), expected(nz)
      rad = Radtran(fixture_file('settings_radtran.yaml'), fixture_file('sun.txt'), &
                    2, 0.15_dp, nz, data_dir, err)
      call check_error(err)
      T = 250.0_dp
      P = [1.0_dp, 0.1_dp, 0.01_dp, 0.001_dp]
      dz = 1.0e5_dp
      allocate(densities(nz,rad%ng), pdensities(nz,rad%np), radii(nz,rad%np))
      do i = 1,rad%ng
        densities(:,i) = P*1.0e6_dp/(k_boltz*T*rad%ng)
      enddo
      pdensities = 0.0_dp
      radii = 1.0e-5_dp
      dtau = 1.0e-7_dp
      w0 = 0.5_dp
      g0 = 0.3_dp
      call rad%set_custom_optical_properties([1.0_dp, 2.0_dp], [1.0_dp, 0.1_dp], dtau, w0, g0, &
                                            tau_mode=TAU_MODE_DZ, err=err)
      call check_error(err)
      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, err=err)
      call check_error(err)
      flux = rad%f_total
      call rad%op%custom_optical_properties(log10(P*1.0e6_dp), dz, 1, tau, w, g)
      if (any(abs(tau - 1.0e-7_dp*dz) > 1.0e-12_dp)) error stop 'Incorrect per-cm optical depth'

      dtau = 1.0e-7_dp*1.0e6_dp*log(10.0_dp)
      call rad%set_custom_optical_properties([1.0_dp, 2.0_dp], [1.0_dp, 0.1_dp], &
                                            dtau, w0, g0, tau_mode=TAU_MODE_DLOG10P, err=err)
      call check_error(err)
      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, err=err)
      if (.not.allocated(err)) error stop 'Missing scale height was accepted'
      height = 1.0e6_dp
      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, scale_height=height, err=err)
      call check_error(err)
      if (any(abs(rad%f_total - flux) > 1.0e-10_dp*max(1.0_dp, abs(flux)))) then
        error stop 'Optical-depth modes give different fluxes'
      endif
      ! Check varying scale heights and both dimensionless properties.
      height = [(1.0e6_dp + i*1.0e4_dp, i=1,nz)]
      expected = dtau(1,1)*dz/(height*log(10.0_dp))
      call rad%op%custom_optical_properties(log10(P*1.0e6_dp), dz, 1, tau, w, g, height)
      if (any(abs(tau - expected) > 1.0e-12_dp)) error stop 'Incorrect per-decade optical depth'
      if (any(abs(w - 0.5_dp) > 1.0e-12_dp)) error stop 'Incorrect custom albedo'
      if (any(abs(g - 0.3_dp) > 1.0e-12_dp)) error stop 'Incorrect custom asymmetry'

      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, scale_height=height(:nz-1), err=err)
      if (.not.allocated(err)) error stop 'Wrong scale-height shape was accepted'
      height(1) = 0.0_dp
      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, scale_height=height, err=err)
      if (.not.allocated(err)) error stop 'Zero scale height was accepted'
      height(1) = ieee_value(1.0_dp, ieee_quiet_nan)
      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, scale_height=height, err=err)
      if (.not.allocated(err)) error stop 'Nonfinite scale height was accepted'
      call rad%set_custom_optical_properties([1.0_dp, 2.0_dp], [1.0_dp, 0.1_dp], &
                                            dtau, w0, g0, tau_mode=-1, err=err)
      if (.not.allocated(err)) error stop 'Invalid tau mode was accepted'

      ! Switching back to per-cm mode must remove the scale-height requirement.
      call rad%set_custom_optical_properties([1.0_dp, 2.0_dp], [1.0_dp, 0.1_dp], &
                                            dtau, w0, g0, tau_mode=TAU_MODE_DZ, err=err)
      call check_error(err)
      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, err=err)
      call check_error(err)
      call rad%unset_custom_optical_properties()
      call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, err=err)
      call check_error(err)
    end block

    pressure = [1.0e6_dp, 1.0e3_dp, 1.0_dp]
    temperature = [300.0_dp, 250.0_dp, 200.0_dp]
    dtau = 0.01_dp
    w0 = 0.5_dp
    g0 = 0.3_dp
    do mode = 1,2
      c = AdiabatClimate(fixture_file('species.yaml'), fixture_file('settings_adiabat.yaml'), &
                        fixture_file('sun.txt'), data_dir, double_radiative_grid=(mode == 2), err=err)
      call check_error(err)
      allocate(mix(3,c%sp%ng), expected(c%nz), tau(c%nz_r), w(c%nz_r), g(c%nz_r))
      mix = 1.0_dp/c%sp%ng
      call c%make_profile_dry(pressure, temperature, mix, err)
      call check_error(err)
      expected = 0.0_dp
      do j = 1,c%sp%ng
        expected = expected + c%f_i(:,j)*c%sp%g(j)%mass
      enddo
      expected = Rgas*c%T/(expected*c%gravity)
      if (any(abs(c%scale_height/expected - 1.0_dp) > 1.0e-12_dp)) error stop 'Incorrect physical scale height'
      call c%rad%set_custom_optical_properties([1.0_dp, 2.0_dp], [1.0e6_dp, 1.0_dp], &
                                              dtau, w0, g0, tau_mode=TAU_MODE_DLOG10P, err=err)
      call check_error(err)
      ! Exercises TOA_fluxes forwarding of scale height and copying to the RT grid.
      call c%TOA_fluxes_dry(pressure, temperature, mix, isr, olr, err)
      call check_error(err)
      if (mode == 1) then
        if (any(c%scale_height_r /= c%scale_height)) error stop 'Incorrect single RT grid scale height'
      else
        do i = 1,c%nz
          if (any(c%scale_height_r(2*i-1:2*i) /= c%scale_height(i))) error stop 'Incorrect doubled RT grid scale height'
        enddo
        if (any(c%scale_height_r(2*c%nz+1:) /= c%scale_height(c%nz))) error stop 'Incorrect ghost scale height'
      endif
      call c%rad%op%custom_optical_properties(log10(c%P_r), c%dz_r, 1, tau, w, g, c%scale_height_r)
      if (mode == 2) then
        do i = 1,c%nz
          if (abs(sum(tau(2*i-1:2*i)) - 0.01_dp*c%dz(i)/(c%scale_height(i)*log(10.0_dp))) &
              > 1.0e-12_dp) error stop 'Splitting a layer changed optical depth'
        enddo
      endif
      deallocate(mix, expected, tau, w, g)
    enddo
  end subroutine

end program
