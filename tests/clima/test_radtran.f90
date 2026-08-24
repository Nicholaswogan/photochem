
program test
  implicit none

  call test_Radtran()

contains

  subroutine test_Radtran()
    use iso_fortran_env, only: output_unit
    use clima, only: Radtran
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
      call rad%set_custom_optical_properties([1.0_dp, 2.0_dp],[1.0_dp, 0.1_dp], dtau_dz, w0, g0, err)
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
  
end program
