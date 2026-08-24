program test
  use clima, only: dp
  use clima, only: AdiabatClimate
  use clima, only: ocean_solubility_fcn
  use clima_test_paths, only: fixture_file, data_dir, output_file
  implicit none
  
  type(AdiabatClimate) :: c
  character(:), allocatable :: err
  real(dp) :: T, OLR, ISR, OLR1, ISR1
  real(dp), allocatable :: P_i_surf(:), N_i_surf(:), f_i(:,:), eddy(:)
  integer :: i
  logical :: converged
  procedure(ocean_solubility_fcn), pointer :: ocean_fcn_ptr
  
  c = AdiabatClimate(fixture_file('species.yaml'), &
                     fixture_file('settings_adiabat.yaml'), &
                     fixture_file('sun.txt'), &
                     data_dir, &
                     double_radiative_grid=.false., &
                     err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  ! Test tidally locked dayside
  c%reference_pressure = 1.0e3_dp
  c%solve_for_T_trop = .true.
  c%tidally_locked_dayside = .true.
  P_i_surf = [270.0_dp, 400e-6_dp, 1.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  P_i_surf = P_i_surf*1.0e6_dp
  T = c%surface_temperature( &
      P_i_surf, &
      T_guess = 380.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,T, c%T_trop

  allocate(eddy(size(c%P)))
  eddy = 1.0e6_dp
  call c%out2atmosphere_txt(output_file('test.txt'), eddy, 5, .true., .true., err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  deallocate(eddy)

  c = AdiabatClimate(fixture_file('species.yaml'), &
                     fixture_file('settings_adiabat.yaml'), &
                     fixture_file('sun.txt'), &
                     data_dir, &
                     err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  ! Test surface_temperature_column
  c%solve_for_T_trop = .true.
  c%tidally_locked_dayside = .false.
  N_i_surf = [15.0e3_dp, 400e-6_dp*23.0_dp, 1.0*36.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  T = c%surface_temperature_column( &
      N_i_surf, &
      T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,T, c%T_trop

  ! Test surface_temperature_bg_gas
  P_i_surf = [270.0_dp, 400e-6_dp, 1.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  P_i_surf = P_i_surf*1.0e6_dp
  T = c%surface_temperature_bg_gas( &
      P_i_surf, &
      P_surf = 1.0e6_dp, bg_gas='N2', T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,T, c%T_trop

  ! Test to_regular_grid
  call c%to_regular_grid(err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  ! Test a case where CO2 begins to condense
  c%T_trop = 150.0_dp
  P_i_surf = [270.0_dp, 10.0_dp, 1.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  P_i_surf = P_i_surf*1.0e6_dp
  call c%TOA_fluxes(280.0_dp, &
      P_i_surf, ISR, OLR, &
      err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  ! Test make_profile_dry
  allocate(f_i(c%nz+1,c%sp%ng))
  f_i(1,:) = c%f_i(1,:)
  do i = 1,c%nz
    f_i(i+1,:) = c%f_i(i,:)
  enddo
  call c%TOA_fluxes_dry( &
    [c%P_surf,c%P], &
    [c%T_surf,c%T], &
    f_i, &
    ISR1, OLR1, &
    err &
  )
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,ISR1/ISR, OLR1/OLR
  deallocate(f_i)

  ! Test ocean solubility functionality
  ocean_fcn_ptr => ocean_fcn
  call c%set_ocean_solubility_fcn('H2O', ocean_fcn_ptr, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  c%use_make_column_P_guess = .true.
  N_i_surf = [15.0e4_dp, 400e-6_dp*23.0_dp, 1.0*36.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  call c%make_column(280.0_dp, &
      N_i_surf, &
      err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  print*,c%N_surface+c%N_atmos+c%N_ocean(:,1)
  print*,N_i_surf
  
  call test_RCE()

  ! deallocate
  deallocate(P_i_surf, N_i_surf)
contains

  subroutine ocean_fcn(T_surf, ng, P_i, m_i, args_p) 
    use iso_c_binding, only: c_double, c_int, c_ptr
    real(c_double), value, intent(in) :: T_surf !! K
    integer(c_int), value, intent(in) :: ng
    real(c_double), intent(in) :: P_i(ng) !! bars
    real(c_double), intent(out) :: m_i(ng) !! mol/kg
    type(c_ptr), value, intent(in) :: args_p
    m_i(:) = 0.0_dp
    m_i(2) = P_i(2)*3.4e-2_dp
    m_i(3) = P_i(3)*6.1e-4_dp
  end subroutine

  subroutine test_RCE()

    c = AdiabatClimate(fixture_file('species.yaml'), &
                     fixture_file('settings_rce.yaml'), &
                     fixture_file('sun.txt'), &
                     data_dir, &
                     err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    c%rce_solve_strategy = 3
    c%convective_max_boundary_shift = 1
    c%P_top = 1.0e2_dp
    c%T_trop = 200.0_dp
    P_i_surf = [270.0_dp, 400e-6_dp, 1.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
    P_i_surf = P_i_surf*1.0e6_dp
    T = c%surface_temperature(P_i_surf, T_guess=280.0_dp, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    c%max_rc_iters = 1
    converged = c%RCE( &
        P_i_surf, &
        c%T_surf, &
        c%T, &
        c%convecting_with_below, &
        err=err&
    )
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    ! Test RCE path 2
    c%rce_solve_strategy = 2
    converged = c%RCE( &
        P_i_surf, &
        c%T_surf + 0.1_dp, &
        c%T, &
        c%convecting_with_below, &
        err=err &
    )
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    ! Test setting particle densities
    block
      real(dp) :: pdensities(1,1), pradii(1,1)
      pdensities(1,1) = 1.0e-100_dp
      pradii(1,1) = 1.0e-4_dp
      ! Test particle density
      call c%set_particle_density_and_radii([1.0_dp],pdensities, pradii, err)
      if (allocated(err)) then
        print*,err
        stop 1
      endif
    end block

    ! Test custom mixing ratio
    block
      character(5), allocatable :: sp_custom(:)
      real(dp), allocatable :: P_custom(:), mix_custom(:,:)
      sp_custom = ['N2']
      P_custom = [1.0e6_dp, 1.0e-7_dp]
      allocate(mix_custom(2,1))
      mix_custom(:,1) = [1.0_dp, 1.0_dp]

      c%rce_solve_strategy = 3
      converged = c%RCE( &
          P_i_surf, &
          c%T_surf, &
          c%T, &
          c%convecting_with_below, &
          sp_custom, &
          P_custom, &
          mix_custom, &
          err=err&
      )
      if (allocated(err)) then
        print*,err
        stop 1
      endif
    end block

  end subroutine

end program
