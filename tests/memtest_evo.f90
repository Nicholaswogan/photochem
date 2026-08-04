
program memtest_evo
  use photochem, only: EvoAtmosphere, dp
  implicit none

  call test()

contains

  subroutine test()

    call test_methods('../data/reaction_mechanisms/zahnle_earth.yaml')
    call test_methods('../tests/no_particle_test.yaml')

  end subroutine

  subroutine test_methods(filename)
    character(*), intent(in) :: filename
    type(EvoAtmosphere) :: pcs

    call test_twoinitializations(pcs, filename) ! Test initialization
    
    ! Test all functions in code. Comments indicates if function is already
    ! tested in another test.

    ! set_trop_ind : test_step
    ! prep_atm_evo_gas : test_production_and_loss
    ! prep_atmosphere : test_production_and_loss
    ! right_hand_side_chem : gas_fluxes
    call test_production_and_loss(pcs)
    ! right_hand_side : test_step
    ! jacobian : test_step
    call test_evolve(pcs)
    ! check_for_convergence : test_step
    ! initialize_stepper : test_step
    call test_step(pcs)
    ! destroy_stepper : test_step
    call test_robust_step(pcs)
    call test_out2atmosphere(pcs)
    call test_gas_fluxes(pcs) 
    call test_set_lower_bc(pcs)
    call test_set_upper_bc(pcs)
    ! set_rate_fcn : NOT TESTED
    call test_set_temperature(pcs)
    call test_set_press_temp_edd(pcs)
    call test_update_vertical_grid(pcs)
    ! rebin_update_vertical_grid : NOT TESTED
    ! regrid_prep_atmosphere : NOT TESTED
    
  end subroutine

  subroutine test_twoinitializations(pc, filename)
    type(EvoAtmosphere), intent(inout) :: pc
    character(*), intent(in) :: filename
    
    character(:), allocatable :: err
    real(dp) :: tn
    
    pc = EvoAtmosphere(filename, &
                       "../examples/ModernEarth/settings.yaml", &
                       "../examples/ModernEarth/Sun_now.txt", &
                       "../examples/ModernEarth/atmosphere.txt", &
                       "../data", &
                       err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    call pc%initialize_stepper(pc%var%usol_init, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    tn = pc%step(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    pc = EvoAtmosphere(filename, &
                       "../examples/ModernEarth/settings.yaml", &
                       "../examples/ModernEarth/Sun_now.txt", &
                       "../examples/ModernEarth/atmosphere.txt", &
                       "../data", &
                       err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    call pc%initialize_stepper(pc%var%usol_init, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    tn = pc%step(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
  end subroutine

  subroutine test_out2atmosphere(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    call pc%out2atmosphere_txt("test.txt", 4, .true., .false., err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine

  subroutine test_gas_fluxes(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    real(dp), allocatable :: surf_fluxes(:)
    real(dp), allocatable :: top_fluxes(:)
    
    allocate(surf_fluxes(pc%dat%nq))
    allocate(top_fluxes(pc%dat%nq))
    
    call pc%gas_fluxes(surf_fluxes, top_fluxes, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine

  subroutine test_set_lower_bc(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    call pc%set_lower_bc("HCN","vdep",vdep=1.0e-3_dp, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine
  
  subroutine test_set_upper_bc(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    call pc%set_upper_bc("HCN","veff",veff=0.0_dp, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine

  subroutine test_production_and_loss(pc)
    use photochem, only: ProductionLoss
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    type(ProductionLoss) :: pl
    
    call pc%production_and_loss("HCN", pc%wrk%usol, pl, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine

  subroutine test_set_temperature(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    real(dp), allocatable :: temperature(:)
    
    allocate(temperature(pc%var%nz))
    temperature = 300.0_dp
    
    call pc%set_temperature(temperature, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine

  subroutine test_update_vertical_grid(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err

    call pc%update_vertical_grid(TOA_pressure=0.01_dp, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

  end subroutine

  subroutine test_set_press_temp_edd(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    real(dp) :: P(2), T(2), edd(2)
    real(dp) :: T_expected(pc%var%nz), log10edd_expected(pc%var%nz)
    real(dp) :: fraction
    integer :: i

    ! Two input points are sufficient because interpolation is linear in log10(P).
    ! First test the default hydrostatic-pressure mode with profiles that differ
    ! substantially from the current atmosphere.
    P = [2.0_dp*pc%var%surface_pressure*1.0e6_dp, &
         0.5_dp*pc%wrk%pressure_hydro(pc%var%nz)]
    T = [310.0_dp, 160.0_dp]
    edd = [1.0e7_dp, 2.0e5_dp]

    call pc%set_press_temp_edd(P, T, edd, pc%wrk%pressure_hydro(10), &
                               hydro_pressure=.true., err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    do i = 1,pc%var%nz
      fraction = (log10(pc%wrk%pressure_hydro(i))-log10(P(2)))/ &
                 (log10(P(1))-log10(P(2)))
      fraction = min(max(fraction,0.0_dp),1.0_dp)
      T_expected(i) = T(2) + fraction*(T(1)-T(2))
      log10edd_expected(i) = log10(edd(2)) + fraction*(log10(edd(1))-log10(edd(2)))
    enddo
    if (maxval(abs(pc%var%temperature-T_expected)) > 1.0e-6_dp) then
      print*,'set_press_temp_edd did not reproduce the hydrostatic P-T profile'
      stop 1
    endif
    if (maxval(abs(log10(pc%var%edd)-log10edd_expected)) > 1.0e-8_dp) then
      print*,'set_press_temp_edd did not reproduce the hydrostatic P-Kzz profile'
      stop 1
    endif

    ! Then test actual-pressure mode. This pressure is density*k*T and each
    ! layer can therefore be solved independently.
    P = [2.0_dp*maxval(pc%wrk%pressure), 0.5_dp*minval(pc%wrk%pressure)]
    T = [280.0_dp, 140.0_dp]
    edd = [2.0e8_dp, 8.0e4_dp]

    call pc%set_press_temp_edd(P, T, edd, pc%wrk%pressure(10), &
                               hydro_pressure=.false., err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    do i = 1,pc%var%nz
      fraction = (log10(pc%wrk%pressure(i))-log10(P(2)))/ &
                 (log10(P(1))-log10(P(2)))
      fraction = min(max(fraction,0.0_dp),1.0_dp)
      T_expected(i) = T(2) + fraction*(T(1)-T(2))
      log10edd_expected(i) = log10(edd(2)) + fraction*(log10(edd(1))-log10(edd(2)))
    enddo
    if (maxval(abs(pc%var%temperature-T_expected)) > 1.0e-6_dp) then
      print*,'set_press_temp_edd did not reproduce the actual P-T profile'
      stop 1
    endif
    if (maxval(abs(log10(pc%var%edd)-log10edd_expected)) > 1.0e-8_dp) then
      print*,'set_press_temp_edd did not reproduce the actual P-Kzz profile'
      stop 1
    endif

  end subroutine

  subroutine test_step(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    
    character(:), allocatable :: err
    real(dp) :: tn
    logical :: converged

    pc%var%autodiff = .true.
    pc%var%atol = 1.0e-20_dp
    
    call pc%initialize_stepper(pc%var%usol_init, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    tn = pc%step(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    converged = pc%check_for_convergence(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    pc%var%autodiff = .false.
    
    call pc%initialize_stepper(pc%var%usol_init, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    tn = pc%step(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    converged = pc%check_for_convergence(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
  end subroutine

  subroutine test_robust_step(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    
    character(:), allocatable :: err
    logical :: give_up, converged
    
    call pc%initialize_robust_stepper(pc%var%usol_init, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    call pc%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    call pc%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
  end subroutine

  subroutine test_evolve(pc)
    use futils, only: linspace
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    logical :: success
    real(dp) :: tstart
    real(dp), allocatable :: t_eval(:)

    pc%var%max_error_reinit_attempts = 0
    pc%var%mxsteps = 3

    allocate(t_eval(100))
    call linspace(5.0_dp, 17.0_dp, t_eval)
    t_eval = 10.0_dp**t_eval
    tstart = 0.0_dp
    success = pc%evolve('tmp.dat',tstart, pc%wrk%usol, t_eval, overwrite=.true., err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

  end subroutine

end program

