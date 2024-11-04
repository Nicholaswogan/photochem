
program memtest_evo
  use photochem, only: EvoAtmosphere, dp
  implicit none

  call test()

contains

  subroutine test()

    call test_climate() ! Test climate version
    call test_methods('../data/reaction_mechanisms/zahnle_earth.yaml')
    call test_methods('../tests/no_particle_test.yaml')

  end subroutine

  subroutine test_methods(filename)
    character(*), intent(in) :: filename
    type(EvoAtmosphere) :: pcs

    call test_twoinitializations(pcs, filename) ! Test initialization
    
    ! Test all functions in code. Comments indicates if function is already
    ! tested in another test.

    ! set_trop_ind : test_climate
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

  subroutine make_new_mechanism(err)
    use fortran_yaml_c, only: YamlFile, type_dictionary, type_list, type_error, type_list_item
    character(:), allocatable, intent(out) :: err
    type(YamlFile) :: f
    type(type_list), pointer :: particles
    type(type_list_item), pointer :: list_item
    type (type_error), allocatable :: io_err

    call f%parse('../data/reaction_mechanisms/zahnle_earth.yaml', err)
    if (allocated(err)) return

    select type (root => f%root)
      class is (type_dictionary)
      particles => root%get_list('particles',.true.,error=io_err)
      if (allocated(io_err)) then; err = io_err%message; return; endif

      ! Save the first list item
      list_item => particles%first

      ! Skip the first (H2O)
      particles%first => particles%first%next

      open(unit=2,file='tmp.yaml',status='replace')
      call f%dump(2, 0)
      close(2)

      ! Replace the first list item with H2O
      particles%first => list_item

    end select

  end subroutine

  subroutine test_climate()
    use futils, only: linspace
    character(:), allocatable :: err
    type(EvoAtmosphere) :: pc
    logical :: success
    real(dp) :: tstart
    real(dp), allocatable :: t_eval(:)

    ! Make new mechanism that deletes H2Oaer
    call make_new_mechanism(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    pc = EvoAtmosphere("tmp.yaml", &
                      "../tests/test_settings1.yaml", &
                      "../examples/ModernEarth/Sun_now.txt", &
                      "../examples/ModernEarth/atmosphere.txt", &
                      "../data", &
                      err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    ! Just take a few steps
    pc%var%max_error_reinit_attempts = 0
    pc%var%mxsteps = 3

    allocate(t_eval(100))
    call linspace(5.0_dp, 17.0_dp, t_eval)
    t_eval = 10.0_dp**t_eval
    tstart = 0.0_dp

    success = pc%evolve('testevo.dat', tstart, pc%var%usol_init, t_eval, overwrite=.true., restart_from_file=.false., err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

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

    call pc%set_press_temp_edd(pc%wrk%pressure, pc%var%temperature, pc%var%edd, 0.22_dp*1e6_dp, err=err)
    if (allocated(err)) then
      print*,trim(err)
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



