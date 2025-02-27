
program memtest
  use photochem, only: Atmosphere, dp
  implicit none

  call test()
  
contains

  subroutine test()
    call test_methods('../data/reaction_mechanisms/zahnle_earth.yaml')
    call test_methods('../tests/no_particle_test.yaml')
  end subroutine

  subroutine test_methods(filename)
    character(*), intent(in) :: filename
    type(Atmosphere) :: pcs

    call test_badfile(pcs, filename)
    call test_twoinitializations(pcs, filename)

    ! prep_atm_background_gas : test_gas_fluxes
    ! prep_atmosphere : test_gas_fluxes
    ! right_hand_side_chem : test_gas_fluxes
    call test_production_and_loss(pcs)
    ! right_hand_side : test_photochemical_equilibrium
    ! jacobian : test_photochemical_equilibrium
    call test_evolve(pcs)
    ! check_for_convergence : test_photochemical_equilibrium
    call test_photochemical_equilibrium(pcs)
    ! initialize_stepper : test_photochemical_equilibrium
    ! step : test_photochemical_equilibrium
    ! destory_stepper : test_photochemical_equilibrium
    call test_out2atmosphere(pcs)
    ! out2in : NOT TESTED
    call test_gas_fluxes(pcs)
    call test_atom_conservation(pcs)
    call test_redox_conservation(pcs)
    call test_set_lower_bc(pcs)
    call test_set_upper_bc(pcs)
    call test_set_temperature(pcs)
    call test_set_press_temp_edd(pcs)
    ! set_rate_fcn :: NOT TESTED
    call test_update_vertical_grid(pcs)

    ! Other tests
    call test_custom_binary_diffusion(pcs)
    call test_autodiff(pcs)

  end subroutine
  
  subroutine test_badfile(pc, filename)
    type(Atmosphere), intent(inout) :: pc
    character(*), intent(in) :: filename
    character(:), allocatable :: err

    pc = Atmosphere(filename, &
                    "../examples/ModernEarth/settings_Atmosphere.yaml", &
                    "../examples/ModernEarth/Sun_now.txt", &
                    "../bad/path.txt", &
                    "../data", &
                    err)
    if (.not. allocated(err)) then
      stop 1
    endif
    print*,trim(err)
    
  end subroutine
  
  subroutine test_twoinitializations(pc, filename)
    type(Atmosphere), intent(inout) :: pc
    character(*), intent(in) :: filename
    
    character(:), allocatable :: err
    real(dp) :: tn
    
    pc = Atmosphere(filename, &
                    "../examples/ModernEarth/settings_Atmosphere.yaml", &
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
    
    pc = Atmosphere(filename, &
                    "../examples/ModernEarth/settings_Atmosphere.yaml", &
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
  
  subroutine test_photochemical_equilibrium(pc)
    type(Atmosphere), intent(inout) :: pc
    
    character(:), allocatable :: err
    logical :: success
    
    pc%var%mxsteps = 2
    call pc%photochemical_equilibrium(success, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    pc%var%mxsteps = 1000000
  
  end subroutine
  
  subroutine test_out2atmosphere(pc)
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    call pc%out2atmosphere_txt("test.txt", 4, .true., .false., err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine
  
  subroutine test_gas_fluxes(pc)
    type(Atmosphere), intent(inout) :: pc
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
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    call pc%set_lower_bc("HCN","vdep",vdep=1.0e-3_dp, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine
  
  subroutine test_set_upper_bc(pc)
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    call pc%set_upper_bc("HCN","veff",veff=0.0_dp, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine
  
  subroutine test_production_and_loss(pc)
    use photochem, only: ProductionLoss
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    type(ProductionLoss) :: pl
    
    call pc%production_and_loss("HCN", pc%wrk%usol, pl, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine
  
  subroutine test_redox_conservation(pc)
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    real(dp) :: redox_factor
    
    redox_factor = pc%redox_conservation(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine
  
  subroutine test_atom_conservation(pc)
    use photochem, only: AtomConservation
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    type(AtomConservation) :: con
    
    con = pc%atom_conservation("H",err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine
  
  subroutine test_evolve(pc)
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    
    logical :: success
    
    success = pc%evolve(filename="test.dat", &
                        tstart=0.0_dp, &
                        usol_start=pc%wrk%usol, &
                        t_eval=[1.0e-7_dp], &
                        overwrite=.true., &
                        err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
  
  end subroutine
  
  subroutine test_set_temperature(pc)
    type(Atmosphere), intent(inout) :: pc
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

  subroutine test_custom_binary_diffusion(pc)
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    real(dp) :: tn

    pc%var%custom_binary_diffusion_fcn => custom_binary_diffusion_param

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

    pc%var%custom_binary_diffusion_fcn => null()

  end subroutine

  function custom_binary_diffusion_param(mu_i, mubar, T) result(b)
    use iso_c_binding, only: c_double
    real(c_double), value, intent(in) :: mu_i, mubar, T
    real(c_double) :: b
    b = 1.52e18_dp*((1.0_dp/mu_i+1.0_dp/mubar)**0.5e0_dp)*(T**0.5e0_dp)
  end function

  subroutine test_update_vertical_grid(pc)
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err

    call pc%update_vertical_grid(TOA_pressure=0.01_dp, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

  end subroutine

  subroutine test_set_press_temp_edd(pc)
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err

    call pc%set_press_temp_edd(pc%wrk%pressure, pc%var%temperature, pc%var%edd, 0.22_dp*1e6_dp, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

  end subroutine

  subroutine test_autodiff(pc)
    type(Atmosphere), intent(inout) :: pc
    
    character(:), allocatable :: err
    real(dp) :: tn

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

    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
  end subroutine
  
end program