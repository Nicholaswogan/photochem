
program memtest_evo
  use photochem, only: EvoAtmosphere, dp
  implicit none

  type(EvoAtmosphere) :: pcs

  call test_twoinitializations(pcs)
  call test_out2atmosphere(pcs)
  call test_gas_fluxes(pcs)
  call test_set_lower_bc(pcs)
  call test_set_upper_bc(pcs)
  call test_production_and_loss(pcs)
  call test_set_temperature(pcs)
  call test_update_vertical_grid(pcs)
  call test_press_temp_edd(pcs)
  call test_autodiff(pcs)

contains

  subroutine test_twoinitializations(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    
    character(:), allocatable :: err
    real(dp) :: tn
    
    pc = EvoAtmosphere("../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
                    "..//tests/testevo_settings2.yaml", &
                    "../templates/ModernEarth/Sun_now.txt", &
                    "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
                    "../photochem/data", &
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
    
    pc = EvoAtmosphere("../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
                    "../tests/testevo_settings2.yaml", &
                    "../templates/ModernEarth/Sun_now.txt", &
                    "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
                    "../photochem/data", &
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

  subroutine test_press_temp_edd(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err

    call pc%set_press_temp_edd(pc%wrk%pressure, pc%var%temperature, pc%var%edd, 0.22_dp*1e6_dp, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

  end subroutine

  subroutine test_autodiff(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    
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



