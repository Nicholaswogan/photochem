
program memtest
  use photochem, only: Atmosphere, dp
  implicit none
  
  type(Atmosphere) :: pcs

  call test_badfile(pcs)
  call test_twoinitializations(pcs)
  call test_photochemical_equilibrium(pcs)
  call test_out2atmosphere(pcs)
  call test_gas_fluxes(pcs)
  call test_set_lower_bc(pcs)
  call test_set_upper_bc(pcs)
  call test_production_and_loss(pcs)
  call test_redox_conservation(pcs)
  call test_atom_conservation(pcs)
  call test_evolve(pcs)
  call test_set_temperature(pcs)
  
contains
  
  subroutine test_badfile(pc)
    type(Atmosphere), intent(inout) :: pc
    character(:), allocatable :: err

    pc = Atmosphere("../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
                    "../templates/ModernEarth/settings_ModernEarth.yaml", &
                    "../templates/ModernEarth/Sun_now.txt", &
                    "../bad/path.txt", &
                    "../photochem/data", &
                    err)
    if (.not. allocated(err)) then
      stop 1
    endif
    print*,trim(err)
    
  end subroutine
  
  subroutine test_twoinitializations(pc)
    type(Atmosphere), intent(inout) :: pc
    
    character(:), allocatable :: err
    real(dp) :: tn
    
    pc = Atmosphere("../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
                    "../templates/Titan/settings_Titan.yaml", &
                    "../templates/ModernEarth/Sun_now.txt", &
                    "../templates/Titan/atmosphere_Titan.txt", &
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
    
    pc = Atmosphere("../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
                    "../templates/ModernEarth/settings_ModernEarth.yaml", &
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
    
    call pc%out2atmosphere_txt("test.txt", .true., .false., err)
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
  
end program