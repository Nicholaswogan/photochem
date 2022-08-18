program testevo
  use futils, only: linspace
  use photochem, only: EvoAtmosphere, version, dp, ProductionLoss
  implicit none
  character(:), allocatable :: err
  type(EvoAtmosphere) :: pc
  logical :: success
  integer :: i, j
  real(dp) :: tstart
  real(dp), allocatable :: t_eval(:)
  type(ProductionLoss) :: pl
  
  print*,'photochem version == ',trim(version)

  call pc%init("../photochem/data", &
               "../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../tests/testevo_settings.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
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

  deallocate(t_eval)

  call pc%production_and_loss('CH4', pc%var%usol_init, pc%var%top_atmos, pl, err)  
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif

end program