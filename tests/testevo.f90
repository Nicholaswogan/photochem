program testevo
  use photochem, only: EvoAtmosphere, version, dp
  implicit none
  character(:), allocatable :: err
  type(EvoAtmosphere) :: pc
  logical :: success
  integer :: i, j
  real(dp), allocatable :: t_eval(:)
  real(dp) :: from, to
  
  print*,'photochem version == ',trim(version)

  call pc%init("../photochem/data", &
               "../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/Hescape/settings.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/Hescape/atmosphere.txt", &
               err)
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif

  pc%var%atol = 1.0e-25_dp
  pc%P_top_min = 1.0e-7_dp
  pc%P_top_max = 1.0e100_dp
  pc%top_atmos_adjust_frac = 0.05

  allocate(t_eval(400))
  from = 5
  to = 15
  do i = 1,size(t_eval)
    t_eval(i) = from + (to - from)*(i-1)/(size(t_eval)-1)
  enddo
  t_eval = 10.0_dp**t_eval

  success = pc%evolve('../test16.dat', 0.0_dp, pc%var%usol_init, t_eval, .true., err)
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif

end program