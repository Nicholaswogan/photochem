
program main
  use photochem, only: Atmosphere, err_len, real_kind
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  logical :: success
  real(real_kind) :: tn
  integer :: i, j

  call pc%init("../data", &
               "../data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif

  call pc%initialize_stepper(pc%var%usol_init, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  tn = pc%step(err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif

end program