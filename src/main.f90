
program main
  use Photochem, only: Atmosphere, err_len
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  logical :: success
  integer :: i
  
  err = ""
  call pc%init("../data", &
               "../data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  call pc%photochemical_equilibrium(success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

  


end program
