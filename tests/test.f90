
program main
  use photochem, only: Atmosphere, err_len
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  logical :: success
  integer :: i, j

  call pc%init("../Photochem/data", &
               "../Photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif

  call pc%photochemical_equilibrium(success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  if (.not. success) then
    print*,'Did not successfully reach equilibrium'
    stop 1
  endif

end program