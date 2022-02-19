
program main
  use Photochem, only: Atmosphere
  use Photochem, only: dp => dp
  implicit none
  character(:), allocatable :: err
  type(Atmosphere) :: pc
  real(dp), allocatable :: temperature(:)
  logical :: success
  
  call pc%init("../Photochem/data", &
               "../Photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif

  
  call pc%photochemical_equilibrium(success, err)
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif
  
end program
