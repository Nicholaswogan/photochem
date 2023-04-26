
program main
  use photochem, only: Atmosphere, version
  implicit none
  character(:), allocatable :: err
  type(Atmosphere) :: pc
  logical :: success
  integer :: i, j
  
  print*,'photochem version == ',trim(version)

  call pc%init("../photochem/data", &
               "../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/Saturn/settings_Saturn.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/Saturn/atmosphere_Saturn.txt", &
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
  
  if (.not. success) then
    print*,'Did not successfully reach equilibrium'
    stop 1
  endif

end program