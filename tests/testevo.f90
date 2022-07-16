program testevo
  use photochem, only: EvoAtmosphere, version
  implicit none
  character(:), allocatable :: err
  type(EvoAtmosphere) :: pc
  logical :: success
  integer :: i, j
  
  print*,'photochem version == ',trim(version)

  call pc%init("../photochem/data", &
               "../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif

  call pc%prep_atmosphere(pc%var%dsol_init, err)
  if (allocated(err)) then
    print*,trim(err)
    stop 1
  endif

end program