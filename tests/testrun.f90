
program main
  use photochem, only: Atmosphere, version
  implicit none
  character(:), allocatable :: err
  type(Atmosphere) :: pc
  logical :: success
  integer :: i, j
  
  print*,'photochem version == ',trim(version)

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
  pc%var%atol = 1.0e-21

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