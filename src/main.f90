
program main
  use Photochem, only: Atmosphere, err_len
  use Photochem, only: dp => real_kind
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  real(dp), allocatable :: temperature(:)
  logical :: success
  
  err = ""
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
  
  allocate(temperature(size(pc%var%temperature)))
  
  temperature = pc%var%temperature + 1
  
  call pc%set_temperature(temperature, -1.d0, err=err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    ! stop 1
  endif
  
  call pc%photochemical_equilibrium(success, err)
  
end program
