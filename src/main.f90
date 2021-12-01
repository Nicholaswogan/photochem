
program main
  use Photochem, only: Atmosphere, err_len, real_kind, ProductionLoss, AtomConservation
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  logical :: success
  
  err = ""
  call pc%init("../Photochem/data", &
               "../test.yaml", &
               "../settings_test.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../atmosphere_test.txt", &
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

end program
