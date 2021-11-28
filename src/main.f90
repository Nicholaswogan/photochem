
program main
  use Photochem, only: Atmosphere, err_len, real_kind, ProductionLoss, AtomConservation
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
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

  call pc%evolve("../test.dat", 0.d0, pc%var%usol_init, [1.d0,1.d17], success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif

end program
