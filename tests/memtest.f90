
program main
  use photochem, only: Atmosphere, err_len, dp
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  real(dp) :: tn

  call pc%init("../photochem/data", &
               "../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../bad/path.txt", &
               err)
  if (len_trim(err) == 0) then
    stop 1
  endif
  print*,trim(err)
  err = ""
  
  call pc%init("../photochem/data", &
               "../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len_trim(err) > 0) then
    print*,trim(err)
    stop 1
  endif

  call pc%initialize_stepper(pc%var%usol_init, err)
  if (len_trim(err) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  tn = pc%step(err)
  if (len_trim(err) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  call pc%init("../photochem/data", &
               "../photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/Titan/settings_Titan.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/Titan/atmosphere_Titan.txt", &
               err)
  if (len_trim(err) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  call pc%initialize_stepper(pc%var%usol_init, err)
  if (len_trim(err) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  tn = pc%step(err)
  if (len_trim(err) > 0) then
    print*,trim(err)
    stop 1
  endif

end program