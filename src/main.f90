
program main
  use photochem_setup, only: setup, out2atmosphere_txt
  use photochem_vars, only: data_dir, max_order, initial_dt, equilibrium_time, verbose
  use photochem, only: photo_equilibrium
  implicit none
  character(len=1024) :: err
  real(8) :: rtol, atol
  logical :: success

  data_dir = "../data"

  call setup("../data/reaction_mechanisms/zahnle_earth.yaml", &
             "../templates/ModernEarth/settings_ModernEarth.yaml", &
             "../templates/ModernEarth/Sun_now.txt", &
             "../templates/ModernEarth/atmosphere_ModernEarth.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

  equilibrium_time = 1.d17
  max_order = 5
  rtol = 1.d-3
  atol = 1.d-25
  initial_dt = 1.d-9
  verbose = 1
  call photo_equilibrium(100000, rtol, atol, success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  ! call out2atmosphere_txt("../atmosphere_ModernEarth_8.txt",.true.,.false.,err)
  ! if (len(trim(err)) > 0) then
  !   print*,trim(err)
  !   stop
  ! endif

end program
