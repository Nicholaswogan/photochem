
program main
  use photochem_setup, only: setup, out2atmosphere_txt
  use photochem_vars, only: data_dir, max_order, initial_dt, equilibrium_time
  use photochem, only: photo_equilibrium
  implicit none
  character(len=1024) :: err
  real(8) :: rtol, atol
  logical :: success

  data_dir = "../data"

  call setup("../zahnle_earth.yaml", "../settings_Hadean.yaml", "../Sun_4.0Ga.txt", "../atmosphere10_12.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

  equilibrium_time = 1.d17
  max_order = 5
  rtol = 1.d-3
  atol = 1.d-27
  initial_dt = 1.d-9
  call photo_equilibrium(100000, rtol, atol, success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  call out2atmosphere_txt("../atmosphere10_13.txt",.true.,.false.,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

end program
