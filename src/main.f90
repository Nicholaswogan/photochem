
program main
  use photochem_setup, only: setup, out2atmosphere_txt
  use photochem_vars, only: data_dir, use_fast_jacobian, max_order, initial_dt
  use photochem, only: photo_equilibrium
  implicit none
  character(len=1024) :: err
  character(len=:), allocatable :: rxstring
  real(8) :: rtol, atol
  logical :: success
  integer i
  
  data_dir = "../data"
  
  call setup("../zahnle_earth.yaml", "../settings.yaml", "../Sun_4.0Ga.txt", "../atmosphere_nothing.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  use_fast_jacobian = .true. 
  max_order = 5
  rtol = 1.d-4
  atol = 1.d-35
  initial_dt = 1.d-15
  call photo_equilibrium(1000000000,rtol, atol, success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  call out2atmosphere_txt("../myatmosphere.txt",.true.,.true.,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

end program
