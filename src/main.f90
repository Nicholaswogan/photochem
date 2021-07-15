
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

  call setup("../zahnle_titan.yaml", "../settings_titan.yaml", "../Sun_4.0Ga.txt", "../myatmosphere2.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

  use_fast_jacobian = .true.
  max_order = 5
  rtol = 1.d-3
  atol = 1.d-30
  initial_dt = 1.d-6
  call photo_equilibrium(100000, rtol, atol, success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

  call out2atmosphere_txt("../myatmosphere3.txt",.true.,.false.,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

end program
