
program main
  ! use photochem_input, only: read_all_files, reaction_string
  use photochem_setup, only: setup
  use photochem_vars, only: data_dir, neqs, usol_init
  use photochem, only: rhs_background_gas, print_reaction_string, photo_equilibrium
  use photochem_data, only: nrT
  implicit none
  character(len=1024) :: err
  character(len=:), allocatable :: rxstring
  integer i
  real(8), pointer :: usol_flat(:)
  real(8), allocatable :: rhs(:)
  integer cr, cm, c1, c2
  
  data_dir = "../data"
  call system_clock(count = c1, count_rate = cr, count_max = cm)
  
  call setup("../zahnle.yaml", "../settings.yaml", "../Sun_4.0Ga.txt", "../atmosphere.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  usol_flat(1:neqs) => usol_init
  allocate(rhs(neqs))
  
  call system_clock(count = c2)
  print*,(c2-c1)/real(cr)

  call rhs_background_gas(neqs, usol_flat, rhs, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  call photo_equilibrium(err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  ! do i=1,nrT
  !   call print_reaction_string(i)
  ! enddo
  ! deallocate(rxstring)
  
end program
