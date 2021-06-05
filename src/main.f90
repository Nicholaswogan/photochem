
program main
  ! use photochem_input, only: read_all_files, reaction_string
  use photochem_setup, only: setup
  use photochem_vars, only: data_dir, neqs, usol_init, nz
  use photochem, only: rhs_background_gas, jac_background_gas, print_reaction_string, photo_equilibrium
  use photochem_data, only: nrT, nq, lda
  implicit none
  character(len=1024) :: err
  character(len=:), allocatable :: rxstring
  integer i, lda_neqs
  real(8), pointer :: usol_flat(:)
  real(8), allocatable :: rhs(:), jac(:)
  real(8), allocatable :: usol_out(:,:)
  integer cr, cm, c1, c2
  
  data_dir = "../data"
  call system_clock(count = c1, count_rate = cr, count_max = cm)
  
  call setup("../zahnle_titan.yaml", "../settings.yaml", "../Sun_4.0Ga.txt", "../atmosphere_nothing.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  call system_clock(count = c2)
  ! print*,(c2-c1)/real(cr)
  
  
  usol_flat(1:neqs) => usol_init
  allocate(rhs(neqs))
  lda_neqs = lda*neqs
  allocate(jac(lda_neqs))
  
  
  call jac_background_gas(lda_neqs, neqs, usol_flat, jac, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  ! print*,nq
  ! do i= 1,lda
  !   print*,jac(i),i
  ! enddo
  ! stop

  call rhs_background_gas(neqs, usol_flat, rhs, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  allocate(usol_out(nq,nz))
  
  call photo_equilibrium(nq,nz,usol_out,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  do i = 1,nq
    print*,i,usol_out(i,200)
  enddo
  

end program
