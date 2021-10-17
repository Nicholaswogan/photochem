
program main
  use Atmos, only: Photochem, err_len
  implicit none
  character(len=err_len) :: err
  type(Photochem) :: pc
  
  integer :: i, j, k
  real(8), allocatable :: usol_flat(:), rhs(:)
  
  err = ""
  call pc%init("../data", &
               "../data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  allocate(usol_flat(pc%var%neqs))
  allocate(rhs(pc%var%neqs))
  
  do j=1,pc%var%nz
    do i=1,pc%dat%nq
      k = i + (j-1)*pc%dat%nq
      usol_flat(k) = pc%var%usol_init(i,j)
    enddo
  enddo
  
  call pc%rhs(pc%var%neqs, usol_flat, rhs, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  do i = 1,pc%dat%nq
    print*,rhs(i)
  enddo
  
  

end program
