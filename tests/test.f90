
program main
  use photochem, only: Atmosphere, err_len
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  logical :: success
  integer :: i, j

  call pc%init("../data", &
               "../data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif

  call pc%photochemical_equilibrium(success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  if (.not. success) then
    print*,'Did not successfully reach equilibrium'
    stop 1
  endif
  
  ! require the solution doesn't change to within 10%
  do j = 1, pc%var%nz
    do i = 1,pc%dat%nq
      if (pc%var%usol_out(i,j) > pc%var%atol) then
        if (pc%var%usol_out(i,j) > pc%var%usol_init(i,j)*(1.d0 + 1.d-1) &
            .or. pc%var%usol_out(i,j) < pc%var%usol_init(i,j)*(1.d0 - 1.d-1)) then
          print*,'The solution changed during integration. It should not.'
          print*,i,pc%var%usol_out(i,j),pc%var%usol_init(i,j)
          stop 1
        endif
      endif
    enddo
  enddo

end program