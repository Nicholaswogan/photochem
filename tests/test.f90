
program main
  ! import functions
  use photochem, only: photo_equilibrium
  use photochem_setup, only: setup
  ! import varibles
  use photochem_data, only: nq
  use photochem_vars, only: data_dir, usol_init, usol_out, nz
  use photochem_const, only: err_len, real_kind
  implicit none
  
  character(len=err_len) :: err
  real(real_kind) :: rtol, atol
  logical :: success
  integer :: i, j

  data_dir = "../data/"

  call setup("../data/reaction_mechanisms/zahnle_earth.yaml", &
             "../templates/ModernEarth/settings_ModernEarth.yaml", &
             "../templates/ModernEarth/Sun_now.txt", &
             "../templates/ModernEarth/atmosphere_ModernEarth.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif

  rtol = 1.d-3
  atol = 1.d-25
  call photo_equilibrium(100000, rtol, atol, success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  if (.not. success) then
    print*,'Did not successfully reach equilibrium'
    stop 1
  endif
  
  ! require the solution doesn't change to within 10%
  do j = 1, nz
    do i = 1,nq
      if (usol_out(i,j) > atol) then
        if (usol_out(i,j) > usol_init(i,j)*(1.d0 + 1.d-1) &
            .or. usol_out(i,j) < usol_init(i,j)*(1.d0 - 1.d-1)) then
          print*,'The solution changed during integration. It should not.'
          print*,i,usol_out(i,j),usol_init(i,j)
          stop 1
        endif
      endif
    enddo
  enddo

end program