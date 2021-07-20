
program main
  use photochem_setup, only: setup, out2atmosphere_txt
  use photochem_data, only: nq
  use photochem_vars, only: data_dir, max_order, initial_dt, equilibrium_time, &
                            usol_init, usol_out, nz
  use photochem, only: photo_equilibrium
  implicit none
  character(len=1024) :: err
  real(8) :: rtol, atol
  logical :: success
  integer :: i, j

  data_dir = "data/"

  call setup("data/reaction_mechanisms/zahnle_earth.yaml", &
             "templates/ModernEarth/settings_ModernEarth.yaml", &
             "templates/ModernEarth/Sun_now.txt", &
             "templates/ModernEarth/atmosphere_ModernEarth.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif

  equilibrium_time = 1.d17
  max_order = 5
  rtol = 1.d-3
  atol = 1.d-30
  initial_dt = 1.d-9
  call photo_equilibrium(100000, rtol, atol, success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  if (.not. success) then
    stop 1
  endif
  
  ! require the solution doesn't change to within 1%
  do j = 1, nz
    do i = 1,nq
      if (usol_out(i,j) > atol) then
        if (usol_out(i,j) > usol_init(i,j)*(1.d0 + 1.d-1) &
            .or. usol_out(i,j) < usol_init(i,j)*(1.d0 - 1.d-1)) then
          print*,i,usol_out(i,j),usol_init(i,j)
          stop 1
        endif
      endif
    enddo
  enddo

end program