program main
  use photochem_mie
  implicit none
  
  integer, parameter :: real_kind = kind(1.0d0)
  integer, parameter :: err_len = 1024
  
  character(len=100) :: filename
  integer, parameter :: nw = 3
  real(real_kind) :: wavl(nw+1)
  
  integer :: nrad
  real(real_kind), allocatable :: radii(:)
  real(real_kind), allocatable :: w0_file(:,:), qext_file(:,:), g_file(:,:)
  character(len=err_len) :: err
  
  integer :: i
  
  
  filename = "../data/aerosol_xsections/khare1984/mie_khare1984.dat"
  err = ""
  wavl = [100.d0, 101.d0, 300.d0, 400.d0]
  
  
  
  call read_mie_data_file(filename, nw, wavl, &
                          nrad, radii, w0_file, qext_file, g_file, err)
  if (len_trim(err) /= 0) then
    print*,trim(err)
    return
  endif                    
  
  
  
  do i =1,nrad
    print*,radii(i),qext_file(1,i)
  enddo
  
  
  
  
end program