program main
  use radtran, only: two_stream
  implicit none
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: nz = 200
  real(real_kind) :: tau(nz)
  real(real_kind) :: w0(nz), u0, Rsfc, amean(nz+1), surface_radiance
  integer :: ierr, i

  ! u0 = dcos(3.14159d0/3.d0)
  u0 = 0.6427876096865394d0
  ! tau = 0.1d0
  ! w0  = 0.5d0
  Rsfc = 0.25d0
  open(2,file='tau.txt',status='old')
  
  do i = 1,nz
    read(2,*) tau(i),w0(i)
  enddo
  close(2)
    

  call two_stream(nz, tau, w0, u0, Rsfc, amean, surface_radiance,ierr)
  ! print*,amean

  ! print*,surface_radiance
  print*,''

  call TWOSTR(nz, tau, w0, u0, Rsfc, amean)

  ! print*,amean

  ! print*,surface_radiance

end program