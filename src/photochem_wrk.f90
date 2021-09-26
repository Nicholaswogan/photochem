
module photochem_wrk
  use, intrinsic :: iso_c_binding, only: c_long, c_ptr
  use photochem_types, only: WrkBackgroundAtm
  implicit none

  integer, parameter :: real_kind = kind(1.0d0)
  
  ! also some pre-allocated work array
  integer(c_long) :: nsteps_previous = -10
  type(c_ptr)    :: cvode_mem  ! CVODE memory
  real(real_kind) :: atol_global
  
  type(WrkBackgroundAtm) :: wrk_out ! output
  
end module