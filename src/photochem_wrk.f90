
module photochem_wrk
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: real_kind = kind(1.0d0)
  ! variables that need to be carried from
  ! one iteration to the next of the photochemical model
  ! e.g.: initial conditions for nonlinear solves
  ! public ...
  
  ! also some pre-allocated work arrays

  real(real_kind), allocatable :: real_nz_nsp(:,:)
  real(real_kind), allocatable :: real_nz(:)
  integer, allocatable :: int_nz(:)
  integer(c_long) :: nsteps_previous = -10
  type(c_ptr)    :: cvode_mem  ! CVODE memory
  
  real(real_kind), allocatable :: usol_wrk(:,:)
  
  real(real_kind) :: wrk_real
  real(real_kind) :: atol_global
  
  
end module