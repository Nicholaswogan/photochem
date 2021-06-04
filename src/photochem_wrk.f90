
module photochem_wrk
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
  real(real_kind) :: time_previous = -tiny(1.d0)
  integer :: step_counter = 1
  
end module