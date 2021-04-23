
module photochem_vars
  implicit none
  private
  integer, parameter :: real_kind = kind(1.0d0)
  ! variables that need to be carried from
  ! one iteration to the next of the photochemical model
  ! e.g.: initial conditions for nonlinear solves
  public ...
  
  ! also some pre-allocated work arrays
  public real_nz_nsp, real_nz, int_nz
  real(real_kind), allocatable :: real_nz_nsp(:,:)
  real(real_kind), allocatable :: real_nz(:)
  integer, allocatable :: int_nz(:)
end module