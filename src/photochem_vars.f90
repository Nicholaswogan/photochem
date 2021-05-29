

module photochem_vars
  implicit none
  public
  integer, private, parameter :: real_kind = kind(1.0d0)
  integer, private, parameter :: str_len = 1024

  ! boundary conditions
  integer, allocatable :: lowerboundcond(:) ! 0, 1, 2 or 3
  real(real_kind), allocatable :: lower_vdep(:)
  real(real_kind), allocatable :: lower_flux(:)
  real(real_kind), allocatable :: lower_dist_height(:)
  real(real_kind), allocatable :: lower_fix_mr(:)
  integer, allocatable :: upperboundcond(:) ! 0 or 1
  real(real_kind), allocatable :: upper_veff(:)
  real(real_kind), allocatable :: upper_flux(:)
  
  ! Atmospheres structure
  real(real_kind) :: bottom_atmos
  real(real_kind) :: top_atmos 
  integer :: nz
  real(real_kind) :: surface_pressure 
  real(real_kind) :: surface_albedo 
  real(real_kind) :: trop_alt 
  
  ! Radiative tranfer
  real(real_kind), allocatable :: photon_flux(:) ! (nw) photonz
  
  
  
  
  ! need to work where to allocate these below
  ! T, z, edd, grav
  real(real_kind), allocatable :: temperature(:)
  real(real_kind), allocatable :: z(:)
  real(real_kind), allocatable :: edd(:)
  real(real_kind), allocatable :: grav(:)
  
  
end module