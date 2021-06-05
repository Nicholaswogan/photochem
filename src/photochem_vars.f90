

module photochem_vars
  implicit none
  public
  integer, private, parameter :: real_kind = kind(1.0d0)
  integer, private, parameter :: str_len = 1024

  ! where the photochem data is
  character(len=str_len) :: data_dir = "data"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! set DURING file read-in !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  real(real_kind) :: diurnal_fac
  real(real_kind) :: solar_zenith
  real(real_kind) :: trop_alt 
  integer :: trop_ind
  
  ! Radiative tranfer
  real(real_kind), allocatable :: photon_flux(:) ! (nw) photonz


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! set AFTER file read-in !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: neqs
  real(real_kind), allocatable :: temperature(:)
  real(real_kind), allocatable :: z(:)
  real(real_kind), allocatable :: dz(:)
  real(real_kind), allocatable :: edd(:)
  real(real_kind), allocatable :: grav(:)
  real(real_kind), allocatable, target :: usol_init(:,:)
  
  real(real_kind), allocatable :: xs_x_qy(:,:,:)
  
  ! other
  real(real_kind) :: epsj = 1.d-9
  logical :: verbose = .true.
  
end module