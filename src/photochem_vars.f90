

module photochem_vars
  use, intrinsic :: iso_c_binding
  implicit none
  ! public
  integer, private, parameter :: real_kind = kind(1.0d0)
  integer, private, parameter :: str_len = 1024

  ! where the photochem data is
  character(len=str_len) :: data_dir = "data"
  ! the name of the xsections folder
  character(len=str_len) :: xs_folder_name = "xsections"

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
  logical :: use_manabe ! use manabe formula
  real(real_kind) :: relative_humidity ! relative humidity if no manabe
  
  ! Radiative tranfer
  real(real_kind), allocatable :: photon_flux(:) ! (nw) photonz
  real(real_kind) :: photon_scale_factor
  
  ! particles
  ! condensation rate of particles
  real(real_kind), allocatable :: condensation_rate(:)
  
  ! switch for dealing with H2O if not read in.
  logical :: no_water_profile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! set AFTER file read-in !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: neqs
  ! integer :: nqL
  real(real_kind), allocatable :: temperature(:)
  real(real_kind), allocatable :: z(:)
  real(real_kind), allocatable :: dz(:)
  real(real_kind), allocatable :: edd(:)
  real(real_kind), allocatable :: grav(:)
  real(real_kind), allocatable, target :: usol_init(:,:)
  real(real_kind), allocatable :: particle_radius(:,:)
  real(real_kind), allocatable :: xs_x_qy(:,:,:)
  real(real_kind), allocatable :: w0_particles(:,:,:)
  real(real_kind), allocatable :: qext_particles(:,:,:)
  real(real_kind), allocatable :: gt_particles(:,:,:)
  
  ! output
  logical :: at_photo_equilibrium = .false.
  real(real_kind), allocatable, target :: usol_out(:,:)
  
  ! other
  real(real_kind) :: equilibrium_time = 1.d17
  real(c_double) :: initial_dt = 1.d-6
  integer(c_int) :: max_err_test_failures = 15
  integer(c_int) :: max_order = 5
  logical :: use_fast_jacobian = .true.
  real(real_kind) :: epsj = 1.d-9
  integer :: verbose = 1
  
end module