
module photochem_types ! make a giant IO object
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long, c_ptr
  use photochem_const, only: real_kind, str_len
  implicit none
  private
  
  public :: PhotochemData, PhotochemVars, PhotochemWrk
  
  type :: XsectionData
    integer :: n_temps
    real(real_kind), allocatable :: xs(:,:) ! (n_temps, nw)
    real(real_kind), allocatable :: xs_temps(:) ! (n_temps)
  end type
  
  type :: PhotochemData
    ! molecules
    integer :: nq
    integer :: ng_1
    integer :: nll
    integer :: nsl
    integer :: ng
    integer :: nsp
    integer :: natoms
    integer :: kd, kl, ku ! not read in. It is nq + nq + 1 (diagonal width of jacobian)
    integer :: lda ! not read in 
    character(len=8), allocatable :: atoms_names(:) 
    real(real_kind), allocatable :: atoms_mass(:) 
    character(len=20), allocatable :: SL_names(:) ! IN PHOTOSETTINGS
    character(len=15), allocatable :: species_names(:)
    integer, allocatable :: species_composition(:,:)
    real(real_kind), allocatable :: species_mass(:) 
    real(real_kind), allocatable :: thermo_data(:,:,:)
    real(real_kind), allocatable :: thermo_temps(:,:)
    real(real_kind), allocatable :: henry_data(:,:)
    
    ! particles
    logical :: there_are_particles
    integer :: np ! number of particles
    integer :: npq ! number of particle equations. for now nq = npq.
    character(len=15), allocatable :: particle_names(:) ! IN PHOTOMECHANISM
    integer, allocatable :: particle_formation_method(:) ! np
    real(real_kind), allocatable :: particle_density(:) ! np
    real(real_kind), allocatable :: particle_sat_params(:,:) ! 3, np
    character(len=15), allocatable :: particle_gas_phase(:) ! IN PHOTOMECHANISM
    integer, allocatable :: particle_gas_phase_ind(:) ! np
    character(len=50), allocatable :: particle_optical_prop(:) ! IN PHOTOMECHANISM
    integer, allocatable :: particle_optical_type(:) ! IN PHOTOMECHANISM
    
    ! reactions
    logical :: reverse
    integer :: nrF ! number of forward reactions
    integer :: nrR ! number of reverse reactions
    integer :: nrT ! number of total reactions
    integer :: max_num_reactants
    integer :: max_num_products
    character(len=8), allocatable :: reactants_names(:,:)
    character(len=8), allocatable :: products_names(:,:)
    integer, allocatable :: reactants_sp_inds(:,:) ! for getting species nums in reactions
    integer, allocatable :: products_sp_inds(:,:)
    integer, allocatable :: nreactants(:) ! number of reactants
    integer, allocatable :: nproducts(:) ! number of products
    integer, allocatable :: reverse_info(:) ! indexs between forward and reverse reactions
    integer, allocatable :: rxtypes(:) ! 0 is photolysis, 1 is elementary, 2 is three-body, 3 is falloff
    real(real_kind), allocatable :: rateparams(:,:) ! (10, nrF)
    real(real_kind), allocatable :: efficiencies(:,:) ! (maxval(num_efficient), nrF)
    integer, allocatable :: eff_sp_inds(:,:) ! (maxval(num_efficient), nrF)
    integer, allocatable :: num_efficient(:) ! number of efficiencies for each reaction
    real(real_kind), allocatable :: def_eff(:) ! default efficiency
    integer, allocatable :: falloff_type(:) ! type of falloff function (0 = none, 1 = Troe without T2,..)
    integer, allocatable :: nump(:) ! number of production mechanisms (rxns) for each sp
    integer, allocatable :: numl(:) ! number of loss mechanisms (rxns) for each sp
    integer, allocatable :: iprod(:,:) ! (nmax,nsp) returns reaction # of production mechanism for sp
    integer, allocatable :: iloss(:,:) ! (nmax,nsp) returns reaction # of loss mechanism for sp
    integer :: kj ! number of photolysis reactions
    integer, allocatable :: photonums(:) ! (kj) the reaction number of each photolysis reaction

    ! raditative transfer
    integer :: nw
    real(real_kind), allocatable :: wavl(:) ! (nw+1)
    type(XsectionData), allocatable :: xs_data(:) ! (kj)
    integer :: nray
    real(real_kind), allocatable :: sigray(:,:) ! (len(raynums), nw)
    integer, allocatable :: raynums(:) ! species number of rayleigh species
    
    ! particle radiative transfer
    integer :: nrad_file
    real(real_kind), allocatable  :: radii_file(:,:) 
    real(real_kind), allocatable  :: w0_file(:,:,:)
    real(real_kind), allocatable  :: qext_file(:,:,:) 
    real(real_kind), allocatable  :: g_file(:,:,:) 
    
    ! initial conditions  
    integer :: nzf
    real(real_kind), allocatable :: z_file(:)
    real(real_kind), allocatable :: T_file(:)
    real(real_kind), allocatable :: edd_file(:)
    real(real_kind), allocatable :: usol_file(:,:)
    real(real_kind), allocatable :: particle_radius_file(:,:)
    
    ! settings
    logical :: regular_grid ! for wavelength
    real(real_kind) :: lower_wavelength
    real(real_kind) :: upper_wavelength
    character(len=str_len) :: grid_file
    logical :: back_gas
    character(len=str_len) :: back_gas_name
    real(real_kind) :: back_gas_mu
    integer :: back_gas_ind
    real(real_kind) :: planet_mass
    real(real_kind) :: planet_radius
    logical :: fix_water_in_trop
    integer :: LH2O
    logical :: stratospheric_cond
    logical :: diff_H_escape
    integer :: LH2
    integer :: LH
    
  end type
  
  type :: PhotochemVars
    ! where the photochem data is
    character(len=str_len) :: data_dir = "../data"
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
    logical :: gas_rainout
    real(real_kind) :: relative_humidity_cold_trap
    real(real_kind) :: H2O_condensation_rate(2)
    
    ! Radiative tranfer
    real(real_kind), allocatable :: photon_flux(:) ! (nw) photonz
    real(real_kind) :: photon_scale_factor
    
    ! particles
    ! condensation rate of particles
    real(real_kind), allocatable :: condensation_rate(:,:)
    
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
    real(real_kind), allocatable :: usol_init(:,:)
    real(real_kind), allocatable :: particle_radius(:,:)
    real(real_kind), allocatable :: xs_x_qy(:,:,:)
    real(real_kind), allocatable :: w0_particles(:,:,:)
    real(real_kind), allocatable :: qext_particles(:,:,:)
    real(real_kind), allocatable :: gt_particles(:,:,:)
    
    ! output
    logical :: at_photo_equilibrium = .false.
    real(real_kind), allocatable :: usol_out(:,:)
    
    ! other 
    real(c_double) :: rtol = 1.d-3
    real(c_double) :: atol = 1.d-25
    integer :: mxsteps = 10000
    real(real_kind) :: equilibrium_time = 1.d17
    real(c_double) :: initial_dt = 1.d-6
    integer(c_int) :: max_err_test_failures = 15
    integer(c_int) :: max_order = 5
    logical :: use_fast_jacobian = .true.
    real(real_kind) :: epsj = 1.d-9
    integer :: verbose = 1
  end type
  
  type :: PhotochemWrk
    
    ! used in cvode
    integer(c_long) :: nsteps_previous = -10
    type(c_ptr) :: cvode_mem
    
    ! Used in prep_all_background_gas
    ! work arrays
    real(real_kind), allocatable :: usol(:,:) ! (nq,nz)
    real(real_kind), allocatable :: densities(:,:) ! (nsp+1,nz)
    real(real_kind), allocatable :: density(:) ! (nz)
    real(real_kind), allocatable :: rx_rates(:,:) ! (nz,nrT)
    real(real_kind), allocatable :: mubar(:) ! (nz)
    real(real_kind), allocatable :: pressure(:) ! (nz)
    real(real_kind), allocatable :: fH2O(:) ! (nz)
    real(real_kind), allocatable :: H2O_sat_mix(:) ! (nz)
    real(real_kind), allocatable :: prates(:,:) ! (nz,kj)
    real(real_kind), allocatable :: surf_radiance(:) ! (nw)
    real(real_kind), allocatable :: upper_veff_copy(:) ! (nq)
    real(real_kind), allocatable :: lower_vdep_copy(:) ! (nq)
    real(real_kind), allocatable :: xp(:) ! (nz)
    real(real_kind), allocatable :: xl(:) ! (nz)
    ! diffusion and H escape
    real(real_kind), allocatable :: DU(:,:) ! (nq,nz)
    real(real_kind), allocatable :: DD(:,:) ! (nq,nz)
    real(real_kind), allocatable :: DL(:,:) ! (nq,nz)
    real(real_kind), allocatable :: ADU(:,:) ! (nq,nz)
    real(real_kind), allocatable :: ADL(:,:) ! (nq,nz)
    real(real_kind) :: VH2_esc
    real(real_kind) :: VH_esc
    ! other
    real(real_kind), allocatable :: sum_usol(:) ! (nz)
    real(real_kind) :: surface_scale_height
    real(real_kind), allocatable :: wfall(:,:)
    real(real_kind), allocatable :: gas_sat_den(:,:)
    real(real_kind), allocatable :: molecules_per_particle(:,:)
    real(real_kind), allocatable :: rainout_rates(:,:)
    ! end used in prep_all_background_gas
    
    
  contains
    procedure :: init => init_PhotochemWrk
  end type
  
contains
 
  subroutine init_PhotochemWrk(self, nsp, np, nq, nz, nrT, kj, nw, trop_ind)
    class(PhotochemWrk), intent(inout) :: self
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw, trop_ind
    
    if (allocated(self%usol)) then
      deallocate(self%usol)
      deallocate(self%mubar)
      deallocate(self%pressure)
      deallocate(self%density)
      deallocate(self%fH2O)
      deallocate(self%H2O_sat_mix)
      deallocate(self%densities)
      deallocate(self%rx_rates)
      deallocate(self%prates)
      deallocate(self%surf_radiance)
      deallocate(self%xp)
      deallocate(self%xl)
      deallocate(self%DU)
      deallocate(self%DD)
      deallocate(self%DL)
      deallocate(self%ADU)
      deallocate(self%ADL)
      deallocate(self%upper_veff_copy)
      deallocate(self%lower_vdep_copy)
      deallocate(self%sum_usol)
      deallocate(self%wfall)
      deallocate(self%gas_sat_den)
      deallocate(self%molecules_per_particle)
      deallocate(self%rainout_rates)
    endif
    
    allocate(self%usol(nq,nz))
    allocate(self%mubar(nz))
    allocate(self%pressure(nz))
    allocate(self%density(nz))
    allocate(self%fH2O(trop_ind))
    allocate(self%H2O_sat_mix(nz))
    allocate(self%densities(nsp+1,nz))
    allocate(self%rx_rates(nz,nrT))
    allocate(self%prates(nz,kj))
    allocate(self%surf_radiance(nw))
    allocate(self%xp(nz))
    allocate(self%xl(nz))
    allocate(self%DU(nq,nz))
    allocate(self%DD(nq,nz))
    allocate(self%DL(nq,nz))
    allocate(self%ADU(nq,nz))
    allocate(self%ADL(nq,nz))
    allocate(self%upper_veff_copy(nq))
    allocate(self%lower_vdep_copy(nq))
    allocate(self%sum_usol(nz))
    allocate(self%wfall(np,nz))
    allocate(self%gas_sat_den(np,nz))
    allocate(self%molecules_per_particle(np,nz))
    allocate(self%rainout_rates(nq,trop_ind))
  end subroutine
  
end module








