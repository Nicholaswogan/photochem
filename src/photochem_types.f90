
module photochem_types ! make a giant IO object
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long, c_ptr, c_null_ptr
  use photochem_const, only: real_kind, str_len, s_str_len, m_str_len
  
  use linear_interpolation_module, only: linear_interp_2d
  use fsundials_nvector_mod, only: N_Vector
  use fsundials_matrix_mod, only: SUNMatrix
  use fsundials_linearsolver_mod, only: SUNLinearSolver
  implicit none
  private
  
  public :: PhotochemData, PhotochemVars, PhotochemWrk
  public :: ProductionLoss, AtomConservation, ThermodynamicData
  
  type :: ProductionLoss
    real(real_kind), allocatable :: production(:,:)
    real(real_kind), allocatable :: loss(:,:)
    real(real_kind), allocatable :: integrated_production(:)
    real(real_kind), allocatable :: integrated_loss(:)
    character(len=m_str_len), allocatable :: production_rx(:)
    character(len=m_str_len), allocatable :: loss_rx(:)
  end type
  
  type :: XsectionData
    integer :: n_temps
    real(real_kind), allocatable :: xs(:,:) ! (n_temps, nw)
    real(real_kind), allocatable :: xs_temps(:) ! (n_temps)
  end type
  
  type :: ParticleXsections
    logical :: ThereIsData
    real(real_kind), allocatable :: w0(:,:) ! (nz,nw) or (nrad_file, nw)
    real(real_kind), allocatable :: qext(:,:)
    real(real_kind), allocatable :: gt(:,:)
  end type
  
  type :: ThermodynamicData
    integer :: dtype ! shomate = 1
    integer :: ntemps
    real(real_kind), allocatable :: temps(:)
    real(real_kind), allocatable :: data(:,:)
  end type
  
  type :: AtomConservation
    real(real_kind) :: in_surf
    real(real_kind) :: in_top
    real(real_kind) :: in_dist
    real(real_kind) :: out_surf
    real(real_kind) :: out_top
    real(real_kind) :: out_rain
    real(real_kind) :: out_other
    real(real_kind) :: net
    real(real_kind) :: factor
  end type
  
  type :: PhotochemData
    ! PhotochemData contains information that is never changed
    ! after file read-in
    
    ! molecules
    integer :: nq ! number of gases + particles which evolve over time from integration
    integer :: ng_1 ! index of first gas
    integer :: nll ! number of long-lived gas molecules
    integer :: nsl ! number of short-lived gas moleules. Short lived abundances are calculated
    ! assuming chemical equilibrium
    integer :: ng  ! number of gases
    integer :: nsp ! total number of species (nq + nsl + 1)
    integer :: natoms ! number of atoms
    integer :: kd, kl, ku ! not read in. It is nq + nq + 1 (diagonal width of jacobian)
    integer :: lda ! not read in. It is nq + nq + nq + 1. leading dimension of array which stores jacobian
    character(len=s_str_len), allocatable :: atoms_names(:) ! (natoms)
    real(real_kind), allocatable :: atoms_mass(:) ! g/mol (natoms)
    real(real_kind), allocatable :: atoms_redox(:) !  (natoms)
    character(len=s_str_len), allocatable :: SL_names(:) ! (nsl)
    character(len=s_str_len), allocatable :: species_names(:) ! (nsp+2) + 2 for hv and M
    integer, allocatable :: species_composition(:,:) ! (natoms, nsp+2)
    real(real_kind), allocatable :: species_mass(:) ! (nsp)
    real(real_kind), allocatable :: species_redox(:) ! (nsp)
    type(ThermodynamicData), allocatable :: thermo_data(:) ! (ng)
    real(real_kind), allocatable :: henry_data(:,:) ! (2, nsp).
    ! henry_data(:,i) = [A, B], and [mol/(kg * Pa)] = A*exp(B*(1.d0/298.15d0 - 1.d0/T))
    
    ! particles
    logical :: there_are_particles
    integer :: np ! number of particles
    integer :: npq ! number of particle equations. for now nq = npq.
    character(len=s_str_len), allocatable :: particle_names(:) ! np
    integer, allocatable :: particle_formation_method(:) ! np. 1 == saturation, 2 == reaction
    real(real_kind), allocatable :: particle_density(:) ! np (g/cm3)
    type(linear_interp_2d) :: H2SO4_sat ! interpolator for H2SO4 saturation, which depends on T and H2O.
    integer, allocatable :: particle_sat_type(:) ! np, 1 == arrhenius, 2 == H2SO4
    real(real_kind), allocatable :: particle_sat_params(:,:) ! (3,np)
    character(len=s_str_len), allocatable :: particle_gas_phase(:) ! (np). gas phase of particle. 
    ! Only for saturation particles
    integer, allocatable :: particle_gas_phase_ind(:) ! np. index of gas phase of particle
    character(len=s_str_len), allocatable :: particle_optical_prop(:) ! (np)
    integer, allocatable :: particle_optical_type(:) ! (np) 1 == mie, 2 == fractal
    
    ! reactions
    logical :: reverse ! True if there are reverse reactions
    integer :: nrF ! number of forward reactions
    integer :: nrR ! number of reverse reactions
    integer :: nrT ! number of total reactions
    integer :: max_num_reactants
    integer :: max_num_products
    character(len=m_str_len), allocatable :: reaction_equations(:)
    character(len=s_str_len), allocatable :: reactants_names(:,:)
    character(len=s_str_len), allocatable :: products_names(:,:)
    ! reactants_sp_inds(1:nreactants(i),i) are the species indicies for reaction i
    integer, allocatable :: reactants_sp_inds(:,:) ! (max_num_reactants, nrT)
    ! same idea as reactants_sp_inds
    integer, allocatable :: products_sp_inds(:,:) ! (max_num_products, nrT)
    integer, allocatable :: nreactants(:) ! (nrT) number of reactants
    integer, allocatable :: nproducts(:) ! (nrT) number of products
    integer, allocatable :: reverse_info(:) ! (nrT) indexs between forward and reverse reactions
    integer, allocatable :: rxtypes(:) ! (nrT) 0 is photolysis, 1 is elementary, 2 is three-body, 3 is falloff
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
    integer :: nw ! number of wavelength bins
    real(real_kind), allocatable :: wavl(:) ! (nw+1) wavelength bins in nm
    type(XsectionData), allocatable :: xs_data(:) ! (kj)
    integer :: nray ! number of species with rayleigh scattering
    real(real_kind), allocatable :: sigray(:,:) ! (len(raynums), nw)
    integer, allocatable :: raynums(:) ! species number of rayleigh species
    
    ! particle radiative transfer
    integer :: nrad_file
    real(real_kind), allocatable  :: radii_file(:,:) ! particle radii in optical data files
    ! We use array of types for particle xs because we want the option
    ! to exclude optical properties, but not take up a ton of useless memory.
    ! So some elements of this array have nothing in it.
    type(ParticleXsections), allocatable :: part_xs_file(:) ! np in length
    
    ! initial conditions  
    integer :: nzf ! number of atmospheric layers in file
    real(real_kind), allocatable :: z_file(:) ! (nzf) cm
    real(real_kind), allocatable :: T_file(:) ! (nzf) K
    real(real_kind), allocatable :: edd_file(:) ! (nzf) cm2/s
    real(real_kind), allocatable :: usol_file(:,:) ! (nq,nzf) mixing ratios
    real(real_kind), allocatable :: particle_radius_file(:,:) ! (np,nzf) cm
    
    ! settings
    logical :: regular_grid ! True of wavelength grid is evenly spaced
    real(real_kind) :: lower_wavelength ! nm
    real(real_kind) :: upper_wavelength ! nm
    character(len=str_len) :: grid_file ! filename of grid file. Only if regular_grid == False
    logical :: back_gas ! True if background gas is used
    character(len=str_len) :: back_gas_name ! Normally N2, but can be most any gas.
    real(real_kind) :: back_gas_mu ! g/mol
    integer :: back_gas_ind
    real(real_kind) :: planet_mass ! grams
    real(real_kind) :: planet_radius ! cm
    logical :: fix_water_in_trop ! True if fixing water in troposphere
    integer :: LH2O ! index of H2O
    logical :: stratospheric_cond ! True if water should condense out of the atmosphere
    logical :: diff_H_escape ! True of diffusion limited H escape
    integer :: LH2 ! H2 index
    integer :: LH ! H index
    
  end type
  
  type :: PhotochemVars
    ! PhotochemVars contains information that can change between
    ! different photochemical integrations, without reading in new
    ! files.
    
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
    integer, allocatable :: upperboundcond(:) ! 0 or 2
    real(real_kind), allocatable :: upper_veff(:)
    real(real_kind), allocatable :: upper_flux(:)
    
    ! Atmospheres structure
    real(real_kind) :: bottom_atmos ! cm
    real(real_kind) :: top_atmos ! cm
    integer :: nz ! number of vertical layers
    real(real_kind) :: surface_pressure ! bars
    real(real_kind) :: surface_albedo
    real(real_kind) :: diurnal_fac ! normally 0.5 cuz planets spin around.
    real(real_kind) :: solar_zenith 
    real(real_kind) :: trop_alt ! cm (only for fix_water_in_trop == true)
    integer :: trop_ind ! index of troposphere (only for fix_water_in_trop == true)
    logical :: use_manabe ! use manabe formula
    real(real_kind) :: relative_humidity ! relative humidity if no manabe
    logical :: gas_rainout ! True if gas rains out
    real(real_kind) :: H2O_condensation_rate(3) 
    
    ! Radiative tranfer
    real(real_kind), allocatable :: photon_flux(:) ! (nw) photonz
    ! for scaling photon flux for different planets in a solar system
    real(real_kind) :: photon_scale_factor 
    
    ! particles
    ! condensation rate of particles
    real(real_kind), allocatable :: condensation_rate(:,:)
    
    ! switch for dealing with H2O if not read in. in some cases
    logical :: no_water_profile

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set AFTER file read-in !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: neqs ! number of equations nq*nz
    real(real_kind), allocatable :: temperature(:) ! (nz) K
    real(real_kind), allocatable :: z(:) ! (nz) cm
    real(real_kind), allocatable :: dz(:) ! (nz) cm
    real(real_kind), allocatable :: edd(:) ! (nz) cm2/s
    real(real_kind), allocatable :: grav(:) ! (nz) cm/s2
    real(real_kind), allocatable :: usol_init(:,:) ! (nq,nz) mixing ratio
    real(real_kind), allocatable :: particle_radius(:,:) ! (np,nz) cm
    real(real_kind), allocatable :: xs_x_qy(:,:,:) ! (nz,kj,nw) cm2/molecule
    type(ParticleXsections), allocatable :: particle_xs(:) ! (np) cm2/molecule
    real(real_kind), allocatable :: gibbs_energy(:,:) ! (nz,ng) Joules/mol
    
    ! output
    logical :: at_photo_equilibrium = .false.
    real(real_kind), allocatable :: usol_out(:,:)
    
    ! other 
    real(c_double) :: rtol = 1.d-3 ! integration relative tolerance
    real(c_double) :: atol = 1.d-25 ! integration absolute tolerance
    integer :: mxsteps = 10000 ! max number of steps before integrator will give up.
    ! seconds. atomsphere considered in equilibrium if integrations reaches this time.
    real(real_kind) :: equilibrium_time = 1.d17 
    real(c_double) :: initial_dt = 1.d-6 ! intial timestep size (seconds)
    integer(c_int) :: max_err_test_failures = 15 
    integer(c_int) :: max_order = 5
    real(real_kind) :: epsj = 1.d-9 ! perturbation for jacobian calculation
    integer :: verbose = 1 ! 0 == no printing. 1 == some printing. 2 == bunch of printing.
  end type
  
  type :: PhotochemWrk
    ! PhotochemWrk are work variables that change
    ! through the course of an integration
    
    ! used in cvode
    integer(c_long) :: nsteps_previous = -10
    type(c_ptr) :: cvode_mem = c_null_ptr
    real(c_double) :: tcur(1)
    type(N_Vector), pointer :: sunvec_y ! sundials vector
    type(SUNMatrix), pointer :: sunmat
    type(SUNLinearSolver), pointer :: sunlin
    real(c_double), allocatable :: yvec(:)
    
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
    real(real_kind), allocatable :: amean_grd(:,:) ! (nz,nw)
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
      deallocate(self%amean_grd)
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
    allocate(self%amean_grd(nz,nw))
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


