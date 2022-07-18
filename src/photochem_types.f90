
module photochem_types ! make a giant IO object
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long, c_ptr, c_null_ptr
  use photochem_const, only: dp, str_len, s_str_len, m_str_len
  
  use linear_interpolation_module, only: linear_interp_2d
  use fsundials_nvector_mod, only: N_Vector
  use fsundials_matrix_mod, only: SUNMatrix
  use fsundials_linearsolver_mod, only: SUNLinearSolver
  implicit none
  private
  
  public :: PhotoSettings, SettingsBC
  public :: PhotochemData, PhotochemVars, PhotochemWrk, PhotochemWrkEvo
  public :: ProductionLoss, AtomConservation, ThermodynamicData
  public :: Reaction, Efficiencies, BaseRate, PhotolysisRate
  public :: ElementaryRate, ThreeBodyRate, FalloffRate, ProdLoss
  
  !!!!!!!!!!!!!!!!
  !!! Settings !!!
  !!!!!!!!!!!!!!!!
  
  type :: SettingsCondensationRate
    real(dp) :: A
    real(dp) :: rhc
    real(dp) :: rh0
  end type
  
  type :: SettingsBC
    integer :: bc_type
    real(dp) :: vel
    real(dp) :: mix
    real(dp) :: flux
    real(dp) :: height
    real(dp) :: den
  end type
  
  type :: PhotoSettings
  
    ! atmosphere-grid
    real(dp) :: bottom
    real(dp) :: top
    integer :: nz
  
    ! photolysis-grid
    logical :: regular_grid
    real(dp) :: lower_wv
    real(dp) :: upper_wv
    integer :: nw
    character(:), allocatable :: grid_file
    real(dp) :: photon_scale_factor 
  
    ! planet
    character(:), allocatable :: back_gas_name
    real(dp), allocatable :: P_surf
    real(dp) :: planet_mass
    real(dp) :: planet_radius
    real(dp) :: surface_albedo
    real(dp) :: diurnal_fac
    real(dp) :: solar_zenith
    logical :: diff_H_escape
    integer :: default_lowerboundcond
    ! water
    logical :: fix_water_in_trop
    logical :: water_cond
    logical :: gas_rainout
    character(:), allocatable :: relative_humidity
    real(dp) :: rainfall_rate
    character(s_str_len), allocatable :: rainout_species(:)
    real(dp) :: trop_alt
    real(dp) :: H2O_condensation_rate(3)
  
    ! particles
    character(s_str_len), allocatable :: con_names(:)
    type(SettingsCondensationRate), allocatable :: con(:)
  
    ! boundary-conditions
    type(SettingsBC), allocatable :: ubcs(:)
    type(SettingsBC), allocatable :: lbcs(:)
    character(s_str_len), allocatable :: sp_names(:)
    character(s_str_len), allocatable :: sp_types(:)
    logical, allocatable :: only_eddy(:) 
    
    integer :: nsl
    character(s_str_len), allocatable :: SL_names(:)
    
  end type
  
  !!!!!!!!!!!!!!!!!
  !!! Utilities !!!
  !!!!!!!!!!!!!!!!!
  
  type :: AtomConservation
    real(dp) :: in_surf
    real(dp) :: in_top
    real(dp) :: in_dist
    real(dp) :: in_other
    real(dp) :: out_surf
    real(dp) :: out_top
    real(dp) :: out_rain
    real(dp) :: out_other
    real(dp) :: net
    real(dp) :: factor
  end type
  
  type :: ProductionLoss
    real(dp), allocatable :: production(:,:)
    real(dp), allocatable :: loss(:,:)
    real(dp), allocatable :: integrated_production(:)
    real(dp), allocatable :: integrated_loss(:)
    character(len=m_str_len), allocatable :: production_rx(:)
    character(len=m_str_len), allocatable :: loss_rx(:)
  end type
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! XS and thermodynamic data !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type :: XsectionData
    integer :: n_temps
    real(dp), allocatable :: xs(:,:) ! (n_temps, nw)
    real(dp), allocatable :: xs_temps(:) ! (n_temps)
  end type
  
  type :: ParticleXsections
    logical :: ThereIsData
    real(dp), allocatable :: w0(:,:) ! (nz,nw) or (nrad_file, nw)
    real(dp), allocatable :: qext(:,:)
    real(dp), allocatable :: gt(:,:)
  end type
  
  type :: ThermodynamicData
    integer :: dtype ! shomate = 1
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type

  !!!!!!!!!!!!!!!!!
  !!! Reactions !!!
  !!!!!!!!!!!!!!!!!
  
  type :: Efficiencies
    integer :: n_eff ! number of efficiencies
    real(dp) :: def_eff ! default efficiency
    real(dp), allocatable :: efficiencies(:) ! 3-body efficiencies
    integer, allocatable :: eff_sp_inds(:) ! species indices for each efficiency
  end type
  
  type, abstract :: BaseRate
    integer :: rxtype 
    ! 0 is PhotolysisRate, 1 is ElementaryRate, 2 is ThreeBodyRate, 3 is FalloffRate
  end type
  
  type, extends(BaseRate) :: PhotolysisRate 
    ! no rate parameters. calculated via radiative transfer
  end type
  
  type, extends(BaseRate) :: ElementaryRate
    ! rate = A*T^b*exp(-Ea/T)
    real(dp) :: A
    real(dp) :: b
    real(dp) :: Ea
  end type
  
  type, extends(ElementaryRate) :: ThreeBodyRate
    type(Efficiencies) :: eff
  end type
  
  type, extends(BaseRate) :: FalloffRate
    real(dp) :: A0
    real(dp) :: b0
    real(dp) :: Ea0
    real(dp) :: Ainf
    real(dp) :: binf
    real(dp) :: Eainf
    
    integer :: falloff_type
    real(dp), allocatable :: A_T
    real(dp), allocatable :: T1
    real(dp), allocatable :: T2
    real(dp), allocatable :: T3
    
    type(Efficiencies) :: eff
  end type
  
  type :: Reaction
    integer :: nreact ! number of reactants
    integer :: nprod ! number of projects
    integer, allocatable :: react_sp_inds(:) ! (nreact) species indexes for reactants
    integer, allocatable :: prod_sp_inds(:) ! (nprod) species indexes for products
    integer, allocatable :: reverse_info ! if a reversed reaction, 
                                         ! then it is the index of the forward
    class(BaseRate), allocatable :: rp ! rate parameters
  end type
  
  type :: ProdLoss
    integer :: nump ! size(iprod)
    integer :: numl ! size(iloss)
    integer, allocatable :: iprod(:) ! reaction #s of production mechanism for sp
    integer, allocatable :: iloss(:) ! reaction #s of loss mechanism for sp
  end type
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Data, Vars, and Wrk !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: PhotochemData
    ! PhotochemData contains information that is never changed
    ! after file read-in
    
    integer :: natoms ! number of atoms
    character(len=s_str_len), allocatable :: atoms_names(:) ! (natoms)
    real(dp), allocatable :: atoms_mass(:) ! g/mol (natoms)
    real(dp), allocatable :: atoms_redox(:) !  (natoms)
    
    ! species
    ! Organization is as follows
    ! [       nsp       ]
    ! [   nq   + nsl + 1]
    ! [np + ng + nsl + 1]
    ! |_______|
    !     |
    ! Only np + ng = nq evolve through time. nsl are assumed to be in equilibrium. +1 is background gas.
    integer :: nq ! number of gases + particles which evolve over time from integration
    integer :: ng_1 ! index of first gas
    integer :: nll ! number of long-lived gas molecules
    integer :: nsl ! number of short-lived gas molecules
    integer :: ng  ! number of gases
    integer :: nsp ! total number of species (nq + nsl + 1)
    integer :: kd, kl, ku ! not read in. It is nq + nq + 1 (diagonal width of jacobian)
    integer :: lda ! not read in. It is nq + nq + nq + 1. leading dimension of array which stores jacobian
    character(len=s_str_len), allocatable :: SL_names(:) ! (nsl)
    character(len=s_str_len), allocatable :: species_names(:) ! (nsp+2) + 2 for hv and M
    integer, allocatable :: species_composition(:,:) ! (natoms, nsp+2)
    real(dp), allocatable :: species_mass(:) ! (nsp)
    real(dp), allocatable :: species_redox(:) ! (nsp)
    type(ThermodynamicData), allocatable :: thermo_data(:) ! (ng)
    real(dp), allocatable :: henry_data(:,:) ! (2, nsp).
    ! henry_data(:,i) = [A, B], and [mol/(kg * Pa)] = A*exp(B*(1.0_dp/298.15e0_dp - 1.0_dp/T))
    
    ! particles
    logical :: there_are_particles
    integer :: np ! number of particles
    integer :: npq ! number of particle equations. for now nq = npq.
    character(len=s_str_len), allocatable :: particle_names(:) ! np
    integer, allocatable :: particle_formation_method(:) ! np. 1 == saturation, 2 == reaction
    real(dp), allocatable :: particle_density(:) ! np (g/cm3)
    type(linear_interp_2d) :: H2SO4_sat ! interpolator for H2SO4 saturation, which depends on T and H2O.
    integer, allocatable :: particle_sat_type(:) ! np, 1 == arrhenius, 2 == H2SO4
    real(dp), allocatable :: particle_sat_params(:,:) ! (3,np)
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
    type(Reaction), allocatable :: rx(:) ! (nrT) array of reaction objects
    character(len=m_str_len), allocatable :: reaction_equations(:) ! (nrT)
    type(ProdLoss), allocatable :: pl(:) ! (nsp) reactions producing and destroying each species
    integer :: kj ! number of photolysis reactions
    integer, allocatable :: photonums(:) ! (kj) the reaction number of each photolysis reaction

    ! raditative transfer
    integer :: nw ! number of wavelength bins
    real(dp), allocatable :: wavl(:) ! (nw+1) wavelength bins in nm
    type(XsectionData), allocatable :: xs_data(:) ! (kj)
    integer :: nray ! number of species with rayleigh scattering
    real(dp), allocatable :: sigray(:,:) ! (len(raynums), nw)
    integer, allocatable :: raynums(:) ! species number of rayleigh species
    
    ! particle radiative transfer
    integer :: nrad_file
    real(dp), allocatable  :: radii_file(:,:) ! particle radii in optical data files
    ! We use array of types for particle xs because we want the option
    ! to exclude optical properties, but not take up a ton of useless memory.
    ! So some elements of this array have nothing in it.
    type(ParticleXsections), allocatable :: part_xs_file(:) ! np in length
    
    ! initial conditions  
    integer :: nzf ! number of atmospheric layers in file
    real(dp), allocatable :: z_file(:) ! (nzf) cm
    real(dp), allocatable :: T_file(:) ! (nzf) K
    real(dp), allocatable :: edd_file(:) ! (nzf) cm2/s
    real(dp), allocatable :: den_file(:) ! (nzf) molecules/cm2
    real(dp), allocatable :: mix_file(:,:) ! (nq,nzf) mixing ratios
    real(dp), allocatable :: particle_radius_file(:,:) ! (np,nzf) cm
    
    ! settings
    logical :: regular_grid ! True of wavelength grid is evenly spaced
    real(dp) :: lower_wavelength ! nm
    real(dp) :: upper_wavelength ! nm
    character(:), allocatable :: grid_file ! filename of grid file. Only if regular_grid == False
    logical :: back_gas ! True if background gas is used
    character(:), allocatable :: back_gas_name ! Normally N2, but can be most any gas.
    real(dp), allocatable :: back_gas_mu ! g/mol
    integer, allocatable :: back_gas_ind
    real(dp) :: planet_mass ! grams
    real(dp) :: planet_radius ! cm
    logical :: fix_water_in_trop ! True if fixing water in troposphere
    integer :: LH2O ! index of H2O
    logical :: water_cond ! True if water should condense out of the atmosphere
    logical :: gas_rainout ! True if gas rains out
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
    real(dp), allocatable :: lower_vdep(:)
    real(dp), allocatable :: lower_flux(:)
    real(dp), allocatable :: lower_dist_height(:)
    real(dp), allocatable :: lower_fix_mr(:)
    real(dp), allocatable :: lower_fix_den(:)
    integer, allocatable :: upperboundcond(:) ! 0 or 2
    real(dp), allocatable :: upper_veff(:)
    real(dp), allocatable :: upper_flux(:)
    logical, allocatable :: only_eddy(:) ! True if only use eddy
    
    ! Atmospheres structure
    real(dp) :: bottom_atmos ! cm
    real(dp) :: top_atmos ! cm
    integer :: nz ! number of vertical layers
    real(dp), allocatable :: surface_pressure ! bars
    real(dp) :: surface_albedo
    real(dp) :: diurnal_fac ! normally 0.5 cuz planets spin around.
    real(dp) :: solar_zenith 
    real(dp) :: trop_alt ! cm (only for fix_water_in_trop == true or gas_rainout == true)
    real(dp) :: rainfall_rate ! relative to modern Earth's average rainfall rate of 1.1e17 molecules/cm2/s
    integer :: trop_ind ! index of troposphere (only for fix_water_in_trop == true or gas_rainout == true)
    logical :: use_manabe ! use manabe formula
    real(dp) :: relative_humidity ! relative humidity if no manabe
    real(dp) :: H2O_condensation_rate(3) 
    
    ! Radiative tranfer
    real(dp), allocatable :: photon_flux(:) ! (nw) photonz
    ! for scaling photon flux for different planets in a solar system
    real(dp) :: photon_scale_factor 
    
    ! particles
    ! condensation rate of particles
    real(dp), allocatable :: condensation_rate(:,:) ! (3,np)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set AFTER file read-in !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: neqs ! number of equations nq*nz
    real(dp), allocatable :: temperature(:) ! (nz) K
    real(dp), allocatable :: z(:) ! (nz) cm
    real(dp), allocatable :: dz(:) ! (nz) cm
    real(dp), allocatable :: edd(:) ! (nz) cm2/s
    real(dp), allocatable :: grav(:) ! (nz) cm/s2
    real(dp), allocatable :: usol_init(:,:) ! (nq,nz) mixing ratio (EvoAtmosphere) or densities (EvoAtmosphere).
    real(dp), allocatable :: particle_radius(:,:) ! (np,nz) cm
    real(dp), allocatable :: xs_x_qy(:,:,:) ! (nz,kj,nw) cm2/molecule
    type(ParticleXsections), allocatable :: particle_xs(:) ! (np) cm2/molecule
    real(dp), allocatable :: gibbs_energy(:,:) ! (nz,ng) Joules/mol
    
    ! output
    logical :: at_photo_equilibrium = .false.
    real(dp), allocatable :: usol_out(:,:)
    
    ! other 
    real(c_double) :: rtol = 1.d-3 ! integration relative tolerance
    real(c_double) :: atol = 1.d-27 ! integration absolute tolerance
    integer :: mxsteps = 10000 ! max number of steps before integrator will give up.
    ! seconds. atomsphere considered in equilibrium if integrations reaches this time.
    real(dp) :: equilibrium_time = 1.e17_dp 
    real(c_double) :: initial_dt = 1.d-6 ! intial timestep size (seconds)
    integer(c_int) :: max_err_test_failures = 15 
    integer(c_int) :: max_order = 5
    real(dp) :: epsj = 1.d-9 ! perturbation for jacobian calculation
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
    real(dp), allocatable :: usol(:,:) ! (nq,nz)
    real(dp), allocatable :: densities(:,:) ! (nsp+1,nz)
    real(dp), allocatable :: density(:) ! (nz)
    real(dp), allocatable :: rx_rates(:,:) ! (nz,nrT)
    real(dp), allocatable :: mubar(:) ! (nz)
    real(dp), allocatable :: pressure(:) ! (nz)
    real(dp), allocatable :: H2O_rh(:) ! (nz)
    real(dp), allocatable :: H2O_sat_mix(:) ! (nz)
    real(dp), allocatable :: prates(:,:) ! (nz,kj)
    real(dp), allocatable :: surf_radiance(:) ! (nw)
    real(dp), allocatable :: amean_grd(:,:) ! (nz,nw)
    real(dp), allocatable :: optical_depth(:,:) ! (nz,nw)
    real(dp), allocatable :: upper_veff_copy(:) ! (nq)
    real(dp), allocatable :: lower_vdep_copy(:) ! (nq)
    real(dp), allocatable :: xp(:) ! (nz)
    real(dp), allocatable :: xl(:) ! (nz)
    ! diffusion and H escape
    real(dp), allocatable :: DU(:,:) ! (nq,nz)
    real(dp), allocatable :: DD(:,:) ! (nq,nz)
    real(dp), allocatable :: DL(:,:) ! (nq,nz)
    real(dp), allocatable :: ADU(:,:) ! (nq,nz)
    real(dp), allocatable :: ADL(:,:) ! (nq,nz)
    real(dp), allocatable :: ADD(:,:) ! (nq,nz)
    real(dp) :: VH2_esc
    real(dp) :: VH_esc
    ! other
    real(dp), allocatable :: sum_usol(:) ! (nz)
    real(dp) :: surface_scale_height
    real(dp), allocatable :: wfall(:,:)
    real(dp), allocatable :: gas_sat_den(:,:)
    real(dp), allocatable :: molecules_per_particle(:,:)
    real(dp), allocatable :: rainout_rates(:,:)
    ! end used in prep_all_background_gas
    
    
  contains
    procedure :: init => init_PhotochemWrk
  end type

  type, extends(PhotochemWrk) :: PhotochemWrkEvo
    real(dp), allocatable :: mix(:,:) ! (nq,nz) mixing ratio.

  contains
    procedure :: init => init_PhotochemWrkEvo

  end type
  
contains

  subroutine init_PhotochemWrkEvo(self, nsp, np, nq, nz, nrT, kj, nw, trop_ind)
    class(PhotochemWrkEvo), intent(inout) :: self
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw, trop_ind

    call init_PhotochemWrk(self, nsp, np, nq, nz, nrT, kj, nw, trop_ind)

    if (allocated(self%mix)) then
      deallocate(self%mix)
    endif

    allocate(self%mix(nq,nz))

  end subroutine
 
  subroutine init_PhotochemWrk(self, nsp, np, nq, nz, nrT, kj, nw, trop_ind)
    class(PhotochemWrk), intent(inout) :: self
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw, trop_ind
    
    if (allocated(self%usol)) then
      deallocate(self%usol)
      deallocate(self%mubar)
      deallocate(self%pressure)
      deallocate(self%density)
      deallocate(self%H2O_rh)
      deallocate(self%H2O_sat_mix)
      deallocate(self%densities)
      deallocate(self%rx_rates)
      deallocate(self%prates)
      deallocate(self%surf_radiance)
      deallocate(self%amean_grd)
      deallocate(self%optical_depth)
      deallocate(self%xp)
      deallocate(self%xl)
      deallocate(self%DU)
      deallocate(self%DD)
      deallocate(self%DL)
      deallocate(self%ADU)
      deallocate(self%ADL)
      deallocate(self%ADD)
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
    allocate(self%H2O_rh(trop_ind))
    allocate(self%H2O_sat_mix(nz))
    allocate(self%densities(nsp+1,nz))
    allocate(self%rx_rates(nz,nrT))
    allocate(self%prates(nz,kj))
    allocate(self%surf_radiance(nw))
    allocate(self%amean_grd(nz,nw))
    allocate(self%optical_depth(nz,nw))
    allocate(self%xp(nz))
    allocate(self%xl(nz))
    allocate(self%DU(nq,nz))
    allocate(self%DD(nq,nz))
    allocate(self%DL(nq,nz))
    allocate(self%ADU(nq,nz))
    allocate(self%ADL(nq,nz))
    allocate(self%ADD(nq,nz))
    allocate(self%upper_veff_copy(nq))
    allocate(self%lower_vdep_copy(nq))
    allocate(self%sum_usol(nz))
    allocate(self%wfall(np,nz))
    allocate(self%gas_sat_den(np,nz))
    allocate(self%molecules_per_particle(np,nz))
    allocate(self%rainout_rates(nq,trop_ind))
  end subroutine
  
end module


