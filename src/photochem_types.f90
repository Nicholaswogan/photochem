
module photochem_types ! make a giant IO object
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long, c_ptr, c_null_ptr
  use photochem_const, only: dp, str_len, s_str_len, m_str_len
  use clima_saturationdata, only: SaturationData
  use fsundials_nvector_mod, only: N_Vector
  use fsundials_matrix_mod, only: SUNMatrix
  use fsundials_linearsolver_mod, only: SUNLinearSolver
  implicit none
  private
  
  public :: PhotoSettings, SettingsBC
  public :: XsectionData
  public :: PhotochemData, PhotochemVars, PhotochemWrk, PhotochemWrkEvo
  public :: ProductionLoss, ThermodynamicData, CondensationParameters
  public :: Reaction, Efficiencies, BaseRate, PhotolysisRate, PressDependentRate, MultiArrheniusRate
  public :: ElementaryRate, ThreeBodyRate, FalloffRate, ProdLoss, ReverseRate
  public :: SundialsDataFinalizer
  public :: time_dependent_flux_fcn, time_dependent_rate_fcn, binary_diffusion_fcn
  
  !!!!!!!!!!!!!!!!
  !!! Settings !!!
  !!!!!!!!!!!!!!!!
  
  type :: SettingsBC
    integer :: bc_type
    real(dp) :: vel
    real(dp) :: mix
    real(dp) :: flux
    real(dp) :: height
    real(dp) :: den
    real(dp) :: press
  end type

  !> Condensation parameters
  type :: CondensationParameters
    real(dp) :: k_cond = 100.0_dp !! rate coefficient for condensation
    real(dp) :: k_evap = 10.0_dp !! rate coefficient for evaporation
    real(dp) :: RHc = 1.0_dp !! RH where condensation occurs
    real(dp) :: smooth_factor = 0.2_dp !! A factor that smooths condensation/evaporation 
                                       !! rate to prevents stiffness
  end type
  
  type :: SettingsParticle
    character(:), allocatable :: name
    type(CondensationParameters) :: params
  endtype

  type :: PhotoSettings
    character(:), allocatable :: filename
  
    ! atmosphere-grid
    real(dp) :: bottom
    character(:), allocatable :: top
    integer :: nz
  
    ! planet
    real(dp) :: planet_mass
    real(dp) :: planet_radius
    real(dp) :: surface_albedo
    real(dp) :: photon_scale_factor 
    real(dp) :: solar_zenith
    integer :: H_escape_type
    real(dp), allocatable :: H_escape_S1
    integer :: default_lowerboundcond
    ! initialization
    logical :: conserving_init
    ! climate
    logical :: evolve_climate
    ! water
    logical :: fix_water_in_trop
    logical :: water_cond
    logical :: gas_rainout
    character(:), allocatable :: relative_humidity
    real(dp) :: rainfall_rate
    character(s_str_len), allocatable :: rainout_species(:)
    real(dp) :: trop_alt

    ! particles
    type(SettingsParticle), allocatable :: particles(:)
  
    ! boundary-conditions
    type(SettingsBC), allocatable :: ubcs(:)
    type(SettingsBC), allocatable :: lbcs(:)
    character(s_str_len), allocatable :: sp_names(:)
    character(s_str_len), allocatable :: sp_types(:)
    logical, allocatable :: only_eddy(:) 
    
    integer :: nsl
    character(s_str_len), allocatable :: SL_names(:)
    
  end type

  interface
    module function create_PhotoSettings(filename, err) result(s)
      character(*), intent(in) :: filename
      character(:), allocatable, intent(out) :: err
      type(PhotoSettings) :: s
    end function
  end interface
  interface PhotoSettings
    module procedure :: create_PhotoSettings
  end interface
  
  !!!!!!!!!!!!!!!!!
  !!! Utilities !!!
  !!!!!!!!!!!!!!!!!

  abstract interface
    !> Custom binary diffusion function
    function binary_diffusion_fcn(mu_i, mubar, T) result(b)
      use iso_c_binding, only: c_double
      real(c_double), value, intent(in) :: mu_i !! molar weight of species i (g/mol)
      real(c_double), value, intent(in) :: mubar !! molar weight of background gas (g/mol)
      real(c_double), value, intent(in) :: T !! Temperature (K)
      real(c_double) :: b !! binary diffusion parameter of species i
                          !! with respect to the background gas (molecules cm^-1 s^1)
    end function
      
    !> Sets the time-dependent photon flux
    subroutine time_dependent_flux_fcn(tn,  nw, photon_flux)
      use iso_c_binding, only: c_double, c_int
      real(c_double), value, intent(in) :: tn
      integer(c_int), value, intent(in) :: nw
      real(c_double), intent(out) :: photon_flux(nw)
    end subroutine

    !> Sets a production or destruction rate for a molecule
    subroutine time_dependent_rate_fcn(tn, nz, rate)
      use iso_c_binding, only: c_double, c_int
      real(c_double), value, intent(in) :: tn !! time (s)
      integer(c_int), value, intent(in) :: nz !! number of atmospheric layers
      real(c_double), intent(out) :: rate(nz) !! molecules/cm^3/s (can be positive or negative)
    end subroutine
  end interface

  !> Container to make an array of functions for time-dependent production rates.
  type :: time_dependent_rate_fcns
    procedure(time_dependent_rate_fcn), nopass, pointer :: fcn => null()
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
    integer :: sp_ind !! Species index
    real(dp), allocatable :: xs(:) !! The cross section in cm^2/molecule (nw)
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
    integer :: n_eff !! number of efficiencies
    real(dp) :: def_eff !! default efficiency
    real(dp), allocatable :: efficiencies(:) !! 3-body efficiencies
    integer, allocatable :: eff_sp_inds(:) !! species indices for each efficiency
  end type
  
  type, abstract :: BaseRate
    integer :: rxtype !! Specifies the types of rate (see enums)
  end type

  type, extends(BaseRate) :: ReverseRate
    ! no rate parameters.
  end type
  
  type, extends(BaseRate) :: PhotolysisRate 
    ! no rate parameters. calculated via radiative transfer
  end type
  
  type, extends(BaseRate) :: ElementaryRate
    ! rate = A*T^b*exp(-Ea/T)
    real(dp) :: A !! pre-exponential factor (various units)
    real(dp) :: b !! temperature exponent (unitless)
    real(dp) :: Ea !! Activate energy (T)
  end type
  
  type, extends(BaseRate) :: ThreeBodyRate
    real(dp) :: A !! pre-exponential factor (various units)
    real(dp) :: b !! temperature exponent (unitless)
    real(dp) :: Ea !! Activate energy (T)
    type(Efficiencies) :: eff
  end type
  
  type, extends(BaseRate) :: FalloffRate
    real(dp) :: A0 !! For low-P rate constant
    real(dp) :: b0 !! For low-P rate constant
    real(dp) :: Ea0 !! For low-P rate constant
    real(dp) :: Ainf !! For high-P rate constant
    real(dp) :: binf !! For high-P rate constant
    real(dp) :: Eainf !! For high-P rate constant
    
    integer :: falloff_type !! Type of falloff
    real(dp), allocatable :: A_T !! Troe falloff parameter 
    real(dp), allocatable :: T1 !! Troe falloff parameter
    real(dp), allocatable :: T2 !! Troe falloff parameter
    real(dp), allocatable :: T3 !! Troe falloff parameter
    
    type(Efficiencies) :: eff
  end type
  
  type :: MultiArrheniusRate
    real(dp), allocatable :: A(:)
    real(dp), allocatable :: b(:)
    real(dp), allocatable :: Ea(:)
  end type

  type, extends(BaseRate) :: PressDependentRate
    real(dp), allocatable :: logP(:)
    type(MultiArrheniusRate), allocatable :: rate(:)
  end type
  
  type :: Reaction
    integer :: nreact !! number of reactants
    integer :: nprod !! number of products
    integer, allocatable :: react_sp_inds(:) !! (nreact) species indexes for reactants
    integer, allocatable :: prod_sp_inds(:) !! (nprod) species indexes for products
    integer, allocatable :: reverse_info !! if a reversed reaction, 
                                         !! then it is the index of the forward
    class(BaseRate), allocatable :: rp !! rate parameters
  end type
  
  type :: ProdLoss
    integer :: nump ! (iprod) number off production processes
    integer :: numl ! (iloss) number off loss processes
    integer, allocatable :: iprod(:) ! reaction #s of production mechanism for sp
    integer, allocatable :: iloss(:) ! reaction #s of loss mechanism for sp
  end type
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Data, Vars, and Wrk !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: PhotochemData
    ! PhotochemData contains information that is never changed
    ! after file read-in
    
    integer :: natoms !! number of atoms
    character(len=s_str_len), allocatable :: atoms_names(:) !! (natoms)
    real(dp), allocatable :: atoms_mass(:) !! g/mol (natoms)
    real(dp), allocatable :: atoms_redox(:) !!  (natoms)
    
    ! species
    ! Organization is as follows
    ! [     nsp     ]
    ! [   nq   + nsl]
    ! [np + ng + nsl]
    ! |_______|
    !     |
    ! Only np + ng = nq evolve through time. nsl are assumed to be in equilibrium.
    integer :: nq !! number of gases + particles which evolve over time from integration
    integer :: ng_1 !! index of first gas
    integer :: nll !! number of long-lived gas molecules
    integer :: nsl !! number of short-lived gas molecules
    integer :: ng  !! number of gases
    integer :: nsp !! total number of species (nq + nsl + 1)
    integer :: kd, kl, ku !! not read in. It is nq + nq + 1 (diagonal width of jacobian)
    integer :: lda !! not read in. It is nq + nq + nq + 1. leading dimension of array which stores jacobian
    character(len=s_str_len), allocatable :: SL_names(:) !! (nsl)
    character(len=s_str_len), allocatable :: species_names(:) !! (nsp+2) + 2 for hv and M
    integer, allocatable :: species_composition(:,:) !! (natoms, nsp+2)
    real(dp), allocatable :: species_mass(:) !! (nsp)
    real(dp), allocatable :: species_redox(:) !! (nsp)
    type(ThermodynamicData), allocatable :: thermo_data(:) !! (ng)
    real(dp), allocatable :: henry_data(:,:) !! (2, nsp).
    ! henry_data(:,i) = [A, B], and [mol/(kg * Pa)] = A*exp(B*(1.0_dp/298.15e0_dp - 1.0_dp/T))
    
    ! particles
    logical :: there_are_particles
    integer :: np !! number of particles
    integer :: npq !! number of particle equations. for now nq = npq.
    character(len=s_str_len), allocatable :: particle_names(:) !! np
    integer, allocatable :: particle_formation_method(:) !! np. 1 == saturation, 2 == reaction
    real(dp), allocatable :: particle_density(:) !! np (g/cm3)
    type(SaturationData), allocatable :: particle_sat(:) !! (np)
    character(len=s_str_len), allocatable :: particle_gas_phase(:) !! (np). gas phase of particle. 
    ! Only for saturation particles
    integer, allocatable :: particle_gas_phase_ind(:) !! np. index of gas phase of particle
    integer, allocatable :: gas_particle_ind(:) !! (nq). Index of particle phase of gas
    character(len=s_str_len), allocatable :: particle_optical_prop(:) !! (np)
    integer, allocatable :: particle_optical_type(:) !! (np) 1 == mie, 2 == fractal
    
    ! reactions
    logical :: reverse !! True if there are reverse reactions
    integer :: nrF !! number of forward reactions
    integer :: nrR !! number of reverse reactions
    integer :: nrT !! number of total reactions
    type(Reaction), allocatable :: rx(:) !! (nrT) array of reaction objects
    character(len=m_str_len), allocatable :: reaction_equations(:) !! (nrT)
    type(ProdLoss), allocatable :: pl(:) !! (nsp) reactions producing and destroying each species
    integer :: kj !! number of photolysis reactions
    integer, allocatable :: photonums(:) !! (kj) the reaction number of each photolysis reaction

    ! raditative transfer
    integer :: nw !! number of wavelength bins
    real(dp), allocatable :: wavl(:) !! (nw+1) wavelength bins in nm
    type(XsectionData), allocatable :: absorp_xs(:) !! Contains photoabsorption species
    real(dp), allocatable :: photolysis_xs(:,:) !! (kj,nw) The xs time qy for each photolysis reaction
    integer :: nray !! number of species with rayleigh scattering
    real(dp), allocatable :: sigray(:,:) !! (len(raynums), nw)
    integer, allocatable :: raynums(:) !! species number of rayleigh species
    
    ! particle radiative transfer
    integer :: nrad_file
    real(dp), allocatable  :: radii_file(:,:) !! particle radii in optical data files
    ! We use array of types for particle xs because we want the option
    ! to exclude optical properties, but not take up a ton of useless memory.
    ! So some elements of this array have nothing in it.
    type(ParticleXsections), allocatable :: part_xs_file(:) !! np in length
    
    ! initial conditions  
    !> Only used in EVOATMOSPHERE.
    !> If True, then during initialization, molecules in the atmosphere will
    !> be conserved.
    logical :: conserving_init = .false.
    integer :: nzf !! number of atmospheric layers in file
    real(dp), allocatable :: z_file(:) !! (nzf) cm
    real(dp), allocatable :: T_file(:) !! (nzf) K
    real(dp), allocatable :: edd_file(:) !! (nzf) cm2/s
    real(dp), allocatable :: den_file(:) !! (nzf) molecules/cm2
    real(dp), allocatable :: mix_file(:,:) !! (nq,nzf) mixing ratios
    real(dp), allocatable :: particle_radius_file(:,:) !! (np,nzf) cm
    
    ! settings
    real(dp) :: planet_mass !! grams
    real(dp) :: planet_radius !! cm
    logical :: fix_water_in_trop !! True if fixing water in troposphere
    integer :: LH2O !! index of H2O
    logical :: water_cond !! True if water should condense out of the atmosphere
    logical :: gas_rainout !! True if gas rains out
    integer :: H_escape_type !! Diffusion-limited, Zahnle, or None
    real(dp), allocatable :: H_escape_coeff ! Coefficient for zahnle hydrogen escape
    integer :: LH2 !! H2 index
    integer :: LH !! H index
    
  end type
  
  type :: PhotochemVars
    ! PhotochemVars contains information that can change between
    ! different photochemical integrations, without reading in new
    ! files.
    
    !> where the photochem data is
    character(:), allocatable :: data_dir

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set DURING file read-in !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! boundary conditions
    !> (nq) Type of lower bc. See photochem_enum.f90 for different types. 
    integer, allocatable :: lowerboundcond(:)
    real(dp), allocatable :: lower_vdep(:)
    real(dp), allocatable :: lower_flux(:)
    real(dp), allocatable :: lower_dist_height(:)
    real(dp), allocatable :: lower_fix_mr(:)
    real(dp), allocatable :: lower_fix_den(:)
    real(dp), allocatable :: lower_fix_press(:)
    integer, allocatable :: upperboundcond(:) ! 0 or 2
    real(dp), allocatable :: upper_veff(:)
    real(dp), allocatable :: upper_flux(:)
    !> True if only use eddy diffusion.
    logical, allocatable :: only_eddy(:) 
    !> Functions for settings time-dependent production rates (nq)
    type(time_dependent_rate_fcns), allocatable :: rate_fcns(:)
    
    ! Atmospheres structure
    real(dp) :: bottom_atmos !! cm
    real(dp) :: top_atmos !! cm
    integer :: nz !! number of vertical layers
    real(dp) :: surface_pressure !! bars
    real(dp) :: surface_albedo
    real(dp) :: diurnal_fac = 0.5_dp !! Default is 0.5, to account for half planet facing the sun.
    real(dp) :: solar_zenith !! degrees
    real(dp) :: trop_alt !! cm (only for fix_water_in_trop == true or gas_rainout == true)
    real(dp) :: rainfall_rate !! relative to modern Earth's average rainfall rate of 1.1e17 molecules/cm2/s
    integer :: trop_ind !! index of troposphere (only for fix_water_in_trop == true or gas_rainout == true)
    logical :: use_manabe !! use manabe formula
    real(dp) :: relative_humidity !! relative humidity if no manabe
    type(CondensationParameters) :: H2O_cond_params !! H2O condensation rate parameters
    
    ! radiative transfer
    real(dp), allocatable :: photon_flux(:) !! (nw) photon/cm^2/s in each wavelength bin hitting planet.
    !> for scaling photon flux for different planets in a solar system
    real(dp) :: photon_scale_factor 
    !> A function for altering the photon flux over time
    procedure(time_dependent_flux_fcn), nopass, pointer :: photon_flux_fcn => null()
    
    ! particles
    !> Parameters describing condensation and evaporation rates and
    !> the RH needed for condensation
    type(CondensationParameters), allocatable :: cond_params(:) ! (np)
    logical :: evaporation = .true. !! If true, then evaporation occurs.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set AFTER file read-in !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: neqs !! number of equations nq*nz
    real(dp), allocatable :: temperature(:) !! (nz) K
    real(dp), allocatable :: z(:) !! (nz) cm
    real(dp), allocatable :: dz(:) !! (nz) cm
    real(dp), allocatable :: edd(:) !! (nz) cm2/s
    !> A function for specifying a custom binary diffusion parameter (b_ij)
    procedure(binary_diffusion_fcn), nopass, pointer :: custom_binary_diffusion_fcn => null()
    real(dp), allocatable :: grav(:) !! (nz) cm/s2
    real(dp), allocatable :: usol_init(:,:) !! (nq,nz) molecules/cm^3.
    real(dp), allocatable :: particle_radius(:,:) !! (np,nz) cm
    real(dp), allocatable :: xs_x_qy(:,:,:) !! (nz,kj,nw) photolysis cross sections times quantum yields (cm2/molecule)
    type(ParticleXsections), allocatable :: particle_xs(:) !! (np)
    real(dp), allocatable :: gibbs_energy(:,:) !! (nz,ng) Joules/mol

    ! Custom optical properties
    real(dp), allocatable :: tauc(:,:) !! (nz,nw) Custom optical depth in each layer
    real(dp), allocatable :: w0c(:,:) !! (nz,nw) Custom single scattering albedo
    real(dp), allocatable :: g0c(:,:) !! (nz,nw) Custom asymetry parameter
    
    ! output
    logical :: at_photo_equilibrium = .false.
    real(dp), allocatable :: usol_out(:,:)
    
    ! other 
    !> number of times we initialize CVODE when it returns
    !> a potentially recoverable error. ONLY USED IN EVOATMOSPHERE (NOT ATMOSPHERE)
    !> in the `evolve` method.
    integer :: max_error_reinit_attempts = 2 
    real(c_double) :: rtol = 1.0e-3_dp !! integration relative tolerance
    !> Integration absolute tolerance. If autodiff == .true., then the model
    !> works better when atol is smaller (e.g., atol = ~1.0e-18).
    real(c_double) :: atol = 1.0e-23_dp 
    integer :: mxsteps = 100000 !! max number of steps before integrator will give up.
    !> seconds. atomsphere considered in equilibrium if integrations reaches this time.
    real(dp) :: equilibrium_time = 1.0e17_dp
    !> For convergence checking. Considers mixing ratio change between t_now and time 
    !> t = t_now*conv_hist_factor to see if atmosphere is changing.
    real(dp) :: conv_hist_factor = 0.5_dp
    !> Minimum mixing ratio considered in convergence checking.
    real(dp) :: conv_min_mix = 1.0e-20_dp
    !> Threshold normalized change in mixing ratios for converchecking check.
    !> A reasonable value is ~1.0e-2.
    real(dp) :: conv_longdy = 0.0_dp
    !> Threshold normalized change in mixing ratios per time change for
    !> convergence checking.
    real(dp) :: conv_longdydt = 1.0e-6_dp
    real(c_double) :: initial_dt = 1.0e-6_dp !! intial timestep size (seconds)
    !> Maximum time step size (seconds).
    real(c_double) :: max_dt = sqrt(huge(1.0_dp))
    integer(c_int) :: max_err_test_failures = 15 !! CVODE max error test failures
    integer(c_int) :: max_order = 5 !! CVODE max order for BDF method.
    !> If .true., then the chemistry terms of the Jacobian are computed uses 
    !> foward mode automatic differentiation.
    logical :: autodiff = .true.
    !> Perturbation for finite difference Jacobian calculation, when autodiff == .false.
    real(dp) :: epsj = 1.0e-4_dp 
    integer :: verbose = 1 !! 0 == no printing. 1 == some printing. 2 == bunch of printing.
    !> Arbitrary rate that is fast (1/s). Used for keeping H2O at saturation in troposphere
    real(dp) :: fast_arbitrary_rate = 1.0e-2_dp 
    !> If True, then the code uses a 1st order upwind method for the advective molecular
    !> diffusion terms instead of a centered scheme. This permits stability (at the cost 
    !> of accuracy) for atmospheres with strong molcular advection in the upper atmosphere.
    logical :: upwind_molec_diff = .false.

    ! Settings for the `robust_stepper` and `find_steady_state` methods
    !> Number of integration errors before giving up completely
    integer :: nerrors_before_giveup = 10
    !> Number of steps to take before checking for convergence
    integer :: nsteps_before_conv_check = 300
    !> Number of steps before reinitializing the integration
    integer :: nsteps_before_reinit = 1000
    !> Number of total steps to take before giving up.
    integer :: nsteps_before_giveup = 100000
    !> When the integrator reinitializes, this is the minimum
    !> density allowed (molecules/cm^3)
    real(dp) :: reinit_min_density = 1.0e-40_dp
    ! End settings for `robust_stepper` and `find_steady_state`

  end type

  type :: SundialsData
    !> cvode memory
    type(c_ptr) :: cvode_mem = c_null_ptr
    ! solution vector
    real(c_double), allocatable :: yvec(:)
    type(N_Vector), pointer :: sunvec_y => NULL()
    ! absolute tolerance (Used in EvoAtmosphere only)
    real(c_double), allocatable :: abstol(:)
    type(N_Vector), pointer :: abstol_nvec => NULL()
    ! matrix and linear solver
    type(SUNMatrix), pointer :: sunmat => NULL()
    type(SUNLinearSolver), pointer :: sunlin => NULL()
  contains
    procedure :: finalize => SundialsData_finalize
    final :: SundialsData_final
  end type

  type :: SundialsDataFinalizer
    type(SundialsData), pointer :: sun => NULL()
  contains
    final :: SundialsDataFinalizer_final
  end type
  
  type :: PhotochemWrk
    ! PhotochemWrk are work variables that change
    ! through the course of an integration

    ! Total step counter for robust_step method
    !> Total number of steps in a robust integration.
    integer :: nsteps_total = -1
    !> Total number of errors experienced in the robust integration.
    integer :: nerrors_total = -1
    
    ! used in cvode
    integer(c_long) :: nsteps_previous = -10 !! For printing
    type(SundialsData) :: sun !! CVODE data

    ! All for determining convergence
    integer :: nsteps = 0 !! Number of integration steps excuted. Updated
                          !! after every successful step.
    !> History of times at previous integration steps. Index 1 is current, 
    !> while index 2, 3, 4 are previous steps. Updated after every successful step.
    real(dp), allocatable :: t_history(:)
    !> History of mixing ratios at previous integration steps. Index 1 is 
    !> current, while index 2, 3, 4 are previous steps. Updated after 
    !> every successful step.
    real(dp), allocatable :: mix_history(:,:,:)
    !> Change in mixing ratio over some number of integrations steps. Updated
    !> during convergence checking.
    real(dp), allocatable :: dmix(:,:)
    !> Normalized change in mixing ratios over some number of integrations steps
    real(dp) :: longdy = 0.0_dp
    !> Normalized change in mixing ratios divided by change in time
    !> over some number of integrations steps.
    real(dp) :: longdydt = 0.0_dp
    ! end stuff for determining convergence

    !> The current time (seconds). The is updated with each call to the
    !> right hand side, and jacobian. It is only important if
    !> var%photon_flux_fcn is set.
    real(dp) :: tn = 0.0_dp
    
    ! Used in prep_all_background_gas
    ! work arrays
    real(dp), allocatable :: usol(:,:) !! (nq,nz)
    real(dp), allocatable :: densities(:,:) !! (nsp+1,nz)
    real(dp), allocatable :: density(:) !! (nz)
    real(dp), allocatable :: rx_rates(:,:) !! (nz,nrT)
    real(dp), allocatable :: mubar(:) !! (nz)
    real(dp), allocatable :: pressure(:) !! (nz)
    real(dp), allocatable :: H2O_rh(:) !! (nz)
    real(dp), allocatable :: H2O_sat_mix(:) !! (nz)
    real(dp), allocatable :: prates(:,:) !! (nz,kj)
    real(dp), allocatable :: surf_radiance(:) !! (nw)
    real(dp), allocatable :: amean_grd(:,:) !! (nz,nw)
    real(dp), allocatable :: optical_depth(:,:) !! (nz,nw)
    real(dp), allocatable :: upper_veff_copy(:) !! (nq)
    real(dp), allocatable :: lower_vdep_copy(:) !! (nq)
    real(dp), allocatable :: xp(:) !! (nz)
    real(dp), allocatable :: xl(:) !! (nz)
    ! diffusion and H escape
    real(dp), allocatable :: DU(:,:) !! (nq,nz)
    real(dp), allocatable :: DD(:,:) !! (nq,nz)
    real(dp), allocatable :: DL(:,:) !! (nq,nz)
    real(dp), allocatable :: ADU(:,:) !! (nq,nz)
    real(dp), allocatable :: ADL(:,:) !! (nq,nz)
    real(dp), allocatable :: ADD(:,:) !! (nq,nz)
    real(dp) :: VH2_esc
    real(dp) :: VH_esc
    ! other
    real(dp), allocatable :: sum_usol(:) !! (nz)
    real(dp), allocatable :: scale_height(:)
    real(dp), allocatable :: wfall(:,:)
    real(dp), allocatable :: gas_sat_den(:,:)
    real(dp), allocatable :: molecules_per_particle(:,:)
    real(dp), allocatable :: rainout_rates(:,:)
    ! end used in prep_all_background_gas

    ! Work space for autodiff jacobian
    !> A sparse representation of the block diagonal Jacobian.
    !> chemistry terms only.
    real(dp), allocatable :: djac_chem(:,:)
    ! end work space for autodiff jacobian
    
  contains
    procedure :: init => init_PhotochemWrk
  end type

  type, extends(PhotochemWrk) :: PhotochemWrkEvo
    real(dp), allocatable :: mix(:,:) !! (nq,nz) mixing ratio.
    real(dp), allocatable :: pressure_hydro(:) !! (nz)
    real(dp), allocatable :: density_hydro(:) !! (nz)

  contains
    procedure :: init => init_PhotochemWrkEvo

  end type
  
contains

  subroutine init_PhotochemWrkEvo(self, nsp, np, nq, nz, nrT, kj, nw)
    class(PhotochemWrkEvo), intent(inout) :: self
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw

    call init_PhotochemWrk(self, nsp, np, nq, nz, nrT, kj, nw)

    if (allocated(self%mix)) then
      deallocate(self%mix)
      deallocate(self%pressure_hydro)
      deallocate(self%density_hydro)
    endif

    allocate(self%mix(nq,nz))
    allocate(self%pressure_hydro(nz))
    allocate(self%density_hydro(nz))

  end subroutine
 
  subroutine init_PhotochemWrk(self, nsp, np, nq, nz, nrT, kj, nw)
    use photochem_const, only: nsteps_save
    class(PhotochemWrk), intent(inout) :: self
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw
    
    if (allocated(self%usol)) then
      deallocate(self%t_history)
      deallocate(self%mix_history)
      deallocate(self%dmix)
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
      deallocate(self%scale_height)
      deallocate(self%wfall)
      deallocate(self%gas_sat_den)
      deallocate(self%molecules_per_particle)
      deallocate(self%rainout_rates)
      deallocate(self%djac_chem)
    endif
    
    allocate(self%t_history(nsteps_save))
    allocate(self%mix_history(nq,nz,nsteps_save))
    allocate(self%dmix(nq,nz))
    allocate(self%usol(nq,nz))
    allocate(self%mubar(nz))
    allocate(self%pressure(nz))
    allocate(self%density(nz))
    allocate(self%H2O_rh(nz))
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
    allocate(self%scale_height(nz))
    allocate(self%wfall(np,nz))
    allocate(self%gas_sat_den(np,nz))
    allocate(self%molecules_per_particle(np,nz))
    allocate(self%rainout_rates(nq,nz))
    allocate(self%djac_chem(nq,nz*nq))
  end subroutine

  subroutine SundialsData_finalize(self, err)
    use iso_c_binding, only: c_associated, c_null_ptr, c_int
    use fcvode_mod, only: FCVodeFree
    use fsundials_nvector_mod, only: FN_VDestroy
    use fsundials_matrix_mod, only: FSUNMatDestroy
    use fsundials_linearsolver_mod, only: FSUNLinSolFree
    class(SundialsData), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    integer(c_int) :: ierr

    if (allocated(self%yvec)) then
      deallocate(self%yvec)
    endif
    if (associated(self%sunvec_y)) then
      call FN_VDestroy(self%sunvec_y)
      nullify(self%sunvec_y)
    endif

    if (allocated(self%abstol)) then
      deallocate(self%abstol)
    endif
    if (associated(self%abstol_nvec)) then
      call FN_VDestroy(self%abstol_nvec)
      nullify(self%abstol_nvec)
    endif

    if (c_associated(self%cvode_mem)) then
      call FCVodeFree(self%cvode_mem)
      self%cvode_mem = c_null_ptr
    endif

    if (associated(self%sunlin)) then
      ierr = FSUNLinSolFree(self%sunlin)
      if (ierr /= 0) then
        err = "Sundials failed to deallocated linear solver"
      end if
      nullify(self%sunlin)
    endif

    if (associated(self%sunmat)) then
      call FSUNMatDestroy(self%sunmat)
      nullify(self%sunmat)
    endif

  end subroutine

  subroutine SundialsData_final(self)
    type(SundialsData), intent(inout) :: self
    character(:), allocatable :: err
    call SundialsData_finalize(self, err)
  end subroutine

  subroutine SundialsDataFinalizer_final(self)
    type(SundialsDataFinalizer), intent(inout) :: self
    character(:), allocatable :: err
    if (associated(self%sun)) then
      call self%sun%finalize(err)
    endif
  end subroutine
  
end module


