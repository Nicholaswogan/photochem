module photochem_data
  use photochem_const, only: dp, s_str_len, m_str_len
  use clima_saturationdata, only: SaturationData
  implicit none
  private

  public :: PhotochemData
  public :: XsectionData, ParticleXsections, ThermodynamicData
  public :: Reaction, Efficiencies, BaseRate, ReverseRate, PhotolysisRate
  public :: ElementaryRate, ThreeBodyRate, FalloffRate, PressDependentRate
  public :: MultiArrheniusRate, ProdLoss

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

    ! radiative transfer
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

    ! settings
    real(dp) :: planet_mass !! grams
    real(dp) :: planet_radius !! cm
    integer :: LH2O !! index of H2O; used by gas rainout
    logical :: gas_rainout !! True if gas rains out
    integer :: H_escape_type !! Diffusion-limited, Zahnle, or None
    real(dp) :: H_escape_coeff ! Coefficient for zahnle hydrogen escape
    integer :: LH2 !! H2 index
    integer :: LH !! H index
  end type

end module
