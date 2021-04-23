
module photochem_types ! make a giant IO object
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  
  public PhotoMechanism, PhotoSettings, PhotoPlanet
  
  type :: PhotoSettings
    real(real_kind) :: bottom_atmosphere
    real(real_kind) :: top_atmosphere 
    integer :: nz
    real(real_kind) :: lower_wavelength
    real(real_kind) :: upper_wavelength
    integer :: nw !number of bins
  end type
  
  type :: PhotoPlanet
    real(real_kind) :: gravity
    real(real_kind) :: surface_pressure
    real(real_kind) :: planet_radius
    real(real_kind) :: surface_albedo
    logical :: water_sat_trop
    real(real_kind) :: trop_alt
    logical :: lightning
    real(real_kind) :: lightning_NO_production
    logical :: rainout
    real(real_kind) :: rainout_multiplier
  end type
  
  type :: PhotoMechanism
    ! settings
    type(PhotoSettings) :: settings
    
    ! planet information
    type(PhotoPlanet) :: planet
    
    ! molecule data
    integer :: nsp
    integer :: natoms
    character(len=8), allocatable :: atoms_names(:) 
    real(real_kind), allocatable :: atoms_mass(:) 
    character(len=8), allocatable :: species_names(:)
    integer, allocatable :: species_composition(:,:)
    real(real_kind), allocatable :: species_mass(:) 
    integer, allocatable :: lowerboundcond(:)
    real(real_kind), allocatable :: lower_vdep(:)
    real(real_kind), allocatable :: lower_flux(:)
    real(real_kind), allocatable :: lower_distributed_height(:)
    real(real_kind), allocatable :: lower_fixed_mr(:)
    integer, allocatable :: upperboundcond(:)
    real(real_kind), allocatable :: upper_veff(:)
    real(real_kind), allocatable :: upper_flux(:)
    real(real_kind), allocatable :: thermo_data(:,:,:)
    real(real_kind), allocatable :: thermo_temps(:,:)
    
    ! reaction
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
    character(len=15), allocatable :: rxtypes(:)
    real(real_kind), allocatable :: rateparams(:,:)
    
    ! to do
    integer, allocatable :: nump(:) ! 
    integer, allocatable :: numl(:)
    integer, allocatable :: iprod(:,:)
    integer, allocatable :: iloss(:,:)
    
  end type
  
end module








