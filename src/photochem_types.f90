
module photochem_types ! make a giant IO object
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  
  public PhotoMechanism, PhotoSettings
  
  type :: PhotoSettings
    
    real(real_kind) :: gravity
    real(real_kind) :: surface_pressure
    real(real_kind) :: planet_radius
    real(real_kind) :: surface_albedo
    logical :: water_sat_trop
    real(real_kind) :: trop_alt
    
    real(real_kind) :: bottom_atmosphere
    real(real_kind) :: top_atmosphere 
    integer :: nz
    real(real_kind) :: lower_wavelength
    real(real_kind) :: upper_wavelength
    integer :: nw !number of bins

    integer, allocatable :: lowerboundcond(:)
    real(real_kind), allocatable :: lower_vdep(:)
    real(real_kind), allocatable :: lower_flux(:)
    real(real_kind), allocatable :: lower_distributed_height(:)
    real(real_kind), allocatable :: lower_fixed_mr(:)
    integer, allocatable :: upperboundcond(:)
    real(real_kind), allocatable :: upper_veff(:)
    real(real_kind), allocatable :: upper_flux(:)
  end type
  
  type :: PhotoMechanism
    
    ! molecules
    integer :: nsp
    integer :: natoms
    character(len=8), allocatable :: atoms_names(:) 
    real(real_kind), allocatable :: atoms_mass(:) 
    character(len=8), allocatable :: species_names(:)
    integer, allocatable :: species_composition(:,:)
    real(real_kind), allocatable :: species_mass(:) 
    real(real_kind), allocatable :: thermo_data(:,:,:)
    real(real_kind), allocatable :: thermo_temps(:,:)
    
    ! reactions
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
    integer, allocatable :: nump(:) ! number of production mechanisms (rxns) for each sp
    integer, allocatable :: numl(:) ! number of loss mechanisms (rxns) for each sp
    integer, allocatable :: iprod(:,:) ! (nmax,nsp) returns reaction # of production mechanism for sp
    integer, allocatable :: iloss(:,:) ! (nmax,nsp) returns reaction # of loss mechanism for sp
    integer :: kj ! number of photolysis reactions
    integer, allocatable :: photonums(:) ! the reaction number of each photolysis reaction
    
    
    ! needs to go somewhere else
    ! raditative properties
    real(real_kind), allocatable :: xs_x_qy(:,:,:) ! (kj,nz,nw) cross section * quantum yield
    real(real_kind), allocatable :: photon_flux(:) ! (nw) photonzzz
    
  end type
  
end module








