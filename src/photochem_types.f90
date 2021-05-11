
module photochem_types ! make a giant IO object
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: str_len = 1000
  
  public PhotoMechanism, PhotoSettings, PhotoRadTran
  
  type :: PhotoSettings
    real(real_kind) :: bottom_atmos
    real(real_kind) :: top_atmos 
    integer :: nz
    real(real_kind), allocatable :: z(:), dz(:)
    
    logical :: regular_grid
    real(real_kind) :: lower_wavelength
    real(real_kind) :: upper_wavelength
    integer :: nw
    character(len=str_len) :: grid_file

    logical :: back_gas
    integer :: back_gas_ind
    real(real_kind) :: surface_pressure
    real(real_kind) :: planet_mass
    real(real_kind) :: planet_radius
    real(real_kind), allocatable :: grav(:) ! compute from planet mass and radius
    real(real_kind) :: surface_albedo
    logical :: water_sat_trop
    real(real_kind) :: trop_alt
    
    ! species types/bookkeeping
    integer, allocatable :: species_type(:)
    integer :: nq ! nubmer of long lived
    integer :: nsl ! number of short lived
    integer, allocatable :: LL_inds(:)
    integer, allocatable :: SL_inds(:)
    
    integer, allocatable :: lowerboundcond(:) ! 0, 1, 2 or 3
    real(real_kind), allocatable :: lower_vdep(:)
    real(real_kind), allocatable :: lower_flux(:)
    real(real_kind), allocatable :: lower_dist_height(:)
    real(real_kind), allocatable :: lower_fix_mr(:)
    integer, allocatable :: upperboundcond(:) ! 0 or 1
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
    character(len=15), allocatable :: rxtypes(:)
    real(real_kind), allocatable :: rateparams(:,:)
    integer, allocatable :: nump(:) ! number of production mechanisms (rxns) for each sp
    integer, allocatable :: numl(:) ! number of loss mechanisms (rxns) for each sp
    integer, allocatable :: iprod(:,:) ! (nmax,nsp) returns reaction # of production mechanism for sp
    integer, allocatable :: iloss(:,:) ! (nmax,nsp) returns reaction # of loss mechanism for sp
    
    integer :: kj ! number of photolysis reactions
    integer, allocatable :: photonums(:) ! (kj) the reaction number of each photolysis reaction

  end type
  
  type :: PhotoRadTran
    
    integer :: nw
    real(real_kind), allocatable :: wavl(:) ! (nw+1)

    integer, allocatable :: num_temp_cols(:) ! (kj)
    integer, allocatable :: sum_temp_cols(:) ! (kj)
    ! All data for every reaction in single vector to save memory
    real(real_kind), allocatable :: xs_data(:) ! (sum(num_temp_cols)*nw) 
    real(real_kind), allocatable :: xs_data_temps(:,:) ! (maxval(num_temp_cols), kj)

    real(real_kind), allocatable :: sigray(:,:) ! (len(raynums), nw)
    integer, allocatable :: raynums(:) ! species number of rayleigh species
    
    ! need some photons
    real(real_kind), allocatable :: photon_flux(:) ! (nw) photonz
    
    real(real_kind), allocatable :: xs_x_qy(:,:,:) ! (kj,nz,nw) cross section * quantum yield

  end type
  
end module








