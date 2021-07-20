
module photochem_types ! make a giant IO object
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: str_len = 1024
  
  public PhotoMechanism, PhotoSettings, PhotoRadTran, PhotoInitAtm
  
  type :: PhotoSettings
    real(real_kind) :: bottom_atmos
    real(real_kind) :: top_atmos 
    integer :: nz
    
    logical :: regular_grid
    real(real_kind) :: lower_wavelength
    real(real_kind) :: upper_wavelength
    integer :: nw
    character(len=str_len) :: grid_file

    logical :: back_gas
    character(len=str_len) :: back_gas_name
    real(real_kind) :: back_gas_mu
    integer :: back_gas_ind
    real(real_kind) :: surface_pressure ! this can be changed
    real(real_kind) :: planet_mass
    real(real_kind) :: planet_radius
    real(real_kind) :: surface_albedo ! this can be changed
    real(real_kind) :: diurnal_fac
    real(real_kind) :: solar_zenith
    logical :: water_sat_trop
    integer :: LH2O
    logical :: diff_H_escape
    integer :: LH2
    integer :: LH
    real(real_kind) :: trop_alt ! this can be changed
    
    integer :: nq ! nubmer of long lived
    integer :: nsl ! number of short lived
    character(len=20), allocatable :: SL_names(:)
    
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
    integer :: nq
    integer :: nsl
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

  end type
  
  type :: PhotoRadTran
    
    integer :: nw
    real(real_kind), allocatable :: wavl(:) ! (nw+1)

    integer, allocatable :: num_temp_cols(:) ! (kj)
    integer, allocatable :: sum_temp_cols(:) ! (kj)
    ! All data for every reaction in single vector to save memory
    real(real_kind), allocatable :: xs_data(:) ! (sum(num_temp_cols)*nw) 
    real(real_kind), allocatable :: xs_data_temps(:,:) ! (maxval(num_temp_cols), kj)

    integer :: nray
    real(real_kind), allocatable :: sigray(:,:) ! (len(raynums), nw)
    integer, allocatable :: raynums(:) ! species number of rayleigh species
    
    
    ! need some photons
    real(real_kind), allocatable :: photon_flux(:) ! (nw) photonz

  end type
  
  type :: PhotoInitAtm
    integer :: nzf
    real(real_kind), allocatable :: z_file(:)
    real(real_kind), allocatable :: T_file(:)
    real(real_kind), allocatable :: edd_file(:)
    real(real_kind), allocatable :: usol_file(:,:)
  end type
  
end module








