
module photochem_types ! make a giant IO object
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: str_len = 1024
  
  public PhotoMechanism, PhotoSettings, PhotoRadTran, PhotoInitAtm
  public WrkBackgroundAtm, WrkTwoStream
  
  type :: PhotoSettings
    real(real_kind) :: bottom_atmos
    real(real_kind) :: top_atmos 
    integer :: nz
    
    logical :: regular_grid
    real(real_kind) :: lower_wavelength
    real(real_kind) :: upper_wavelength
    integer :: nw
    character(len=str_len) :: grid_file
    real(real_kind) :: photon_scale_factor

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
    logical :: use_manabe ! use manabe formula
    real(real_kind) :: relative_humidity ! relative humidity if no manabe
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
    integer :: nq ! number of PDEs
    integer :: nll ! number of long lived gas molecules
    integer :: nsl ! number of short lived gas molecules
    integer :: ng ! nll + nsl, total number of gas-phase species
    integer :: nsp ! nll + nsl + np
    integer :: natoms
    character(len=8), allocatable :: atoms_names(:) 
    real(real_kind), allocatable :: atoms_mass(:) 
    character(len=15), allocatable :: species_names(:)
    integer, allocatable :: species_composition(:,:)
    real(real_kind), allocatable :: species_mass(:) 
    real(real_kind), allocatable :: thermo_data(:,:,:)
    real(real_kind), allocatable :: thermo_temps(:,:)
    
    ! particles
    logical :: there_are_particles
    integer :: np
    character(len=15), allocatable :: particle_names(:)
    integer, allocatable :: particle_formation_method(:)
    real(real_kind), allocatable :: particle_density(:)
    real(real_kind), allocatable :: particle_monomer_radius(:)
    integer, allocatable :: particle_nr(:)
    integer, allocatable :: ind_smaller_polymer_maker(:,:) ! (particle_nr(), np)
    real(real_kind), allocatable :: particle_eta(:,:)
    
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
    logical :: no_water_profile
  end type
  
  type :: WrkBackgroundAtm
    
    ! dimensions
    integer :: nsp
    integer :: nq
    integer :: nz
    integer :: nrT
    integer :: kj
    integer :: nw
    integer :: trop_ind
    
    ! work arrays
    real(real_kind), allocatable :: usol(:,:) ! (nq,nz)
    real(real_kind), allocatable :: mubar(:) ! (nz)
    real(real_kind), allocatable :: pressure(:) ! (nz)
    real(real_kind), allocatable :: density(:) ! (nz)
    real(real_kind), allocatable :: fH2O(:) ! (nz)
    real(real_kind), allocatable :: densities(:,:) ! (nsp+1,nz)
    real(real_kind), allocatable :: rx_rates(:,:) ! (nz,nrT)
    real(real_kind), allocatable :: prates(:,:) ! (nz,kj)
    real(real_kind), allocatable :: surf_radiance(:) ! (nw)
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
    ! boundary conditions copy
    real(real_kind), allocatable :: upper_veff_copy(:) ! (nq)
    
  contains
    procedure :: init => init_WrkBackgroundAtm
  end type
  
  type WrkTwoStream
    real(real_kind), allocatable :: gt(:)
    real(real_kind), allocatable :: gam1(:), gam2(:), gam3(:), gam4(:)
    real(real_kind), allocatable :: lambda(:), cap_gam(:)
    real(real_kind), allocatable :: e1(:), e2(:), e3(:), e4(:)
    real(real_kind), allocatable :: tauc(:), direct(:)
    real(real_kind), allocatable :: cp0(:), cpb(:), cm0(:), cmb(:)
    real(real_kind), allocatable :: A(:), B(:), D(:), E(:)
    real(real_kind), allocatable :: y1(:), y2(:)
  contains
    procedure :: init => init_WrkTwoStream
  end type
  
contains
  
  subroutine init_WrkTwoStream(self, nz)
    class(WrkTwoStream), intent(inout) :: self
    integer, intent(in) :: nz
    
    allocate(self%gt(nz))
    allocate(self%gam1(nz))
    allocate(self%gam2(nz))
    allocate(self%gam3(nz))
    allocate(self%gam4(nz))
    allocate(self%lambda(nz))
    allocate(self%cap_gam(nz))
    allocate(self%e1(nz))
    allocate(self%e2(nz))
    allocate(self%e3(nz))
    allocate(self%e4(nz))
    allocate(self%tauc(nz+1))
    allocate(self%direct(nz+1))
    allocate(self%cp0(nz))
    allocate(self%cpb(nz))
    allocate(self%cm0(nz))
    allocate(self%cmb(nz))
    allocate(self%A(nz*2))
    allocate(self%B(nz*2))
    allocate(self%D(nz*2))
    allocate(self%E(nz*2))
    allocate(self%y1(nz))
    allocate(self%y2(nz))

  end subroutine
  
  
  subroutine init_WrkBackgroundAtm(self, nsp, nq, nz, nrT, kj, nw, trop_ind)
    class(WrkBackgroundAtm), intent(inout) :: self
    integer, intent(in) :: nsp, nq, nz, nrT, kj, nw, trop_ind
    
    self%nsp = nsp
    self%nq = nq
    self%nz = nz
    self%nrT = nrT
    self%kj = kj
    self%nw = nw
    self%trop_ind = trop_ind
    
    allocate(self%usol(nq,nz))
    allocate(self%mubar(nz))
    allocate(self%pressure(nz))
    allocate(self%density(nz))
    allocate(self%fH2O(nz))
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
  end subroutine
  
end module








