
module photochem_types ! make a giant IO object
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: str_len = 1024
  
  public PhotoMechanism, PhotoSettings, PhotoRadTran, PhotoInitAtm
  public WrkBackgroundAtm
  
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
    logical :: fix_water_in_trop
    integer :: LH2O
    logical :: use_manabe ! use manabe formula
    real(real_kind) :: relative_humidity ! relative humidity if no manabe
    logical :: gas_rainout
    logical :: stratospheric_cond
    real(real_kind) :: relative_humidity_cold_trap
    real(real_kind) :: H2O_condensation_rate(2)
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
    
    ! condensation rate of particles
    real(real_kind), allocatable :: condensation_rate(:,:)
    
  end type
  
  type :: PhotoMechanism
    
    ! molecules
    integer :: nq ! nqp + nll, total number of PDEs
    
    integer :: ng_1
    integer :: nll ! number of long lived gas molecules (number of gas-phase PDEs)
    integer :: nsl ! number of short lived gas molecules (their index is between nq + 1 and nsp - 1)
    integer :: ng ! nll + nsl + 1, total number of gas-phase species, including background

    integer :: nsp ! nll + nsl + npq - 1, background gas index
    ! nsp + 1 is hv index
    ! nsp + 2 is M index
    
    integer :: natoms
    character(len=8), allocatable :: atoms_names(:) 
    real(real_kind), allocatable :: atoms_mass(:) 
    character(len=15), allocatable :: species_names(:) ! nsp (particle names + gas names)
    integer, allocatable :: species_composition(:,:) ! natoms, nsp
    real(real_kind), allocatable :: species_mass(:) ! nsp
    real(real_kind), allocatable :: thermo_data(:,:,:) ! ng
    real(real_kind), allocatable :: thermo_temps(:,:) ! ng
    real(real_kind), allocatable :: henry_data(:,:)
    
    ! particles
    logical :: there_are_particles
    integer :: np ! number of particles
    integer :: npq ! number of particle equations. for now nq = npq.
    character(len=15), allocatable :: particle_names(:) ! nqp
    integer, allocatable :: particle_formation_method(:) ! np
    real(real_kind), allocatable :: particle_density(:) ! np
    real(real_kind), allocatable :: particle_sat_params(:,:) ! 3, np
    character(len=15), allocatable :: particle_gas_phase(:) ! np
    integer, allocatable :: particle_gas_phase_ind(:) ! np
    character(len=50), allocatable :: particle_optical_prop(:) ! np
    integer, allocatable :: particle_optical_type(:) ! np
    
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
    
    ! particles
    integer :: nrad_file
    real(real_kind), allocatable :: radii_file(:,:) ! (nrad, np)
    real(real_kind), allocatable :: w0_file(:,:,:) ! (nw, nrad, np)
    real(real_kind), allocatable :: qext_file(:,:,:) ! (nw, nrad, np)
    real(real_kind), allocatable :: g_file(:,:,:) ! (nw, nrad, np)

  end type
  
  type :: PhotoInitAtm
    integer :: nzf
    real(real_kind), allocatable :: z_file(:)
    real(real_kind), allocatable :: T_file(:)
    real(real_kind), allocatable :: edd_file(:)
    real(real_kind), allocatable :: usol_file(:,:)
    real(real_kind), allocatable :: particle_radius_file(:,:)
    logical :: no_water_profile
  end type
  
  type :: WrkBackgroundAtm
    
    ! Used in prep_all_background_gas
    ! dimensions
    integer :: nsp
    integer :: np
    integer :: nq
    integer :: nz
    integer :: nrT
    integer :: kj
    integer :: nw
    integer :: trop_ind
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
    procedure :: init => init_WrkBackgroundAtm
  end type
  
contains
 
  subroutine init_WrkBackgroundAtm(self, nsp, np, nq, nz, nrT, kj, nw, trop_ind)
    class(WrkBackgroundAtm), intent(inout) :: self
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw, trop_ind
    
    self%nsp = nsp
    self%np = np
    self%nq = nq
    self%nz = nz
    self%nrT = nrT
    self%kj = kj
    self%nw = nw
    self%trop_ind = trop_ind
    
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
    allocate(self%rainout_rates(self%nq,trop_ind))
  end subroutine
  
end module








