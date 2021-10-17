
module photochem_data
  use photochem_const, only: real_kind, err_len, str_len
  use photochem_types, only: XsectionData
  implicit none  
  ! public ! public but protected
  ! Data that doesn't after file read-in

  ! molecules
  integer, protected :: nq
  integer, protected :: ng_1
  integer, protected :: nll
  integer, protected :: nsl
  integer, protected :: ng
  integer, protected :: nsp
  integer, protected :: natoms
  integer, protected :: kd, kl, ku ! not read in. It is nq + nq + 1 (diagonal width of jacobian)
  integer, protected :: lda ! not read in 
  character(len=8), allocatable :: atoms_names(:) 
  real(real_kind), allocatable :: atoms_mass(:) 
  character(len=15), allocatable :: species_names(:)
  integer, allocatable :: species_composition(:,:)
  real(real_kind), allocatable :: species_mass(:) 
  real(real_kind), allocatable :: thermo_data(:,:,:)
  real(real_kind), allocatable :: thermo_temps(:,:)
  real(real_kind), allocatable :: henry_data(:,:)
  
  ! particles
  logical, protected :: there_are_particles
  integer, protected :: np ! number of particles
  integer, protected :: npq ! number of particle equations. for now nq = npq.
  integer, allocatable :: particle_formation_method(:) ! np
  real(real_kind), allocatable :: particle_density(:) ! np
  real(real_kind), allocatable :: particle_sat_params(:,:) ! 3, np
  integer, allocatable :: particle_gas_phase_ind(:) ! np
  
  ! reactions
  logical, protected :: reverse
  integer, protected :: nrF ! number of forward reactions
  integer, protected :: nrR ! number of reverse reactions
  integer, protected :: nrT ! number of total reactions
  integer, protected :: max_num_reactants
  integer, protected :: max_num_products
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
  integer, protected :: kj ! number of photolysis reactions
  integer, allocatable :: photonums(:) ! (kj) the reaction number of each photolysis reaction

  ! raditative transfer
  integer, protected :: nw
  real(real_kind), allocatable :: wavl(:) ! (nw+1)
  !f2py integer(8), allocatable :: xs_data(:)
  type(XsectionData), allocatable :: xs_data(:) ! (kj)
  integer, protected :: nray
  real(real_kind), allocatable :: sigray(:,:) ! (len(raynums), nw)
  integer, allocatable :: raynums(:) ! species number of rayleigh species
  
  ! particle radiative transfer
  integer, protected :: nrad_file
  real(real_kind), allocatable  :: radii_file(:,:) 
  real(real_kind), allocatable  :: w0_file(:,:,:)
  real(real_kind), allocatable  :: qext_file(:,:,:) 
  real(real_kind), allocatable  :: g_file(:,:,:) 
  
  ! initial conditions  
  integer, protected :: nzf
  real(real_kind), allocatable :: z_file(:)
  real(real_kind), allocatable :: T_file(:)
  real(real_kind), allocatable :: edd_file(:)
  real(real_kind), allocatable :: usol_file(:,:)
  real(real_kind), allocatable :: particle_radius_file(:,:)
  
  ! settings
  logical, protected :: back_gas
  real(real_kind), protected :: back_gas_mu
  character(len=str_len), protected :: back_gas_name
  integer, protected :: back_gas_ind
  real(real_kind), protected :: planet_mass
  real(real_kind), protected :: planet_radius
  logical, protected :: fix_water_in_trop
  integer, protected :: LH2O
  logical, protected :: stratospheric_cond
  logical, protected :: diff_H_escape
  integer, protected :: LH2
  integer, protected :: LH
  
  
contains
  
  subroutine setup_files(mechanism_file, settings_file, flux_file, atmosphere_txt, err)
    use photochem_types, only: PhotoMechanism, PhotoSettings, PhotoRadTran, PhotoInitAtm
    use photochem_input, only: read_all_files
    use photochem_vars
    
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(len=err_len), intent(out) :: err
    type(PhotoMechanism) :: photomech
    type(PhotoSettings) :: photoset
    type(PhotoRadTran) :: photorad
    type(PhotoInitAtm) :: photoinit
    
    call read_all_files(mechanism_file, settings_file, flux_file, atmosphere_txt, &
                        photomech, photoset, photorad, photoinit, err)
    if (len_trim(err) /= 0) return
        
    !!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_data !!!
    !!!!!!!!!!!!!!!!!!!!!!
    
    if (allocated(atoms_names)) then
      deallocate(atoms_names)
      deallocate(atoms_mass)
      deallocate(species_names)
      deallocate(species_composition)
      deallocate(species_mass)
      if (allocated(thermo_data)) then
        deallocate(thermo_data)
        deallocate(thermo_temps)
      endif
      deallocate(henry_data)
      
      if (allocated(particle_formation_method)) then
        deallocate(particle_formation_method)
        deallocate(particle_density) 
        deallocate(particle_sat_params)
        deallocate(particle_gas_phase_ind)
      endif
      
      deallocate(reactants_names)
      deallocate(products_names)
      deallocate(reactants_sp_inds)
      deallocate(products_sp_inds)
      deallocate(nreactants)
      deallocate(nproducts)
      if (allocated(reverse_info)) then
        deallocate(reverse_info)
      endif
      deallocate(rxtypes)
      deallocate(rateparams)
      deallocate(efficiencies)
      deallocate(eff_sp_inds)
      deallocate(num_efficient)
      deallocate(def_eff)
      deallocate(falloff_type)
      deallocate(nump)
      deallocate(numl)
      deallocate(iprod)
      deallocate(iloss)
      deallocate(photonums)
      
      deallocate(wavl)
      deallocate(xs_data)
      deallocate(sigray)
      deallocate(raynums)
      
      if (allocated(radii_file)) then 
        deallocate(radii_file)
        deallocate(w0_file)
        deallocate(qext_file)
        deallocate(g_file)
      endif
        
      deallocate(z_file)
      deallocate(T_file)
      deallocate(edd_file)
      deallocate(usol_file)
      if (allocated(particle_radius_file)) then
        deallocate(particle_radius_file)
      endif
    
    endif
    
    ! species
    nq = photomech%nq
    ng_1 = photomech%ng_1
    nll = photomech%nll
    nsl = photomech%nsl
    ng = photomech%ng
    nsp = photomech%nsp
    natoms = photomech%natoms
    kd = 2*nq + 1
    kl = kd + nq
    ku = kd - nq
    lda = 3*nq + 1
    allocate(atoms_names(natoms))
    atoms_names = photomech%atoms_names
    allocate(atoms_mass(natoms))
    atoms_mass = photomech%atoms_mass
    allocate(species_names(nsp+2))
    species_names = photomech%species_names
    allocate(species_composition(natoms,nsp+2))
    species_composition = photomech%species_composition
    allocate(species_mass(nsp))
    species_mass = photomech%species_mass
    if (photomech%reverse) then
      allocate(thermo_data(7,2,ng))
      thermo_data = photomech%thermo_data
      allocate(thermo_temps(3,ng))
      thermo_temps = photomech%thermo_temps
    endif
    allocate(henry_data(2,nsp))
    henry_data = photomech%henry_data
    
    ! particles
    there_are_particles = photomech%there_are_particles
    np = photomech%np
    npq = photomech%npq
    if (there_are_particles) then
      allocate(particle_formation_method(np))
      particle_formation_method = photomech%particle_formation_method
      allocate(particle_density(np))
      particle_density = photomech%particle_density
      allocate(particle_sat_params(3,np))
      particle_sat_params = photomech%particle_sat_params
      allocate(particle_gas_phase_ind(np))
      particle_gas_phase_ind = photomech%particle_gas_phase_ind
    endif
    
    ! reactions
    reverse = photomech%reverse
    nrF = photomech%nrF
    nrR = photomech%nrR
    nrT = photomech%nrT
    max_num_reactants = photomech%max_num_reactants
    max_num_products = photomech%max_num_products
    allocate(reactants_names(max_num_reactants, nrF))
    reactants_names = photomech%reactants_names
    allocate(products_names(max_num_products, nrF))
    products_names = photomech%products_names
    allocate(reactants_sp_inds(max_num_reactants, nrT))
    reactants_sp_inds = photomech%reactants_sp_inds
    allocate(products_sp_inds(max_num_products, nrT))
    products_sp_inds = photomech%products_sp_inds
    allocate(nreactants(nrT))
    nreactants = photomech%nreactants
    allocate(nproducts(nrT))
    nproducts = photomech%nproducts
    if (reverse) then
      allocate(reverse_info(nrT))
      reverse_info = photomech%reverse_info
    endif
    allocate(rxtypes(nrF))
    rxtypes = photomech%rxtypes
    allocate(rateparams(10,nrF))
    rateparams = photomech%rateparams
    allocate(efficiencies(maxval(photomech%num_efficient), nrF))
    efficiencies = photomech%efficiencies
    allocate(eff_sp_inds(maxval(photomech%num_efficient), nrF))
    eff_sp_inds = photomech%eff_sp_inds
    allocate(num_efficient(nrF))
    num_efficient = photomech%num_efficient
    allocate(def_eff(nrF))
    def_eff = photomech%def_eff
    allocate(falloff_type(nrF))
    falloff_type = photomech%falloff_type
    allocate(nump(nsp))
    nump = photomech%nump
    allocate(numl(nsp))
    numl = photomech%numl
    allocate(iprod(maxval(nump),nsp))
    iprod = photomech%iprod
    allocate(iloss(maxval(numl),nsp))
    iloss = photomech%iloss
    kj = photomech%kj
    allocate(photonums(kj))
    photonums = photomech%photonums
    
    ! raditative transfer
    nw = photorad%nw
    allocate(wavl(nw))
    wavl = photorad%wavl
    allocate(xs_data(kj))
    xs_data = photorad%xs_data    
    nray = photorad%nray
    allocate(sigray(size(photorad%sigray,1),nw))
    sigray = photorad%sigray
    allocate(raynums(size(photorad%raynums)))
    raynums = photorad%raynums
    
    ! particle raditative transfer
    if (there_are_particles) then
      nrad_file = photorad%nrad_file
      allocate(radii_file(np, nrad_file))
      radii_file = photorad%radii_file
      allocate(w0_file(np,nrad_file,nw))
      w0_file = photorad%w0_file
      allocate(qext_file(np,nrad_file,nw))
      qext_file = photorad%qext_file
      allocate(g_file(np,nrad_file,nw))
      g_file = photorad%g_file
    endif
    
    ! initial conditions  
    nzf = photoinit%nzf
    allocate(z_file(nzf))
    z_file = photoinit%z_file
    allocate(T_file(nzf))
    T_file = photoinit%T_file
    allocate(edd_file(nzf))
    edd_file = photoinit%edd_file
    allocate(usol_file(nq,nzf))
    usol_file = photoinit%usol_file
    if (there_are_particles) then
      allocate(particle_radius_file(npq,nzf))
      particle_radius_file = photoinit%particle_radius_file
    endif
    no_water_profile = photoinit%no_water_profile
    
    ! settings
    back_gas = photoset%back_gas
    back_gas_mu = photoset%back_gas_mu
    back_gas_name = photoset%back_gas_name
    planet_mass = photoset%planet_mass
    planet_radius = photoset%planet_radius
    fix_water_in_trop = photoset%fix_water_in_trop
    if (fix_water_in_trop) then
      stratospheric_cond = photoset%stratospheric_cond
      LH2O = photoset%LH2O
    endif
    LH2O = photoset%LH2O
    diff_H_escape = photoset%diff_H_escape
    LH2 = photoset%LH2
    LH = photoset%LH
    
    !!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_vars !!!
    !!!!!!!!!!!!!!!!!!!!!!
    
    if (allocated(lowerboundcond)) then
      deallocate(lowerboundcond)
      deallocate(lower_vdep)
      deallocate(lower_flux)
      deallocate(lower_dist_height)
      deallocate(lower_fix_mr)
      deallocate(upperboundcond)
      deallocate(upper_veff)
      deallocate(upper_flux)
      deallocate(photon_flux)
      if (allocated(condensation_rate)) then
        deallocate(condensation_rate)
      endif
    endif
    
    ! boundary conditions
    allocate(lowerboundcond(nq))
    lowerboundcond = photoset%lowerboundcond
    allocate(lower_vdep(nq))
    lower_vdep = photoset%lower_vdep
    allocate(lower_flux(nq))
    lower_flux = photoset%lower_flux
    allocate(lower_dist_height(nq))
    lower_dist_height = photoset%lower_dist_height
    allocate(lower_fix_mr(nq))
    lower_fix_mr = photoset%lower_fix_mr
    allocate(upperboundcond(nq))
    upperboundcond = photoset%upperboundcond
    allocate(upper_veff(nq))
    upper_veff= photoset%upper_veff
    allocate(upper_flux(nq))
    upper_flux = photoset%upper_flux
    
    ! Atmospheres structure
    bottom_atmos = photoset%bottom_atmos
    top_atmos = photoset%top_atmos
    nz = photoset%nz
    surface_pressure = photoset%surface_pressure
    surface_albedo = photoset%surface_albedo 
    solar_zenith = photoset%solar_zenith
    diurnal_fac = photoset%diurnal_fac
    if (fix_water_in_trop) then
      trop_alt = photoset%trop_alt
    else
      trop_alt = -1.d0
    endif
    
    allocate(photon_flux(nw))
    photon_flux = photorad%photon_flux
    photon_scale_factor = photoset%photon_scale_factor  
    
    ! particles
    if (photomech%there_are_particles) then
      allocate(condensation_rate(2,np))
      condensation_rate = photoset%condensation_rate
    endif
    
    ! settings
    if (fix_water_in_trop) then
      use_manabe = photoset%use_manabe ! use manabe formula
      relative_humidity = photoset%relative_humidity ! relative humidity if no manabe
      gas_rainout = photoset%gas_rainout
      if (stratospheric_cond) then
        relative_humidity_cold_trap = photoset%relative_humidity_cold_trap
        H2O_condensation_rate = photoset%H2O_condensation_rate
      endif
    else
      gas_rainout = .false.
    endif

  end subroutine

  
end module