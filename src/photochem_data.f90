
module photochem_data
  implicit none  
  ! public ! public but protected
  integer, private, parameter :: real_kind = kind(1.0d0)
  integer, private, parameter :: str_len = 1024
  integer, private, parameter :: err_len = 1024
  ! Data that doesn't after file read-in

  ! molecules
  integer, protected :: nsp
  integer, protected :: nq
  integer, protected :: kd, kl, ku ! not read in. It is nq + nq + 1 (diagonal width of jacobian)
  integer, protected :: lda ! not read in 
  integer, protected :: nsl
  integer, protected :: natoms
  character(len=8), allocatable, protected :: atoms_names(:) 
  real(real_kind), allocatable, protected :: atoms_mass(:) 
  character(len=8), allocatable, protected :: species_names(:)
  integer, allocatable, protected :: species_composition(:,:)
  real(real_kind), allocatable, protected :: species_mass(:) 
  real(real_kind), allocatable, protected :: thermo_data(:,:,:)
  real(real_kind), allocatable, protected :: thermo_temps(:,:)
  
  ! reactions
  logical, protected :: reverse
  integer, protected :: nrF ! number of forward reactions
  integer, protected :: nrR ! number of reverse reactions
  integer, protected :: nrT ! number of total reactions
  integer, protected :: max_num_reactants
  integer, protected :: max_num_products
  character(len=8), allocatable, protected :: reactants_names(:,:)
  character(len=8), allocatable, protected :: products_names(:,:)
  integer, allocatable, protected :: reactants_sp_inds(:,:) ! for getting species nums in reactions
  integer, allocatable, protected :: products_sp_inds(:,:)
  integer, allocatable, protected :: nreactants(:) ! number of reactants
  integer, allocatable, protected :: nproducts(:) ! number of products
  integer, allocatable, protected :: reverse_info(:) ! indexs between forward and reverse reactions
  integer, allocatable, protected :: rxtypes(:) ! 0 is photolysis, 1 is elementary, 2 is three-body, 3 is falloff
  real(real_kind), allocatable, protected :: rateparams(:,:) ! (10, nrF)
  real(real_kind), allocatable, protected :: efficiencies(:,:) ! (maxval(num_efficient), nrF)
  integer, allocatable, protected :: eff_sp_inds(:,:) ! (maxval(num_efficient), nrF)
  integer, allocatable, protected :: num_efficient(:) ! number of efficiencies for each reaction
  real(real_kind), allocatable, protected :: def_eff(:) ! default efficiency
  integer, allocatable, protected :: falloff_type(:) ! type of falloff function (0 = none, 1 = Troe without T2,..)
  integer, allocatable, protected :: nump(:) ! number of production mechanisms (rxns) for each sp
  integer, allocatable, protected :: numl(:) ! number of loss mechanisms (rxns) for each sp
  integer, allocatable, protected :: iprod(:,:) ! (nmax,nsp) returns reaction # of production mechanism for sp
  integer, allocatable, protected :: iloss(:,:) ! (nmax,nsp) returns reaction # of loss mechanism for sp
  integer, protected :: kj ! number of photolysis reactions
  integer, allocatable, protected :: photonums(:) ! (kj) the reaction number of each photolysis reaction

  ! raditative transfer
  integer, protected :: nw
  real(real_kind), allocatable, protected :: wavl(:) ! (nw+1)
  integer, allocatable, protected :: num_temp_cols(:) ! (kj)
  integer, allocatable, protected :: sum_temp_cols(:) ! (kj)
  ! All data for every reaction in single vector to save memory
  real(real_kind), allocatable, protected :: xs_data(:) ! (sum(num_temp_cols)*nw) 
  real(real_kind), allocatable, protected :: xs_data_temps(:,:) ! (maxval(num_temp_cols), kj)
  integer, protected :: nray
  real(real_kind), allocatable, protected :: sigray(:,:) ! (len(raynums), nw)
  integer, allocatable, protected :: raynums(:) ! species number of rayleigh species
    
  ! initial conditions  
  integer, protected :: nzf
  real(real_kind), allocatable, protected :: z_file(:)
  real(real_kind), allocatable, protected :: T_file(:)
  real(real_kind), allocatable, protected :: edd_file(:)
  real(real_kind), allocatable, protected :: usol_file(:,:)
  
  ! settings
  logical, protected :: back_gas
  real(real_kind), protected :: back_gas_mu
  character(len=str_len), protected :: back_gas_name
  integer, protected :: back_gas_ind
  real(real_kind), protected :: planet_mass
  real(real_kind), protected :: planet_radius
  logical, protected :: water_sat_trop
  integer, protected :: LH2O
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
      deallocate(thermo_data)
      deallocate(thermo_temps)
      
      deallocate(reactants_names)
      deallocate(products_names)
      deallocate(reactants_sp_inds)
      deallocate(products_sp_inds)
      deallocate(nreactants)
      deallocate(nproducts)
      deallocate(reverse_info)
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
      deallocate(num_temp_cols)
      deallocate(sum_temp_cols)
      deallocate(xs_data)
      deallocate(xs_data_temps)
      deallocate(sigray)
      deallocate(raynums)
      
      deallocate(z_file)
      deallocate(T_file)
      deallocate(edd_file)
      deallocate(usol_file)
    
    endif
    
    ! species
    nsp = photomech%nsp
    nq = photomech%nq
    kd = 2*nq + 1
    kl = kd + nq
    ku = kd - nq
    lda = 3*nq + 1
    nsl = photomech%nsl
    natoms = photomech%natoms
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
      allocate(thermo_data(7,2,nsp))
      thermo_data = photomech%thermo_data
      allocate(thermo_temps(3,nsp))
      thermo_temps = photomech%thermo_temps
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
    allocate(num_temp_cols(kj))
    num_temp_cols = photorad%num_temp_cols
    allocate(sum_temp_cols(kj))
    sum_temp_cols = photorad%sum_temp_cols
    allocate(xs_data(sum(num_temp_cols)*nw))
    xs_data = photorad%xs_data
    allocate(xs_data_temps(maxval(num_temp_cols), kj))
    xs_data_temps = photorad%xs_data_temps
    nray = photorad%nray
    allocate(sigray(size(photorad%sigray,1),nw))
    sigray = photorad%sigray
    allocate(raynums(size(photorad%raynums)))
    raynums = photorad%raynums
    
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
    no_water_profile = photoinit%no_water_profile
    
    ! settings
    back_gas = photoset%back_gas
    back_gas_mu = photoset%back_gas_mu
    back_gas_name = photoset%back_gas_name
    planet_mass = photoset%planet_mass
    planet_radius = photoset%planet_radius
    water_sat_trop = photoset%water_sat_trop
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
    if (water_sat_trop) then
      trop_alt = photoset%trop_alt
    else
      trop_alt = -1.d0
    endif
    
    allocate(photon_flux(nw))
    photon_flux = photorad%photon_flux
    photon_scale_factor = photoset%photon_scale_factor  
    
    ! settings
    use_manabe = photoset%use_manabe ! use manabe formula
    relative_humidity = photoset%relative_humidity ! relative humidity if no manabe

  end subroutine

  
end module