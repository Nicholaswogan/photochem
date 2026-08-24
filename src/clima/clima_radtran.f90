module clima_radtran
  use clima_const, only: dp, s_str_len
  use clima_radtran_types, only: OpticalProperties, OpticalPropertiesWork, OpticalPropertiesResult
  use clima_radtran_types, only: RTChannel, RadiateWork
  implicit none
  private

  public :: ClimaRadtranWrk ! Work type that holds results
  public :: Radtran ! IR and solar radiative transfer
  
  type :: ClimaRadtranWrk
    
    !> (nz+1,nw) mW/m2/Hz in each wavelength bin
    !> at the edges of the vertical grid
    real(dp), allocatable :: fup_a(:,:), fdn_a(:,:)
    !> (nz+1) mW/m2 at the edges of the vertical grid 
    !> (integral of fup_a and fdn_a over wavelength grid)
    real(dp), allocatable :: fup_n(:), fdn_n(:)
    !> (nz+1,nw) Mean intensity in photons/cm^2/s. Only used
    !> in solar radiative transfer.
    real(dp), allocatable :: amean(:,:)
    !> Band optical thickness (nz,nw)
    real(dp), allocatable :: tau_band(:,:)
    
  end type
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR and Solar Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type :: Radtran

    integer :: ng
    character(s_str_len), allocatable :: species_names(:) ! (ng) copy of species names

    integer :: np
    character(s_str_len), allocatable :: particle_names(:)
  
    ! number of layers
    integer :: nz
  
    !!! Optical properties !!!
    type(OpticalProperties) :: op
    type(OpticalPropertiesWork) :: opw
    type(OpticalPropertiesResult) :: opr

    type(RTChannel) :: ir
    type(RTChannel) :: sol
    type(RadiateWork) :: rw
    
    real(dp) :: diurnal_fac = 0.5_dp
    real(dp), allocatable :: zenith_u(:) !! cosine of the zenith angle in radians
    real(dp), allocatable :: zenith_weights(:)
    !> surface albedo in each solar wavelength bin (sol%nw) 
    real(dp), allocatable :: surface_albedo(:) 
    !> surface emissivity in each IR wavelength bin (ir%nw) 
    real(dp), allocatable :: surface_emissivity(:) 
    !> If true, use a hard-surface lower thermal boundary; if false, use
    !> a gas-giant style no-hard-surface diffusion boundary.
    logical :: has_hard_surface = .true.
    !> Minimum optical depth used in IR two-stream thin-layer guard.
    real(dp) :: ir_tau_min = 1.0e-6_dp
    !> (nw) mW/m2/Hz in each bin from the input star file. This is 
    !> later scaled by the variable `photon_scale_factor`.
    real(dp), allocatable :: photons_sol(:)
    !> A scale factor that is applied to `photons_sol` so that
    !> bolometric luminosity can be easily changed.
    real(dp) :: photon_scale_factor = 1.0_dp 
  
    type(ClimaRadtranWrk) :: wrk_ir
    type(ClimaRadtranWrk) :: wrk_sol
    real(dp), allocatable :: f_total(:)
    
  contains
    procedure :: radiate => Radtran_radiate
    procedure :: TOA_fluxes => Radtran_TOA_fluxes
    procedure :: set_bolometric_flux => Radtran_set_bolometric_flux
    procedure :: bolometric_flux => Radtran_bolometric_flux
    procedure :: skin_temperature => Radtran_skin_temperature
    procedure :: equilibrium_temperature => Radtran_equilibrium_temperature
    procedure :: opacities2yaml => Radtran_opacities2yaml
    procedure :: apply_radiation_enhancement => Radtran_apply_radiation_enhancement
    procedure :: set_custom_optical_properties => Radtran_set_custom_optical_properties
    procedure :: unset_custom_optical_properties => Radtran_unset_custom_optical_properties
  end type

  interface Radtran
    module procedure :: create_Radtran_1
    module procedure :: create_Radtran_2
  end interface
  
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR and Solar Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function create_Radtran_1(settings_f, star_f, num_zenith_angles, surface_albedo, nz, datadir, err) result(rad)
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    integer, intent(in) :: num_zenith_angles
    real(dp), intent(in) :: surface_albedo
    integer, intent(in) :: nz
    character(*), intent(in) :: datadir
    character(:), allocatable, intent(out) :: err
    
    type(Radtran) :: rad
    
    type(ClimaSettings) :: s
    
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    if (.not. allocated(s%gases)) then
      err = '"'//settings_f//'/optical-properties/gases" does not exist'
      return
    endif

    if (.not. allocated(s%particles)) allocate(s%particles(0))
    
    rad = create_Radtran_2(s%gases, s%particles, s, star_f, num_zenith_angles, surface_albedo, nz, datadir, err)
    if (allocated(err)) return
    
  end function
  
  function create_Radtran_2(species_names, particle_names, s, star_f, &
                            num_zenith_angles, surface_albedo, nz, datadir, err) result(rad)
    use clima_radtran_types, only: SolarChannel, IRChannel, read_stellar_flux
    use clima_types, only: ClimaSettings
    use clima_eqns, only: zenith_angles_and_weights
    use clima_const, only: pi
    
    character(*), intent(in) :: species_names(:)
    character(*), intent(in) :: particle_names(:)
    type(ClimaSettings), intent(in) :: s
    character(*), intent(in) :: star_f
    integer, intent(in) :: num_zenith_angles
    real(dp), intent(in) :: surface_albedo
    integer, intent(in) :: nz
    character(*), intent(in) :: datadir
    character(:), allocatable, intent(out) :: err
    
    type(Radtran) :: rad

    real(dp) :: photon_scale_factor
    
    if (nz < 1) then
      err = '"nz" can not be less than 1.'
      return
    endif
    
    rad%ng = size(species_names)
    rad%species_names = species_names
    rad%np = size(particle_names)
    if (rad%np > 0) then
      rad%particle_names = particle_names
    endif
    rad%nz = nz

    allocate(rad%zenith_u(num_zenith_angles))
    allocate(rad%zenith_weights(num_zenith_angles))
    call zenith_angles_and_weights(num_zenith_angles, rad%zenith_u, rad%zenith_weights)
    rad%zenith_u = cos(rad%zenith_u*pi/180.0_dp)

    if (.not. allocated(s%op)) then
      err = '"'//s%filename//'/optical-properties" does not contain opacity information.'
      return
    endif
    rad%op = OpticalProperties(datadir, species_names, particle_names, s%op, err)
    if (allocated(err)) return
    rad%opw = OpticalPropertiesWork(rad%op, rad%nz)
    rad%opr = OpticalPropertiesResult(rad%op, rad%nz)
    
    rad%ir = RTChannel(datadir, IRChannel, s%wavelength_bins_file, rad%op, err)
    if (allocated(err)) return
    rad%sol = RTChannel(datadir, SolarChannel, s%wavelength_bins_file, rad%op, err)
    if (allocated(err)) return
    rad%rw = RadiateWork(rad%op%nw, rad%nz)

    allocate(rad%surface_albedo(rad%sol%nw))
    rad%surface_albedo(:) = surface_albedo

    allocate(rad%surface_emissivity(rad%ir%nw))
    rad%surface_emissivity(:) = 1.0_dp

    if (s%planet_is_present) then
      photon_scale_factor = s%photon_scale_factor
    else
      photon_scale_factor = 1.0_dp
    endif
    rad%photon_scale_factor = photon_scale_factor
    ! photons hitting the planet
    allocate(rad%photons_sol(rad%sol%nw))
    call read_stellar_flux(star_f, rad%sol%nw, rad%sol%wavl, rad%photons_sol, err)
    if (allocated(err)) return

    ! IR work arrays
    allocate(rad%wrk_ir%fup_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fdn_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fup_n(nz+1))
    allocate(rad%wrk_ir%fdn_n(nz+1))
    allocate(rad%wrk_ir%amean(nz+1, rad%ir%nw))
    rad%wrk_ir%amean = 0.0_dp ! not used
    allocate(rad%wrk_ir%tau_band(nz,rad%ir%nw))

    ! Solar work arrays
    allocate(rad%wrk_sol%fup_a(nz+1, rad%sol%nw))
    allocate(rad%wrk_sol%fdn_a(nz+1, rad%sol%nw))
    allocate(rad%wrk_sol%fup_n(nz+1))
    allocate(rad%wrk_sol%fdn_n(nz+1))
    allocate(rad%wrk_sol%amean(nz+1, rad%sol%nw))
    allocate(rad%wrk_sol%tau_band(nz,rad%sol%nw))

    ! total flux
    allocate(rad%f_total(nz+1))

  end function

  subroutine Radtran_radiate(self, T_surface, T, P, densities, dz, pdensities, radii, compute_solar, compute_opacity, err)
    use clima_radtran_radiate, only: radiate
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), optional, target, intent(in) :: pdensities(:,:), radii(:,:)
    logical, optional, intent(in) :: compute_solar
    logical, optional, intent(in) :: compute_opacity
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err

    logical :: compute_solar_, compute_opacity_

    type(ClimaRadtranWrk), pointer :: wrk_ir
    type(ClimaRadtranWrk), pointer :: wrk_sol
    
    wrk_ir => self%wrk_ir
    wrk_sol => self%wrk_sol                                    
    
    call check_inputs(self%nz, self%ng, self%np, T, P, densities, dz, pdensities, radii, err)
    if (allocated(err)) return

    compute_solar_ = .true.
    if (present(compute_solar)) then
      compute_solar_ = compute_solar
    endif
    compute_opacity_ = .true.
    if (present(compute_opacity)) then
      compute_opacity_ = compute_opacity
    endif

    if (compute_opacity_) then
      ! First compute opacity
      call self%op%compute_opacity(P, T, densities, dz, pdensities, radii, self%opw, self%opr, err)
      if (allocated(err)) return
    endif

    ! IR radiative transfer                                     
    call radiate( &
      rtc=self%ir, &
      op=self%op, &
      opr=self%opr, &
      rw=self%rw, &
      surface_emissivity=self%surface_emissivity, &
      surface_albedo=self%surface_albedo, &
      has_hard_surface=self%has_hard_surface, &
      ir_tau_min=self%ir_tau_min, &
      diurnal_fac=0.0_dp, &
      photons_sol=[0.0_dp], &
      zenith_u=[0.0_dp], &
      zenith_weights=[1.0_dp], &
      T_surface=T_surface, &
      T=T, &
      fup_a=wrk_ir%fup_a, &
      fdn_a=wrk_ir%fdn_a, &
      fup_n=wrk_ir%fup_n, &
      fdn_n=wrk_ir%fdn_n, &
      amean=wrk_ir%amean, &
      tau_band=wrk_ir%tau_band &
    )

    ! If we don't want to compute solar RT, then we compute f_total and return
    if (.not.compute_solar_) then
      self%f_total = (wrk_sol%fdn_n - wrk_sol%fup_n) + (wrk_ir%fdn_n - wrk_ir%fup_n)
      return
    endif

    ! Solar radiative transfer                                    
    call radiate( &
      rtc=self%sol, &
      op=self%op, &
      opr=self%opr, &
      rw=self%rw, &
      surface_emissivity=self%surface_emissivity, &
      surface_albedo=self%surface_albedo, &
      has_hard_surface=self%has_hard_surface, &
      ir_tau_min=self%ir_tau_min, &
      diurnal_fac=self%diurnal_fac, &
      photons_sol=self%photons_sol*self%photon_scale_factor, &
      zenith_u=self%zenith_u, &
      zenith_weights=self%zenith_weights, &
      T_surface=T_surface, &
      T=T, &
      fup_a=wrk_sol%fup_a, &
      fdn_a=wrk_sol%fdn_a, &
      fup_n=wrk_sol%fup_n, &
      fdn_n=wrk_sol%fdn_n, &
      amean=wrk_sol%amean, &
      tau_band=wrk_sol%tau_band &
    )

    ! Total flux
    self%f_total = (wrk_sol%fdn_n - wrk_sol%fup_n) + (wrk_ir%fdn_n - wrk_ir%fup_n)

  end subroutine

  subroutine Radtran_TOA_fluxes(self, T_surface, T, P, densities, dz, pdensities, radii, &
      compute_solar, compute_opacity, ISR, OLR, err)
    use clima_radtran_radiate, only: radiate
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    real(dp), optional, target, intent(in) :: pdensities(:,:), radii(:,:) !! (nz,np)
    logical, optional, intent(in) :: compute_solar
    logical, optional, intent(in) :: compute_opacity
    real(dp), intent(out) :: ISR, OLR
    character(:), allocatable, intent(out) :: err
      
    call self%radiate(T_surface, T, P, densities, dz, pdensities, radii, compute_solar, compute_opacity, err)
    if (allocated(err)) return
    
    ISR = (self%wrk_sol%fdn_n(self%nz+1) - self%wrk_sol%fup_n(self%nz+1)) ! Incoming short wave
    OLR = - (self%wrk_ir%fdn_n(self%nz+1) - self%wrk_ir%fup_n(self%nz+1)) ! Outgoing long wave

  end subroutine

  !> Sets the bolometric stellar flux by adjusting the `photon_scale_factor`.
  subroutine Radtran_set_bolometric_flux(self, flux)
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: flux !! Bolometric flux (W/m^2)
    self%photon_scale_factor = 1.0_dp
    self%photon_scale_factor = flux/self%bolometric_flux()
  end subroutine

  !> The bolometric stellar flux at the planet in W/m^2
  function Radtran_bolometric_flux(self) result(flux)
    class(Radtran), target, intent(inout) :: self
    real(dp) :: flux !! Bolometric flux (W/m^2)
    integer :: i

    flux = 0.0_dp
    do i = 1,self%sol%nw
      flux = flux + self%photons_sol(i)*(self%sol%freq(i) - self%sol%freq(i+1))
    enddo
    flux = self%photon_scale_factor*flux/1.0e3_dp

  end function

  !> The skin temperature
  function Radtran_skin_temperature(self, bond_albedo) result(T_skin)
    use clima_eqns, only: skin_temperature
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: bond_albedo !! The bond albedo of a planet
    real(dp) :: T_skin
    T_skin = skin_temperature(self%bolometric_flux(), bond_albedo)
  end function

  !> The equilibrium temperature
  function Radtran_equilibrium_temperature(self, bond_albedo) result(T_eq)
    use clima_eqns, only: equilibrium_temperature
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: bond_albedo !! The bond albedo of a planet
    real(dp) :: T_eq
    T_eq = equilibrium_temperature(self%bolometric_flux(), bond_albedo)
  end function

  !> Returns a yaml string representing all opacities in the model.
  function Radtran_opacities2yaml(self) result(out)
    class(Radtran), target, intent(inout) :: self
    character(:), allocatable :: out

    character(:), allocatable :: line

    out = ''
    
    line = ''
    line = line//'optical-properties:'
    out = out//line

    out = out//new_line('(a)')
    out = out//self%op%opacities2yaml()

  end function

  subroutine Radtran_apply_radiation_enhancement(self, rad_enhancement)
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: rad_enhancement
    self%wrk_sol%fdn_n = self%wrk_sol%fdn_n*rad_enhancement
    self%wrk_sol%fdn_a = self%wrk_sol%fdn_a*rad_enhancement
    self%wrk_sol%fup_n = self%wrk_sol%fup_n*rad_enhancement
    self%wrk_sol%fup_a = self%wrk_sol%fup_a*rad_enhancement
    self%f_total = (self%wrk_sol%fdn_n - self%wrk_sol%fup_n) + &
                   (self%wrk_ir%fdn_n - self%wrk_ir%fup_n)
  end subroutine

  !!!!!!!!!!!!!!!!!
  !!! Utilities !!!
  !!!!!!!!!!!!!!!!!

  subroutine check_inputs(nz, ng, np, T, P, densities, dz, pdensities, radii, err)
    integer, intent(in) :: nz, ng, np
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: densities(:,:)
    real(dp), intent(in) :: dz(:)
    real(dp), optional, target, intent(in) :: pdensities(:,:), radii(:,:)
    character(:), allocatable, intent(out) :: err

    if ((present(pdensities) .and. .not. present(radii)) &
      .or.(present(radii) .and. .not. present(pdensities))) then
      err = 'Both pdensities and radii must be arguments.'
      return
    endif
    if (np > 0) then
      if (.not. present(radii)) then
        err = 'The model contains particles but "pdensities" and "radii" are not arguments.'
        return
      endif
    endif
    
    call check_dimensions(nz, ng, T, P, densities, dz, err)
    if (allocated(err)) return
    if (present(radii)) then
      call check_dimensions_p(nz, np, pdensities, radii, err)
      if (allocated(err)) return
    endif

  end subroutine

  subroutine check_dimensions_p(nz, np, pdensities, radii, err)
    integer, intent(in) :: nz, np
    real(dp), intent(in) :: pdensities(:,:)
    real(dp), intent(in) :: radii(:,:)
    character(:), allocatable, intent(out) :: err

    if (size(pdensities,1) /= nz .or. size(pdensities,2) /= np) then
      err = '"pdensities" has the wrong input dimension.'
      return
    endif

    if (size(radii,1) /= nz .or. size(radii,2) /= np) then
      err = '"radii" has the wrong input dimension.'
      return
    endif

  end subroutine
  
  subroutine check_dimensions(nz, ng, T, P, densities, dz, err)
    integer, intent(in) :: nz, ng
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err
    
    if (size(T) /= nz) then
      err = '"T" has the wrong input dimension.'
      return
    endif
    if (size(P) /= nz) then
      err = '"P" has the wrong input dimension.'
      return
    endif
    if (size(densities,1) /= nz .or. size(densities,2) /= ng) then
      err = '"densities" has the wrong input dimension.'
      return
    endif
    if (size(dz) /= nz) then
      err = '"dz" has the wrong input dimension.'
      return
    endif
    
  end subroutine

  !> Sets custom optical properties
  subroutine Radtran_set_custom_optical_properties(self, wv, P, dtau_dz, w0, g0, err)
    class(Radtran), intent(inout) :: self
    real(dp), intent(in) :: wv(:) !! Array of of wavelengths in nm
    real(dp), intent(in) :: P(:) !! Array of pressures in dynes/cm^2. Must be decreasing.
    real(dp), intent(in) :: dtau_dz(:,:) !! (size(P),size(wv)), Optical depth per altitude (1/cm).
    real(dp), intent(in) :: w0(:,:) !! (size(P),size(wv)), Single scattering albedo
    real(dp), intent(in) :: g0(:,:) !! (size(P),size(wv)), Asymetry parameter
    character(:), allocatable, intent(out) :: err

    call self%op%set_custom_optical_properties(wv, P, dtau_dz, w0, g0, err)
    if (allocated(err)) return

  end subroutine

  !> Unsets custom optical properties set with `set_custom_optical_properties`.
  subroutine Radtran_unset_custom_optical_properties(self)
    class(Radtran), intent(inout) :: self
    call self%op%unset_custom_optical_properties()
  end subroutine
  
end module
  
