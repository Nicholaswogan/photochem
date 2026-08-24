module clima_adiabat
  use clima_const, only: dp, s_str_len
  use clima_types, only: Species
  use clima_radtran, only: Radtran
  use linear_interpolation_module, only: linear_interp_1d
  use clima_eqns, only: ocean_solubility_fcn, temp_dependent_albedo_fcn
  use clima_adiabat_general, only: OceanFunction
  use iso_c_binding, only: c_ptr, c_null_ptr
  implicit none
  private

  public :: AdiabatClimate
  public :: RCE_SOLVE_HYBRJ_ONLY, RCE_SOLVE_PTC_THEN_HYBRJ, RCE_SOLVE_HYBRJ_THEN_PTC_THEN_HYBRJ

  integer, parameter :: RCE_SOLVE_HYBRJ_ONLY = 1
  integer, parameter :: RCE_SOLVE_PTC_THEN_HYBRJ = 2
  integer, parameter :: RCE_SOLVE_HYBRJ_THEN_PTC_THEN_HYBRJ = 3

  type :: AdiabatClimate

    ! settings and free parameters
    integer :: nz
    real(dp) :: P_top = 1.0_dp !! (dynes/cm2)
    real(dp) :: T_trop = 180.0_dp !! (T)
    real(dp), allocatable :: RH(:) !! relative humidity (ng)
    !> If .true., then any function that calls `make_column` will
    !> use the initial guess in `self%make_column_P_guess`
    logical :: use_make_column_P_guess = .true.
    !> Initial guess for surface pressure of all gases for `make_column`
    !> routine (ng).
    real(dp), allocatable :: make_column_P_guess(:)

    !> If .true., then Tropopause temperature is non-linearly solved for such that
    !> it matches the skin temperature. The initial guess will always be self%T_trop.
    logical :: solve_for_T_trop = .false.
    !> Callback that sets the surface albedo based on the surface temperature.
    !> This can be used to parameterize the ice-albedo feedback.
    procedure(temp_dependent_albedo_fcn), nopass, pointer :: albedo_fcn => null()

    !> Function describing how gases disolve in oceans. This allows for multiple oceans, each
    !> made of different condensed volatiles.
    type(OceanFunction), allocatable :: ocean_fcns(:)
    !> A c-ptr which will be passed to each ocean solubility function.
    type(c_ptr) :: ocean_args_p = c_null_ptr

    ! Heat re-distribution terms for Equation (10) in Koll (2020), ApJ when
    ! considering tidally locked planets. Default values are from the paper.
    !> If true, then will attempt to compute the climate corresponding to the
    !> observed dayside temperature of a tidally locked planet.
    logical :: tidally_locked_dayside = .false.
    real(dp) :: L !! = planet radius. Circulation’s horizontal scale (cm)
    real(dp) :: chi = 0.2_dp !! Heat engine efficiency term (no units)
    real(dp) :: n_LW = 2.0_dp !! = 1 or 2 (no units)
    real(dp) :: Cd = 1.9e-3_dp !! Drag coefficient (no units)

    !> Heat flow from the surface into the atmosphere (mW/m^2)
    real(dp) :: surface_heat_flow = 0.0_dp
    
    ! planet properties
    real(dp) :: planet_mass !! (g)
    real(dp) :: planet_radius !! (cm)
    !> Reference pressure for which `planet_radius` is defined (dynes/cm^2).
    !> If <= 0, then `planet_radius` is defined at the surface pressure.
    real(dp) :: reference_pressure = -1.0_dp
    
    ! species in the model
    character(s_str_len), allocatable :: species_names(:) !! copy of species names
    character(s_str_len), allocatable :: particle_names(:) !! copy of particle names
    type(Species) :: sp
    
    ! Radiative transfer
    type(Radtran) :: rad
    logical :: double_radiative_grid
    integer :: nz_r
    ! Work space for radiative transfer
    real(dp), allocatable :: T_r(:)
    real(dp), allocatable :: P_r(:)
    real(dp), allocatable :: densities_r(:,:)
    real(dp), allocatable :: dz_r(:)
    real(dp), allocatable :: pdensities_r(:,:)
    real(dp), allocatable :: pradii_r(:,:)

    ! For custom mixing ratios
    !> Length ng. If true, then the species has a custom dry mixing ratio
    logical, private, allocatable :: sp_custom(:)
    !> Length ng. Contains 1-D interpolators for log10P (dynes/cm^2) and 
    !> log10mix (vmr) for each custom species, if specified
    type(linear_interp_1d), private, allocatable :: mix_custom(:)

    ! Information about convection
    integer :: n_convecting_zones !! number of convecting zones
    !> Describes the lower index of each convecting zone. Note index 1
    !> is the ground layer.
    integer, private, allocatable :: ind_conv_lower(:) 
    !> Describes the upper index of each convecting zone.
    integer, private, allocatable :: ind_conv_upper(:)
    integer, private, allocatable :: ind_conv_lower_x(:)
    !> Indicates if a layer is super-saturated
    logical, private, allocatable :: super_saturated(:)
    !> Another representation of where convection is occuring. If True,
    !> then the layer below is convecting with the current layer. Index 1 determines
    !> if the first atomspheric layer is convecting with the ground.
    logical, allocatable :: convecting_with_below(:)
    integer, allocatable :: inds_Tx(:)
    real(dp), allocatable :: lapse_rate(:) !! The true lapse rate (dlnT/dlnP)
    real(dp), allocatable :: lapse_rate_intended(:) !! The computed lapse rate (dlnT/dlnP)
    !> Fraction of the Newton step used in convective classification (0..1).
    real(dp) :: convective_newton_step_size = 1.0e-1_dp
    ! Hysteresis parameters for updating convective mask in RCE.
    ! A radiative layer becomes convective only if superadiabaticity exceeds
    ! `max(convective_hysteresis_min, convective_hysteresis_frac_on*|lapse_rate_intended|)`.
    ! A convective layer becomes radiative only if it is stable by more than
    ! `max(convective_hysteresis_min, convective_hysteresis_frac_off*|lapse_rate_intended|)`.
    real(dp) :: convective_hysteresis_frac_on = 2.0e-2_dp
    real(dp) :: convective_hysteresis_frac_off = 2.0e-2_dp
    real(dp) :: convective_hysteresis_min = 1.0e-3_dp
    !> Boundary-motion limiter for convective mask updates. If < 0, no limiter
    !> is applied. Set to 1 to limit expansion of convective zones to 1 layer at
    !> a time
    integer :: convective_max_boundary_shift = -1
    !> If true, shrink convective tops when a strong inversion exists just above.
    logical :: prevent_overconvection = .true.
    !> If true, require passing through mode 2 when mode 1 converges
    !> before optional mode 3 polishing.
    logical :: require_mode2 = .true.
    !> Per-layer lockout counter used by prevent_overconvection polishing.
    !> A positive value temporarily prevents immediate re-demotion of the same
    !> convective top boundary to avoid ABAB oscillations.
    integer, allocatable :: prevent_overconvection_lock(:)

    ! tolerances
    !> Relative tolerance of integration
    real(dp) :: rtol = 1.0e-9_dp
    !> Absolute tolerance of integration
    real(dp) :: atol = 1.0e-12_dp
    !> Tolerance for nonlinear solve in make_column
    real(dp) :: tol_make_column = 1.0e-8_dp
    !> Perturbation for the jacobian
    real(dp) :: epsj = 1.0e-2_dp
    !> xtol for RC equilibrium
    real(dp) :: xtol_rc = 1.0e-5_dp
    !> Multiplicative growth factor for PTC timestep updates.
    !> Values > 1 increase `dt` after successful steps.
    real(dp) :: dt_increment = 1.5_dp
    !> Max number of iterations in the RCE routine
    integer :: max_rc_iters = 30
    !> Max number of iterations for which convective layers can
    !> be converged to radiative layers in the RCE routine
    integer :: max_rc_iters_convection = 5
    !> If False, then the jacobian calculation in RCE does not recompute
    !> solar radiative transfer.
    logical :: compute_solar_in_jac = .false.
    !> Strategy for nonlinear solve inside each RCE outer iteration.
    !> 1 = HYBRJ only.
    !> 2 = PTC then HYBRJ tighten if no convergence.
    !> 3 = HYBRJ first, then fallback to PTC then HYBRJ tighten if no convergence.
    integer :: rce_solve_strategy = RCE_SOLVE_HYBRJ_THEN_PTC_THEN_HYBRJ
    logical :: verbose = .true. !! verbosity
    
    ! State of the atmosphere
    real(dp), allocatable :: f_i_surf(:) !! Surface mixing ratios (ng)
    real(dp) :: P_surf !! Surface pressure (dynes/cm^2)
    real(dp) :: P_trop !! Tropopause pressure (dynes/cm^2)
    real(dp), allocatable :: P(:) !! pressure in each grid cell, dynes/cm^2 (nz)
    real(dp) :: T_surf !! Surface Temperature (K)
    real(dp), allocatable :: T(:) !! Temperature in each grid cell, K (nz) 
    real(dp), allocatable :: f_i(:,:) !! mixing ratios of species in each grid cell (nz,ng)
    real(dp), allocatable :: z(:) !! Altitude at the center of the grid cell, cm (nz)
    real(dp), allocatable :: dz(:) !! Thickness of each grid cell, cm (nz)
    real(dp) :: gravity_surf !! Surface gravity (cm/s^2)
    real(dp), allocatable :: gravity(:) !! Gravity at each layer center (cm/s^2) (nz)
    real(dp), allocatable :: densities(:,:) !! densities in each grid cell, molecules/cm^3 (nz,ng)
    real(dp), allocatable :: N_atmos(:) !! reservoir of gas in atmosphere mol/cm^2 (ng)
    !> reservoir of gas on surface mol/cm^2 (ng).
    !> Strictly consistent when `reference_pressure <= 0` (i.e., `planet_radius` at `P_surf`).
    real(dp), allocatable :: N_surface(:)
    !> reservoir of gas dissolved in oceans in mol/cm^2 (ng, ng). There can be multiple oceans.
    !> The gases dissolved in ocean made of species 1 is given by `N_ocean(:,1)`.
    !> Strictly consistent when `reference_pressure <= 0` (i.e., `planet_radius` at `P_surf`).
    real(dp), allocatable :: N_ocean(:,:)
    !> particle densities in particles/cm^3. (nz,np)
    real(dp), allocatable :: pdensities(:,:)
    !> For interpolating particle densities
    type(linear_interp_1d), private, allocatable :: pdensities_interp(:)
    !> particle radii in cm. (nz,np)
    real(dp), allocatable :: pradii(:,:)
    !> For interpolating particle radii
    type(linear_interp_1d), private, allocatable :: pradii_interp(:)

  contains
    ! Constructs atmospheres
    procedure :: make_profile => AdiabatClimate_make_profile
    procedure :: make_column => AdiabatClimate_make_column
    procedure :: make_profile_bg_gas => AdiabatClimate_make_profile_bg_gas
    procedure :: make_profile_dry => AdiabatClimate_make_profile_dry
    ! Constructs atmosphere and does radiative transfer
    procedure :: TOA_fluxes => AdiabatClimate_TOA_fluxes
    procedure :: TOA_fluxes_column => AdiabatClimate_TOA_fluxes_column
    procedure :: TOA_fluxes_bg_gas => AdiabatClimate_TOA_fluxes_bg_gas
    procedure :: TOA_fluxes_dry => AdiabatClimate_TOA_fluxes_dry
    ! Non-linear solves for equilibrium climate
    procedure :: surface_temperature => AdiabatClimate_surface_temperature
    procedure :: surface_temperature_column => AdiabatClimate_surface_temperature_column
    procedure :: surface_temperature_bg_gas => AdiabatClimate_surface_temperature_bg_gas
    ! Routines for full radiative convective equilibrium
    procedure, private :: make_profile_rc => AdiabatClimate_make_profile_rc
    procedure :: RCE => AdiabatClimate_RCE
    ! Utilities
    procedure, private :: interpolate_particles => AdiabatClimate_interpolate_particles
    procedure, private :: compute_altitude => AdiabatClimate_compute_altitude
    procedure :: set_particle_density_and_radii => AdiabatClimate_set_particle_density_and_radii
    procedure :: set_ocean_solubility_fcn => AdiabatClimate_set_ocean_solubility_fcn
    procedure :: to_regular_grid => AdiabatClimate_to_regular_grid
    procedure :: out2atmosphere_txt => AdiabatClimate_out2atmosphere_txt
    ! For tidally locked planets
    procedure :: heat_redistribution_parameters => AdiabatClimate_heat_redistribution_parameters

    ! Private methods
    procedure, private :: copy_atm_to_radiative_grid => AdiabatClimate_copy_atm_to_radiative_grid
  end type
  
  interface AdiabatClimate
    module procedure :: create_AdiabatClimate
  end interface

  interface
    module subroutine AdiabatClimate_make_profile_rc(self, P_i_surf, T_in, err)
      class(AdiabatClimate), intent(inout) :: self
      real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
      real(dp), intent(in) :: T_in(:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine AdiabatClimate_compute_altitude(self, err)
      class(AdiabatClimate), intent(inout) :: self
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Compute full radiative-convective equilibrium.
    module function AdiabatClimate_RCE(self, P_i_surf, T_surf_guess, T_guess, convecting_with_below, &
      sp_custom, P_custom, mix_custom, err) result(converged)
      class(AdiabatClimate), intent(inout) :: self
      !> Array of surface pressures of each species (dynes/cm^2)
      real(dp), intent(in) :: P_i_surf(:)
      !> A guess for the surface temperature (K)
      real(dp), intent(in) :: T_surf_guess
      !> A guess for the temperature in each atmospheric layer (K)
      real(dp), intent(in) :: T_guess(:)
      !> An array describing a guess for the radiative vs. convective 
      !> regions of the atmosphere
      logical, optional, intent(in) :: convecting_with_below(:)
      character(*), optional, intent(in) :: sp_custom(:)
      real(dp), optional, intent(in) :: P_custom(:)
      real(dp), optional, intent(in) :: mix_custom(:,:)
      character(:), allocatable, intent(out) :: err
      logical :: converged !! Whether the routine converged or not.
    end function
  end interface

  abstract interface
    subroutine toa_fluxes_fcn(T_surf, T_trop, ASR, OLR, err)
      import :: dp
      real(dp), intent(in) :: T_surf, T_trop
      real(dp), intent(out) :: ASR, OLR
      character(:), allocatable, intent(out) :: err
    end subroutine
  end interface

contains
  
  function create_AdiabatClimate(species_f, settings_f, star_f, data_dir, double_radiative_grid, err) result(c)
    use futils, only: linspace
    use clima_types, only: ClimaSettings 
    character(*), intent(in) :: species_f !! Species yaml file
    character(*), intent(in) :: settings_f !! Settings yaml file
    character(*), intent(in) :: star_f !! Star text file
    character(*), intent(in) :: data_dir !! Directory with radiative transfer data
    !> If true, radiative transfer is computed on a refined vertical grid.
    !> Each physical layer is split into two RT layers, and two ghost RT layers
    !> are added above the physical top to stabilize TOA behavior.
    logical, optional, intent(in) :: double_radiative_grid
    character(:), allocatable, intent(out) :: err
    
    type(AdiabatClimate) :: c
    
    type(ClimaSettings) :: s
    integer :: i

    ! species
    c%sp = Species(species_f, err)
    if (allocated(err)) return
    
    ! unpack species
    allocate(c%species_names(c%sp%ng))
    do i = 1,c%sp%ng
      c%species_names(i) = c%sp%g(i)%name
    enddo

    ! unpack particles
    allocate(c%particle_names(c%sp%np))
    do i = 1,c%sp%np
      c%particle_names(i) = c%sp%p(i)%name
    enddo

    ! default relative humidty is 1
    allocate(c%RH(c%sp%ng))
    c%RH(:) = 1.0_dp

    allocate(c%make_column_P_guess(c%sp%ng))
    c%make_column_P_guess(:) = 1.0_dp
    
    ! There must be more than 1 species
    if (c%sp%ng == 1) then 
      err = 'There must be more than 1 species in "'//species_f//'"'
      return
    endif
    
    ! settings
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    ! unpack setting
    if (s%atmos_grid_is_present) then
      c%nz = s%nz
    else
      err = '"atmosphere-grid" is missing from file "'//settings_f//'"'
      return
    endif
    if (.not.s%planet_is_present) then
      err = '"planet" is missing from file "'//settings_f//'"'
      return
    endif
    c%planet_mass = s%planet_mass
    c%planet_radius = s%planet_radius
    if (.not.allocated(s%number_of_zenith_angles)) then
      err = '"number-of-zenith-angles" is missing from file "'//settings_f//'"'
      return
    endif
    if (.not.allocated(s%surface_albedo)) then
      err = '"surface-albedo" is missing from file "'//settings_f//'"'
      return
    endif

    if (present(double_radiative_grid)) then
      c%double_radiative_grid = double_radiative_grid
    else
      c%double_radiative_grid = .true.
    endif
    if (c%double_radiative_grid) then
      ! Doubled RT grid plus two ghost layers above the physical top layer.
      c%nz_r = c%nz*2 + 2
    else
      c%nz_r = c%nz
    endif
    allocate(c%T_r(c%nz_r),c%P_r(c%nz_r),c%densities_r(c%nz_r,c%sp%ng),c%dz_r(c%nz_r))
    allocate(c%pdensities_r(c%nz_r,c%sp%np),c%pradii_r(c%nz_r,c%sp%np))
    ! Make radiative transfer with list of species, from species file
    ! and the optical-properties from the settings file
    c%rad = Radtran(c%species_names, c%particle_names, s, star_f, s%number_of_zenith_angles, &
                    s%surface_albedo, c%nz_r, data_dir, err)
    if (allocated(err)) return

    ! Custom mixing ratios
    allocate(c%sp_custom(c%sp%ng))
    allocate(c%mix_custom(c%sp%ng))

    ! Convection
    allocate(c%super_saturated(c%nz))
    allocate(c%convecting_with_below(c%nz))
    allocate(c%lapse_rate(c%nz),c%lapse_rate_intended(c%nz))
    allocate(c%prevent_overconvection_lock(c%nz))

    ! allocate ocean functions
    allocate(c%ocean_fcns(c%sp%ng))

    ! allocate work variables
    allocate(c%f_i_surf(c%sp%ng), c%P(c%nz), c%T(c%nz), c%f_i(c%nz,c%sp%ng), c%z(c%nz), c%dz(c%nz), &
             c%gravity(c%nz))
    allocate(c%densities(c%nz,c%sp%ng))
    allocate(c%N_atmos(c%sp%ng),c%N_surface(c%sp%ng),c%N_ocean(c%sp%ng,c%sp%ng))
    allocate(c%pdensities(c%nz,c%sp%np),c%pradii(c%nz,c%sp%np))
    allocate(c%pdensities_interp(c%sp%np),c%pradii_interp(c%sp%np))

    ! Initialize particle density interpolators with no particles
    call linspace(0.0_dp, -5.0_dp, c%P)
    c%P = 10.0**c%P
    c%pdensities = 0.0_dp
    c%pradii = 1.0e-4_dp ! 1 micron
    call c%set_particle_density_and_radii(c%P, c%pdensities, c%pradii, err)
    if (allocated(err)) return
    c%P = 0.0_dp

    ! Heat redistribution parameter
    c%L = c%planet_radius

  end function
  
  !> Constructs an atmosphere using a multispecies pseudoadiabat (Eq. 1 in Graham et al. 2021, PSJ)
  !> troposphere connected to an isothermal stratosphere. Input are surface partial pressures
  !> of each gas.
  subroutine AdiabatClimate_make_profile(self, T_surf, P_i_surf, err)
    use clima_adiabat_general, only: make_profile
    use clima_const, only: k_boltz, N_avo
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: P_e(:), z_e(:), T_e(:), f_i_e(:,:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    if (size(P_i_surf) /= self%sp%ng) then
      err = "P_i_surf has the wrong dimension"
      return
    endif
    
    allocate(P_e(2*self%nz+1),  z_e(2*self%nz+1), T_e(2*self%nz+1))
    allocate(f_i_e(2*self%nz+1,self%sp%ng))
    allocate(density(self%nz))

    call make_profile(T_surf, P_i_surf, &
                      self%sp, self%nz, self%planet_mass, &
                      self%planet_radius, self%P_top, self%T_trop, self%RH, &
                      self%rtol, self%atol, &
                      self%ocean_fcns, self%ocean_args_p, &
                      P_e, z_e, T_e, f_i_e, self%P_trop, &
                      self%N_surface, self%N_ocean, &
                      err)
    if (allocated(err)) return

    self%f_i_surf = f_i_e(1,:)
    self%T_surf = T_surf
    self%P_surf = P_e(1)
    do i = 1,self%nz
      self%P(i) = P_e(2*i)
      self%T(i) = T_e(2*i)
      do j =1,self%sp%ng
        self%f_i(i,j) = f_i_e(2*i,j)
      enddo
    enddo

    call self%compute_altitude(err)
    if (allocated(err)) return
    
    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    call self%interpolate_particles(self%P, err)
    if (allocated(err)) return

    do i = 1,self%sp%ng
      ! mol/cm^2 in atmosphere
      self%N_atmos(i) = sum(density*self%f_i(:,i)*self%dz)/N_avo
    enddo

    do i = 1,self%nz
      if (self%P(i) > self%P_trop) then
        self%convecting_with_below(i) = .true.
      else
        self%convecting_with_below(i) = .false.
      endif
    enddo

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo
    
  end subroutine

  !> Similar to `make_profile`, but instead the input is column reservoirs 
  !> of each gas (mol/cm^2). 
  subroutine AdiabatClimate_make_column(self, T_surf, N_i_surf, err)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    use clima_eqns, only: gravity
    use ieee_arithmetic, only: ieee_positive_inf, ieee_value
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: N_i_surf(:) !! mole/cm^2
    character(:), allocatable, intent(out) :: err
    
    type(MinpackHybrd1Vars) :: mv
    real(dp), parameter :: scale_factors(*) = [1.0_dp, 0.5_dp, 2.0_dp, 0.1_dp, 5.0_dp, 0.01_dp]
    real(dp), allocatable :: P_i(:), N_i(:)
    real(dp) :: grav
    integer :: i, ii

    if (size(N_i_surf) /= self%sp%ng) then
      err = "N_i_surf has the wrong dimension"
      return
    endif

    ! Allocate memory for minpack
    mv = MinpackHybrd1Vars(n=self%sp%ng, tol=self%tol_make_column)
    mv%info = 0

    ! Allocate some work memory
    allocate(P_i(self%sp%ng), N_i(self%sp%ng))

    ! Gravity at the surface
    grav = gravity(self%planet_radius, self%planet_mass, 0.0_dp)

    ! First, if a initial guess is provided, then we try it
    if (self%use_make_column_P_guess) then
      do i = 1,mv%n
        mv%x(i) = log10(max(self%make_column_P_guess(i), sqrt(tiny(1.0_dp))))
      enddo
      ! Attempt the root solve
      call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    endif

    ! If the initial guess fails, or if an initial guess is not provided
    ! then we try to solve anyway using a variety of initial guesses.
    if (mv%info /= 1) then
      do ii = 1,size(scale_factors)
        ! Initial guess will be crude conversion of moles/cm2 (column) to 
        ! to dynes/cm2 (pressure). I try a few different `scale_factors` 
        ! to this guess.
        do i = 1,mv%n
          mv%x(i) = log10(max(N_i_surf(i)*self%sp%g(i)%mass*grav*scale_factors(ii), sqrt(tiny(1.0_dp))))
        enddo
        ! Attempt the root solve
        call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
        if (mv%info == 1) exit ! Success, so we exit.
      enddo
    endif

    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in make_column.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in make_column: '//err
      return
    elseif (mv%info == 1) then
      ! Converged, ensure there no allocated error
      if (allocated(err)) deallocate(err)
    endif

    ! Call one more time with solution
    call fcn(mv%n, mv%x, mv%fvec, mv%info)
    if (allocated(err)) return
    self%make_column_P_guess = 10.0_dp**mv%x

  contains

    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_

      P_i(:) = 10.0_dp**x_(:)
      if (any(P_i == ieee_value(1.0_dp, ieee_positive_inf))) then
        err = 'infinity values were encountered.'
        iflag_ = -1
        return
      endif

      call self%make_profile(T_surf, P_i, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif

      ! mol/cm^2 in atmosphere, on surface, and in ocean.
      N_i(:) = self%N_atmos(:)
      N_i(:) = N_i(:) + self%N_surface(:)
      do i = 1,self%sp%ng
        N_i(:) = N_i(:) + self%N_ocean(:,i)
      enddo

      ! Residual
      fvec_(:) = N_i - N_i_surf

    end subroutine

  end subroutine

  !> Similar to `make_profile`, but instead imposes a background gas and fixed
  !> surface pressure. We do a non-linear solve for the background gas pressure
  !> so that the desired total surface pressure is satisfied.
  subroutine AdiabatClimate_make_profile_bg_gas(self, T_surf, P_i_surf, P_surf, bg_gas, err)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: P_surf !! dynes/cm^2
    character(*), intent(in) :: bg_gas !! background gas
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: P_i_surf_copy(:)
    integer :: i, ind
    type(MinpackHybrd1Vars) :: mv
    real(dp), parameter :: scale_factors(*) = [1.0_dp, 0.1_dp]

    if (P_surf <= 0.0_dp) then
      err = 'P_surf must be greater than zero.'
      return
    endif

    ind = findloc(self%species_names, bg_gas, 1)
    if (ind == 0) then
      err = 'Gas "'//bg_gas//'" is not in the list of species'
      return
    endif

    allocate(P_i_surf_copy(size(P_i_surf)))
    P_i_surf_copy = P_i_surf

    mv = MinpackHybrd1Vars(1)
  
    do i = 1,size(scale_factors)
      mv%x(1) = log10(P_surf*scale_factors(i))
      call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
      if (mv%info == 1) then
        exit
      endif
    enddo
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in make_profile_bg_gas.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in make_profile_bg_gas: '//err
      return
    endif

    call fcn(mv%n, mv%x, mv%fvec, mv%info)

  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_

      P_i_surf_copy(ind) = 10.0_dp**x_(1)

      call self%make_profile(T_surf, P_i_surf_copy, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif

      fvec_(1) = self%P_surf - P_surf
    end subroutine
  end subroutine

  !> Given a P, T and mixing ratios, this function will update all atmospheric variables
  !> (except self%P_trop and self%convecting_with_below) to reflect these inputs. 
  !> The atmosphere is assumed to be dry (no condensibles). Any gas exceeding saturation 
  !> will not be altered.
  subroutine AdiabatClimate_make_profile_dry(self, P, T, f_i, err)
    use clima_adiabat_dry, only: make_profile_dry
    use clima_const, only: k_boltz, N_avo
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P(:) !! Pressure (dynes/cm^2). The first element is the surface.
    real(dp), intent(in) :: T(:) !! Temperature (K) defined on `P`.
    real(dp), intent(in) :: f_i(:,:) !! Mixing ratios defined on `P` of shape (size(P),ng).
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: P_e(:), z_e(:), T_e(:), f_i_e(:,:), lapse_rate_e(:)
    real(dp), allocatable :: density(:)
    integer :: i, j
    
    allocate(P_e(2*self%nz+1),z_e(2*self%nz+1),T_e(2*self%nz+1))
    allocate(f_i_e(2*self%nz+1,self%sp%ng),lapse_rate_e(2*self%nz+1))
    allocate(density(self%nz))

    call make_profile_dry(P, T, f_i, &
                          self%sp, self%nz, &
                          self%planet_mass, self%planet_radius, self%P_top, &
                          self%rtol, self%atol, &
                          P_e, z_e, T_e, f_i_e, lapse_rate_e, &
                          err)
    if (allocated(err)) return

    ! We assume the atmosphere is dry
    self%N_surface = 0.0_dp
    self%N_ocean = 0.0_dp

    self%f_i_surf = f_i_e(1,:)
    self%T_surf = T_e(1)
    self%P_surf = P_e(1)
    do i = 1,self%nz
      self%P(i) = P_e(2*i)
      self%T(i) = T_e(2*i)
      do j =1,self%sp%ng
        self%f_i(i,j) = f_i_e(2*i,j)
      enddo
    enddo

    call self%compute_altitude(err)
    if (allocated(err)) return

    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    call self%interpolate_particles(self%P, err)
    if (allocated(err)) return

    do i = 1,self%sp%ng
      ! mol/cm^2 in atmosphere
      self%N_atmos(i) = sum(density*self%f_i(:,i)*self%dz)/N_avo
    enddo

    ! lapse rates
    self%lapse_rate_intended(1) = lapse_rate_e(1)
    do i = 2,self%nz
      self%lapse_rate_intended(i) = lapse_rate_e(2*i-2)
    enddo

    ! Convecting with below is not changed. We don't get convecting information.

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo

  end subroutine

  !> Copies the atmosphere to the data used for radiative transfer
  subroutine AdiabatClimate_copy_atm_to_radiative_grid(self)
    class(AdiabatClimate), intent(inout) :: self
    integer :: i

    if (self%double_radiative_grid) then
      ! In doubled-grid mode, each physical layer i maps to RT layers 2*i-1 and 2*i.
      ! We then append two ghost RT layers that copy top-layer properties.
      do i = 1,self%nz
        self%T_r(2*(i-1)+1) = self%T(i)
        self%T_r(2*(i-1)+2) = self%T(i)

        self%P_r(2*(i-1)+1) = self%P(i)
        self%P_r(2*(i-1)+2) = self%P(i)

        self%densities_r(2*(i-1)+1,:) = self%densities(i,:)
        self%densities_r(2*(i-1)+2,:) = self%densities(i,:)

        self%pdensities_r(2*(i-1)+1,:) = self%pdensities(i,:)
        self%pdensities_r(2*(i-1)+2,:) = self%pdensities(i,:)

        self%pradii_r(2*(i-1)+1,:) = self%pradii(i,:)
        self%pradii_r(2*(i-1)+2,:) = self%pradii(i,:)

        self%dz_r(2*(i-1)+1) = 0.5_dp*self%dz(i)
        self%dz_r(2*(i-1)+2) = 0.5_dp*self%dz(i)
      enddo
      ! Fill in ghost layers
      do i = 1,2
        self%T_r(2*self%nz+i) = self%T_r(2*self%nz)
        self%P_r(2*self%nz+i) = self%P_r(2*self%nz)
        self%densities_r(2*self%nz+i,:) = self%densities_r(2*self%nz,:)
        self%pdensities_r(2*self%nz+i,:) = self%pdensities_r(2*self%nz,:)
        self%pradii_r(2*self%nz+i,:) = self%pradii_r(2*self%nz,:)
        self%dz_r(2*self%nz+i) = self%dz_r(2*self%nz)
      enddo
    else
      self%T_r = self%T
      self%P_r = self%P
      self%densities_r = self%densities
      self%pdensities_r = self%pdensities
      self%pradii_r = self%pradii
      self%dz_r = self%dz
    endif

  end subroutine
  
  !> Calls `make_profile`, then does radiative transfer on the constructed atmosphere
  subroutine AdiabatClimate_TOA_fluxes(self, T_surf, P_i_surf, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(out) :: ISR !! Top-of-atmosphere incoming solar radiation (mW/m^2)
    real(dp), intent(out) :: OLR !! Top-of-atmosphere outgoing longwave radiation (mW/m^2)
    character(:), allocatable, intent(out) :: err
    
    ! make atmosphere profile
    call self%make_profile(T_surf, P_i_surf, err)
    if (allocated(err)) return

    ! set albedo
    if (associated(self%albedo_fcn)) then
      self%rad%surface_albedo = self%albedo_fcn(T_surf)
    endif
    
    ! Do radiative transfer
    call self%copy_atm_to_radiative_grid()
    call self%rad%TOA_fluxes(T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                             self%pdensities_r, self%pradii_r, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return
    
  end subroutine

  !> Calls `make_column`, then does radiative transfer on the constructed atmosphere
  subroutine AdiabatClimate_TOA_fluxes_column(self, T_surf, N_i_surf, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: N_i_surf(:) !! mole/cm^2
    real(dp), intent(out) :: ISR !! Top-of-atmosphere incoming solar radiation (mW/m^2)
    real(dp), intent(out) :: OLR !! Top-of-atmosphere outgoing longwave radiation (mW/m^2)
    character(:), allocatable, intent(out) :: err    
    
    ! make atmosphere profile
    call self%make_column(T_surf, N_i_surf, err)
    if (allocated(err)) return

    ! set albedo
    if (associated(self%albedo_fcn)) then
      self%rad%surface_albedo = self%albedo_fcn(T_surf)
    endif
    
    ! Do radiative transfer
    call self%copy_atm_to_radiative_grid()
    call self%rad%TOA_fluxes(T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                             self%pdensities_r, self%pradii_r, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return

  end subroutine

  !> Calls `make_profile_bg_gas`, then does radiative transfer on the constructed atmosphere
  subroutine AdiabatClimate_TOA_fluxes_bg_gas(self, T_surf, P_i_surf, P_surf, bg_gas, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: P_surf !! dynes/cm^2
    character(*), intent(in) :: bg_gas !! background gas
    real(dp), intent(out) :: ISR !! Top-of-atmosphere incoming solar radiation (mW/m^2)
    real(dp), intent(out) :: OLR !! Top-of-atmosphere outgoing longwave radiation (mW/m^2)
    character(:), allocatable, intent(out) :: err
    
    ! make atmosphere profile
    call self%make_profile_bg_gas(T_surf, P_i_surf, P_surf, bg_gas, err)
    if (allocated(err)) return

    ! set albedo
    if (associated(self%albedo_fcn)) then
      self%rad%surface_albedo = self%albedo_fcn(T_surf)
    endif
    
    ! Do radiative transfer
    call self%copy_atm_to_radiative_grid()
    call self%rad%TOA_fluxes(T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                             self%pdensities_r, self%pradii_r, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return

  end subroutine

  !> Calls `make_profile_dry`, then does radiative transfer on the constructed atmosphere
  subroutine AdiabatClimate_TOA_fluxes_dry(self, P, T, f_i, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P(:) !! Pressure (dynes/cm^2). The first element is the surface.
    real(dp), intent(in) :: T(:) !! Temperature (K) defined on `P`.
    real(dp), intent(in) :: f_i(:,:) !! Mixing ratios defined on `P` of shape (size(P),ng).
    real(dp), intent(out) :: ISR !! Top-of-atmosphere incoming solar radiation (mW/m^2)
    real(dp), intent(out) :: OLR !! Top-of-atmosphere outgoing longwave radiation (mW/m^2)
    character(:), allocatable, intent(out) :: err
    
    ! make atmosphere profile
    call self%make_profile_dry(P, T, f_i, err)
    if (allocated(err)) return

    ! set albedo
    if (associated(self%albedo_fcn)) then
      self%rad%surface_albedo = self%albedo_fcn(self%T_surf)
    endif
    
    ! Do radiative transfer
    call self%copy_atm_to_radiative_grid()
    call self%rad%TOA_fluxes(self%T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                             self%pdensities_r, self%pradii_r, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return
    
  end subroutine

  function AdiabatClimate_simple_solver(self, fcn, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(AdiabatClimate), intent(inout) :: self
    procedure(toa_fluxes_fcn) :: fcn
    real(dp), optional, intent(in) :: T_guess !! K
    character(:), allocatable, intent(out) :: err
    real(dp) :: T_surf

    real(dp) :: T_guess_
    type(MinpackHybrd1Vars) :: mv
    
    if (present(T_guess)) then
      T_guess_ = T_guess
    else
      T_guess_ = 280.0_dp
    endif
    
    if (self%solve_for_T_trop) then
      mv = MinpackHybrd1Vars(2)
      mv%x(1) = log10(T_guess_)
      mv%x(2) = log10(self%T_trop)
    else
      mv = MinpackHybrd1Vars(1)
      mv%x(1) = log10(T_guess_)
    endif
    call hybrd1(fcn_hybrj, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed: '//err
      return
    endif
    
    T_surf = 10.0_dp**mv%x(1)
    call fcn_hybrj(mv%n, mv%x, mv%fvec, mv%info)

  contains
    subroutine fcn_hybrj(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_
      real(dp) :: T, T_trop, ISR, OLR, rad_enhancement
      T = 10.0_dp**x_(1)
      if (self%solve_for_T_trop) then
        T_trop = 10.0_dp**x_(2)
      else
        T_trop = self%T_trop
      endif
      call fcn(T, T_trop, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      rad_enhancement = 1.0_dp
      if (self%tidally_locked_dayside) then; block
        real(dp) :: tau_LW, k_term, f_term
        call self%heat_redistribution_parameters(tau_LW, k_term, f_term, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
        ! Increase the stellar flux, because we are computing the climate of
        ! observed dayside.
        rad_enhancement = 4.0_dp*f_term
        call self%rad%apply_radiation_enhancement(rad_enhancement)
      endblock; endif
      fvec_(1) = ISR*rad_enhancement - OLR + self%surface_heat_flow

      if (self%solve_for_T_trop) then; block
        use clima_eqns, only: skin_temperature
        real(dp) :: bond_albedo, stellar_radiation
        bond_albedo = self%rad%wrk_sol%fup_n(self%nz_r+1)/self%rad%wrk_sol%fdn_n(self%nz_r+1)
        stellar_radiation = self%rad%bolometric_flux()
        fvec_(2) = skin_temperature(stellar_radiation*rad_enhancement, bond_albedo) - T_trop
      endblock; endif
    end subroutine    
  end function
  
  !> Does a non-linear solve for the surface temperature that balances incoming solar
  !> and outgoing longwave radiation. Uses `make_profile`.
  function AdiabatClimate_surface_temperature(self, P_i_surf, T_guess, err) result(T_surf)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), optional, intent(in) :: T_guess !! K
    character(:), allocatable, intent(out) :: err
    real(dp) :: T_surf
    T_surf = AdiabatClimate_simple_solver(self, fcn, T_guess, err)
  contains
    subroutine fcn(T_surf_, T_trop_, ASR_, OLR_, err_)
      real(dp), intent(in) :: T_surf_, T_trop_
      real(dp), intent(out) :: ASR_, OLR_
      character(:), allocatable, intent(out) :: err_
      self%T_trop = T_trop_
      call self%TOA_fluxes(T_surf_, P_i_surf, ASR_, OLR_, err_)
    end subroutine
  end function

  function AdiabatClimate_surface_temperature_column(self, N_i_surf, T_guess, err) result(T_surf)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: N_i_surf(:)
    real(dp), optional, intent(in) :: T_guess
    character(:), allocatable, intent(out) :: err
    real(dp) :: T_surf
    T_surf = AdiabatClimate_simple_solver(self, fcn, T_guess, err)
  contains
    subroutine fcn(T_surf_, T_trop_, ASR_, OLR_, err_)
      real(dp), intent(in) :: T_surf_, T_trop_
      real(dp), intent(out) :: ASR_, OLR_
      character(:), allocatable, intent(out) :: err_
      self%T_trop = T_trop_
      call self%TOA_fluxes_column(T_surf_, N_i_surf, ASR_, OLR_, err_)
    end subroutine
  end function

  !> Similar to surface_temperature. The difference is that this function imposes
  !> a background gas and fixed surface pressure.
  function AdiabatClimate_surface_temperature_bg_gas(self, P_i_surf, P_surf, bg_gas, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: P_surf !! dynes/cm^2
    character(*), intent(in) :: bg_gas !! background gas
    real(dp), optional, intent(in) :: T_guess
    character(:), allocatable, intent(out) :: err
    real(dp) :: T_surf
    T_surf = AdiabatClimate_simple_solver(self, fcn, T_guess, err)
  contains
    subroutine fcn(T_surf_, T_trop_, ASR_, OLR_, err_)
      real(dp), intent(in) :: T_surf_, T_trop_
      real(dp), intent(out) :: ASR_, OLR_
      character(:), allocatable, intent(out) :: err_
      self%T_trop = T_trop_
      call self%TOA_fluxes_bg_gas(T_surf_, P_i_surf, P_surf, bg_gas, ASR_, OLR_, err_)
    end subroutine
  end function

  subroutine AdiabatClimate_interpolate_particles(self, P, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P(:)
    character(:), allocatable, intent(out) :: err

    integer :: i, j

    ! Check inputs
    if (size(P) /= self%nz) then
      err = '`P` has the wrong shape'
      return
    endif

    do i = 1,self%sp%np
      do j = 1,self%nz
        call self%pdensities_interp(i)%evaluate(log10(P(j)), self%pdensities(j,i))
        self%pdensities(j,i) = 10.0_dp**self%pdensities(j,i)
        call self%pradii_interp(i)%evaluate(log10(P(j)), self%pradii(j,i))
        self%pradii(j,i) = 10.0_dp**self%pradii(j,i)
      enddo
    enddo

  end subroutine

  !> Sets particle densities and radii
  subroutine AdiabatClimate_set_particle_density_and_radii(self, P, pdensities, pradii, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P(:) !! Array of pressures in dynes/cm^2
    !> Particle densities in particles/cm^3 at each pressure and for each particle
    !> in the model. Shape (nz, np).
    real(dp), intent(in) :: pdensities(:,:) 
    !> Particle radii in cm at each pressure and for each particle
    !> in the model. Shape (nz, np).
    real(dp), intent(in) :: pradii(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: log10P(:), tmp(:)
    integer :: i, istat

    ! Check inputs
    if (size(P) < 1) then
      err = '`P` must have a length greater than zero'
      return
    endif
    if (size(P) /= size(pdensities,1)) then
      err = '`P` and `pdensities` have incompatible shapes'
      return
    endif
    if (size(P) /= size(pradii,1)) then
      err = '`P` and `pradii` have incompatible shapes'
      return
    endif
    if (size(pdensities,2) /= self%sp%np) then
      err = 'The second dimension of `pdensities` does not match the number of particles'
      return
    endif
    if (size(pradii,2) /= self%sp%np) then
      err = 'The second dimension of `pradii` does not match the number of particles'
      return
    endif
    if (any(P <= 0.0_dp)) then
      err = 'All elements of `P` must be larger than zero'
      return
    endif
    if (any(pdensities < 0.0_dp)) then
      err = 'All elements of `pdensities` must be larger than zero'
      return
    endif
    if (any(pradii < 0.0_dp)) then
      err = 'All elements of `pradii` must be larger than zero'
      return
    endif

    log10P = [[huge(1.0_dp),P], tiny(1.0_dp)]
    log10P = log10(log10P(size(log10P):1:-1))
    allocate(tmp(size(log10P)))

    do i = 1,self%sp%np
      tmp(2:size(tmp)-1) = pdensities(:,i)
      tmp(1) = pdensities(1,i)
      tmp(size(tmp)) = pdensities(size(pdensities,1),i)
      tmp = log10(max(tmp(size(tmp):1:-1),tiny(1.0_dp)))

      call self%pdensities_interp(i)%initialize(log10P, tmp, istat)
      if (istat /= 0) then
        err = 'Interpolation initialization in `set_particle_density_and_radii` failed'
        return
      endif

      tmp(2:size(tmp)-1) = pradii(:,i)
      tmp(1) = pradii(1,i)
      tmp(size(tmp)) = pradii(size(pradii,1),i)
      tmp = log10(max(tmp(size(tmp):1:-1),tiny(1.0_dp)))

      call self%pradii_interp(i)%initialize(log10P, tmp, istat)
      if (istat /= 0) then
        err = 'Interpolation initialization in `set_particle_density_and_radii` failed'
        return
      endif
    enddo

  end subroutine

  !> Sets a function for describing how gases dissolve in a liquid ocean.
  subroutine AdiabatClimate_set_ocean_solubility_fcn(self, sp, fcn, err)
    class(AdiabatClimate), intent(inout) :: self
    character(*), intent(in) :: sp !! name of species that makes the ocean
    !> Function describing solubility of other gases in ocean
    procedure(ocean_solubility_fcn), pointer, intent(in) :: fcn 
    character(:), allocatable, intent(out) :: err

    integer :: ind

    ind = findloc(self%species_names, sp, 1)
    if (ind == 0) then
      err = 'Gas "'//sp//'" is not in the list of species'
      return
    endif

    self%ocean_fcns(ind)%fcn => fcn

  end subroutine

  !> Re-grids atmosphere so that each grid cell is equally spaced in altitude.
  subroutine AdiabatClimate_to_regular_grid(self, err)
    use futils, only: rebin, interp
    use clima_eqns, only: vertical_grid
    use clima_const, only: k_boltz
    class(AdiabatClimate), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: ze(:), ze_new(:)
    real(dp), allocatable :: z_new(:), dz_new(:)
    real(dp), allocatable :: densities_new(:,:)
    real(dp), allocatable :: f_i_new(:,:)
    real(dp), allocatable :: T_new(:)
    real(dp), allocatable :: density_new(:)
    real(dp), allocatable :: P_new(:)

    integer :: i, j, ierr

    allocate(z_new(self%nz), dz_new(self%nz))

    ! compute the new grid
    call vertical_grid(0.0_dp, self%z(self%nz)+0.5_dp*self%dz(self%nz), &
                       self%nz, z_new, dz_new)

    ! rebin of the densities
    allocate(ze(self%nz+1), ze_new(self%nz+1))
    ze_new(1) = z_new(1) - 0.5_dp*dz_new(1)
    do i = 1,self%nz
      ze_new(i+1) = z_new(i) + 0.5_dp*dz_new(i)
    enddo
    ze(1) = self%z(1) - 0.5_dp*self%dz(1)
    do i = 1,self%nz
      ze(i+1) = self%z(i) + 0.5_dp*self%dz(i)
    enddo

    allocate(densities_new(self%nz, self%sp%ng))
    do i = 1,self%sp%ng
      call rebin(ze, self%densities(:,i), ze_new, densities_new(:,i), ierr)
      if (ierr /= 0) then
        err = 'subroutine `rebin` returned an error'
        return
      endif
    enddo

    allocate(T_new(self%nz))
    call interp(self%nz, self%nz, z_new, self%z, self%T, T_new, ierr)
    if (ierr /= 0) then
      err = 'Subroutine `interp` returned an error.'
      return
    endif

    allocate(density_new(self%nz))
    allocate(f_i_new(self%nz,self%sp%ng))
    do i = 1,self%nz
      density_new(i) = sum(densities_new(i,:))
      do j = 1,self%sp%ng
        f_i_new(i,j) = densities_new(i,j)/density_new(i)
      enddo
    enddo
    allocate(P_new(self%nz))
    P_new = density_new(:)*k_boltz*T_new(:)

    self%P(:) = P_new(:)
    self%T(:) = T_new(:)
    self%f_i(:,:) = f_i_new(:,:)
    self%z(:) = z_new(:)
    self%dz(:) = dz_new(:)
    self%densities(:,:) = densities_new(:,:)

  end subroutine

  subroutine AdiabatClimate_out2atmosphere_txt(self, filename, eddy, number_of_decimals, overwrite, clip, err)
    use futils, only: FileCloser
    class(AdiabatClimate), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: eddy(:)
    integer, intent(in) :: number_of_decimals
    logical, intent(in) :: overwrite, clip
    character(:), allocatable, intent(out) :: err
    
    character(len=100) :: tmp
    integer :: io, i, j
    integer :: number_of_spaces
    character(len=10) :: number_of_decimals_str, number_of_spaces_str
    character(:), allocatable :: fmt_label, fmt_number
    real(dp) :: clip_value
    type(FileCloser) :: file

    call self%to_regular_grid(err)
    if (allocated(err)) return

    if (size(eddy,1) /= self%nz) then
      err = '"eddy" has the wrong size'
      return
    endif

    if (clip) then
      clip_value = 1.0e-40_dp
    else
      clip_value = - huge(1.0_dp)
    endif
    
    if (overwrite) then
      open(2, file=filename, form='formatted', status='replace', iostat=io)
      file%unit = 2
      if (io /= 0) then
        err = "Unable to overwrite file "//trim(filename)
        return
      endif
    else
      open(2, file=filename, form='formatted', status='new', iostat=io)
      file%unit = 2
      if (io /= 0) then
        err = "Unable to create file "//trim(filename)//" because it already exists"
        return
      endif
    endif

    ! number of decimals must be reasonable
    if (number_of_decimals < 2 .or. number_of_decimals > 17) then
      err = '"number_of_decimals" should be between 1 and 17.'
      return
    endif
    number_of_spaces = number_of_decimals + 9
    ! make sure number of spaces works with the length of species names
    do i = 1,self%sp%ng
      number_of_spaces = max(number_of_spaces,len_trim(self%species_names(i)) + 3)
    enddo
    write(number_of_decimals_str,'(i10)') number_of_decimals
    write(number_of_spaces_str,'(i10)') number_of_spaces

    fmt_label = "(a"//trim(adjustl(number_of_spaces_str))//")"
    fmt_number = "(es"//trim(adjustl(number_of_spaces_str))//"."//trim(adjustl(number_of_decimals_str))//"e3)"
    
    tmp = 'alt'
    write(unit=2,fmt=fmt_label,advance='no') tmp
    tmp = 'press'
    write(unit=2,fmt=fmt_label,advance='no') tmp
    tmp = 'den'
    write(unit=2,fmt=fmt_label,advance='no') tmp
    tmp = 'temp'
    write(unit=2,fmt=fmt_label,advance='no') tmp
    tmp = 'eddy'
    write(unit=2,fmt=fmt_label,advance='no') tmp
    do j = 1,self%sp%ng
      tmp = self%species_names(j)
      write(unit=2,fmt=fmt_label,advance='no') tmp
    enddo
    
    do i = 1,self%nz
      write(2,*)
      write(tmp,fmt=fmt_number) self%z(i)/1.e5_dp
      write(unit=2,fmt=fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt=fmt_number) self%P(i)/1.e6_dp
      write(unit=2,fmt=fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt=fmt_number) sum(self%densities(i,:))
      write(unit=2,fmt=fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt=fmt_number) self%T(i)
      write(unit=2,fmt=fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt=fmt_number) eddy(i)
      write(unit=2,fmt=fmt_label,advance='no') adjustl(tmp)

      do j = 1,self%sp%ng
        write(tmp,fmt=fmt_number) max(self%f_i(i,j), clip_value)
        write(unit=2,fmt=fmt_label,advance='no') adjustl(tmp)
      enddo
    enddo
    
  end subroutine

  !> For considering a tidally locked planet. This function computes key parameters for
  !> Equation (10) in Koll (2022), ApJ. The function must be called after calling a function
  !> `self%TOA_fluxes`, because it uses atmosphere properties, and radiative properties.
  subroutine AdiabatClimate_heat_redistribution_parameters(self, tau_LW, k_term, f_term, err)
    use clima_const, only: c_light
    use clima_eqns, only: k_term_heat_redistribution, f_heat_redistribution, &
                          gravity, heat_capacity_eval, planck_fcn
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(out) :: tau_LW !! Long wavelength optical depth
    real(dp), intent(out) :: k_term !! k term in Equation (10)
    real(dp), intent(out) :: f_term !! The heat redistribution parameter, f in Equation (10)
    character(:), allocatable, intent(out) :: err

    integer :: i

    real(dp) :: bond_albedo, Teq
    real(dp) :: grav
    real(dp) :: mubar
    logical :: found
    real(dp) ::  cp_tmp, cp

    ! equilibrium temperature
    bond_albedo = self%rad%wrk_sol%fup_n(self%nz_r+1)/self%rad%wrk_sol%fdn_n(self%nz_r+1)
    Teq = self%rad%equilibrium_temperature(bond_albedo)
    
    ! gravity
    grav = gravity(self%planet_radius, self%planet_mass, 0.0_dp)

    ! mean molecular weight in first layer
    mubar = 0.0_dp
    do i = 1,self%sp%ng
      mubar = mubar + self%f_i(1,i)*self%sp%g(i)%mass
    enddo

    ! heat capacity at surface
    cp = tiny(0.0_dp)
    do i = 1,self%sp%ng
      call heat_capacity_eval(self%sp%g(i)%thermo, self%T_surf, found, cp_tmp) ! J/(mol*K)
      if (.not. found) then
        err = "Failed to compute heat capacity"
        return
      endif
      ! J/(mol*K)
      cp = cp + cp_tmp*self%f_i(1,i) ! J/(mol*K)
    enddo
    ! J/(mol*K) * (mol/kg) = J/(kg*K)
    cp = cp*(1.0_dp/(mubar*1.0e-3_dp))
    ! J/(kg*K) * (erg/J) * (kg/g) = erg/(g*K)
    cp = cp*1.0e4_dp

    ! tau_LW. Computed with Equation (13) in Koll (2020), ApJ
    block
      real(dp) :: numerator, denominator, dlambda, tau_lambda, bplank, avg_freq, avg_lam
      numerator = 0.0_dp
      denominator = 0.0_dp
      do i = 1,self%rad%ir%nw
        dlambda = self%rad%ir%wavl(i+1) - self%rad%ir%wavl(i) ! Width of a wavelength bin (nm)
        tau_lambda = sum(self%rad%wrk_ir%tau_band(:,i)) ! IR optical depth

        avg_freq = 0.5_dp*(self%rad%ir%freq(i) + self%rad%ir%freq(i+1))
        avg_lam = (c_light*1.0e9_dp/avg_freq) ! (m/s) * (nm/m) * (s) = nm

        bplank = planck_fcn(avg_freq, self%T_surf) ! (mW sr^−1 m^−2 Hz^-1)
        bplank = bplank * (avg_freq/avg_lam) ! (mW sr^−1 m^−2 Hz^-1) * (Hz/nm) = (mW sr^−1 m^−2 nm^-1)

        ! integrate the numerator
        numerator = numerator + exp(-tau_lambda)*bplank*dlambda
        ! integrate the denominator
        denominator = denominator + bplank*dlambda
      enddo
      tau_LW = - log(numerator/denominator)
    endblock

    k_term = k_term_heat_redistribution(self%L, grav, self%chi, mubar, cp, self%n_LW, self%Cd)
    f_term = f_heat_redistribution(tau_LW, self%P_surf, Teq, k_term)

  end subroutine

end module
