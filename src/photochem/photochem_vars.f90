module photochem_vars
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use photochem_const, only: dp
  use photochem_data, only: PhotochemData, ParticleXsections
  use photochem_enum, only: AnalyticalJacobian, ContinuousPressTempEdd, TAU_MODE_DZ, TAU_MODE_DLOG10P
  use photochem_settings, only: PhotoSettings, CondensationParameters
  implicit none
  private

  public :: CustomOpticalProperties
  public :: PhotochemVars
  public :: PressureTempEddProfile, TOAPressureMaintenance
  public :: time_dependent_flux_fcn, time_dependent_rate_fcn
  public :: binary_diffusion_fcn
  public :: refresh_temperature_dependent_vars, interp2particlexsdata
  public :: tropopause_layer_index

  !> Pressure-based external aerosol optics, independent of the atmospheric grid.
  type :: CustomOpticalProperties
    logical :: enabled = .false.
    integer :: tau_mode = TAU_MODE_DZ
    real(dp), allocatable :: wavelength(:) !! nm, increasing
    real(dp), allocatable :: pressure(:) !! dyn/cm^2, decreasing
    real(dp), allocatable :: dtau(:,:), w0(:,:), g0(:,:) !! (pressure,wavelength)
    ! Derived lazily for a model wavelength grid; never serialized.
    real(dp), private, allocatable :: wavelength_edges_cache(:), log10_pressure_cache(:)
    real(dp), private, allocatable :: dtau_cache(:,:), w0_cache(:,:), g0_cache(:,:)
  contains
    procedure :: set_custom_optical_properties => CustomOpticalProperties_set_custom_optical_properties
    procedure :: unset_custom_optical_properties => CustomOpticalProperties_unset_custom_optical_properties
    procedure :: validate => CustomOpticalProperties_validate
    procedure :: evaluate => CustomOpticalProperties_evaluate
  end type

  !> Persistent temperature and eddy-diffusion profiles defined on pressure.
  type :: PressureTempEddProfile
    logical :: enabled = .false.
    !> Continuous mapping is the backward-compatible default. Periodic mapping
    !! is performed by the robust stepper at resynchronization points.
    integer :: mode = ContinuousPressTempEdd
    logical :: hydro_pressure = .true.
    real(dp) :: trop_p = -1.0_dp !! Non-positive values disable the tropopause pressure.
    real(dp), allocatable :: pressure(:)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: edd(:)
    ! Relative mismatch tolerances and multiplier for classifying extreme drift.
    real(dp) :: temperature_tol = 0.005_dp
    real(dp) :: edd_tol = 0.01_dp
    real(dp) :: extreme_factor = 10.0_dp
  contains
    procedure :: validate => PressureTempEddProfile_validate
  end type

  !> Settings for optional robust-stepper maintenance of model-top pressure.
  !! This requires a persistent pressure-based temperature and eddy-diffusion profile.
  type :: TOAPressureMaintenance
    logical :: enabled = .false. !! Enable automatic TOA-pressure maintenance.
    real(dp) :: target_pressure = 0.0_dp !! Target pressure (dyn/cm^2).
    real(dp) :: pressure_factor = 3.0_dp !! Multiplicative acceptable pressure factor
    real(dp) :: extreme_pressure_factor = 10.0_dp
  contains
    procedure :: validate => TOAPressureMaintenance_validate
  end type

  abstract interface
    !> Compute a custom binary diffusion parameter.
    function binary_diffusion_fcn(mu_i, mubar, T) result(b)
      use iso_c_binding, only: c_double
      real(c_double), value, intent(in) :: mu_i !! molar weight of species i (g/mol)
      real(c_double), value, intent(in) :: mubar !! molar weight of background gas (g/mol)
      real(c_double), value, intent(in) :: T !! Temperature (K)
      real(c_double) :: b !! Binary diffusion parameter relative to the background gas
                          !! (molecules cm^-1 s^-1)
    end function

    !> Set the incident photon flux at a specified model time.
    subroutine time_dependent_flux_fcn(tn, nw, photon_flux)
      use iso_c_binding, only: c_double, c_int
      real(c_double), value, intent(in) :: tn !! Time (s)
      integer(c_int), value, intent(in) :: nw !! Number of wavelength bins
      real(c_double), intent(out) :: photon_flux(nw) !! photons cm^-2 s^-1
    end subroutine

    !> Set a time-dependent production or destruction rate for a species.
    subroutine time_dependent_rate_fcn(tn, nz, rate)
      use iso_c_binding, only: c_double, c_int
      real(c_double), value, intent(in) :: tn !! time (s)
      integer(c_int), value, intent(in) :: nz !! number of atmospheric layers
      real(c_double), intent(out) :: rate(nz) !! molecules/cm^3/s (can be positive or negative)
    end subroutine
  end interface

  !> Container to make an array of functions for time-dependent production rates.
  type :: time_dependent_rate_fcns
    procedure(time_dependent_rate_fcn), nopass, pointer :: fcn => null()
  end type

  !> Configured model variables and persistent prepared atmospheric state.
  !! The constructor allocates atmosphere-dependent arrays without initializing
  !! their physical values; those become valid when the owning atmosphere is initialized.
  type :: PhotochemVars

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set DURING file read-in !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! boundary conditions
    !> (nq) Type of lower bc. See photochem_enum.f90 for different types.
    integer, allocatable :: lowerboundcond(:)
    real(dp), allocatable :: lower_vdep(:)
    real(dp), allocatable :: lower_flux(:)
    real(dp), allocatable :: lower_dist_height(:)
    real(dp), allocatable :: lower_fix_den(:)
    real(dp), allocatable :: lower_fix_press(:)
    integer, allocatable :: upperboundcond(:) ! 0 or 2
    real(dp), allocatable :: upper_veff(:)
    real(dp), allocatable :: upper_flux(:)
    !> Functions for settings time-dependent production rates (nq)
    type(time_dependent_rate_fcns), allocatable :: rate_fcns(:)

    ! Atmospheres structure
    real(dp) :: bottom_atmos = 0.0_dp !! Bottom altitude of the model domain (cm).
    real(dp) :: top_atmos = 0.0_dp !! Top altitude of the model domain (cm).
    integer :: nz !! Number of vertical layers.
    real(dp) :: surface_albedo
    !> Factor multiplying the incident photon flux during photolysis-rate calculations.
    !! The default 0.5 accounts for half the planet facing the star.
    real(dp) :: diurnal_fac = 0.5_dp
    real(dp) :: solar_zenith !! degrees
    real(dp) :: trop_alt = 0.0_dp !! Tropopause altitude (cm); used for gas rainout.
    real(dp) :: rainfall_rate !! Relative to Earth's average rainfall of 1.1e17 molecules/cm^2/s.
    integer :: trop_ind !! Internal Fortran tropopause index; used for gas rainout.

    ! radiative transfer
    real(dp), allocatable :: photon_flux(:) !! (nw) photon/cm^2/s in each wavelength bin hitting planet.
    !> for scaling photon flux for different planets in a solar system
    real(dp) :: photon_scale_factor
    !> Callback that supplies the photon flux as a function of model time.
    procedure(time_dependent_flux_fcn), nopass, pointer :: photon_flux_fcn => null()

    ! particles
    !> Condensation and evaporation parameters for each particle species (np).
    type(CondensationParameters), allocatable :: cond_params(:) ! (np)
    logical :: evaporation = .true. !! If true, then evaporation occurs.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! set AFTER file read-in !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Most fields in this section change only through public model-mutation
    ! routines. The prepared subset temperature, edd, trop_alt, trop_ind,
    ! xs_x_qy, and gibbs_energy may also be updated during integration when a
    ! persistent pressure-temperature-Kzz profile is enabled. These remain the
    ! single authoritative prepared values; successful integration routines
    ! leave them consistent with the committed wrk%usol before returning.
    integer :: neqs !! number of equations nq*nz
    real(dp), allocatable :: temperature(:) !! Layer temperatures (K), shape (nz).
    real(dp), allocatable :: z(:) !! Layer-center altitudes (cm), shape (nz).
    real(dp), allocatable :: dz(:) !! Layer thicknesses (cm), shape (nz).
    real(dp), allocatable :: edd(:) !! Eddy diffusion coefficients (cm^2/s), shape (nz).
    !> State for a persistent pressure-based temperature and Kzz profile.
    type(PressureTempEddProfile) :: press_temp_edd_profile
    !> Callback for specifying a custom binary diffusion parameter (b_ij).
    procedure(binary_diffusion_fcn), nopass, pointer :: custom_binary_diffusion_fcn => null()
    real(dp), allocatable :: grav(:) !! Layer-center gravity (cm/s^2), shape (nz).
    real(dp), allocatable :: particle_radius(:,:) !! Particle radii (cm), shape (np,nz).
    real(dp), allocatable :: xs_x_qy(:,:,:) !! (nz,kj,nw) photolysis cross sections times quantum yields (cm2/molecule)
    type(ParticleXsections), allocatable :: particle_xs(:) !! (np)
    real(dp), allocatable :: gibbs_energy(:,:) !! (nz,ng) Joules/mol

    ! Persistent custom optical prescription and its wavelength interpolation cache.
    type(CustomOpticalProperties) :: custom_optical_properties

    ! other
    !> Number of times to reinitialize CVODE when it returns
    !> a potentially recoverable error. ONLY USED IN EVOATMOSPHERE (NOT ATMOSPHERE)
    !> in the `evolve` method.
    integer :: max_error_reinit_attempts = 2
    real(c_double) :: rtol = 1.0e-3_dp !! integration relative tolerance
    !> Dimensionless absolute tolerance in mixing ratio. This is multiplied by
    !> each layer's hydrostatic number density for CVODE.
    real(c_double) :: atol = 1.0e-23_dp
    integer :: mxsteps = 100000 !! max number of steps before integrator will give up.
    !> Time in seconds at which an integration is considered converged.
    real(dp) :: equilibrium_time = 1.0e17_dp
    !> For convergence checking. Considers mixing ratio change between t_now and time
    !> t = t_now*conv_hist_factor to see if atmosphere is changing.
    real(dp) :: conv_hist_factor = 0.5_dp
    !> Minimum mixing ratio considered in convergence checking.
    real(dp) :: conv_min_mix = 1.0e-20_dp
    !> Threshold normalized change in mixing ratios for convergence.
    !> A reasonable value is ~1.0e-2.
    real(dp) :: conv_longdy = 0.0_dp
    !> Threshold normalized mixing-ratio change rate for convergence (s^-1).
    real(dp) :: conv_longdydt = 1.0e-6_dp
    real(c_double) :: initial_dt = 1.0e-6_dp !! Initial CVODE time step (seconds).
    !> Maximum time step size (seconds).
    real(c_double) :: max_dt = sqrt(huge(1.0_dp))
    integer(c_int) :: max_err_test_failures = 15 !! CVODE max error test failures
    integer(c_int) :: max_order = 5 !! CVODE max order for BDF method.
    !> Method used to compute the chemistry Jacobian. `1` selects autodiff,
    !> `2` finite differences, and `3` the default analytical implementation.
    !> See photochem_enum.f90.
    integer(c_int) :: jacobian_method = AnalyticalJacobian
    !> Relative perturbation used by the finite-difference Jacobian method.
    real(dp) :: epsj = 1.0e-4_dp
    integer :: verbose = 1 !! 0 == no printing. 1 == some printing. 2 == bunch of printing.
    !> If True, then the code uses a 1st order upwind method for the advective molecular
    !> diffusion terms instead of a centered scheme. This permits stability (at the cost
    !> of accuracy) for atmospheres with strong molecular advection in the upper atmosphere.
    logical :: upwind_molec_diff = .false.

    ! Settings for the `robust_stepper` and `find_steady_state` methods
    !> Number of failed-step recovery restarts allowed. The next integration
    !> error ends the robust session.
    integer :: nerrors_before_giveup = 10
    !> Recoverable grid or periodic profile resynchronization failures allowed
    !! before robust integration stops. The failure counter resets after a successful resynchronization.
    integer :: max_resync_failures = 0
    !> Limit for chemistry-converged states that still require a maintenance
    !! resynchronization of TOA pressure or a periodic P-T-Kzz profile.
    integer :: nconverged_but_restarted_limit = 7
    !> Number of accepted steps after initialization or restart to take before
    !> checking the non-time convergence criteria.
    integer :: nsteps_before_conv_check = 300
    !> Accepted steps per integration segment. At this count CVODE is restarted
    !> and the segment-local convergence history is discarded.
    integer :: nsteps_before_reinit = 1000
    !> Maximum total accepted steps before `robust_step` returns give_up.
    integer :: nsteps_before_giveup = 100000
    !> When the integrator reinitializes, this is the minimum
    !> density allowed (molecules/cm^3)
    real(dp) :: reinit_min_density = 1.0e-40_dp
    !> Optional robust-stepper maintenance of the model-top pressure. This mode
    !> requires press_temp_edd_profile%enabled to preserve the physical T-Kzz
    !> structure when the vertical grid changes.
    type(TOAPressureMaintenance) :: toa_pressure_maintenance
    ! End settings for `robust_stepper` and `find_steady_state`
  end type

  interface PhotochemVars
    module procedure :: create_PhotochemVars
  end interface

contains

  !> Validate and copy a pressure-based optical prescription, preserving self on error.
  subroutine CustomOpticalProperties_set_custom_optical_properties(self, wv, P, dtau, w0, g0, tau_mode, err)
    class(CustomOpticalProperties), intent(inout) :: self
    real(dp), intent(in) :: wv(:) !! Increasing wavelength (nm).
    real(dp), intent(in) :: P(:) !! Decreasing hydrostatic pressure (dyn/cm^2).
    real(dp), intent(in) :: dtau(:,:) !! Extinction per cm or pressure decade, shape (size(P),size(wv)).
    real(dp), intent(in) :: w0(:,:), g0(:,:) !! Dimensionless albedo and asymmetry, shape as dtau.
    integer, intent(in) :: tau_mode !! TAU_MODE_DZ or TAU_MODE_DLOG10P.
    character(:), allocatable, intent(out) :: err
    type(CustomOpticalProperties) :: candidate

    candidate%enabled = .true.
    candidate%tau_mode = tau_mode
    candidate%wavelength = wv
    candidate%pressure = P
    candidate%dtau = dtau
    candidate%w0 = w0
    candidate%g0 = g0
    call candidate%validate(err)
    if (allocated(err)) return

    call self%unset_custom_optical_properties()
    self%enabled = candidate%enabled
    self%tau_mode = candidate%tau_mode
    self%wavelength = candidate%wavelength
    self%pressure = candidate%pressure
    self%dtau = candidate%dtau
    self%w0 = candidate%w0
    self%g0 = candidate%g0
  end subroutine

  !> Disable custom optics and release both the prescription and derived cache.
  subroutine CustomOpticalProperties_unset_custom_optical_properties(self)
    class(CustomOpticalProperties), intent(inout) :: self
    if (allocated(self%wavelength_edges_cache)) deallocate(self%wavelength_edges_cache)
    if (allocated(self%log10_pressure_cache)) deallocate(self%log10_pressure_cache)
    if (allocated(self%dtau_cache)) deallocate(self%dtau_cache)
    if (allocated(self%w0_cache)) deallocate(self%w0_cache)
    if (allocated(self%g0_cache)) deallocate(self%g0_cache)
    if (allocated(self%wavelength)) deallocate(self%wavelength)
    if (allocated(self%pressure)) deallocate(self%pressure)
    if (allocated(self%dtau)) deallocate(self%dtau)
    if (allocated(self%w0)) deallocate(self%w0)
    if (allocated(self%g0)) deallocate(self%g0)
    self%tau_mode = TAU_MODE_DZ
    self%enabled = .false.
  end subroutine

  subroutine CustomOpticalProperties_validate(self, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(CustomOpticalProperties), intent(in) :: self
    character(:), allocatable, intent(out) :: err
    integer :: n, m

    if (self%tau_mode /= TAU_MODE_DZ .and. self%tau_mode /= TAU_MODE_DLOG10P) then
      err = 'Invalid custom optical properties tau_mode.'
      return
    endif
    if (.not.self%enabled) return
    if (.not.allocated(self%wavelength) .or. .not.allocated(self%pressure) .or. &
        .not.allocated(self%dtau) .or. .not.allocated(self%w0) .or. .not.allocated(self%g0)) then
      err = 'Enabled custom optical properties require all input arrays.'
      return
    endif
    n = size(self%pressure)
    m = size(self%wavelength)
    if (n < 2 .or. m < 2) then
      err = 'Custom optics require at least two pressures and wavelengths.'
      return
    endif
    if (.not.all(ieee_is_finite(self%pressure)) .or. .not.all(ieee_is_finite(self%wavelength))) then
      err = 'Custom optics coordinates must be finite.'
      return
    endif
    if (any(self%pressure <= 0.0_dp) .or. any(self%wavelength <= 0.0_dp) .or. &
        any(self%pressure(2:) >= self%pressure(:n-1)) .or. &
        any(self%wavelength(2:) <= self%wavelength(:m-1))) then
      err = 'Custom optics require positive decreasing pressure and increasing wavelength.'
      return
    endif
    if (any(shape(self%dtau) /= [n,m]) .or. any(shape(self%w0) /= [n,m]) .or. &
        any(shape(self%g0) /= [n,m])) then
      err = 'Custom optical arrays must have shape (n_pressure,n_wavelength).'
      return
    endif
    if (.not.all(ieee_is_finite(self%dtau)) .or. .not.all(ieee_is_finite(self%w0)) .or. &
        .not.all(ieee_is_finite(self%g0))) then
      err = 'Custom optical properties must be finite.'
      return
    endif
    if (any(self%dtau < 0.0_dp) .or. any(self%w0 < 0.0_dp) .or. any(self%w0 > 1.0_dp) .or. &
        any(abs(self%g0) > 1.0_dp)) then
      err = 'Custom optics require dtau >= 0, 0 <= w0 <= 1, and -1 <= g0 <= 1.'
    endif
  end subroutine

  !> Lazily map wavelengths, then interpolate all optics with shared pressure weights.
  !! The first evaluation caches optical arrays at the target wavelength-bin
  !! centers and a single log10-pressure grid. Later calls only map pressure.
  !! A changed target wavelength grid rebuilds the cache. Both coordinates use
  !! constant extrapolation. The per-decade mode converts at layer centers using
  !! dtau*dz/(H*ln(10)). Original arrays remain the authoritative prescription.
  subroutine CustomOpticalProperties_evaluate(self, wavelength_edges, pressure, scale_height, dz, tauc, w0c, g0c, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(CustomOpticalProperties), intent(inout) :: self
    real(dp), intent(in) :: wavelength_edges(:) !! Increasing target bin edges (nm).
    real(dp), intent(in) :: pressure(:) !! Target layer-center hydrostatic pressures (dyn/cm^2).
    real(dp), intent(in) :: scale_height(:), dz(:) !! Layer-center scale heights and thicknesses (cm).
    real(dp), intent(out) :: tauc(:,:), w0c(:,:), g0c(:,:) !! Shape (size(pressure),size(wavelength_edges)-1).
    character(:), allocatable, intent(out) :: err
    integer :: ip(size(pressure))
    real(dp) :: fp(size(pressure)), width(size(pressure))
    real(dp) :: wavelength, fraction, logP
    integer :: i, j, k, lo, hi, mid, nz, nw, np_source, nw_source
    logical :: rebuild_cache

    tauc = 0.0_dp
    w0c = 0.0_dp
    g0c = 0.0_dp
    if (.not.self%enabled) return
    nz = size(pressure)
    nw = size(wavelength_edges)-1
    if (nz < 1 .or. nw < 1 .or. size(scale_height) /= nz .or. size(dz) /= nz) then
      err = 'Custom optics require target layers, wavelength bins, and matching atmospheric array sizes.'
      return
    endif
    if (any(shape(tauc) /= [nz,nw]) .or. any(shape(w0c) /= [nz,nw]) .or. any(shape(g0c) /= [nz,nw])) then
      err = 'Custom optics outputs must have shape (n_layers,n_wavelength_bins).'
      return
    endif
    if (any(.not.ieee_is_finite(pressure)) .or. any(pressure <= 0.0_dp) .or. &
        any(.not.ieee_is_finite(scale_height)) .or. any(scale_height <= 0.0_dp) .or. &
        any(.not.ieee_is_finite(dz)) .or. any(dz <= 0.0_dp)) then
      err = 'Custom optics mapping requires finite positive pressure, scale height, and layer thickness.'
      return
    endif
    if (any(.not.ieee_is_finite(wavelength_edges)) .or. any(wavelength_edges <= 0.0_dp) .or. &
        any(wavelength_edges(2:) <= wavelength_edges(:nw))) then
      err = 'Custom optics target wavelength edges must be finite, positive, and increasing.'
      return
    endif

    rebuild_cache = .true.
    if (allocated(self%wavelength_edges_cache)) then
      if (size(self%wavelength_edges_cache) == size(wavelength_edges)) then
        rebuild_cache = any(self%wavelength_edges_cache /= wavelength_edges)
      endif
    endif
    if (rebuild_cache) then
      call self%validate(err)
      if (allocated(err)) return
      np_source = size(self%pressure)
      nw_source = size(self%wavelength)
      if (allocated(self%dtau_cache)) deallocate(self%dtau_cache, self%w0_cache, self%g0_cache)
      allocate(self%dtau_cache(np_source,nw), self%w0_cache(np_source,nw), self%g0_cache(np_source,nw))
      self%log10_pressure_cache = log10(self%pressure)
      do j = 1,nw
        wavelength = 0.5_dp*(wavelength_edges(j)+wavelength_edges(j+1))
        lo = 1
        hi = nw_source
        do while (hi-lo > 1)
          mid = (lo+hi)/2
          if (wavelength < self%wavelength(mid)) then
            hi = mid
          else
            lo = mid
          endif
        enddo
        wavelength = max(self%wavelength(1), min(self%wavelength(nw_source), wavelength))
        fraction = min(1.0_dp, max(0.0_dp, &
          (wavelength-self%wavelength(lo))/(self%wavelength(lo+1)-self%wavelength(lo))))
        ! All source pressures and all properties share this wavelength bracket.
        self%dtau_cache(:,j) = (1.0_dp-fraction)*self%dtau(:,lo)+fraction*self%dtau(:,lo+1)
        self%w0_cache(:,j) = (1.0_dp-fraction)*self%w0(:,lo)+fraction*self%w0(:,lo+1)
        self%g0_cache(:,j) = (1.0_dp-fraction)*self%g0(:,lo)+fraction*self%g0(:,lo+1)
      enddo
      self%wavelength_edges_cache = wavelength_edges
    endif

    np_source = size(self%log10_pressure_cache)
    ! Locate each layer's pressure bracket once, shared by every wavelength/property.
    do i = 1,nz
      lo = 1
      hi = np_source
      do while (hi-lo > 1)
        mid = (lo+hi)/2
        if (pressure(i) >= self%pressure(mid)) then
          hi = mid
        else
          lo = mid
        endif
      enddo
      ip(i) = lo
      logP = log10(max(self%pressure(np_source), min(self%pressure(1), pressure(i))))
      fp(i) = min(1.0_dp, max(0.0_dp, (logP-self%log10_pressure_cache(lo))/ &
                                   (self%log10_pressure_cache(lo+1)-self%log10_pressure_cache(lo))))
    enddo
    width = dz
    if (self%tau_mode == TAU_MODE_DLOG10P) width = width/(scale_height*log(10.0_dp))
    do j = 1,nw
      do i = 1,nz
        k = ip(i)
        fraction = fp(i)
        tauc(i,j) = ((1.0_dp-fraction)*self%dtau_cache(k,j)+fraction*self%dtau_cache(k+1,j))*width(i)
        w0c(i,j) = (1.0_dp-fraction)*self%w0_cache(k,j)+fraction*self%w0_cache(k+1,j)
        g0c(i,j) = (1.0_dp-fraction)*self%g0_cache(k,j)+fraction*self%g0_cache(k+1,j)
      enddo
    enddo
  end subroutine

  !> Return the layer index nearest a tropopause altitude.
  pure integer function tropopause_layer_index(z, trop_alt) result(trop_ind)
    real(dp), intent(in) :: z(:), trop_alt

    trop_ind = max(minloc(abs(z - trop_alt), 1) - 1, 1)

  end function

  !> Validate the settings owned by the pressure-temperature-eddy profile.
  subroutine PressureTempEddProfile_validate(self, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_enum, only: PeriodicPressTempEdd, ContinuousPressTempEdd
    class(PressureTempEddProfile), intent(in) :: self
    character(:), allocatable, intent(out) :: err

    if (.not.self%enabled) return

    if (self%mode /= PeriodicPressTempEdd .and. self%mode /= ContinuousPressTempEdd) then
      err = 'Unknown pressure-temperature-eddy profile mode.'
    elseif (self%mode == PeriodicPressTempEdd .and. &
            (.not.ieee_is_finite(self%temperature_tol) .or. self%temperature_tol < 0.0_dp)) then
      err = 'Periodic profile temperature tolerance must be finite and nonnegative.'
    elseif (self%mode == PeriodicPressTempEdd .and. &
            (.not.ieee_is_finite(self%edd_tol) .or. self%edd_tol < 0.0_dp)) then
      err = 'Periodic profile eddy-diffusion tolerance must be finite and nonnegative.'
    elseif (self%mode == PeriodicPressTempEdd .and. &
            (.not.ieee_is_finite(self%extreme_factor) .or. self%extreme_factor <= 1.0_dp)) then
      err = 'Periodic profile extreme factor must be finite and greater than one.'
    endif

  end subroutine

  !> Validate the settings owned by TOA-pressure maintenance.
  subroutine TOAPressureMaintenance_validate(self, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(TOAPressureMaintenance), intent(in) :: self
    character(:), allocatable, intent(out) :: err

    if (.not.self%enabled) return

    if (.not.ieee_is_finite(self%target_pressure) .or. self%target_pressure <= 0.0_dp) then
      err = "`toa_pressure_maintenance%target_pressure` must be finite and positive"
    elseif (.not.ieee_is_finite(self%pressure_factor) .or. self%pressure_factor < 1.0_dp) then
      err = "`toa_pressure_maintenance%pressure_factor` must be finite and at least one"
    elseif (.not.ieee_is_finite(self%extreme_pressure_factor) .or. &
            self%extreme_pressure_factor <= self%pressure_factor) then
      err = "`toa_pressure_maintenance%extreme_pressure_factor` must be finite and greater than `pressure_factor`"
    endif

  end subroutine

  !> Construct configured model variables and allocate atmospheric storage.
  function create_PhotochemVars(dat, settings, flux_file, err) result(var)
    type(PhotochemData), intent(in) :: dat
    type(PhotoSettings), intent(in) :: settings
    character(len=*), intent(in) :: flux_file
    character(:), allocatable, intent(out) :: err
    type(PhotochemVars) :: var

    call apply_vars_settings(dat, settings, var, err)
    if (allocated(err)) return

    allocate(var%photon_flux(dat%nw))
    call read_stellar_flux(flux_file, dat%nw, dat%wavl, var%photon_flux, err)
    if (allocated(err)) return

    call allocate_model_grid(dat, var)
  end function

  subroutine apply_vars_settings(dat, s, var, err)
    use photochem_enum, only: VelocityBC, DiffusionLimHydrogenEscape
    type(PhotochemData), intent(in) :: dat
    type(PhotoSettings), intent(in) :: s
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    integer :: i, j, ind(1)

    var%bottom_atmos = 0.0_dp
    var%top_atmos = 0.0_dp
    var%nz = s%nz

    var%surface_albedo = s%surface_albedo
    var%photon_scale_factor = s%photon_scale_factor
    var%solar_zenith = s%solar_zenith

    if (dat%gas_rainout) then
      var%rainfall_rate = s%rainfall_rate
      var%trop_alt = s%trop_alt
    endif

    ! Condensation rate parameters. Size is zero when there are no particles.
    allocate(var%cond_params(dat%np))
    if (allocated(s%particles) .and. dat%there_are_particles) then
      do i = 1,size(s%particles)
        ind = findloc(dat%species_names(1:dat%np), s%particles(i)%name)
        if (ind(1) == 0) then
          err = 'Particle '//s%particles(i)%name// &
                ' in the settings file is not a particle in the reaction file.'
          return
        endif
        var%cond_params(ind(1)) = s%particles(i)%params
      enddo
    endif

    allocate(var%lowerboundcond(dat%nq))
    allocate(var%lower_vdep(dat%nq))
    allocate(var%lower_flux(dat%nq))
    allocate(var%lower_dist_height(dat%nq))
    allocate(var%lower_fix_den(dat%nq))
    allocate(var%lower_fix_press(dat%nq))
    allocate(var%upperboundcond(dat%nq))
    allocate(var%upper_veff(dat%nq))
    allocate(var%upper_flux(dat%nq))
    allocate(var%rate_fcns(dat%nq))

    var%lowerboundcond(:dat%np) = VelocityBC
    var%lowerboundcond(dat%ng_1:) = s%default_lowerboundcond
    var%lower_vdep = 0.0_dp
    var%upperboundcond = VelocityBC
    var%upper_veff = 0.0_dp

    do j = 1,size(s%ubcs)
      ind = findloc(dat%species_names, s%sp_names(j))
      if (ind(1) == 0) then
        err = 'IOError: Species '//trim(s%sp_names(j))// &
              ' in settings file is not in the reaction mechanism file.'
        return
      endif

      if (s%sp_types(j) == 'long lived') then
        var%lowerboundcond(ind(1)) = s%lbcs(j)%bc_type
        var%lower_vdep(ind(1)) = s%lbcs(j)%vel
        var%lower_flux(ind(1)) = s%lbcs(j)%flux
        var%lower_dist_height(ind(1)) = s%lbcs(j)%height
        var%lower_fix_den(ind(1)) = s%lbcs(j)%den
        var%lower_fix_press(ind(1)) = s%lbcs(j)%press

        var%upperboundcond(ind(1)) = s%ubcs(j)%bc_type
        var%upper_veff(ind(1)) = s%ubcs(j)%vel
        var%upper_flux(ind(1)) = s%ubcs(j)%flux
      endif
    enddo

    if (dat%H_escape_type == DiffusionLimHydrogenEscape) then
      if (var%upperboundcond(dat%LH2) /= VelocityBC) then
        err = 'IOError: H2 must have a have a effusion velocity upper boundary'// &
              ' for diffusion limited hydrogen escape'
        return
      endif
      if (var%upperboundcond(dat%LH) /= VelocityBC) then
        err = 'IOError: H must have a have a effusion velocity upper boundary'// &
              ' for diffusion limited hydrogen escape'
        return
      endif
    endif

    if (dat%there_are_particles) then
      do i = 1,dat%npq
        if (var%lowerboundcond(i) /= VelocityBC) then
          err = 'Particle "'//trim(dat%species_names(i))// &
                '" must have deposition velocity lower boundary condition.'
          return
        endif
      enddo
    endif
  end subroutine

  subroutine allocate_model_grid(dat, var)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var

    integer :: i

    var%neqs = dat%nq*var%nz

    allocate(var%temperature(var%nz))
    allocate(var%z(var%nz))
    allocate(var%dz(var%nz))
    allocate(var%edd(var%nz))
    allocate(var%grav(var%nz))
    allocate(var%particle_radius(dat%npq,var%nz))
    allocate(var%xs_x_qy(var%nz,dat%kj,dat%nw))

    allocate(var%particle_xs(dat%np))
    do i = 1,dat%np
      if (dat%part_xs_file(i)%ThereIsData) then
        var%particle_xs(i)%ThereIsData = .true.
        allocate(var%particle_xs(i)%w0(var%nz,dat%nw))
        allocate(var%particle_xs(i)%qext(var%nz,dat%nw))
        allocate(var%particle_xs(i)%gt(var%nz,dat%nw))
      else
        var%particle_xs(i)%ThereIsData = .false.
      endif
    enddo

    if (dat%reverse) allocate(var%gibbs_energy(var%nz,dat%ng))

  end subroutine

  subroutine read_stellar_flux(star_file, nw, wavl, photon_flux, err)
    use futils, only: inter2, addpnt, FileCloser
    use photochem_const, only: c_light, plank
    character(len=*), intent(in) :: star_file
    integer, intent(in) :: nw
    real(dp), intent(in) :: wavl(nw+1)
    real(dp), intent(out) :: photon_flux(nw)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: file_wav(:), file_flux(:)
    real(dp) :: flux(nw)
    real(dp) :: dum1, dum2
    integer :: io, i, n, ierr
    real(dp), parameter :: rdelta = 1.0e-4_dp
    type(FileCloser) :: file

    open(1,file=star_file,status='old',iostat=io)
    file%unit = 1
    if (io /= 0) then
      err = 'The input file '//star_file//' does not exist.'
      return
    endif

    n = -1
    read(1,*)
    do while (io == 0)
      read(1,*,iostat=io) dum1, dum2
      n = n + 1
    enddo

    allocate(file_wav(n+4), file_flux(n+4))

    rewind(1)
    read(1,*)
    do i = 1,n
      read(1,*,iostat=io) file_wav(i), file_flux(i)
      if (io /= 0) then
        err = 'Problem reading '//star_file
        return
      endif
    enddo

    i = n
    call addpnt(file_wav, file_flux, n+4, i, file_wav(1)*(1.0_dp-rdelta), 0.0_dp, ierr)
    call addpnt(file_wav, file_flux, n+4, i, 0.0_dp, 0.0_dp, ierr)
    call addpnt(file_wav, file_flux, n+4, i, file_wav(i)*(1.0_dp+rdelta), 0.0_dp,ierr)
    call addpnt(file_wav, file_flux, n+4, i, huge(rdelta), 0.0_dp,ierr)
    if (ierr /= 0) then
      err = 'Problem interpolating '//trim(star_file)
      return
    endif

    call inter2(nw+1, wavl, flux, n+4, file_wav, file_flux, ierr)
    if (ierr /= 0) then
      err = 'Problem interpolating '//trim(star_file)
      return
    endif

    do i = 1,nw
      photon_flux(i) = (1/(plank*c_light*1.e16_dp))*flux(i)*(wavl(i+1)-wavl(i))* &
                       ((wavl(i+1)+wavl(i))/2.0_dp)
    enddo
  end subroutine

  !> Recompute persistent variables affected by a temperature-profile change.
  subroutine refresh_temperature_dependent_vars(dat, temperature, z, &
                                                bottom_atmos, top_atmos, &
                                                trop_alt_new, xs_x_qy, &
                                                gibbs_energy, trop_alt, &
                                                trop_ind, err)
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: temperature(:), z(:)
    real(dp), intent(in) :: bottom_atmos, top_atmos
    real(dp), optional, intent(in) :: trop_alt_new
    real(dp), intent(inout) :: xs_x_qy(:,:,:)
    real(dp), allocatable, intent(inout) :: gibbs_energy(:,:)
    real(dp), intent(inout) :: trop_alt
    integer, intent(inout) :: trop_ind
    character(:), allocatable, intent(out) :: err

    call interp2xsdata(dat, xs_x_qy)

    if (dat%reverse) then
      call compute_gibbs_energy(dat, temperature, gibbs_energy, err)
      if (allocated(err)) return
    endif

    call set_tropopause(dat, z, bottom_atmos, top_atmos, trop_alt_new, &
                        trop_alt, trop_ind, err)
  end subroutine

  subroutine set_tropopause(dat, z, bottom_atmos, top_atmos, trop_alt_new, &
                            trop_alt, trop_ind, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: z(:), bottom_atmos, top_atmos
    real(dp), optional, intent(in) :: trop_alt_new
    real(dp), intent(inout) :: trop_alt
    integer, intent(inout) :: trop_ind
    character(:), allocatable, intent(out) :: err

    trop_ind = 1
    if (.not. dat%gas_rainout) return

    if (present(trop_alt_new)) trop_alt = trop_alt_new

    if (.not. ieee_is_finite(trop_alt) .or. &
        trop_alt < bottom_atmos .or. trop_alt > top_atmos) then
      err = 'trop_alt is above or bellow the atmosphere!'
      return
    endif

    trop_ind = tropopause_layer_index(z, trop_alt)
    if (trop_ind < 3) then
      err = 'Tropopause is too low.'
      return
    elseif (trop_ind > size(z) - 2) then
      err = 'Tropopause is too high.'
      return
    endif
  end subroutine

  subroutine compute_gibbs_energy(dat, temperature, gibbs_energy, err)
    use photochem_data, only: gibbs_energy_eval
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: temperature(:)
    real(dp), intent(inout) :: gibbs_energy(:,:)
    character(:), allocatable, intent(out) :: err

    integer :: i, j
    logical :: found

    do i = 1,dat%ng
      do j = 1,size(temperature)
        call gibbs_energy_eval(dat%thermo_data(i), temperature(j), &
                               found, gibbs_energy(j,i))
        if (.not. found) then
          err = 'The temperature is not within the ranges '// &
                'given for the thermodynamic data for '// &
                trim(dat%species_names(i+dat%npq))
          return
        endif
      enddo
    enddo
  end subroutine

  subroutine interp2xsdata(dat, xs_x_qy)
    use photochem_const, only: smaller_real
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(out) :: xs_x_qy(:,:,:)

    integer :: i, k

    ! No temperature dependence yet, so copy the static cross sections to
    ! every atmospheric layer.
    do k = 1,dat%nw
      do i = 1,dat%kj
        xs_x_qy(:,i,k) = abs(dat%photolysis_xs(i,k)) + smaller_real
      enddo
    enddo
  end subroutine

  !> Interpolate particle optical properties onto the model particle radii.
  subroutine interp2particlexsdata(dat, particle_radius, particle_xs, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: particle_radius(:,:)
    type(ParticleXsections), intent(inout) :: particle_xs(:)
    character(:), allocatable, intent(out) :: err

    integer :: j, k, lo, hi, middle, nrad, nz
    real(dp) :: radius, fraction

    if (.not. dat%there_are_particles) return

    nz = size(particle_radius, 2)

    ! Validate all radii before changing any output table.
    do k = 1,dat%np
      if (.not.dat%part_xs_file(k)%ThereIsData) cycle
      nrad = size(dat%part_xs_file(k)%radii)
      do j = 1,nz
        radius = particle_radius(k,j)
        if (.not.ieee_is_finite(radius) .or. &
            radius <= dat%part_xs_file(k)%radii(1) .or. &
            radius >= dat%part_xs_file(k)%radii(nrad)) then
          err = 'There is not any optical data for the particle radii '// &
                'specified in the atmosphere.'
          return
        endif
      enddo
    enddo

    do k = 1,dat%np
      if (.not.dat%part_xs_file(k)%ThereIsData) cycle
      nrad = size(dat%part_xs_file(k)%radii)
      do j = 1,nz
        radius = particle_radius(k,j)
        lo = 1
        hi = nrad
        do while (hi - lo > 1)
          middle = lo + (hi - lo)/2
          if (radius < dat%part_xs_file(k)%radii(middle)) then
            hi = middle
          else
            lo = middle
          endif
        enddo
        fraction = (radius - dat%part_xs_file(k)%radii(lo))/ &
                   (dat%part_xs_file(k)%radii(hi) - dat%part_xs_file(k)%radii(lo))

        particle_xs(k)%w0(j,:) = (1.0_dp - fraction)*dat%part_xs_file(k)%w0(lo,:) + &
                                fraction*dat%part_xs_file(k)%w0(hi,:)
        particle_xs(k)%qext(j,:) = (1.0_dp - fraction)*dat%part_xs_file(k)%qext(lo,:) + &
                                  fraction*dat%part_xs_file(k)%qext(hi,:)
        particle_xs(k)%gt(j,:) = (1.0_dp - fraction)*dat%part_xs_file(k)%gt(lo,:) + &
                                fraction*dat%part_xs_file(k)%gt(hi,:)
      enddo
    enddo
  end subroutine

end module
