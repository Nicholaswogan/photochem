module photochem_vars
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use photochem_const, only: dp
  use photochem_data, only: ParticleXsections
  use photochem_settings, only: CondensationParameters
  implicit none
  private

  public :: PhotochemVars
  public :: PressureTempEddProfile, TOAPressureMaintenance
  public :: time_dependent_flux_fcn, time_dependent_rate_fcn
  public :: binary_diffusion_fcn

  type :: PressureTempEddProfile
    logical :: enabled = .false.
    logical :: hydro_pressure = .true.
    real(dp) :: trop_p = -1.0_dp !! Non-positive values disable the tropopause pressure.
    real(dp), allocatable :: pressure(:)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: edd(:)
  end type

  ! Settings for optional robust-stepper maintenance of the model-top
  ! pressure. The feature is valid only while a persistent pressure-based
  ! temperature and eddy-diffusion profile is enabled.
  type :: TOAPressureMaintenance
    logical :: enabled = .false.
    real(dp) :: target_pressure = 0.0_dp !! Target pressure (dynes/cm^2)
    real(dp) :: pressure_factor = 3.0_dp !! Multiplicative acceptable pressure factor
    integer :: nsteps_between_updates = 100 !! Minimum accepted steps between updates
    integer :: max_failures = 0 !! Failed updates allowed before robust integration stops
  end type

  abstract interface
    !> Custom binary diffusion function
    function binary_diffusion_fcn(mu_i, mubar, T) result(b)
      use iso_c_binding, only: c_double
      real(c_double), value, intent(in) :: mu_i !! molar weight of species i (g/mol)
      real(c_double), value, intent(in) :: mubar !! molar weight of background gas (g/mol)
      real(c_double), value, intent(in) :: T !! Temperature (K)
      real(c_double) :: b !! binary diffusion parameter of species i
                          !! with respect to the background gas (molecules cm^-1 s^1)
    end function

    !> Sets the time-dependent photon flux
    subroutine time_dependent_flux_fcn(tn, nw, photon_flux)
      use iso_c_binding, only: c_double, c_int
      real(c_double), value, intent(in) :: tn
      integer(c_int), value, intent(in) :: nw
      real(c_double), intent(out) :: photon_flux(nw)
    end subroutine

    !> Sets a production or destruction rate for a molecule
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

  type :: PhotochemVars
    ! PhotochemVars contains information that can change between
    ! different photochemical integrations, without reading in new
    ! files.

    !> where the photochem data is
    character(:), allocatable :: data_dir

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
    real(dp) :: bottom_atmos = 0.0_dp !! cm
    real(dp) :: top_atmos = 0.0_dp !! cm; determined during atmosphere initialization
    integer :: nz !! number of vertical layers
    real(dp) :: surface_albedo
    real(dp) :: diurnal_fac = 0.5_dp !! Default is 0.5, to account for half planet facing the sun.
    real(dp) :: solar_zenith !! degrees
    real(dp) :: trop_alt = 0.0_dp !! cm (only for gas_rainout == true)
    real(dp) :: rainfall_rate !! relative to modern Earth's average rainfall rate of 1.1e17 molecules/cm2/s
    integer :: trop_ind !! index of troposphere (only for gas_rainout == true)

    ! radiative transfer
    real(dp), allocatable :: photon_flux(:) !! (nw) photon/cm^2/s in each wavelength bin hitting planet.
    !> for scaling photon flux for different planets in a solar system
    real(dp) :: photon_scale_factor
    !> A function for altering the photon flux over time
    procedure(time_dependent_flux_fcn), nopass, pointer :: photon_flux_fcn => null()

    ! particles
    !> Parameters describing condensation and evaporation rates and
    !> the RH needed for condensation
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
    real(dp), allocatable :: temperature(:) !! (nz) K
    real(dp), allocatable :: z(:) !! (nz) cm
    real(dp), allocatable :: dz(:) !! (nz) cm
    real(dp), allocatable :: edd(:) !! (nz) cm2/s
    !> State for a persistent pressure-based temperature and Kzz profile.
    type(PressureTempEddProfile) :: press_temp_edd_profile
    !> A function for specifying a custom binary diffusion parameter (b_ij)
    procedure(binary_diffusion_fcn), nopass, pointer :: custom_binary_diffusion_fcn => null()
    real(dp), allocatable :: grav(:) !! (nz) cm/s2
    real(dp), allocatable :: particle_radius(:,:) !! (np,nz) cm
    real(dp), allocatable :: xs_x_qy(:,:,:) !! (nz,kj,nw) photolysis cross sections times quantum yields (cm2/molecule)
    type(ParticleXsections), allocatable :: particle_xs(:) !! (np)
    real(dp), allocatable :: gibbs_energy(:,:) !! (nz,ng) Joules/mol

    ! Custom optical properties
    real(dp), allocatable :: tauc(:,:) !! (nz,nw) Custom optical depth in each layer
    real(dp), allocatable :: w0c(:,:) !! (nz,nw) Custom single scattering albedo
    real(dp), allocatable :: g0c(:,:) !! (nz,nw) Custom asymetry parameter

    ! other
    !> number of times we initialize CVODE when it returns
    !> a potentially recoverable error. ONLY USED IN EVOATMOSPHERE (NOT ATMOSPHERE)
    !> in the `evolve` method.
    integer :: max_error_reinit_attempts = 2
    real(c_double) :: rtol = 1.0e-3_dp !! integration relative tolerance
    !> Integration absolute tolerance. If autodiff == .true., then the model
    !> works better when atol is smaller (e.g., atol = ~1.0e-18).
    real(c_double) :: atol = 1.0e-23_dp
    integer :: mxsteps = 100000 !! max number of steps before integrator will give up.
    !> seconds. atomsphere considered in equilibrium if integrations reaches this time.
    real(dp) :: equilibrium_time = 1.0e17_dp
    !> For convergence checking. Considers mixing ratio change between t_now and time
    !> t = t_now*conv_hist_factor to see if atmosphere is changing.
    real(dp) :: conv_hist_factor = 0.5_dp
    !> Minimum mixing ratio considered in convergence checking.
    real(dp) :: conv_min_mix = 1.0e-20_dp
    !> Threshold normalized change in mixing ratios for converchecking check.
    !> A reasonable value is ~1.0e-2.
    real(dp) :: conv_longdy = 0.0_dp
    !> Threshold normalized change in mixing ratios per time change for
    !> convergence checking.
    real(dp) :: conv_longdydt = 1.0e-6_dp
    real(c_double) :: initial_dt = 1.0e-6_dp !! intial timestep size (seconds)
    !> Maximum time step size (seconds).
    real(c_double) :: max_dt = sqrt(huge(1.0_dp))
    integer(c_int) :: max_err_test_failures = 15 !! CVODE max error test failures
    integer(c_int) :: max_order = 5 !! CVODE max order for BDF method.
    !> If .true., then the chemistry terms of the Jacobian are computed uses
    !> foward mode automatic differentiation.
    logical :: autodiff = .true.
    !> Perturbation for finite difference Jacobian calculation, when autodiff == .false.
    real(dp) :: epsj = 1.0e-4_dp
    integer :: verbose = 1 !! 0 == no printing. 1 == some printing. 2 == bunch of printing.
    !> If True, then the code uses a 1st order upwind method for the advective molecular
    !> diffusion terms instead of a centered scheme. This permits stability (at the cost
    !> of accuracy) for atmospheres with strong molcular advection in the upper atmosphere.
    logical :: upwind_molec_diff = .false.

    ! Settings for the `robust_stepper` and `find_steady_state` methods
    !> Number of failed-step recovery restarts allowed. The next integration
    !> error ends the robust session.
    integer :: nerrors_before_giveup = 10
    !> Number of accepted steps after initialization or restart to take before
    !> checking the non-time convergence criteria.
    integer :: nsteps_before_conv_check = 300
    !> Accepted steps per integration segment. At this count CVODE is restarted
    !> and the segment-local convergence history is discarded.
    integer :: nsteps_before_reinit = 1000
    !> Maximum total accepted steps before returning give_up.
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

end module
