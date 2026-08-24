module photochem_evoatmosphere
  use photochem_const, only: dp, m_str_len
  use photochem_data, only: PhotochemData, ParticleXsections
  use photochem_vars, only: PhotochemVars, PressureTempEddProfile, &
                            TOAPressureMaintenance
  use photochem_wrk, only: PhotochemWrk
  implicit none

  private
  public :: EvoAtmosphere, ProductionLoss

  !> Process-resolved production and loss rates for one chemical species.
  !!
  !! Instances are returned by [[EvoAtmosphere:production_and_loss]]. Reaction
  !! and process columns are ordered from largest to smallest vertically
  !! integrated rate.
  type :: ProductionLoss
    real(dp), allocatable :: production(:,:) !! Production rates, shape `(nz,nproduction)` (molecules/cm^3/s).
    real(dp), allocatable :: loss(:,:) !! Loss rates, shape `(nz,nloss)` (molecules/cm^3/s).
    real(dp), allocatable :: integrated_production(:) !! Column production rates (molecules/cm^2/s).
    real(dp), allocatable :: integrated_loss(:) !! Column loss rates (molecules/cm^2/s).
    character(len=m_str_len), allocatable :: production_rx(:) !! Reaction or process labels for `production` columns.
    character(len=m_str_len), allocatable :: loss_rx(:) !! Reaction or process labels for `loss` columns.
    real(dp), allocatable :: total_production(:) !! Total production profile (molecules/cm^3/s).
    real(dp), allocatable :: total_loss(:) !! Total loss profile (molecules/cm^3/s).
    real(dp), allocatable :: net(:) !! Net reported tendency, `total_production-total_loss` (molecules/cm^3/s).
    real(dp) :: integrated_total_production = 0.0_dp !! Total column production (molecules/cm^2/s).
    real(dp) :: integrated_total_loss = 0.0_dp !! Total column loss (molecules/cm^2/s).
    real(dp) :: integrated_net = 0.0_dp !! Net column tendency (molecules/cm^2/s).
  end type

  !> Internal transactional handoff used to validate atmospheric changes
  !! before committing them to the live model.
  type :: AtmosphereState
    real(dp) :: bottom_atmos
    real(dp) :: top_atmos
    real(dp) :: trop_alt
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: dz(:)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: edd(:)
    real(dp), allocatable :: particle_radius(:,:)
    real(dp), allocatable :: usol(:,:)

    ! State derived before commit.
    real(dp), allocatable :: grav(:)
    integer :: trop_ind
    real(dp), allocatable :: xs_x_qy(:,:,:)
    type(ParticleXsections), allocatable :: particle_xs(:)
    real(dp), allocatable :: gibbs_energy(:,:)

    ! Persistent atmospheric-profile configuration associated with the state.
    type(PressureTempEddProfile) :: press_temp_edd_profile
    type(TOAPressureMaintenance) :: toa_pressure_maintenance
  contains
    procedure :: allocate => AtmosphereState_allocate
  end type

  !> A one-dimensional photochemical model.
  !!
  !! Construction loads the static mechanism and model configuration. The
  !! atmosphere may be initialized during construction or later with one of
  !! the explicit atmosphere initialization procedures. Once initialized, the
  !! model can prepare chemical rates and integrate toward steady state.
  type :: EvoAtmosphere
    type(PhotochemData), allocatable :: dat !! Static mechanism, species, reaction, and optical data.
    type(PhotochemVars), allocatable :: var !! Prepared model configuration and atmospheric profiles.
    type(PhotochemWrk), allocatable :: wrk !! Mutable integration and atmospheric work state.

    !> True only after atmosphere-dependent initialization and preparation have
    !! completed successfully.
    logical :: atmosphere_initialized = .false.

  contains

    !~~ photochem_evoatmosphere_init.f90 ~~!
    procedure :: initialize_from_atmosphere_file
    procedure :: initialize_atmosphere_z
    procedure :: initialize_atmosphere_p
    procedure, private :: require_atmosphere_initialized

    !~~ photochem_evoatmosphere_rhs.f90 ~~!
    procedure, private :: apply_lower_boundary_conditions
    procedure, private :: prepare_atmosphere_structure
    procedure, private :: prep_atmosphere_unchecked
    procedure :: prep_atmosphere
    procedure :: chemistry_right_hand_side
    procedure :: production_and_loss
    procedure :: right_hand_side
    procedure :: jacobian

    !~~ photochem_evoatmosphere_integrate.f90 ~~!
    procedure :: evolve
    procedure :: check_for_convergence
    procedure :: initialize_stepper
    procedure, private :: initialize_stepper_at_time
    procedure, private :: prepare_stepper_state
    procedure, private :: configure_stepper
    procedure, private :: restart_robust_stepper
    procedure, private :: validate_robust_stepper_settings
    procedure, private :: maybe_maintain_toa_pressure
    procedure :: step
    procedure :: destroy_stepper
    procedure :: initialize_robust_stepper
    procedure :: robust_step
    procedure :: find_steady_state

    !~~ photochem_evoatmosphere_output.f90 ~~!
    procedure :: out2atmosphere_txt

    !~~ photochem_evoatmosphere_boundary.f90 ~~!
    procedure :: gas_fluxes
    procedure :: set_lower_bc
    procedure :: set_upper_bc
    procedure :: set_rate_fcn

    !~~ photochem_evoatmosphere_profiles.f90 ~~!
    procedure :: set_temperature
    procedure :: set_press_temp_edd
    procedure :: set_press_temp_edd_profile
    procedure :: clear_press_temp_edd_profile

    !~~ photochem_evoatmosphere_grid.f90 ~~!
    procedure :: update_vertical_grid

  end type
  interface EvoAtmosphere
    module procedure :: create_EvoAtmosphere
    module procedure :: create_EvoAtmosphere_static
  end interface

  interface

    !~~ photochem_evoatmosphere_init.f90 ~~!

    module subroutine AtmosphereState_allocate(self, dat, nz)
      class(AtmosphereState), intent(inout) :: self
      type(PhotochemData), intent(in) :: dat
      integer, intent(in) :: nz
    end subroutine

    module subroutine finalize_atmosphere_state(dat, state, err)
      type(PhotochemData), intent(in) :: dat
      type(AtmosphereState), intent(inout) :: state
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Construct and initialize a photochemical model from a legacy atmosphere file.
    !!
    !! This compatibility constructor performs static model setup and then
    !! calls [[initialize_from_atmosphere_file]].
    module function create_EvoAtmosphere(mechanism_file, &
                                         settings_file, flux_file, atmosphere_txt, data_dir, err) result(self)
      character(len=*), intent(in) :: mechanism_file !! Path to the reaction mechanism file (yaml format).
      character(len=*), intent(in) :: settings_file !! Path to the settings file (yaml format).
      character(len=*), intent(in) :: flux_file !! Path to the file describing the stellar flux.
      !> Path to the file containing altitude, total number density, temperature, 
      !> eddy diffusion, initial concentrations of each gas (mixing ratios), 
      !> and particle radii.
      character(len=*), intent(in) :: atmosphere_txt
      !> Path to the data directory containing photolysis cross sections and other data
      !> needed to run the model
      character(len=*), intent(in) :: data_dir
      character(:), allocatable, intent(out) :: err
      type(EvoAtmosphere) :: self
    end function

    !> Construct a photochemical model without initializing an atmosphere.
    !!
    !! This reads the mechanism, settings, stellar spectrum, and static optical
    !! data and allocates arrays whose dimensions are known from those inputs.
    !! On success, `atmosphere_initialized` is false. Static configuration may
    !! be inspected or changed, but atmosphere-dependent operations return an
    !! error until an explicit atmosphere initializer succeeds.
    module function create_EvoAtmosphere_static(mechanism_file, settings_file, &
                                                flux_file, data_dir, err) result(self)
      character(len=*), intent(in) :: mechanism_file !! Path to the reaction mechanism file.
      character(len=*), intent(in) :: settings_file !! Path to the settings file.
      character(len=*), intent(in) :: flux_file !! Path to the stellar-flux file.
      !> Path to the data directory containing photolysis cross sections and other data.
      character(len=*), intent(in) :: data_dir
      character(:), allocatable, intent(out) :: err
      type(EvoAtmosphere) :: self
    end function

    !> Initialize atmosphere-dependent model state from a legacy atmosphere file.
    !!
    !! Static setup and work-array allocation must already be complete. On
    !! success, the vertical grid, atmospheric profiles, derived data, and the
    !! initial prepared RHS state are ready for use, and any active integrator
    !! is destroyed. If reading or preparing the replacement atmosphere fails,
    !! the existing atmospheric state and integrator are retained.
    module subroutine initialize_from_atmosphere_file(self, atmosphere_txt, err)
      class(EvoAtmosphere), intent(inout) :: self
      character(len=*), intent(in) :: atmosphere_txt !! Path to the legacy atmosphere text file.
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine copy_model_to_state(self, state)
      class(EvoAtmosphere), intent(in) :: self
      type(AtmosphereState), intent(inout) :: state
    end subroutine

    module subroutine copy_state_to_model(self, state)
      class(EvoAtmosphere), intent(inout) :: self
      type(AtmosphereState), intent(in) :: state
    end subroutine

    !> Initialize atmosphere-dependent model state from altitude-based profiles.
    !!
    !! The input altitude points define the lower and upper edges of the model
    !! domain and must begin at zero. Temperature is interpolated linearly in
    !! altitude. Eddy diffusion, mixing ratios, and particle radii are
    !! interpolated in log10 space. Gas mixing ratios are normalized at each
    !! input point and again on the model grid. Total gas number density is
    !! derived by hydrostatic integration upward from `surface_pressure`.
    !!
    !! Rows of `mix` must follow evolved-species order (`1:dat%nq`), matching
    !! `usol`: particle species first, followed by evolved gases.
    !! `particle_radius` contains the particle species (`1:dat%npq`). Fixed
    !! density and partial-pressure lower boundary conditions override the
    !! corresponding supplied mixing ratios in the bottom model layer.
    !!
    !! This procedure selects altitude-based temperature and eddy-diffusion
    !! profiles, replacing any persistent pressure-based profile. On success,
    !! any active integrator is destroyed and the replacement atmosphere is
    !! committed. If mapping or preparation fails, the existing atmosphere and
    !! integrator are retained.
    module subroutine initialize_atmosphere_z(self, z, temperature, edd, &
                                              surface_pressure, mix, &
                                              particle_radius, err)
      class(EvoAtmosphere), intent(inout) :: self
      real(dp), intent(in) :: z(:) !! Altitude profile knots (cm), including both domain edges.
      real(dp), intent(in) :: temperature(:) !! Temperature at `z` (K).
      real(dp), intent(in) :: edd(:) !! Eddy diffusion at `z` (cm^2/s).
      real(dp), intent(in) :: surface_pressure !! Pressure at the lower domain edge (dyn/cm^2).
      !> Mixing ratios in evolved-species order at `z`, matching `usol`.
      real(dp), intent(in) :: mix(:,:)
      !> Particle radii in mechanism order at `z` (cm).
      real(dp), intent(in) :: particle_radius(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Initialize atmosphere-dependent model state from pressure-based profiles.
    !!
    !! `pressure` must be strictly decreasing. Its first and last points define
    !! the lower and upper boundaries of the model domain, respectively, and
    !! the configured planet radius is interpreted at the lower boundary.
    !! Altitude is constructed by hydrostatic integration using gas composition
    !! only. Rows of `mix` follow evolved-species order (`1:dat%nq`), matching
    !! `usol`: particle species first, followed by evolved gases.
    !!
    !! By default, pressure is used only to construct the initial altitude
    !! grid. If `persistent` is true, temperature and eddy diffusion are also
    !! retained as functions of hydrostatic pressure and reapplied during every
    !! subsequent atmospheric preparation. Composition always remains part of
    !! the evolving ODE state. When gas rainout and persistence are both
    !! enabled, `trop_p` is required. When `persistent` is enabled, TOA-pressure
    !! maintenance is enabled by default and may be configured with
    !! `maintain_toa_pressure` and `target_pressure`.
    !! The requested initial grid is retained even when its TOA pressure lies
    !! outside the maintenance band; robust-stepper initialization performs
    !! the preflight regrid before CVODE starts.
    !!
    !! On success, any active integrator is destroyed and the replacement
    !! atmosphere is committed. If mapping or preparation fails, the existing
    !! atmosphere, persistent profile, and integrator are retained.
    module subroutine initialize_atmosphere_p(self, pressure, temperature, edd, &
                                              mix, particle_radius, persistent, &
                                              trop_p, maintain_toa_pressure, &
                                              target_pressure, err)
      class(EvoAtmosphere), intent(inout) :: self
      !> Strictly decreasing pressure profile knots (dyn/cm^2), including both domain edges.
      real(dp), intent(in) :: pressure(:)
      real(dp), intent(in) :: temperature(:) !! Temperature at `pressure` (K).
      real(dp), intent(in) :: edd(:) !! Eddy diffusion at `pressure` (cm^2/s).
      !> Mixing ratios in evolved-species order at `pressure`, matching `usol`.
      real(dp), intent(in) :: mix(:,:)
      !> Particle radii in mechanism order at `pressure` (cm).
      real(dp), intent(in) :: particle_radius(:,:)
      !> Retain temperature and eddy diffusion as functions of hydrostatic pressure.
      logical, optional, intent(in) :: persistent
      !> Tropopause pressure (dyn/cm^2). May be supplied whenever gas rainout is
      !! enabled; it overrides the settings-file tropopause for the initial state.
      !! It is required when persistence and gas rainout are both enabled.
      real(dp), optional, intent(in) :: trop_p
      !> Enable approximate TOA-pressure maintenance when `persistent` is true.
      logical, optional, intent(in) :: maintain_toa_pressure
      !> Target TOA pressure for approximate maintenance (dyn/cm^2).
      real(dp), optional, intent(in) :: target_pressure
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Return an error when an operation requires an initialized atmosphere.
    module subroutine require_atmosphere_initialized(self, operation, err)
      class(EvoAtmosphere), intent(in) :: self
      character(len=*), intent(in) :: operation !! Name of the attempted operation.
      character(:), allocatable, intent(out) :: err
    end subroutine

    !~~ photochem_evoatmosphere_rhs.f90 ~~!

    !> Apply fixed-density and fixed-pressure lower boundary conditions to a
    !! bottom-layer number-density vector without modifying any other state.
    module subroutine apply_lower_boundary_conditions(self, temperature, usol_bottom, err)
      class(EvoAtmosphere), target, intent(in) :: self
      real(dp), intent(in) :: temperature
      real(dp), intent(inout) :: usol_bottom(:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Copy a number-density state while applying the signed minimum density.
    module subroutine clip_usol(usol_in, usol_out)
      real(dp), intent(in) :: usol_in(:,:)
      real(dp), intent(out) :: usol_out(:,:)
    end subroutine

    !> Prepare composition, pressure, density, and mean molecular mass from an
    !! evolved-species number-density state.
    module subroutine prepare_atmosphere_structure(self, usol_in, usol, &
                                                   molecules_per_particle, pressure, density, mix, mubar, &
                                                   pressure_hydro, density_hydro, apply_persistent_profile, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:)
      real(dp), intent(out) :: usol(:,:)
      real(dp), intent(out) :: molecules_per_particle(:,:)
      real(dp), intent(out) :: pressure(:), density(:), mix(:,:), mubar(:)
      real(dp), intent(out) :: pressure_hydro(:), density_hydro(:)
      logical, optional, intent(in) :: apply_persistent_profile
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Given `usol`, the densities of each species in the atmosphere,
    !> this subroutine calculates reaction rates, photolysis rates, etc.
    !> and puts this information into self.wrk. self.wrk contains all the
    !> information needed for `dochem` to compute chemistry.
    module subroutine prep_atmosphere_unchecked(self, usol_in, apply_persistent_profile, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:) !! Number densities (molecules/cm^3)
      logical, optional, intent(in) :: apply_persistent_profile
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Prepare all atmosphere-dependent work arrays for a number-density state.
    !!
    !! `usol_in` has shape `(dat%nq,var%nz)` and contains evolved-species
    !! number densities. On success, `wrk%usol` and the chemical, transport,
    !! pressure, density, and radiative quantities in `wrk` correspond to that
    !! state.
    module subroutine prep_atmosphere(self, usol_in, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:) !! Number densities (molecules/cm^3)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Compute the chemistry-only tendency for an evolved-species state.
    !!
    !! This prepares `usol`, evaluates chemical production and loss, and adds
    !! configured custom rates and hydrogen escape. Transport is not included.
    !! The result uses altitude-major flattened storage: species vary fastest
    !! within each atmospheric layer.
    module subroutine chemistry_right_hand_side(self, usol, rhs, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol(:,:) !! Number densities, shape `(dat%nq,var%nz)` (molecules/cm^3).
      real(dp), intent(out) :: rhs(:) !! Chemistry tendency, size `var%neqs` (molecules/cm^3/s).
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Compute the full chemistry-and-transport ODE right-hand side.
    !!
    !! `usol_flat` and `rhs` use altitude-major flattened storage: the
    !! `dat%nq` evolved species vary fastest within each atmospheric layer.
    module subroutine right_hand_side(self, neqs, tn, usol_flat, rhs, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      integer, intent(in) :: neqs !! Number of ODE equations, `dat%nq*var%nz`.
      real(dp), intent(in) :: tn !! Integration time (s).
      real(dp), target, intent(in) :: usol_flat(neqs) !! Flattened number-density state (molecules/cm^3).
      real(dp), intent(out) :: rhs(neqs) !! Flattened tendency (molecules/cm^3/s).
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Compute the banded Jacobian of [[EvoAtmosphere:right_hand_side]].
    !!
    !! The input state uses the same altitude-major flattened storage as the
    !! right-hand side. `jac` uses the extended band-matrix storage expected by
    !! the SUNDIALS band linear solver.
    module subroutine jacobian(self, lda_neqs, neqs, usol_flat, jac, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      integer, intent(in) :: lda_neqs !! Size of the extended band-matrix storage.
      integer, intent(in) :: neqs !! Number of ODE equations, `dat%nq*var%nz`.
      real(dp), target, intent(in) :: usol_flat(neqs) !! Flattened number-density state (molecules/cm^3).
      real(dp), intent(out), target :: jac(lda_neqs) !! Extended band-matrix Jacobian storage.
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Compute process-resolved production and loss rates for one species.
    !!
    !! The supplied atmosphere is prepared before rates are evaluated. The
    !! returned contributions are sorted from largest to smallest vertically
    !! integrated rate. The current implementation includes reactions, gas
    !! rainout, condensation and evaporation, custom rates, and Zahnle hydrogen
    !! escape. Transport and boundary fluxes are not yet included.
    module subroutine production_and_loss(self, species, usol, pl, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: species !! Species name.
      real(dp), intent(in) :: usol(:,:) !! Number densities, shape `(dat%nq,var%nz)` (molecules/cm^3).
      type(ProductionLoss), intent(out) :: pl !! Reaction-resolved production and loss information.
      character(:), allocatable, intent(out) :: err
    end subroutine

    !~~ photochem_evoatmosphere_integrate.f90 ~~!

    !> Run a fixed-grid batch integration and save it in a binary Fortran file.
    !!
    !! This is separate from the basic and robust step-by-step pathways: it
    !! creates, owns, and releases its CVODE session internally. The vertical
    !! grid remains fixed, and TOA-pressure maintenance is not performed. Grid
    !! changes must be requested explicitly with
    !! [[EvoAtmosphere:update_vertical_grid]] outside the integration.
    module function evolve(self, filename, tstart, usol_start, t_eval, overwrite, restart_from_file, err) result(success)
      use, intrinsic :: iso_c_binding
      class(EvoAtmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename !! Filename to save results.
      real(c_double), intent(inout) :: tstart !! start time in seconds
      real(dp), intent(inout) :: usol_start(:,:) !! Initial number densities (molecules/cm^3)
      real(c_double), intent(in) :: t_eval(:) !! times to evaluate the solution
      logical, optional, intent(in) :: overwrite !! If true, then overwrites pre-existing files with `filename`
      logical, optional, intent(in) :: restart_from_file !! If true, then the integration restarts from the input file.
      logical :: success !! If True, then integration was successful.
      character(:), allocatable, intent(out) :: err
    end function

    !> Determine whether an active stepper's integration history satisfies the
    !! configured photochemical steady-state criteria.
    !!
    !! This check is valid for both the basic [[EvoAtmosphere:step]] pathway and
    !! the policy-managed [[EvoAtmosphere:robust_step]] pathway. It does not
    !! advance the integration and requires that either stepper has first been
    !! initialized.
    module function check_for_convergence(self, err) result(converged)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      logical :: converged !! True when all configured convergence criteria are satisfied.
    end function

    !> Initialize the basic CVODE stepping pathway at time zero from
    !! `usol_start`.
    !!
    !! Advance this pathway with [[EvoAtmosphere:step]]. The caller is
    !! responsible for responding to solver failures, checking convergence,
    !! and deciding when to restart or stop. Robust recovery, scheduled
    !! restarts, and TOA-pressure maintenance are provided instead by
    !! [[EvoAtmosphere:initialize_robust_stepper]] and
    !! [[EvoAtmosphere:robust_step]]. Any existing stepper is replaced; call
    !! [[EvoAtmosphere:destroy_stepper]] when finished.
    module subroutine initialize_stepper(self, usol_start, err)      
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_start(:,:) !! Initial number densities (molecules/cm^3)
      character(:), allocatable, intent(out) :: err
    end subroutine

    ! Build a fresh CVODE stepper at an explicitly supplied integration time.
    ! This is private so public initialization continues to begin at time zero.
    module subroutine initialize_stepper_at_time(self, usol_start, tstart, err, initial_step)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_start(:,:) !! Initial number densities (molecules/cm^3)
      real(dp), intent(in) :: tstart !! Initial integration time (seconds).
      character(:), allocatable, intent(out) :: err
      real(dp), optional, intent(in) :: initial_step !! Override for CVODE's initial step.
    end subroutine

    ! Load an integration state and reset restart-sensitive Photochem history.
    module subroutine prepare_stepper_state(self, usol_start, tstart, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_start(:,:)
      real(dp), intent(in) :: tstart
      character(:), allocatable, intent(out) :: err
    end subroutine

    ! Apply mutable CVODE settings after either CVodeInit or CVodeReInit.
    module subroutine configure_stepper(self, initial_step, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: initial_step
      character(:), allocatable, intent(out) :: err
    end subroutine

    ! Restart a robust session in place when possible, rebuilding as fallback.
    module subroutine restart_robust_stepper(self, usol_restart, tstart, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_restart(:,:)
      real(dp), intent(in) :: tstart
      character(:), allocatable, intent(out) :: err
    end subroutine

    ! Validate mutable robust-stepper policy before creating or using a session.
    module subroutine validate_robust_stepper_settings(self, err)
      class(EvoAtmosphere), intent(in) :: self
      character(:), allocatable, intent(out) :: err
    end subroutine

    ! Apply one optional pressure-maintenance update after an accepted step.
    module subroutine maybe_maintain_toa_pressure(self, chemistry_converged, updated, failed, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      logical, intent(in) :: chemistry_converged
      logical, intent(out) :: updated
      logical, intent(out) :: failed
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Advance the basic CVODE pathway by one internal solver step.
    !!
    !! [[EvoAtmosphere:initialize_stepper]] must have been called first. This
    !! low-level pathway does not perform robust failure recovery, scheduled
    !! restarts, convergence management, or TOA-pressure maintenance; the
    !! caller is responsible for that policy. Use [[EvoAtmosphere:robust_step]]
    !! for the policy-managed pathway.
    module function step(self, err) result(tn)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      real(dp) :: tn !! Current time in the integration.
    end function
    
    !> Destroy the active basic or robust CVODE stepper and release its
    !! integration resources. It is safe to call this procedure when no
    !! stepper is active.
    module subroutine destroy_stepper(self, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Initialize the policy-managed integration pathway at time zero from
    !! `usol_start`.
    !!
    !! Advance this pathway with [[EvoAtmosphere:robust_step]]. It uses the
    !! same underlying CVODE stepper as the basic pathway, while adding failure
    !! recovery, scheduled restarts, convergence management, and integration
    !! counters. When TOA-pressure maintenance is enabled, the starting
    !! composition is prepared and the model top is brought inside the
    !! configured pressure band before CVODE is initialized. This preflight
    !! allows pressure-based atmosphere initialization to retain its requested
    !! domain endpoints. Total accepted-step and failed-step counters are reset.
    module subroutine initialize_robust_stepper(self, usol_start, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_start(:,:) !! Initial number densities (molecules/cm^3)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Advance the policy-managed integration pathway by one robust step.
    !!
    !! [[EvoAtmosphere:initialize_robust_stepper]] must have been called first.
    !! Unlike the basic [[EvoAtmosphere:step]] pathway, this routine manages
    !! failure recovery, scheduled restarts, convergence checks, and integration
    !! counters. Failed steps recover from the last committed state without
    !! advancing logical time. Scheduled restarts retain logical time and total
    !! counters but discard segment-local convergence history. When TOA-pressure
    !! maintenance is enabled, an accepted step may trigger a pressure-targeted
    !! vertical-grid update after chemistry has converged. A successful update
    !! restarts CVODE while retaining logical time and total integration counters.
    module subroutine robust_step(self, give_up, converged, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      logical, intent(out) :: give_up !! If .true., then the algorithm thinks it is time to give up.
      logical, intent(out) :: converged !! If .true., then the integration has converged to a steady state.
      character(:), allocatable, intent(out) :: err  
    end subroutine
    
    !> Run the policy-managed pathway until convergence or give-up.
    !!
    !! This convenience routine initializes the robust stepper from the current
    !! `wrk%usol`, then repeatedly calls [[EvoAtmosphere:robust_step]]. It
    !! returns true on convergence and false if the robust policy gives up.
    module function find_steady_state(self, err) result(converged)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      logical :: converged !! If .true., then the integration has converged to a steady state.
    end function
    
    module function right_hand_side_callback(tn, sunvec_y, sunvec_f, user_data) &
                          result(ierr) bind(c, name='right_hand_side_callback')
      use, intrinsic :: iso_c_binding
      use fcvode_mod
      use fsundials_nvector_mod
      real(c_double), value :: tn
      type(N_Vector)        :: sunvec_y
      type(N_Vector)        :: sunvec_f
      type(c_ptr), value    :: user_data
      integer(c_int)        :: ierr
    end function
    
    module function jacobian_callback(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                              tmp1, tmp2, tmp3) &
                          result(ierr) bind(C,name='jacobian_callback')
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      use fnvector_serial_mod
      use fsunmatrix_band_mod
      use fsundials_matrix_mod
      real(c_double), value :: tn
      type(N_Vector)        :: sunvec_y 
      type(N_Vector)        :: sunvec_f
      type(SUNMatrix)       :: sunmat_J 
      type(c_ptr), value    :: user_data 
      type(N_Vector)        :: tmp1, tmp2, tmp3
      integer(c_int)        :: ierr
    end function

    !~~ photochem_evoatmosphere_output.f90 ~~!

    !> Save the prepared atmosphere in the legacy text atmosphere format.
    !! The output uses the concentrations currently stored in `self%wrk%usol`.
    module subroutine out2atmosphere_txt(self, filename, number_of_decimals, overwrite, clip, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename !! Output filename
      integer, intent(in) :: number_of_decimals !! Number of decimals
      logical, intent(in) :: overwrite !! If true, then output file can be overwritten, by default False
      !> If true, mixing ratios are clipped at a small positive value.
      logical, intent(in) :: clip
      character(:), allocatable, intent(out) :: err
    end subroutine

    !~~ photochem_evoatmosphere_boundary.f90 ~~!

    !> Computes gas fluxes at model boundaries in order to maintain
    !> current atmospheric concentrations. Uses the densities stored in
    !> self%wrk%usol.
    module subroutine gas_fluxes(self, surf_fluxes, top_fluxes, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(out) :: surf_fluxes(:) !! surface fluxes (molecules/cm^2/s)
      real(dp), intent(out) :: top_fluxes(:) !! top-of-atmosphere fluxes (molecules/cm^2/s)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Set the lower boundary condition for an evolved species.
    !!
    !! Valid `bc_type` values are `vdep`, `den`, `press`, `flux`,
    !! `vdep + dist flux`, and `Moses`. They require `vdep`, `den`, `press`,
    !! `flux`, (`vdep`, `flux`, and `height`), and no additional argument,
    !! respectively. Particle species only support `vdep`.
    module subroutine set_lower_bc(self, species, bc_type, vdep, den, press, flux, height, err)
      class(EvoAtmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species !! Species to set boundary condition
      character(len=*), intent(in) :: bc_type !! Boundary condition type
      real(dp), optional, intent(in) :: vdep !! Deposition velocity (cm/s)
      real(dp), optional, intent(in) :: den !! density (molecules/cm^3)
      real(dp), optional, intent(in) :: press !! pressure (dynes/cm^2)
      real(dp), optional, intent(in) :: flux !! Flux (molecules/cm^2/s)
      real(dp), optional, intent(in) :: height !! Height in atmosphere (km)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Set the upper boundary condition for an evolved species.
    !!
    !! Valid `bc_type` values are `veff` and `flux`, requiring the matching
    !! optional argument. Boundary conditions controlled by configured
    !! diffusion-limited hydrogen escape cannot be replaced.
    module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
      class(EvoAtmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species !! Species to set boundary condition
      character(len=*), intent(in) :: bc_type !! Boundary condition type
      real(dp), optional, intent(in) :: veff !! effusion velocity (cm/s)
      real(dp), optional, intent(in) :: flux !! Flux (molecules/cm^2/s)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Set or clear a custom time-dependent production or loss rate.
    !!
    !! `fcn` is called as `fcn(time, nz, rate)` and fills `rate` with a signed
    !! tendency in molecules/cm^3/s. A disassociated pointer clears the custom
    !! rate for the species.
    module subroutine set_rate_fcn(self, species, fcn, err)
      use photochem_vars, only: time_dependent_rate_fcn
      class(EvoAtmosphere), target, intent(inout) :: self
      character(*), intent(in) :: species !! Species name
      procedure(time_dependent_rate_fcn), pointer :: fcn
      character(:), allocatable, intent(inout) :: err
    end subroutine

    !~~ photochem_evoatmosphere_profiles.f90 ~~!

    !> Changes the altitude-based temperature profile.
    !!
    !! This procedure cannot be used while a CVODE stepper is initialized or
    !! while a persistent pressure-temperature-eddy profile is enabled. Call
    !! `destroy_stepper` or `clear_press_temp_edd_profile` first, respectively.
    module subroutine set_temperature(self, temperature, trop_alt, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: temperature(:) !! new temperature at each atomspheric layer
      real(dp), optional, intent(in) :: trop_alt !! Tropopause altitude (cm). Only necessary if
                                                 !! gas rainout is enabled.
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Maps pressure-temperature and pressure-eddy-diffusion profiles onto
    !! the model's current altitude grid.
    !!
    !! Temperature is interpolated linearly in log10 pressure, while eddy
    !! diffusion is interpolated linearly in log10 pressure-log10 eddy space.
    !! The input pressure must be strictly decreasing, and all three input
    !! arrays must have the same size with at least two elements. If the input
    !! profile does not reach the model surface, its deepest two points are
    !! used to extrapolate it to the surface. Values above the input profile
    !! are held constant.
    !!
    !! By default, the mapping uses hydrostatic pressure and solves for it
    !! sequentially from the bottom of the atmosphere upward. Alternatively,
    !! the mapping can use the actual gas pressure, `density*k_boltz*T`.
    !! `trop_p` is valid only when gas rainout is enabled; when supplied, the
    !! mapped pressure must decrease strictly with altitude so the tropopause
    !! altitude is unambiguous.
    !!
    !! On success, this procedure updates the model temperature, eddy
    !! diffusion, and derived atmospheric working state. Atmospheric species
    !! number densities are not evolved by this procedure. This procedure
    !! cannot be used while a CVODE stepper is initialized or while a
    !! persistent pressure-temperature-eddy profile is enabled; call
    !! `destroy_stepper` or `clear_press_temp_edd_profile` first, respectively.
    module subroutine set_press_temp_edd(self, P, T, edd, trop_p, hydro_pressure, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: P(:) !! Strictly decreasing pressure profile (dynes/cm^2)
      real(dp), intent(in) :: T(:) !! Temperature corresponding to `P` (K)
      real(dp), intent(in) :: edd(:) !! Eddy diffusion corresponding to `P` (cm^2/s)
      !> Tropopause pressure (dynes/cm^2). Only valid and required when gas
      !! rainout is enabled; omit it otherwise.
      real(dp), optional, intent(in) :: trop_p
      !> If .true., then use hydrostatic pressure. If .false. then use the
      !> actual gas pressure in the atmosphere. Default is .true..
      logical, optional, intent(in) :: hydro_pressure
      !> Allocated with an error message if the profiles are invalid or cannot
      !> be mapped onto the altitude grid.
      character(:), allocatable, intent(out) :: err
    end subroutine

    ! Internal, non-mutating mapping kernel shared by the public setter and
    ! routines that need to map a profile for an arbitrary trial composition.
    module subroutine map_press_temp_edd(self, usol, P, T, edd, trop_p, hydro_pressure, &
                                         grid_z, grid_dz, grid_grav, temperature_reference, &
                                         pressure_reference, &
                                         T_grid, edd_grid, log10P_grid, trop_alt, err)
      class(EvoAtmosphere), target, intent(in) :: self
      real(dp), intent(in) :: usol(:,:)
      real(dp), intent(in) :: P(:)
      real(dp), intent(in) :: T(:)
      real(dp), intent(in) :: edd(:)
      real(dp), intent(in) :: trop_p
      logical, intent(in) :: hydro_pressure
      real(dp), intent(in) :: grid_z(:)
      real(dp), intent(in) :: grid_dz(:)
      real(dp), intent(in) :: grid_grav(:)
      real(dp), intent(in) :: temperature_reference(:)
      real(dp), intent(in) :: pressure_reference(:)
      real(dp), intent(out) :: T_grid(:)
      real(dp), intent(out) :: edd_grid(:)
      real(dp), intent(out) :: log10P_grid(:)
      real(dp), intent(out) :: trop_alt
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Prescribes persistent pressure-temperature and pressure-eddy-diffusion
    !! profiles.
    !!
    !! The profiles are mapped onto the altitude grid immediately and again
    !! during every subsequent atmospheric preparation, including every ODE
    !! right-hand-side evaluation. Thus, temperature and eddy diffusion remain
    !! functions of the trial atmospheric pressure as composition evolves.
    !! Temperature is interpolated linearly in log10 pressure, while eddy
    !! diffusion is interpolated linearly in log10 pressure-log10 eddy space.
    !!
    !! The input pressure must be strictly decreasing, and all three arrays
    !! must have the same size with at least two elements. If the input profile
    !! does not reach the model surface, its deepest two points are extrapolated
    !! to the surface. Values above the input profile are held constant.
    !! `trop_p` is valid only when gas rainout is enabled; when supplied, the
    !! mapped pressure must decrease strictly with altitude so the tropopause
    !! altitude is unambiguous.
    !!
    !! A persistent pressure profile enables TOA-pressure maintenance by
    !! default. Supply `maintain_toa_pressure=.false.` to disable it, or
    !! supply `target_pressure` (dynes/cm^2) to choose an explicit target;
    !! when omitted, a target of 0.1 dynes/cm^2 is used.
    !! This procedure cannot be called while a CVODE stepper is initialized.
    !! Call `destroy_stepper` first. While the persistent mode is enabled,
    !! `set_temperature` and `set_press_temp_edd` cannot be used; call
    !! `clear_press_temp_edd_profile` first. Vertical-grid updates preserve and
    !! remap the persistent profiles.
    module subroutine set_press_temp_edd_profile(self, P, T, edd, trop_p, hydro_pressure, &
                                                 maintain_toa_pressure, target_pressure, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: P(:) !! Strictly decreasing pressure profile (dynes/cm^2)
      real(dp), intent(in) :: T(:) !! Temperature corresponding to `P` (K)
      real(dp), intent(in) :: edd(:) !! Eddy diffusion corresponding to `P` (cm^2/s)
      !> Tropopause pressure (dynes/cm^2). Only valid and required when gas
      !! rainout is enabled; omit it otherwise.
      real(dp), optional, intent(in) :: trop_p
      !> If .true., use hydrostatic pressure. If .false., use actual gas
      !! pressure, `density*k_boltz*T`. Default is .true..
      logical, optional, intent(in) :: hydro_pressure
      !> Enable automatic TOA-pressure maintenance. Default is .true..
      logical, optional, intent(in) :: maintain_toa_pressure
      !> Explicit TOA-pressure maintenance target (dynes/cm^2). If omitted
      !! while maintenance is enabled, 0.1 dynes/cm^2 is used.
      real(dp), optional, intent(in) :: target_pressure
      !> Allocated with an error message if the profiles are invalid, cannot be
      !! mapped, or a CVODE stepper is currently initialized.
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Disables the persistent pressure-temperature-eddy profile.
    !!
    !! The most recently mapped altitude-based temperature and eddy-diffusion
    !! profiles remain in place, and TOA-pressure maintenance is disabled.
    !! This procedure cannot be called while a CVODE stepper is initialized;
    !! call `destroy_stepper` first.
    module subroutine clear_press_temp_edd_profile(self, err)
      class(EvoAtmosphere), intent(inout) :: self
      !> Allocated with an error message if a CVODE stepper is initialized.
      character(:), allocatable, intent(out) :: err
    end subroutine

    ! Reset persistent pressure-temperature-eddy profile state without public
    ! lifecycle or integrator checks.
    module subroutine reset_press_temp_edd_profile(var)
      type(PhotochemVars), intent(inout) :: var
    end subroutine

    module subroutine apply_press_temp_edd_profile(self, usol_in, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !~~ photochem_evoatmosphere_grid.f90 ~~!

    !> Rebuilds the vertical grid for a new `TOA_alt` or `TOA_pressure`.
    !! Inside the old model domain, atmospheric properties are interpolated.
    !! Above it, normalized gas mixing ratios, particle abundances relative to
    !! gas, and particle radii are held constant, while gas density is extended
    !! hydrostatically. Altitude-based temperature and eddy diffusion are held
    !! constant above the old domain. If a persistent pressure-based profile is
    !! enabled, its mapped temperature and eddy diffusion are instead reconciled
    !! with the hydrostatic extension before the new grid is committed. This is
    !! a profile-preserving remap, not a column-conservative remap: gas mixing
    !! ratios are normalized on the new grid, but integrated species columns
    !! are not constrained to equal their values on the old grid.
    !!
    !! The update is failure atomic. An error leaves both the atmosphere and any
    !! active CVODE stepper unchanged. A successful update invalidates an active
    !! stepper, which must be initialized again before integration can continue.
    module subroutine update_vertical_grid(self, TOA_alt, TOA_pressure, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), optional, intent(in) :: TOA_alt !! New top of atmosphere altitude (cm)
      real(dp), optional, intent(in) :: TOA_pressure !! New top of atmosphere pressure (dynes/cm^2)
      character(:), allocatable, intent(out) :: err
    end subroutine

  end interface

end module
