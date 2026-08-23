module photochem_vars
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use photochem_const, only: dp
  use photochem_data, only: PhotochemData, ParticleXsections
  use photochem_enum, only: AutodiffJacobian
  use photochem_settings, only: PhotoSettings, CondensationParameters
  implicit none
  private

  public :: PhotochemVars
  public :: PressureTempEddProfile, TOAPressureMaintenance
  public :: time_dependent_flux_fcn, time_dependent_rate_fcn
  public :: binary_diffusion_fcn
  public :: refresh_temperature_dependent_vars, interp2particlexsdata

  !> Persistent temperature and eddy-diffusion profiles defined on pressure.
  type :: PressureTempEddProfile
    logical :: enabled = .false.
    logical :: hydro_pressure = .true.
    real(dp) :: trop_p = -1.0_dp !! Non-positive values disable the tropopause pressure.
    real(dp), allocatable :: pressure(:)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: edd(:)
  end type

  !> Settings for optional robust-stepper maintenance of model-top pressure.
  !! This requires a persistent pressure-based temperature and eddy-diffusion profile.
  type :: TOAPressureMaintenance
    logical :: enabled = .false.
    real(dp) :: target_pressure = 0.0_dp !! Target pressure (dynes/cm^2)
    real(dp) :: pressure_factor = 3.0_dp !! Multiplicative acceptable pressure factor
    integer :: nsteps_between_updates = 100 !! Minimum accepted steps between updates
    integer :: max_failures = 0 !! Failed updates allowed before robust integration stops
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
    !> Integration absolute tolerance. Automatic differentiation generally
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
    !> Method used to compute the chemistry Jacobian. See photochem_enum.f90.
    integer(c_int) :: jacobian_method = AutodiffJacobian
    !> Relative perturbation used by the finite-difference Jacobian method.
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

  interface PhotochemVars
    module procedure :: create_PhotochemVars
  end interface

contains

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

    allocate(var%tauc(var%nz,dat%nw))
    var%tauc = 0.0_dp
    allocate(var%w0c(var%nz,dat%nw))
    var%w0c = 0.0_dp
    allocate(var%g0c(var%nz,dat%nw))
    var%g0c = 0.0_dp
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

    trop_ind = max(minloc(abs(z - trop_alt), 1) - 1, 1)
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
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: particle_radius(:,:)
    type(ParticleXsections), intent(inout) :: particle_xs(:)
    character(:), allocatable, intent(out) :: err

    integer :: i, j, k, jj, nz
    real(dp) :: dr, slope, intercept

    if (.not. dat%there_are_particles) return

    nz = size(particle_radius, 2)

    do j = 1,nz
      do k = 1,dat%np
        ! If there is optical data, check that it covers the particle radii
        ! in the atmosphere.
        if (dat%part_xs_file(k)%ThereIsData) then
          if (particle_radius(k,j) <= dat%radii_file(1,k) .or. &
              particle_radius(k,j) >= dat%radii_file(dat%nrad_file,k)) then
            err = 'There is not any optical data for the particle radii '// &
                  'specified in the atmosphere.'
            return
          endif
        endif
      enddo
    enddo

    do i = 1,dat%nw
      do j = 1,nz
        do k = 1,dat%np
          if (dat%part_xs_file(k)%ThereIsData) then
            do jj = 1,dat%nrad_file-1
              if (particle_radius(k,j) >= dat%radii_file(jj,k) .and. &
                  particle_radius(k,j) < dat%radii_file(jj+1,k)) then
                dr = dat%radii_file(jj+1,k) - dat%radii_file(jj,k)

                slope = (dat%part_xs_file(k)%w0(jj+1,i) - &
                         dat%part_xs_file(k)%w0(jj,i))/dr
                intercept = dat%part_xs_file(k)%w0(jj,i) - &
                            dat%radii_file(jj,k)*slope
                particle_xs(k)%w0(j,i) = slope*particle_radius(k,j) + intercept

                slope = (dat%part_xs_file(k)%qext(jj+1,i) - &
                         dat%part_xs_file(k)%qext(jj,i))/dr
                intercept = dat%part_xs_file(k)%qext(jj,i) - &
                            dat%radii_file(jj,k)*slope
                particle_xs(k)%qext(j,i) = slope*particle_radius(k,j) + intercept

                slope = (dat%part_xs_file(k)%gt(jj+1,i) - &
                         dat%part_xs_file(k)%gt(jj,i))/dr
                intercept = dat%part_xs_file(k)%gt(jj,i) - &
                            dat%radii_file(jj,k)*slope
                particle_xs(k)%gt(j,i) = slope*particle_radius(k,j) + intercept
              endif
            enddo
          endif
        enddo
      enddo
    enddo
  end subroutine

end module
