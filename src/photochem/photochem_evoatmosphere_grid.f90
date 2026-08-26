
submodule(photochem_evoatmosphere) photochem_evoatmosphere_grid
  implicit none

  ! Scratch storage used only while constructing a vertical-grid state.
  type :: VerticalGridWork
    real(dp), allocatable :: pressure(:)
    real(dp), allocatable :: mix(:,:), mix_new(:,:)
    real(dp), allocatable :: density(:), density_new(:)
    real(dp), allocatable :: temperature_reference(:)
    real(dp), allocatable :: pressure_reference(:)
    real(dp), allocatable :: source_pressure(:)
    real(dp), allocatable :: log10P(:)
  end type

contains

  function TOA_at_pressure(self, TOA_pressure, state, work, err) result(top_atmos)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use futils, only: brent_class
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: TOA_pressure !! Target top-of-atmosphere pressure (dyn/cm^2).
    type(AtmosphereState), intent(inout) :: state
    type(VerticalGridWork), intent(inout) :: work
    character(:), allocatable, intent(out) :: err
    real(dp) :: top_atmos !! cm

    integer, parameter :: max_bracket_expansions = 60
    real(dp), parameter :: pressure_residual_tolerance = 1.0e-12_dp
    real(dp) :: bottom, current_top, minimum_top, span, search_distance
    real(dp) :: lower, upper, flower, fupper, fcurrent, fnew, fzero
    real(dp) :: altitude_tolerance
    integer :: expansion, iflag
    type(brent_class) :: root_solver

    bottom = self%var%bottom_atmos
    current_top = self%var%top_atmos
    span = current_top - bottom
    ! Lowest trial model top allowed during downward bracket search.
    minimum_top = bottom + max(100.0_dp*epsilon(1.0_dp)* &
                               max(abs(bottom),span,1.0_dp), tiny(1.0_dp))
    ! Absolute altitude tolerance passed to the Brent root solver.
    altitude_tolerance = max(1.0e-10_dp*span, &
                             100.0_dp*epsilon(1.0_dp)*max(abs(current_top),1.0_dp))
    ! Initial altitude increment used to expand the pressure bracket.
    search_distance = max(0.05_dp*span, maxval(self%var%dz))

    ! Set the residual function in the root solver.
    call root_solver%set_function(altitude_residual)
    fcurrent = altitude_residual(root_solver, current_top)
    if (allocated(err)) then
      err = 'Could not evaluate TOA pressure at the current model top: '//err
      return
    endif
    if (abs(fcurrent) <= pressure_residual_tolerance) then
      ! We are at the right pressure, so we can just return
      top_atmos = current_top
      return
    endif

    ! Now determine the bracket for which to do the solve
    if (fcurrent > 0.0_dp) then
      ! The requested pressure is lower than the current TOA pressure, so the
      ! model top must move upward.
      lower = current_top
      flower = fcurrent
      fnew = flower
      do expansion = 1,max_bracket_expansions
        upper = current_top+search_distance
        fupper = altitude_residual(root_solver, upper)
        if (allocated(err)) then
          err = 'Could not bracket the requested TOA pressure while raising the model top: '//err
          return
        endif
        if (fupper >= fnew) then
          err = 'TOA pressure did not decrease while raising the model top; '// &
                'a monotonic pressure bracket could not be constructed.'
          return
        endif
        if (fupper <= 0.0_dp) exit
        fnew = fupper
        search_distance = 2.0_dp*search_distance
      enddo
      if (fupper > 0.0_dp) then
        err = 'Could not bracket the requested TOA pressure above the current model top.'
        return
      endif
    else
      ! The requested pressure is higher than the current TOA pressure, so the
      ! model top must move downward without crossing the model bottom.
      upper = current_top
      fupper = fcurrent
      fnew = fupper
      do expansion = 1,max_bracket_expansions
        lower = max(current_top-search_distance,minimum_top)
        flower = altitude_residual(root_solver, lower)
        if (allocated(err)) then
          err = 'Could not bracket the requested TOA pressure while lowering the model top: '//err
          return
        endif
        if (flower <= fnew) then
          err = 'TOA pressure did not increase while lowering the model top; '// &
                'a monotonic pressure bracket could not be constructed.'
          return
        endif
        if (flower >= 0.0_dp) exit
        if (lower == minimum_top) then
          err = 'The requested TOA pressure exceeds the maximum reachable pressure '// &
                'above the model bottom.'
          return
        endif
        fnew = flower
        search_distance = 2.0_dp*search_distance
      enddo
      if (flower < 0.0_dp) then
        err = 'Could not bracket the requested TOA pressure below the current model top.'
        return
      endif
    endif

    ! Call the root solver for the determined bracket
    call root_solver%find_zero(lower, upper, altitude_tolerance, top_atmos, &
                               fzero, iflag, flower, fupper)
    if (allocated(err)) then
      err = 'Evaluating the bracketed TOA-pressure solve failed: '//err
      return
    endif
    if (iflag /= 0 .or. .not.ieee_is_finite(top_atmos) .or. &
        top_atmos <= bottom .or. .not.ieee_is_finite(fzero)) then
      err = 'The bracketed TOA-pressure solve failed.'
      return
    endif

    ! Successful if reached the end of the function.

  contains
    function altitude_residual(me, altitude) result(residual)
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: altitude
      real(dp) :: residual

      real(dp) :: pressure

      call build_vertical_grid_state(self, altitude, state, work, err)
      if (allocated(err)) then
        residual = ieee_value(0.0_dp,ieee_quiet_nan)
        return
      endif
      pressure = work%pressure(size(work%pressure))
      if (.not.ieee_is_finite(pressure) .or. pressure <= 0.0_dp) then
        err = 'Candidate TOA pressure was not finite and positive.'
        residual = ieee_value(0.0_dp,ieee_quiet_nan)
        return
      endif
      residual = log10(pressure)-log10(TOA_pressure)
    end function
  end function

  subroutine build_vertical_grid_state(self, top_atmos_new, state, work, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: top_atmos_new !! cm
    type(AtmosphereState), intent(inout) :: state
    type(VerticalGridWork), intent(inout) :: work
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    dat => self%dat

    if (.not.ieee_is_finite(top_atmos_new) .or. top_atmos_new <= self%var%bottom_atmos) then
      err = 'The candidate TOA altitude must be finite and above the model bottom.'
      return
    endif
    if (size(self%wrk%usol,1) /= dat%nq .or. &
        size(self%wrk%usol,2) /= self%var%nz) then
      err = 'The candidate regrid state has the wrong dimensions.'
      return
    endif
    if (.not.all(ieee_is_finite(self%wrk%usol))) then
      err = 'The candidate regrid state must be finite.'
      return
    endif

    if (self%var%press_temp_edd_profile%enabled) then
      call build_vertical_grid_state_p(self, top_atmos_new, state, work, err)
    else
      call build_vertical_grid_state_z(self, top_atmos_new, state, work, err)
    endif

  end subroutine

  subroutine initialize_candidate_grid(self, top_atmos_new, state, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_eqns, only: vertical_grid, gravity
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: top_atmos_new !! cm
    type(AtmosphereState), intent(inout) :: state
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    dat => self%dat

    state%top_atmos = top_atmos_new
    state%trop_alt = self%var%trop_alt

    call vertical_grid(self%var%bottom_atmos, top_atmos_new, &
                       state%z, state%dz)
    call gravity(dat%planet_radius, dat%planet_mass, &
                 state%z, state%grav)
    if (.not.all(ieee_is_finite(state%z)) .or. &
        .not.all(ieee_is_finite(state%dz)) .or. &
        .not.all(ieee_is_finite(state%grav)) .or. &
        any(state%dz <= 0.0_dp) .or. any(state%grav <= 0.0_dp)) then
      err = 'The candidate TOA altitude produced invalid grid geometry or gravity.'
      return
    endif

  end subroutine

  subroutine seed_candidate_temperature_edd(self, state, err)
    use futils, only: interp
    class(EvoAtmosphere), target, intent(in) :: self
    type(AtmosphereState), intent(inout) :: state
    character(:), allocatable, intent(out) :: err

    integer :: ierr

    ! Map the altitude-based profiles. This is the final profile mapping for
    ! the altitude path and an initial guess for the pressure-profile path.
    call interp(state%z, self%var%z, self%var%temperature, &
                state%temperature, ierr=ierr)
    if (ierr /= 0) then
      err = 'Temperature interpolation failed while constructing a candidate vertical grid.'
      return
    endif
    call interp(state%z, self%var%z, log10(max(self%var%edd,1.0e-40_dp)), &
                state%edd, ierr=ierr)
    if (ierr /= 0) then
      err = 'Eddy-diffusion interpolation failed while constructing a candidate vertical grid.'
      return
    endif
    state%edd = 10.0_dp**state%edd
    call validate_candidate_temperature_edd(state, err)
    if (allocated(err)) return

  end subroutine

  subroutine validate_candidate_temperature_edd(state, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(AtmosphereState), intent(in) :: state
    character(:), allocatable, intent(out) :: err

    if (.not.all(ieee_is_finite(state%temperature)) .or. &
        any(state%temperature <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(state%edd)) .or. &
        any(state%edd <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid temperature or eddy diffusion.'
      return
    endif

  end subroutine

  subroutine prepare_candidate_composition(self, state, work, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_const, only: small_real
    class(EvoAtmosphere), target, intent(in) :: self
    type(AtmosphereState), intent(in) :: state
    type(VerticalGridWork), intent(inout) :: work
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    real(dp) :: gas_mix_total
    integer :: i, j, ierr

    dat => self%dat

    ! Convert the input atmosphere to mixing ratios and total gas density.
    do j = 1,self%var%nz
      work%density(j) = sum(self%wrk%usol(dat%ng_1:,j))
      if (.not.ieee_is_finite(work%density(j)) .or. &
          work%density(j) <= 0.0_dp) then
        err = 'Gas density must be finite and positive to construct a candidate vertical grid.'
        return
      endif
      work%mix(:,j) = self%wrk%usol(:,j)/work%density(j)
    enddo

    ! Interpolate mixing ratios in log-space. Gas mixing ratios are
    ! normalized below; particle abundances remain ratios to total gas.
    do i = 1,dat%nq
      call interp(state%z, self%var%z, &
                  log10(max(work%mix(i,:),small_real)), &
                  work%mix_new(i,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Mixing-ratio interpolation failed while constructing a candidate vertical grid.'
        return
      endif
    enddo
    work%mix_new = 10.0_dp**work%mix_new
    if (.not.all(ieee_is_finite(work%mix_new))) then
      err = 'The candidate vertical grid produced nonfinite mixing ratios.'
      return
    endif

    ! Normalize gas mixing ratios so they sum to one. Particle abundances
    ! remain ratios relative to total gas density.
    do j = 1,size(state%z)
      gas_mix_total = sum(work%mix_new(dat%ng_1:dat%nq,j))
      if (.not.ieee_is_finite(gas_mix_total) .or. gas_mix_total <= 0.0_dp) then
        err = 'The candidate vertical grid produced invalid gas mixing ratios.'
        return
      endif
      work%mix_new(dat%ng_1:dat%nq,j) = &
          work%mix_new(dat%ng_1:dat%nq,j)/gas_mix_total
    enddo

  end subroutine

  subroutine extend_candidate_density(self, state, work, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_const, only: k_boltz, N_avo
    class(EvoAtmosphere), target, intent(in) :: self
    type(AtmosphereState), intent(in) :: state
    type(VerticalGridWork), intent(inout) :: work
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    real(dp) :: pressure_previous, temperature_previous
    real(dp) :: delta_z, mubar
    integer :: j, ierr, first_extended

    dat => self%dat

    ! Interpolate density inside the old grid. Above the old top, replace the
    ! interpolation with a hydrostatic continuation using candidate T and z.
    call interp(state%z, self%var%z, log10(work%density), &
                work%density_new, ierr=ierr)
    if (ierr /= 0) then
      err = 'Gas-density interpolation failed while constructing a candidate vertical grid.'
      return
    endif
    work%density_new = 10.0_dp**work%density_new
    if (.not.all(ieee_is_finite(work%density_new)) .or. &
        any(work%density_new <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid gas density.'
      return
    endif

    ! If the model was extended upward, replace the extrapolation above the
    ! old top with a hydrostatic continuation.
    first_extended = self%var%nz + 1
    do j = 1,size(state%z)
      if (state%z(j) > self%var%z(self%var%nz)) then
        first_extended = j
        exit
      endif
    enddo
    if (first_extended > size(state%z)) return

    pressure_previous = work%density(self%var%nz)*k_boltz* &
                        self%var%temperature(self%var%nz)
    temperature_previous = self%var%temperature(self%var%nz)
    delta_z = state%z(first_extended)-self%var%z(self%var%nz)
    do j = first_extended,size(state%z)
      if (j > first_extended) delta_z = state%z(j)-state%z(j-1)
      mubar = sum(dat%species_mass(dat%ng_1:dat%nq)* &
                  work%mix_new(dat%ng_1:dat%nq,j))
      pressure_previous = pressure_previous*exp( &
          -(mubar*state%grav(j)*delta_z)/ &
           (N_avo*k_boltz*0.5_dp* &
            (temperature_previous+state%temperature(j))))
      work%density_new(j) = pressure_previous/(k_boltz*state%temperature(j))
      temperature_previous = state%temperature(j)
    enddo

  end subroutine

  subroutine map_candidate_particle_radii(self, state, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_const, only: small_real
    class(EvoAtmosphere), target, intent(in) :: self
    type(AtmosphereState), intent(inout) :: state
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    integer :: i, ierr

    dat => self%dat

    if (.not.dat%there_are_particles) return

    do i = 1,dat%npq
      call interp(state%z, self%var%z, &
                  log10(max(self%var%particle_radius(i,:),small_real)), &
                  state%particle_radius(i,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Particle-radius interpolation failed while constructing a candidate vertical grid.'
        return
      endif
    enddo
    state%particle_radius = 10.0_dp**state%particle_radius
    if (.not.all(ieee_is_finite(state%particle_radius)) .or. &
        any(state%particle_radius <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid particle radii.'
      return
    endif

  end subroutine

  subroutine assemble_candidate_atmosphere(self, state, work, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_const, only: k_boltz
    class(EvoAtmosphere), target, intent(in) :: self
    type(AtmosphereState), intent(inout) :: state
    type(VerticalGridWork), intent(inout) :: work
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    integer :: j

    dat => self%dat

    ! Construct the candidate atmospheric state and enforce the lower BCs.
    do j = 1,size(state%z)
      state%usol(:,j) = work%mix_new(:,j)*work%density_new(j)
    enddo
    call self%apply_lower_boundary_conditions(state%temperature(1), &
                                              state%usol(:,1), err)
    if (allocated(err)) return

    do j = 1,size(state%z)
      work%pressure(j) = sum(state%usol(dat%ng_1:,j))* &
                         k_boltz*state%temperature(j)
    enddo
    if (.not.all(ieee_is_finite(state%usol)) .or. &
        .not.all(ieee_is_finite(work%pressure)) .or. &
        any(work%pressure <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid composition or pressure.'
      return
    endif

  end subroutine

  subroutine build_vertical_grid_state_z(self, top_atmos_new, state, work, err)
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: top_atmos_new !! cm
    type(AtmosphereState), intent(inout) :: state
    type(VerticalGridWork), intent(inout) :: work
    character(:), allocatable, intent(out) :: err

    call initialize_candidate_grid(self, top_atmos_new, state, err)
    if (allocated(err)) return

    call seed_candidate_temperature_edd(self, state, err)
    if (allocated(err)) return

    call prepare_candidate_composition(self, state, work, err)
    if (allocated(err)) return

    call extend_candidate_density(self, state, work, err)
    if (allocated(err)) return

    call map_candidate_particle_radii(self, state, err)
    if (allocated(err)) return

    call assemble_candidate_atmosphere(self, state, work, err)

  end subroutine

  subroutine build_vertical_grid_state_p(self, top_atmos_new, state, work, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: top_atmos_new !! cm
    type(AtmosphereState), intent(inout) :: state
    type(VerticalGridWork), intent(inout) :: work
    character(:), allocatable, intent(out) :: err

    real(dp), parameter :: persistent_tolerance = 1.0e-10_dp
    integer, parameter :: persistent_max_iterations = 50
    real(dp) :: temperature_change
    integer :: iteration
    logical :: converged

    ! Build an initial candidate using the shared altitude-based seed. The
    ! pressure-profile mapping replaces that seed during the fixed-point loop.
    call initialize_candidate_grid(self, top_atmos_new, state, err)
    if (allocated(err)) return
    call seed_candidate_temperature_edd(self, state, err)
    if (allocated(err)) return
    call prepare_candidate_composition(self, state, work, err)
    if (allocated(err)) return
    call extend_candidate_density(self, state, work, err)
    if (allocated(err)) return
    call assemble_candidate_atmosphere(self, state, work, err)
    if (allocated(err)) return

    ! Use the current hydrostatic pressure as the initial reference for the
    ! pressure-profile solve. Each later iteration uses the pressure profile
    ! produced by the preceding mapping pass.
    work%pressure_reference = work%source_pressure
    converged = .false.
    do iteration = 1,persistent_max_iterations
      work%temperature_reference = state%temperature
      call map_press_temp_edd( &
        self, &
        state%usol, &
        self%var%press_temp_edd_profile%pressure, &
        self%var%press_temp_edd_profile%temperature, &
        self%var%press_temp_edd_profile%edd, &
        trop_p=self%var%press_temp_edd_profile%trop_p, &
        hydro_pressure=self%var%press_temp_edd_profile%hydro_pressure, &
        grid_z=state%z, &
        grid_dz=state%dz, &
        grid_grav=state%grav, &
        temperature_reference=work%temperature_reference, &
        pressure_reference=work%pressure_reference, &
        T_grid=state%temperature, &
        edd_grid=state%edd, &
        log10P_grid=work%log10P, &
        trop_alt=state%trop_alt, &
        err=err &
      )
      if (allocated(err)) return

      call validate_candidate_temperature_edd(state, err)
      if (allocated(err)) return
      if (.not.all(ieee_is_finite(work%log10P))) then
        err = 'The pressure-profile mapping produced nonfinite log-pressure values.'
        return
      endif

      temperature_change = maxval(abs(state%temperature - &
                                      work%temperature_reference)/ &
                                  max(state%temperature,1.0_dp))
      work%pressure_reference = 10.0_dp**work%log10P
      if (.not.all(ieee_is_finite(work%pressure_reference)) .or. &
          any(work%pressure_reference <= 0.0_dp)) then
        err = 'The pressure-profile mapping produced invalid pressure values.'
        return
      endif

      if (temperature_change <= persistent_tolerance) then
        converged = .true.
        exit
      endif
      if (iteration == persistent_max_iterations) exit

      ! Rebuild the state used by the next mapping pass from the updated
      ! pressure-profile temperature.
      call extend_candidate_density(self, state, work, err)
      if (allocated(err)) return
      call assemble_candidate_atmosphere(self, state, work, err)
      if (allocated(err)) return
    enddo

    if (.not.converged) then
      err = 'The persistent pressure profile did not converge on the candidate vertical grid.'
      return
    endif

    ! Rebuild the final candidate state using the converged temperature.
    call extend_candidate_density(self, state, work, err)
    if (allocated(err)) return
    call map_candidate_particle_radii(self, state, err)
    if (allocated(err)) return
    call assemble_candidate_atmosphere(self, state, work, err)

  end subroutine

  subroutine VerticalGridWork_allocate(work, nq, nz)
    type(VerticalGridWork), intent(inout) :: work
    integer, intent(in) :: nq, nz

    allocate(work%pressure(nz))
    allocate(work%mix(nq, nz))
    allocate(work%mix_new(nq, nz))
    allocate(work%density(nz))
    allocate(work%density_new(nz))
    allocate(work%temperature_reference(nz))
    allocate(work%pressure_reference(nz))
    allocate(work%source_pressure(nz))
    allocate(work%log10P(nz))

  end subroutine

  module subroutine update_vertical_grid(self, TOA_alt, TOA_pressure, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), optional, intent(in) :: TOA_alt !! cm
    real(dp), optional, intent(in) :: TOA_pressure !! Target top-of-atmosphere pressure (dyn/cm^2).
    character(:), allocatable, intent(out) :: err

    real(dp) :: top_atmos_new
    type(AtmosphereState) :: previous_state, state
    type(VerticalGridWork) :: work

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    character(:), allocatable :: original_err

    dat => self%dat
    var => self%var
    wrk => self%wrk

    ! Check inputs
    call self%require_atmosphere_initialized('update_vertical_grid', err)
    if (allocated(err)) return

    if (present(TOA_alt) .and. present(TOA_pressure)) then
      err = 'Both "TOA_alt" and "TOA_pressure" cannot be specified.'
      return
    endif
    if (.not.present(TOA_alt) .and. .not.present(TOA_pressure)) then
      err = 'Either "TOA_alt" and "TOA_pressure" must be specified'
      return
    endif

    if (present(TOA_alt)) then
      if (.not.ieee_is_finite(TOA_alt) .or. TOA_alt <= var%bottom_atmos) then
        err = '"TOA_alt" must be finite and greater than the model bottom.'
        return
      endif
    endif

    if (present(TOA_pressure)) then
      if (.not.ieee_is_finite(TOA_pressure) .or. TOA_pressure <= 0.0_dp) then
        err = '"TOA_pressure" must be finite and positive.'
        return
      endif
    endif

    call previous_state%allocate(dat, var%nz)
    call state%allocate(dat, var%nz)
    call VerticalGridWork_allocate(work, dat%nq, var%nz)

    ! Snapshot the committed state before refreshing any derived work arrays.
    ! The hydrostatic pressure used by candidate construction must be derived
    ! from the canonical atmospheric state, not a possibly stale work cache.
    call copy_model_to_state(self, previous_state)

    ! Refresh all atmospheric work state from the committed composition. When
    ! enabled, the persistent profile is applied before pressure and
    ! hydrostatic quantities are recomputed.
    call self%prepare_atmosphere_structure(previous_state%usol, wrk%usol, &
                               wrk%molecules_per_particle, wrk%pressure, &
                               wrk%density, wrk%mix, wrk%mubar, &
                               wrk%pressure_hydro, wrk%density_hydro, &
                               var%press_temp_edd_profile%enabled, err)
    if (allocated(err)) then
      original_err = err
      call restore_previous_state(original_err, err)
      return
    endif

    ! Candidate construction treats the refreshed live model as a read-only
    ! source. Only the local state and work objects are modified.
    call copy_model_to_state(self, state)
    work%source_pressure = wrk%pressure_hydro

    if (present(TOA_alt)) then
      top_atmos_new = TOA_alt
    else
      top_atmos_new = TOA_at_pressure(self, TOA_pressure, state, work, err)
      if (allocated(err)) then
        original_err = err
        call restore_previous_state(original_err, err)
        return
      endif
    endif

    ! Compute properties associated with new TOA
    call build_vertical_grid_state(self, top_atmos_new, state, work, err)
    if (allocated(err)) then
      original_err = err
      call restore_previous_state(original_err, err)
      return
    endif

    ! Finalize all derived state before touching the live model.
    call finalize_atmosphere_state(dat, state, err)
    if (allocated(err)) then
      original_err = err
      call restore_previous_state(original_err, err)
      return
    endif

    ! Commit and prepare while the old stepper is still intact. If candidate
    ! preparation fails, restore both model state and atmospheric work arrays.
    call copy_state_to_model(self, state)
    call self%prep_atmosphere(state%usol, err)
    if (allocated(err)) then
      original_err = err
      call restore_previous_state(original_err, err)
      return
    endif

    ! A successful grid change invalidates the old CVODE infrastructure.
    call self%destroy_stepper(err)
    if (allocated(err)) then
      original_err = err
      call restore_previous_state(original_err, err)
      return
    endif

  contains

    subroutine restore_previous_state(original_err_, err_)
      character(*), intent(in) :: original_err_
      character(:), allocatable, intent(out) :: err_
      character(:), allocatable :: rollback_err

      call copy_state_to_model(self, previous_state)
      call self%prep_atmosphere(previous_state%usol, rollback_err)
      if (allocated(rollback_err)) then
        err_ = original_err_//' Rollback failed: '//rollback_err
      else
        err_ = original_err_
      endif

    end subroutine

  end subroutine

end submodule
