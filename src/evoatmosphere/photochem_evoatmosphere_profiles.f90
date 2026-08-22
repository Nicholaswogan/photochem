submodule(photochem_evoatmosphere) photochem_evoatmosphere_profiles
  use photochem_vars, only: PressureTempEddProfile
  implicit none

  type :: PressTempEddState
    real(dp) :: trop_alt
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: edd(:)
    integer :: trop_ind ! derived
    real(dp), allocatable :: xs_x_qy(:,:,:) ! derived
    real(dp), allocatable :: gibbs_energy(:,:) ! derived
  end type

contains

  function initialize_PressTempEddState(self) result(state)
    class(EvoAtmosphere), intent(in) :: self
    type(PressTempEddState) :: state

    state%trop_alt = self%var%trop_alt
    state%temperature = self%var%temperature
    state%edd = self%var%edd

    state%trop_ind = self%var%trop_ind
    state%xs_x_qy = self%var%xs_x_qy
    if (self%dat%reverse) state%gibbs_energy = self%var%gibbs_energy

  end function

  subroutine copy_PressTempEddState_to_model(self, state)
    class(EvoAtmosphere), intent(inout) :: self
    type(PressTempEddState), intent(in) :: state

    self%var%trop_alt = state%trop_alt
    self%var%temperature = state%temperature
    self%var%edd = state%edd

    self%var%trop_ind = state%trop_ind
    self%var%xs_x_qy = state%xs_x_qy
    if (self%dat%reverse) self%var%gibbs_energy = state%gibbs_energy

  end subroutine

  module subroutine set_temperature(self, temperature, trop_alt, err)
    use iso_c_binding, only: c_associated
    use photochem_vars, only: refresh_temperature_dependent_vars

    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: temperature(:)
    real(dp), optional, intent(in) :: trop_alt
    character(:), allocatable, intent(out) :: err

    type :: TemperatureState
      real(dp) :: trop_alt
      real(dp), allocatable :: temperature(:)
      integer :: trop_ind ! derived
      real(dp), allocatable :: xs_x_qy(:,:,:) ! derived
      real(dp), allocatable :: gibbs_energy(:,:) ! derived
    end type

    type(TemperatureState) :: state, previous_state
    type(PhotochemVars), pointer :: var
    real(dp), allocatable :: usol_start(:,:)
    character(:), allocatable :: original_err, rollback_err

    var => self%var

    call self%require_atmosphere_initialized('set_temperature', err)
    if (allocated(err)) return
    if (c_associated(self%wrk%sun%cvode_mem)) then
      err = "Can not change the temperature profile while a CVODE stepper "// &
            "is initialized. Call 'destroy_stepper' first."
      return
    endif
    if (var%press_temp_edd_profile%enabled) then
      err = "The persistent pressure-temperature-eddy profile is enabled. "// &
            "Call 'clear_press_temp_edd_profile' before 'set_temperature'."
      return
    endif
    if (size(temperature) /= var%nz) then
      err = "temperature has the wrong input dimension"
      return
    endif

    ! Build the candidate from the current state so fields not changed by
    ! this operation retain their existing values.
    previous_state = initialize_TemperatureState(self)
    state = previous_state

    ! Commit the new temperature
    state%temperature = temperature

    ! Update temperature-dependent state
    call refresh_temperature_dependent_vars( &
      self%dat, state%temperature, var%z, var%bottom_atmos, var%top_atmos, &
      trop_alt, state%xs_x_qy, state%gibbs_energy, state%trop_alt, &
      state%trop_ind, err &
    )
    if (allocated(err)) return

    ! Commit only after the candidate and all derived values are valid.
    call copy_state_to_model(self, state)

    ! Fill wrk with new values. Keep the input separate from wrk%usol because
    ! prep_atmosphere has distinct input and output arrays.
    usol_start = self%wrk%usol
    call self%prep_atmosphere(usol_start, err)
    if (allocated(err)) then
      ! Restore both the configuration and derived work state if preparation
      ! failed after making partial changes.
      original_err = err
      call copy_state_to_model(self, previous_state)
      call self%prep_atmosphere(usol_start, rollback_err)
      if (allocated(rollback_err)) then
        err = original_err//' Rollback failed: '//rollback_err
      else
        err = original_err
      endif
      return
    endif

  contains

    function initialize_TemperatureState(self_) result(state_)
      class(EvoAtmosphere), intent(in) :: self_
      type(TemperatureState) :: state_

      state_%trop_alt = self_%var%trop_alt
      state_%temperature = self_%var%temperature

      state_%trop_ind = self_%var%trop_ind
      state_%xs_x_qy = self_%var%xs_x_qy
      if (self_%dat%reverse) state_%gibbs_energy = self_%var%gibbs_energy

    end function

    subroutine copy_state_to_model(self_, state_)
      class(EvoAtmosphere), intent(inout) :: self_
      type(TemperatureState), intent(in) :: state_

      self_%var%trop_alt = state_%trop_alt
      self_%var%temperature = state_%temperature

      self_%var%trop_ind = state_%trop_ind
      self_%var%xs_x_qy = state_%xs_x_qy
      if (self_%dat%reverse) self_%var%gibbs_energy = state_%gibbs_energy

    end subroutine

  end subroutine

  module subroutine set_press_temp_edd(self, P, T, edd, trop_p, hydro_pressure, err)
    use iso_c_binding, only: c_associated
    use photochem_vars, only: refresh_temperature_dependent_vars
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: edd(:)
    real(dp), optional, intent(in) :: trop_p
    logical, optional, intent(in) :: hydro_pressure
    character(:), allocatable, intent(out) :: err

    type(PressTempEddState) :: state, previous_state
    real(dp), allocatable :: T_grid(:), edd_grid(:), log10P_grid(:)
    real(dp), allocatable :: usol_start(:,:)
    real(dp) :: trop_alt
    real(dp) :: trop_p_
    logical :: has_trop_p, hydro_pressure_
    character(:), allocatable :: original_err, rollback_err

    call self%require_atmosphere_initialized('set_press_temp_edd', err)
    if (allocated(err)) return

    if (c_associated(self%wrk%sun%cvode_mem)) then
      err = "Can not set a pressure-temperature-eddy profile while a CVODE "// &
            "stepper is initialized. Call 'destroy_stepper' first."
      return
    endif

    if (self%var%press_temp_edd_profile%enabled) then
      err = "The persistent pressure-temperature-eddy profile is enabled. "// &
            "Call 'clear_press_temp_edd_profile' before 'set_press_temp_edd'."
      return
    endif

    trop_p_ = -1.0_dp
    if (present(trop_p)) trop_p_ = trop_p
    has_trop_p = trop_p_ > 0.0_dp
    if (present(hydro_pressure)) then
      hydro_pressure_ = hydro_pressure
    else
      hydro_pressure_ = .true.
    endif

    ! Allocate only after the lifecycle checks above. In particular, an
    ! uninitialized object must return the normal API error without using an
    ! undefined `var%nz` as an automatic-array extent.
    allocate(T_grid(self%var%nz), edd_grid(self%var%nz), log10P_grid(self%var%nz))

    ! First compute the mapping without changing model state. This kernel is
    ! also suitable for applying a persistent pressure-based profile to an
    ! arbitrary trial composition during a future RHS evaluation.
    call map_press_temp_edd( &
      self, &
      self%wrk%usol, &
      P, &
      T, &
      edd, &
      trop_p=trop_p_, &
      hydro_pressure=hydro_pressure_, &
      grid_z=self%var%z, &
      grid_dz=self%var%dz, &
      grid_grav=self%var%grav, &
      temperature_reference=self%var%temperature, &
      pressure_reference=self%wrk%pressure_hydro, &
      T_grid=T_grid, &
      edd_grid=edd_grid, &
      log10P_grid=log10P_grid, &
      trop_alt=trop_alt, &
      err=err &
    )
    if (allocated(err)) return

    ! Build and validate the complete candidate without changing the model.
    previous_state = initialize_PressTempEddState(self)
    state = previous_state
    state%temperature = T_grid
    state%edd = edd_grid
    if (.not. has_trop_p) trop_alt = state%trop_alt

    call refresh_temperature_dependent_vars( &
      self%dat, state%temperature, self%var%z, self%var%bottom_atmos, &
      self%var%top_atmos, trop_alt, state%xs_x_qy, state%gibbs_energy, &
      state%trop_alt, state%trop_ind, err &
    )
    if (allocated(err)) return

    usol_start = self%wrk%usol
    call copy_PressTempEddState_to_model(self, state)
    call self%prep_atmosphere(usol_start, err)
    if (allocated(err)) then
      original_err = err
      call copy_PressTempEddState_to_model(self, previous_state)
      call self%prep_atmosphere(usol_start, rollback_err)
      if (allocated(rollback_err)) then
        err = original_err//' Rollback failed: '//rollback_err
      else
        err = original_err
      endif
      return
    endif

  end subroutine

  module subroutine set_press_temp_edd_profile(self, P, T, edd, trop_p, hydro_pressure, &
                                               maintain_toa_pressure, target_pressure, err)
    use iso_c_binding, only: c_associated
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: edd(:)
    real(dp), optional, intent(in) :: trop_p
    logical, optional, intent(in) :: hydro_pressure
    logical, optional, intent(in) :: maintain_toa_pressure
    real(dp), optional, intent(in) :: target_pressure
    character(:), allocatable, intent(out) :: err

    type(PressTempEddState) :: previous_state
    type(PressureTempEddProfile) :: previous_profile
    real(dp), allocatable :: usol_start(:,:)
    logical :: maintain_toa_pressure_
    character(:), allocatable :: original_err, rollback_err

    call self%require_atmosphere_initialized('set_press_temp_edd_profile', err)
    if (allocated(err)) return

    if (c_associated(self%wrk%sun%cvode_mem)) then
      err = "Can not set a persistent pressure-temperature-eddy profile while "// &
            "a stepper is initialized. Call 'destroy_stepper' first."
      return
    endif

    maintain_toa_pressure_ = .true.
    if (present(maintain_toa_pressure)) maintain_toa_pressure_ = maintain_toa_pressure
    if (present(target_pressure)) then
      if (.not.ieee_is_finite(target_pressure) .or. target_pressure <= 0.0_dp) then
        err = '`target_pressure` must be finite and positive.'
        return
      endif
    endif

    previous_state = initialize_PressTempEddState(self)
    previous_profile = self%var%press_temp_edd_profile

    self%var%press_temp_edd_profile%pressure = P
    self%var%press_temp_edd_profile%temperature = T
    self%var%press_temp_edd_profile%edd = edd
    if (present(trop_p)) then
      self%var%press_temp_edd_profile%trop_p = trop_p
    else
      self%var%press_temp_edd_profile%trop_p = -1.0_dp
    endif
    if (present(hydro_pressure)) then
      self%var%press_temp_edd_profile%hydro_pressure = hydro_pressure
    else
      self%var%press_temp_edd_profile%hydro_pressure = .true.
    endif
    self%var%press_temp_edd_profile%enabled = .true.

    ! prep_atmosphere validates and applies the stored profile through the
    ! same path used by RHS and Jacobian evaluations.
    ! Keep the input separate from wrk%usol because preparation writes the
    ! canonical working state through a distinct output argument.
    usol_start = self%wrk%usol
    call self%prep_atmosphere(usol_start, err)
    if (allocated(err)) then
      original_err = err
      self%var%press_temp_edd_profile = previous_profile
      call copy_PressTempEddState_to_model(self, previous_state)
      call self%prep_atmosphere(usol_start, rollback_err)
      if (allocated(rollback_err)) then
        err = original_err//' Rollback failed: '//rollback_err
      else
        err = original_err
      endif
      return
    endif

    ! The profile is now committed. Couple TOA maintenance to the same
    ! successful profile installation, using the requested target or the
    ! standard default when the caller did not provide one.
    self%var%toa_pressure_maintenance%enabled = maintain_toa_pressure_
    if (maintain_toa_pressure_) then
      if (present(target_pressure)) then
        self%var%toa_pressure_maintenance%target_pressure = target_pressure
      else
        self%var%toa_pressure_maintenance%target_pressure = 0.1_dp
      endif
    endif

  end subroutine

  module subroutine clear_press_temp_edd_profile(self, err)
    use iso_c_binding, only: c_associated
    class(EvoAtmosphere), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    if (c_associated(self%wrk%sun%cvode_mem)) then
      err = "Can not clear the persistent pressure-temperature-eddy profile while "// &
            "a stepper is initialized. Call 'destroy_stepper' first."
      return
    endif

    call reset_press_temp_edd_profile(self%var)

  end subroutine

  module subroutine reset_press_temp_edd_profile(var)
    type(PhotochemVars), intent(inout) :: var

    var%press_temp_edd_profile%enabled = .false.
    var%press_temp_edd_profile%hydro_pressure = .true.
    var%press_temp_edd_profile%trop_p = -1.0_dp
    ! A cleared persistent profile cannot support automatic TOA maintenance.
    ! Preserve the maintenance tuning values, but clear its active state and
    ! any target that was coupled to the previous profile.
    var%toa_pressure_maintenance%enabled = .false.
    var%toa_pressure_maintenance%target_pressure = 0.0_dp
    if (allocated(var%press_temp_edd_profile%pressure)) &
      deallocate(var%press_temp_edd_profile%pressure)
    if (allocated(var%press_temp_edd_profile%temperature)) &
      deallocate(var%press_temp_edd_profile%temperature)
    if (allocated(var%press_temp_edd_profile%edd)) &
      deallocate(var%press_temp_edd_profile%edd)

  end subroutine

  module subroutine apply_press_temp_edd_profile(self, usol_in, err)
    use photochem_vars, only: refresh_temperature_dependent_vars
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp) :: T_grid(self%var%nz), edd_grid(self%var%nz)
    real(dp) :: log10P_grid(self%var%nz), trop_alt

    if (.not. self%var%press_temp_edd_profile%enabled) return

    call map_press_temp_edd( &
      self, &
      usol_in, &
      self%var%press_temp_edd_profile%pressure, &
      self%var%press_temp_edd_profile%temperature, &
      self%var%press_temp_edd_profile%edd, &
      trop_p=self%var%press_temp_edd_profile%trop_p, &
      hydro_pressure=self%var%press_temp_edd_profile%hydro_pressure, &
      grid_z=self%var%z, &
      grid_dz=self%var%dz, &
      grid_grav=self%var%grav, &
      temperature_reference=self%var%temperature, &
      pressure_reference=self%wrk%pressure_hydro, &
      T_grid=T_grid, &
      edd_grid=edd_grid, &
      log10P_grid=log10P_grid, &
      trop_alt=trop_alt, &
      err=err &
    )
    if (allocated(err)) return

    ! Commit only after the profile has mapped successfully and the
    ! tropopause has been validated. Do not call set_temperature here: that
    ! routine calls prep_atmosphere and would recurse back into this routine.
    self%var%temperature = T_grid
    self%var%edd = edd_grid
    call refresh_temperature_dependent_vars( &
      self%dat, self%var%temperature, self%var%z, self%var%bottom_atmos, &
      self%var%top_atmos, trop_alt, self%var%xs_x_qy, &
      self%var%gibbs_energy, self%var%trop_alt, self%var%trop_ind, err &
    )

  end subroutine

  module subroutine map_press_temp_edd(self, usol, P, T, edd, trop_p, hydro_pressure, &
                                       grid_z, grid_dz, grid_grav, temperature_reference, &
                                       pressure_reference, &
                                       T_grid, edd_grid, log10P_grid, trop_alt, err)
    use futils, only: interp, brent_class
    use ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use photochem_const, only: small_real, k_boltz, N_avo
    use photochem_enum, only: DensityBC, PressureBC
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: edd(:)
    real(dp), intent(in) :: trop_p
    logical, intent(in) :: hydro_pressure
    real(dp), intent(in) :: grid_z(:), grid_dz(:), grid_grav(:)
    real(dp), intent(in) :: temperature_reference(:), pressure_reference(:)
    real(dp), intent(out) :: T_grid(:)
    real(dp), intent(out) :: edd_grid(:)
    real(dp), intent(out) :: log10P_grid(:)
    real(dp), intent(out) :: trop_alt
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: log10P_in(:), T_in(:), log10edd_in(:)
    real(dp), allocatable :: P_wrk(:), log10edd_new(:)
    real(dp), allocatable :: usol_base(:,:), density_base(:), mubar_base(:)
    real(dp), allocatable :: usol_bottom(:)
    real(dp) :: xzero, Psurf_initial, Psurf_final
    real(dp) :: log10P_previous, temperature_previous
    real(dp) :: trop_alt_array(1)
    logical :: has_trop_p
    integer :: ierr, i, j, residual_layer
    character(32) :: layer_string
    character(:), allocatable :: residual_error

    real(dp), parameter :: log10e = log10(exp(1.0_dp))
    real(dp), parameter :: root_tol = 1.0e-10_dp
    integer, parameter :: max_bracket_expansions = 60

    type(brent_class) :: root_solver

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! Check inputs
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = '"usol" has the wrong dimensions'
      return
    endif
    if (size(T_grid) /= var%nz .or. size(edd_grid) /= var%nz .or. &
        size(log10P_grid) /= var%nz) then
      err = 'Pressure-profile mapping output arrays have the wrong dimensions'
      return
    endif
    if (size(grid_z) /= var%nz .or. size(grid_dz) /= var%nz .or. &
        size(grid_grav) /= var%nz .or. &
        size(temperature_reference) /= var%nz .or. &
        size(pressure_reference) /= var%nz) then
      err = 'Pressure-profile mapping grid arrays have the wrong dimensions'
      return
    endif
    if (.not.all(ieee_is_finite(grid_z)) .or. &
        .not.all(ieee_is_finite(grid_dz)) .or. any(grid_dz <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(grid_grav)) .or. any(grid_grav <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(temperature_reference)) .or. &
        any(temperature_reference <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(pressure_reference)) .or. &
        any(pressure_reference <= 0.0_dp)) then
      err = 'Pressure-profile mapping requires finite, positive grid reference state'
      return
    endif
    if (size(P) /= size(T)) then
      err = '"P" and "T" not the same size'
      return
    endif
    if (size(P) /= size(edd)) then
      err = '"P" and "edd" not the same size'
      return
    endif
    if (size(P) < 2) then
      err = 'size(P) must be >= 2'
      return
    endif
    if (.not. all(ieee_is_finite(P)) .or. .not. all(ieee_is_finite(T)) .or. &
        .not. all(ieee_is_finite(edd))) then
      err = 'All elements of "P", "T" and "edd" must be finite'
      return
    endif
    if (any(P <= 0.0_dp) .or. any(T <= 0.0_dp) .or. any(edd <= 0.0_dp)) then
      err = 'All elements of "P", "T" and "edd" must be positive'
      return
    endif
    if (.not. all(P(2:) < P(:size(P)-1))) then
      err = '"P" must be strictly decreasing'
      return
    endif
    if (.not. ieee_is_finite(trop_p)) then
      err = '"trop_p" must be finite; use a non-positive value to omit it'
      return
    endif
    has_trop_p = trop_p > 0.0_dp
    if (has_trop_p .and. .not.dat%gas_rainout) then
      err = '"trop_p" can only be supplied when gas rainout is enabled.'
      return
    endif
    if (dat%gas_rainout .and. .not.has_trop_p) then
      err = '"trop_p" is a required input.'
      return
    endif

    ! Allocate work space
    allocate(P_wrk(var%nz), log10edd_new(var%nz))
    allocate(usol_base(dat%nq,var%nz), density_base(var%nz), mubar_base(var%nz))
    allocate(usol_bottom(dat%nq))

    ! Copy and clip the current state in the same way as prep_atm_evo_gas.
    call clip_usol(usol, usol_base)

    ! Above the bottom layer, density and mean molecular weight do not
    ! depend on temperature. The bottom layer is handled separately because
    ! fixed-pressure boundary conditions depend on temperature.
    do j = 2,var%nz
      density_base(j) = sum(usol_base(dat%ng_1:,j))
      if (.not. ieee_is_finite(density_base(j)) .or. density_base(j) <= 0.0_dp) then
        write(layer_string,'(i0)') j
        err = 'The gas density is not finite and positive in layer '//trim(layer_string)
        return
      endif
      mubar_base(j) = sum(dat%species_mass(dat%ng_1:dat%nq)* &
                          usol_base(dat%ng_1:dat%nq,j))/density_base(j)
      if (.not. ieee_is_finite(mubar_base(j)) .or. mubar_base(j) <= 0.0_dp) then
        write(layer_string,'(i0)') j
        err = 'The mean molecular weight is not finite and positive in layer '//trim(layer_string)
        return
      endif
    enddo

    call bottom_column_state(temperature_reference(1), density_base(1), &
                             mubar_base(1), Psurf_initial, err)
    if (allocated(err)) return

    ! Copy the input profiles into log-pressure space. The arrays remain at
    ! the supplied pressure points; values below the deepest point are
    ! extrapolated directly from the deepest two points by evaluate_profile.
    allocate(log10P_in(size(P)),T_in(size(P)),log10edd_in(size(P)))
    log10P_in = log10(P)
    T_in = T
    log10edd_in = log10(edd)
    ! Flip order for interpolation purposes
    log10P_in = log10P_in(size(log10P_in):1:-1)
    T_in = T_in(size(log10P_in):1:-1)
    log10edd_in = log10edd_in(size(log10P_in):1:-1)

    call root_solver%set_function(pressure_residual)

    ! Solve one scalar pressure equation per layer. In hydrostatic mode the
    ! recurrence is triangular, so each solved layer provides the lower
    ! pressure and temperature boundary for the next layer.
    do i = 1,var%nz

      ! Perform the solve
      call solve_layer_pressure(i, xzero, err)
      if (allocated(err)) return

      ! Save the result, and evaluate T and Kzz at the root.
      log10P_grid(i) = xzero
      P_wrk(i) = 10.0_dp**xzero
      call evaluate_profile(xzero, log10P_in, T_in, T_grid(i), 'temperature', .true., err)
      if (allocated(err)) return

      call evaluate_profile(xzero, log10P_in, log10edd_in, log10edd_new(i), &
                            'eddy-diffusion', .false., err)
      if (allocated(err)) return

      ! Check for problems
      if (.not. ieee_is_finite(P_wrk(i)) .or. P_wrk(i) <= 0.0_dp) then
        write(layer_string,'(i0)') i
        err = 'The pressure solve produced a non-finite or non-positive pressure in layer '// &
              trim(layer_string)
        return
      endif
      if (i > 1) then
        if (P_wrk(i) >= P_wrk(i-1)) then
          if (hydro_pressure) then
            write(layer_string,'(i0)') i
            err = 'The hydrostatic pressure does not decrease upward at layer '//trim(layer_string)
            return
          elseif (has_trop_p) then
            write(layer_string,'(i0," and ",i0)') i-1, i
            err = 'Cannot determine the tropopause altitude because actual pressure is not '// &
                  'strictly decreasing between layers '//trim(layer_string)
            return
          endif
        endif
      endif

      if (i == 1) then; block
        real(dp) :: density_final, mubar_final

        ! Little sanity check to make sure pressure is decreasing.

        call bottom_column_state(T_grid(1), density_final, mubar_final, Psurf_final, err)
        if (allocated(err)) return
        if (hydro_pressure .and. P_wrk(1) >= Psurf_final) then
          err = 'The solved bottom-layer pressure is not below the surface pressure'
          return
        endif
      end block; endif

      ! Save the previous solve. This is needed for initial guesses and for the residuals on
      ! subsequent layers.
      log10P_previous = xzero
      temperature_previous = T_grid(i)
    enddo

    ! Ensure eddy diffusion is valid
    if (.not. all(ieee_is_finite(log10edd_new))) then
      err = 'The mapped eddy-diffusion profile is not finite'
      return
    endif
    edd_grid = 10.0_dp**log10edd_new
    if (.not. all(ieee_is_finite(edd_grid)) .or. any(edd_grid <= 0.0_dp)) then
      err = 'The mapped eddy-diffusion profile is not finite and positive'
      return
    endif

    ! Compute the tropopause altitude
    trop_alt = 0.0_dp
    if (has_trop_p) then
      call interp([log10(trop_p)], log10P_grid(var%nz:1:-1), &
                  grid_z(var%nz:1:-1), trop_alt_array, ierr=ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        return
      endif
      trop_alt = trop_alt_array(1)
    endif

  contains
    subroutine solve_layer_pressure(layer, xzero_, err_)
      integer, intent(in) :: layer
      real(dp), intent(out) :: xzero_
      character(:), allocatable, intent(out) :: err_

      real(dp) :: fzero, xcenter, xlower, xupper, flower, fupper, width
      real(dp) :: reference_density
      integer :: iflag, nexpand

      ! pressure_residual has the callback signature required by brent_class.
      ! residual_layer selects the layer; all other callback context is
      ! read-only and consists of the input profiles, atmospheric state, and
      ! the previously solved pressure and temperature in hydrostatic mode.
      residual_layer = layer
      if (allocated(residual_error)) deallocate(residual_error)

      ! `xcenter` is the root estimate, while `width` is how far above and below
      ! `xcenter` we will search for a bracket that we can use for the solve.
      if (hydro_pressure) then
        if (layer == 1) then
          ! Existing bottom pressure is the root estimate.
          xcenter = log10(max(pressure_reference(1), small_real))
          width = max(abs(log10(Psurf_initial) - xcenter), 1.0e-3_dp)
        else
          ! Subtract the expected hydrostatic drop from the previous pressure;
          ! use that drop for the initial bracket.
          xcenter = log10P_previous - &
                    (mubar_base(layer)*grid_grav(layer)*grid_dz(layer)*log10e)/ &
                    (N_avo*k_boltz*max(temperature_previous,small_real))
          width = max(abs(log10P_previous - xcenter), 1.0e-3_dp)
        endif
      else
        ! Use the ideal-gas pressure from the current density and temperature.
        reference_density = density_base(layer)
        xcenter = log10(reference_density*k_boltz* &
                        max(temperature_reference(layer),small_real))
        ! Start narrow; the bracket expands if the residual has no sign change.
        width = 1.0e-3_dp
      endif

      ! Now find a bracket that we can solve on.
      xlower = xcenter - width
      xupper = xcenter + width
      flower = pressure_residual(root_solver, xlower)
      if (allocated(residual_error)) then
        err_ = residual_error
        return
      endif
      fupper = pressure_residual(root_solver, xupper)
      if (allocated(residual_error)) then
        err_ = residual_error
        return
      endif
      nexpand = 0
      do while (.not. opposite_signs(flower, fupper))
        width = 2.0_dp*width
        xlower = xcenter - width
        xupper = xcenter + width
        flower = pressure_residual(root_solver, xlower)
        if (allocated(residual_error)) then
          err_ = residual_error
          return
        endif
        fupper = pressure_residual(root_solver, xupper)
        if (allocated(residual_error)) then
          err_ = residual_error
          return
        endif
        nexpand = nexpand + 1
        if (nexpand >= max_bracket_expansions) then
          write(layer_string,'(i0)') layer
          err_ = 'Could not bracket the pressure root in layer '//trim(layer_string)
          return
        endif
      enddo

      ! Perform the root solve.
      call root_solver%find_zero(xlower, xupper, root_tol, xzero_, fzero, iflag, flower, fupper)
      if (allocated(residual_error)) then
        err_ = residual_error
        return
      endif
      if (iflag /= 0 .or. .not. ieee_is_finite(xzero_)) then
        write(layer_string,'(i0)') layer
        err_ = 'The scalar pressure solve failed in layer '//trim(layer_string)
        return
      endif

    end subroutine

    function pressure_residual(me, x) result(residual)
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: x
      real(dp) :: residual

      real(dp) :: temperature_trial, density_trial, mubar_trial, Psurf_trial

      ! Preserve the first invalid evaluation while Brent is unwinding. This
      ! prevents a later callback from replacing the useful diagnostic.
      if (allocated(residual_error)) then
        residual = ieee_value(0.0_dp, ieee_quiet_nan)
        return
      endif

      ! For out input log10 pressure `x`, compute the `temperature_trial`
      call evaluate_profile(x, log10P_in, T_in, temperature_trial, 'temperature', .true., residual_error)
      if (allocated(residual_error)) then
        residual = ieee_value(0.0_dp, ieee_quiet_nan)
        return
      endif

      if (hydro_pressure) then
        ! x is trial log10 pressure. Since the input temperature is T(P),
        ! each trial pressure determines the temperature used in the residual.
        if (residual_layer == 1) then
          ! Recompute the surface pressure because lower boundary conditions
          ! can depend on temperature. Enforce balance to the layer center,
          ! hence the half-layer thickness.
          call bottom_column_state(temperature_trial, density_trial, mubar_trial, &
                                   Psurf_trial, residual_error)
          if (allocated(residual_error)) then
            residual = ieee_value(0.0_dp, ieee_quiet_nan)
            return
          endif
          residual = x - log10(Psurf_trial) + &
                     (mubar_trial*grid_grav(1)*0.5_dp*grid_dz(1)*log10e)/ &
                     (N_avo*k_boltz*temperature_trial)
        else
          ! For higher layers, enforce the hydrostatic pressure drop from the
          ! previous layer using the average of the two temperatures.
          residual = x - log10P_previous + &
                     (mubar_base(residual_layer)*grid_grav(residual_layer)* &
                      grid_dz(residual_layer)*log10e)/ &
                     (N_avo*k_boltz*0.5_dp*(temperature_previous + temperature_trial))
        endif
      else
        ! Without hydrostatic coupling, enforce the local ideal-gas relation
        ! P = density*k_B*T at the trial pressure.
        if (residual_layer == 1) then
          call bottom_column_state(temperature_trial, density_trial, mubar_trial, &
                                   Psurf_trial, residual_error)
          if (allocated(residual_error)) then
            residual = ieee_value(0.0_dp, ieee_quiet_nan)
            return
          endif
        else
          density_trial = density_base(residual_layer)
        endif
        residual = x - log10(density_trial*k_boltz*temperature_trial)
      endif

      if (.not. ieee_is_finite(residual)) then
        residual_error = 'The pressure residual became non-finite during the pressure solve'
        residual = ieee_value(0.0_dp, ieee_quiet_nan)
      endif

    end function

    subroutine bottom_column_state(temperature, density_bottom, mubar_bottom, surface_pressure, err_)
      real(dp), intent(in) :: temperature
      real(dp), intent(out) :: density_bottom, mubar_bottom, surface_pressure
      character(:), allocatable, intent(out) :: err_

      real(dp) :: column_mass
      integer :: gas_ind

      ! Apply boundary conditions to a local bottom-layer copy; self is never
      ! modified by this calculation.
      usol_bottom = usol_base(:,1)
      call self%apply_lower_boundary_conditions(temperature, usol_bottom, err_)
      if (allocated(err_)) return

      density_bottom = sum(usol_bottom(dat%ng_1:))
      if (.not. ieee_is_finite(density_bottom) .or. density_bottom <= 0.0_dp) then
        err_ = 'The bottom-layer number density is not finite and positive'
        return
      endif
      mubar_bottom = sum(dat%species_mass(dat%ng_1:dat%nq)* &
                         usol_bottom(dat%ng_1:dat%nq))/density_bottom
      if (.not. ieee_is_finite(mubar_bottom) .or. mubar_bottom <= 0.0_dp) then
        err_ = 'The bottom-layer mean molecular weight is not finite and positive'
        return
      endif

      column_mass = density_bottom*mubar_bottom*grid_grav(1)*grid_dz(1)
      do gas_ind = 2,var%nz
        column_mass = column_mass + density_base(gas_ind)*mubar_base(gas_ind)* &
                                    grid_grav(gas_ind)*grid_dz(gas_ind)
      enddo
      surface_pressure = column_mass/N_avo
      if (.not. ieee_is_finite(surface_pressure) .or. surface_pressure <= 0.0_dp) then
        err_ = 'The bottom-layer surface pressure is not finite and positive'
        return
      endif

    end subroutine

    subroutine evaluate_profile(x, pressure_grid, values, value, profile_name, require_positive, err_)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: pressure_grid(:)
      real(dp), intent(in) :: values(:)
      real(dp), intent(out) :: value
      character(*), intent(in) :: profile_name
      logical, intent(in) :: require_positive
      character(:), allocatable, intent(out) :: err_

      real(dp) :: fraction
      integer :: low, high, mid

      if (x <= pressure_grid(1)) then
        ! Above the supplied top of the atmosphere, retain the top value.
        value = values(1)
      elseif (x >= pressure_grid(size(pressure_grid))) then
        ! Below the supplied deepest point, continue the deepest two-point
        ! slope in log-pressure rather than attaching a moving surface knot.
        fraction = (x-pressure_grid(size(pressure_grid))) / &
                   (pressure_grid(size(pressure_grid))-pressure_grid(size(pressure_grid)-1))
        value = values(size(values)) + fraction * &
                (values(size(values))-values(size(values)-1))
      else
        ! Otherwise, we perform linear interpolation
        low = 1
        high = size(pressure_grid)
        do while (high-low > 1)
          mid = (low+high)/2
          if (x >= pressure_grid(mid)) then
            low = mid
          else
            high = mid
          endif
        enddo
        fraction = (x-pressure_grid(low))/(pressure_grid(high)-pressure_grid(low))
        value = values(low) + fraction*(values(high)-values(low))
      endif

      if (.not. ieee_is_finite(value) .or. (require_positive .and. value <= 0.0_dp)) then
        if (require_positive) then
          err_ = 'Linear extrapolation of the input '//trim(profile_name)// &
                 ' profile produced a non-finite or non-positive value'
        else
          err_ = 'Linear extrapolation of the input '//trim(profile_name)// &
                 ' profile produced a non-finite value'
        endif
      endif

    end subroutine

    pure function opposite_signs(a, b) result(opposite)
      real(dp), intent(in) :: a, b
      logical :: opposite

      opposite = (a <= 0.0_dp .and. b >= 0.0_dp) .or. &
                 (a >= 0.0_dp .and. b <= 0.0_dp)

    end function
  end subroutine

end submodule

