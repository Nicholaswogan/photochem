submodule(photochem_evoatmosphere) photochem_evoatmosphere_init
  implicit none
  
  ! Contains the Constructor for the EvoAtmosphere derived type.
  
contains
  
  module function create_EvoAtmosphere(mechanism_file, settings_file, flux_file, atmosphere_txt, data_dir, err) result(self)
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphere) :: self

    call setup_evoatmosphere_static(self, mechanism_file, settings_file, &
                                    flux_file, data_dir, err)
    if (allocated(err)) return

    call self%initialize_from_atmosphere_file(atmosphere_txt, err)

  end function

  module function create_EvoAtmosphere_static(mechanism_file, settings_file, &
                                              flux_file, data_dir, err) result(self)
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphere) :: self

    call setup_evoatmosphere_static(self, mechanism_file, settings_file, &
                                    flux_file, data_dir, err)

  end function

  subroutine setup_evoatmosphere_static(self, mechanism_file, settings_file, &
                                        flux_file, data_dir, err)
    use photochem_input, only: setup_static
    use photochem_types, only: PhotoSettings

    type(EvoAtmosphere), intent(out) :: self
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err

    type(PhotoSettings) :: s

    s = PhotoSettings(settings_file, err)
    if (allocated(err)) return

    allocate(self%dat)
    allocate(self%var)
    allocate(self%wrk)

    self%var%data_dir = data_dir
    call setup_static(mechanism_file, s, flux_file, self%dat, self%var, err)
    if (allocated(err)) return

    call self%wrk%init(self%dat%nsp, self%dat%np, self%dat%nq, &
                       self%var%nz, self%dat%nrT, self%dat%kj, &
                       self%dat%nw)
    self%atmosphere_initialized = .false.

  end subroutine

  module subroutine initialize_from_atmosphere_file(self, atmosphere_txt, err)
    use photochem_input, only: finalize_atmosphere_initialization, &
                               setup_atmosphere_from_file

    class(EvoAtmosphere), intent(inout) :: self
    character(len=*), intent(in) :: atmosphere_txt
    character(:), allocatable, intent(out) :: err

    type(EvoAtmosphere) :: candidate
    real(dp), allocatable :: usol_start(:,:)

    call create_atmosphere_candidate(self, candidate, err)
    if (allocated(err)) return

    call reset_press_temp_edd_profile(candidate%var)

    allocate(usol_start(candidate%dat%nq,candidate%var%nz))
    call setup_atmosphere_from_file(atmosphere_txt, candidate%dat, candidate%var, &
                                    usol_start, err)
    if (allocated(err)) return

    call finalize_atmosphere_initialization(candidate%dat, candidate%var, err)
    if (allocated(err)) return

    call prepare_atmosphere_candidate(candidate, usol_start, err)
    if (allocated(err)) return

    call commit_atmosphere_candidate(self, candidate, err)

  end subroutine

  module subroutine initialize_atmosphere_z(self, z, temperature, edd, &
                                            surface_pressure, mix, &
                                            particle_radius, err)
    use photochem_input, only: finalize_atmosphere_state, &
                               map_atmosphere_z_to_grid
    use photochem_types, only: AtmosphereState

    class(EvoAtmosphere), intent(inout) :: self
    real(dp), intent(in) :: z(:), temperature(:), edd(:)
    real(dp), intent(in) :: surface_pressure
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    character(:), allocatable, intent(out) :: err

    type(AtmosphereState) :: state, previous_state
    logical :: was_initialized

    ! Allocate some work space
    was_initialized = self%atmosphere_initialized
    call state%allocate(self%dat, self%var%nz)
    call copy_model_to_state(self, state)

    ! If initialized, then save the previous state
    if (was_initialized) then
      call previous_state%allocate(self%dat, self%var%nz)
      call copy_model_to_state(self, previous_state)
    endif

    ! Build most of the atmospheric state from inputs
    call map_atmosphere_z_to_grid(self%dat, self%var%nz, self%var%trop_alt, &
                                  z, temperature, edd, &
                                  surface_pressure, mix, particle_radius, &
                                  state, err)
    if (allocated(err)) return

    ! finalized atmospheric state
    call finalize_atmosphere_state(self%dat, state, err)
    if (allocated(err)) return

    ! Destroy integrator, if relevant
    call self%destroy_stepper(err)
    if (allocated(err)) return

    ! Copy state to model
    call copy_state_to_model(self, state)

    ! Prepare the atmosphere
    call self%prep_atmosphere_unchecked(state%usol, apply_persistent_profile=.false., err=err)
    if (allocated(err)) then
      call restore_previous_state()
      return
    endif

    ! Destroy persistent P-T-Kzz profile if relevant
    call self%clear_press_temp_edd_profile(err)
    if (allocated(err)) return

    ! Atmosphere is now initialized
    self%atmosphere_initialized = .true.

  contains

    subroutine restore_previous_state()
      if (was_initialized) then
        call copy_state_to_model(self, previous_state)
        self%atmosphere_initialized = .true.
      else
        self%atmosphere_initialized = .false.
      endif
    end subroutine

  end subroutine

  subroutine copy_model_to_state(self, state)
    use photochem_types, only: AtmosphereState
    class(EvoAtmosphere), intent(in) :: self
    type(AtmosphereState), intent(inout) :: state

    state%bottom_atmos = self%var%bottom_atmos
    state%top_atmos = self%var%top_atmos
    state%trop_alt = self%var%trop_alt
    state%z = self%var%z
    state%dz = self%var%dz
    state%temperature = self%var%temperature
    state%edd = self%var%edd
    state%particle_radius = self%var%particle_radius
    state%usol = self%wrk%usol

    state%trop_ind = self%var%trop_ind
    state%grav = self%var%grav
    state%xs_x_qy = self%var%xs_x_qy
    state%particle_xs = self%var%particle_xs
    if (self%dat%reverse) state%gibbs_energy = self%var%gibbs_energy
    state%press_temp_edd_profile = self%var%press_temp_edd_profile
    state%toa_pressure_maintenance = self%var%toa_pressure_maintenance

  end subroutine

  subroutine copy_state_to_model(self, state)
    use photochem_types, only: AtmosphereState
    class(EvoAtmosphere), intent(inout) :: self
    type(AtmosphereState), intent(in) :: state

    self%var%bottom_atmos = state%bottom_atmos
    self%var%top_atmos = state%top_atmos
    self%var%trop_alt = state%trop_alt
    self%var%z = state%z
    self%var%dz = state%dz
    self%var%temperature = state%temperature
    self%var%edd = state%edd
    self%var%particle_radius = state%particle_radius
    self%wrk%usol = state%usol

    self%var%trop_ind = state%trop_ind
    self%var%grav = state%grav
    self%var%xs_x_qy = state%xs_x_qy
    self%var%particle_xs = state%particle_xs
    if (self%dat%reverse) self%var%gibbs_energy = state%gibbs_energy
    self%var%press_temp_edd_profile = state%press_temp_edd_profile
    self%var%toa_pressure_maintenance = state%toa_pressure_maintenance

  end subroutine

  module subroutine initialize_atmosphere_p(self, pressure, temperature, edd, &
                                            mix, particle_radius, persistent, &
                                            trop_p, maintain_toa_pressure, &
                                            target_pressure, err)
    use photochem_input, only: finalize_atmosphere_state, &
                               map_atmosphere_p_to_grid
    use photochem_types, only: AtmosphereState

    class(EvoAtmosphere), intent(inout) :: self
    real(dp), intent(in) :: pressure(:), temperature(:), edd(:)
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    logical, optional, intent(in) :: persistent
    real(dp), optional, intent(in) :: trop_p
    logical, optional, intent(in) :: maintain_toa_pressure
    real(dp), optional, intent(in) :: target_pressure
    character(:), allocatable, intent(out) :: err

    type(AtmosphereState) :: state, previous_state
    logical :: was_initialized
    logical :: persistent_, maintain_toa_pressure_

    persistent_ = .false.
    if (present(persistent)) persistent_ = persistent
    maintain_toa_pressure_ = .true.
    if (present(maintain_toa_pressure)) maintain_toa_pressure_ = maintain_toa_pressure
    if (.not. persistent_ .and. &
        (present(target_pressure) .or. present(maintain_toa_pressure))) then
      err = '"maintain_toa_pressure" and "target_pressure" '// &
            'can only be specified when "persistent" is true.'
      return
    endif
    if (persistent_ .and. .not.maintain_toa_pressure_ .and. present(target_pressure)) then
      err = '"target_pressure" cannot be specified when "maintain_toa_pressure" is false.'
      return
    endif

    was_initialized = self%atmosphere_initialized

    call state%allocate(self%dat, self%var%nz)
    call copy_model_to_state(self, state)
    if (was_initialized) then
      call previous_state%allocate(self%dat, self%var%nz)
      call copy_model_to_state(self, previous_state)
    endif

    call map_atmosphere_p_to_grid(self%dat, self%var%nz, self%var%trop_alt, &
                                  pressure, temperature, edd, mix, particle_radius, &
                                  trop_p=trop_p, state=state, err=err)
    if (allocated(err)) return

    call finalize_atmosphere_state(self%dat, state, err)
    if (allocated(err)) return

    call self%destroy_stepper(err)
    if (allocated(err)) return

    call copy_state_to_model(self, state)
    call self%prep_atmosphere_unchecked(state%usol, apply_persistent_profile=.false., err=err)
    if (allocated(err)) then
      call restore_previous_state()
      return
    endif

    self%atmosphere_initialized = .true.
    if (persistent_) then
      ! Optional dummy arguments may be forwarded directly to matching
      ! optional dummies. If either input was omitted, it remains absent in
      ! set_press_temp_edd_profile and its default behavior applies.
      call self%set_press_temp_edd_profile(pressure, temperature, edd, &
                                           trop_p=trop_p, hydro_pressure=.true., &
                                           maintain_toa_pressure=maintain_toa_pressure_, &
                                           target_pressure=target_pressure, err=err)
      if (allocated(err)) then
        call restore_previous_state()
        return
      endif
    else
      call self%clear_press_temp_edd_profile(err)
      if (allocated(err)) then
        call restore_previous_state()
        return
      endif
    endif

    self%atmosphere_initialized = .true.

  contains

    subroutine restore_previous_state()
      if (was_initialized) then
        call copy_state_to_model(self, previous_state)
        self%atmosphere_initialized = .true.
      else
        self%atmosphere_initialized = .false.
      endif
    end subroutine

  end subroutine

  subroutine create_atmosphere_candidate(self, candidate, err)
    class(EvoAtmosphere), intent(in) :: self
    type(EvoAtmosphere), intent(out) :: candidate
    character(:), allocatable, intent(out) :: err

    if (.not. allocated(self%dat) .or. .not. allocated(self%var) .or. &
        .not. allocated(self%wrk)) then
      err = 'EvoAtmosphere static setup is not complete.'
      return
    endif

    allocate(candidate%dat)
    allocate(candidate%var)
    allocate(candidate%wrk)
    candidate%dat = self%dat
    candidate%var = self%var
    call candidate%wrk%init(candidate%dat%nsp, candidate%dat%np, candidate%dat%nq, &
                            candidate%var%nz, candidate%dat%nrT, candidate%dat%kj, &
                            candidate%dat%nw)

  end subroutine

  subroutine prepare_atmosphere_candidate(candidate, usol_start, err)
    type(EvoAtmosphere), intent(inout) :: candidate
    real(dp), intent(in) :: usol_start(:,:)
    character(:), allocatable, intent(out) :: err

    call candidate%prep_atmosphere_unchecked(usol_start, err=err)
    if (allocated(err)) return

    candidate%atmosphere_initialized = .true.

  end subroutine

  subroutine commit_atmosphere_candidate(self, candidate, err)
    class(EvoAtmosphere), intent(inout) :: self
    type(EvoAtmosphere), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    character(:), allocatable :: destroy_err

    ! Do not disturb the current model or its integrator until the replacement
    ! atmosphere has been prepared successfully.
    call self%destroy_stepper(destroy_err)
    if (allocated(destroy_err)) then
      err = destroy_err
      return
    endif

    call move_alloc(candidate%dat, self%dat)
    call move_alloc(candidate%var, self%var)
    call move_alloc(candidate%wrk, self%wrk)
    self%atmosphere_initialized = .true.

  end subroutine

  module subroutine require_atmosphere_initialized(self, operation, err)
    class(EvoAtmosphere), intent(in) :: self
    character(len=*), intent(in) :: operation
    character(:), allocatable, intent(out) :: err

    if (.not. self%atmosphere_initialized) then
      err = 'EvoAtmosphere atmosphere is not initialized. Initialize the atmosphere before calling "'// &
            trim(operation)//'".'
    endif

  end subroutine
  
end submodule
