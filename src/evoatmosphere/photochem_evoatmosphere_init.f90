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
    use photochem_types, only: AtmosphereState, AtmosphereStateDerived

    class(EvoAtmosphere), intent(inout) :: self
    real(dp), intent(in) :: z(:), temperature(:), edd(:)
    real(dp), intent(in) :: surface_pressure
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    character(:), allocatable, intent(out) :: err

    type(AtmosphereState) :: state
    type(AtmosphereStateDerived) :: derived
    type(AtmosphereState) :: previous_state
    type(AtmosphereStateDerived) :: previous_derived
    logical :: was_initialized

    if (.not. allocated(self%dat) .or. .not. allocated(self%var) .or. &
        .not. allocated(self%wrk)) then
      err = 'EvoAtmosphere static setup is not complete.'
      return
    endif

    was_initialized = self%atmosphere_initialized
    if (was_initialized) then
      call state_from_model(self, previous_state, previous_derived)
    endif

    call map_atmosphere_z_to_grid(self%dat, self%var%nz, self%var%trop_alt, &
                                  z, temperature, edd, &
                                  surface_pressure, mix, particle_radius, &
                                  state, derived, err)
    if (allocated(err)) return

    derived%surface_pressure = surface_pressure/1.0e6_dp

    call finalize_atmosphere_state(self%dat, state, derived, err)
    if (allocated(err)) return

    if (was_initialized) then
      call self%destroy_stepper(err)
      if (allocated(err)) return
    endif

    call copy_state_to_model(self, state, derived)
    call self%prep_atmosphere_unchecked(state%usol, &
                                        apply_persistent_profile=.false., err=err)
    if (allocated(err)) then
      call restore_previous_state()
      return
    endif

    call reset_press_temp_edd_profile(self%var)
    self%atmosphere_initialized = .true.

  contains

    subroutine restore_previous_state()
      if (was_initialized) then
        call copy_state_to_model(self, previous_state, previous_derived)
        self%atmosphere_initialized = .true.
      else
        self%atmosphere_initialized = .false.
      endif
    end subroutine

  end subroutine

  subroutine state_from_model(self, state, derived)
    use photochem_types, only: AtmosphereState, AtmosphereStateDerived
    class(EvoAtmosphere), intent(in) :: self
    type(AtmosphereState), intent(inout) :: state
    type(AtmosphereStateDerived), intent(inout) :: derived

    call state%ensure(self%var%nz, self%dat%nq, self%dat%npq)
    call derived%ensure(self%var%nz, self%dat%np, self%dat%ng, &
                        self%dat%kj, self%dat%nw)

    state%bottom_atmos = self%var%bottom_atmos
    state%top_atmos = self%var%top_atmos
    state%trop_alt = self%var%trop_alt
    state%z = self%var%z
    state%dz = self%var%dz
    state%temperature = self%var%temperature
    state%edd = self%var%edd
    state%particle_radius = self%var%particle_radius
    state%usol = self%wrk%usol

    derived%surface_pressure = self%var%surface_pressure
    derived%trop_ind = self%var%trop_ind
    derived%grav = self%var%grav
    derived%pressure = self%wrk%pressure
    derived%density = self%wrk%density
    derived%mubar = self%wrk%mubar
    derived%xs_x_qy = self%var%xs_x_qy
    derived%particle_xs = self%var%particle_xs
    derived%gibbs_energy = self%var%gibbs_energy

  end subroutine

  subroutine copy_state_to_model(self, state, derived)
    use photochem_types, only: AtmosphereState, AtmosphereStateDerived
    class(EvoAtmosphere), intent(inout) :: self
    type(AtmosphereState), intent(in) :: state
    type(AtmosphereStateDerived), intent(in) :: derived

    self%var%bottom_atmos = state%bottom_atmos
    self%var%top_atmos = state%top_atmos
    self%var%trop_alt = state%trop_alt
    self%var%z = state%z
    self%var%dz = state%dz
    self%var%temperature = state%temperature
    self%var%edd = state%edd
    self%var%particle_radius = state%particle_radius
    self%var%surface_pressure = derived%surface_pressure
    self%var%trop_ind = derived%trop_ind
    self%var%grav = derived%grav
    self%var%xs_x_qy = derived%xs_x_qy
    self%var%particle_xs = derived%particle_xs
    self%var%gibbs_energy = derived%gibbs_energy
    self%wrk%usol = state%usol

  end subroutine

  module subroutine initialize_atmosphere_p(self, pressure, temperature, edd, &
                                            mix, particle_radius, persistent, &
                                            trop_p, maintain_toa_pressure, &
                                            target_pressure, err)
    use photochem_input, only: finalize_atmosphere_initialization, &
                               map_atmosphere_p_to_grid

    class(EvoAtmosphere), intent(inout) :: self
    real(dp), intent(in) :: pressure(:), temperature(:), edd(:)
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    logical, optional, intent(in) :: persistent
    real(dp), optional, intent(in) :: trop_p
    logical, optional, intent(in) :: maintain_toa_pressure
    real(dp), optional, intent(in) :: target_pressure
    character(:), allocatable, intent(out) :: err

    type(EvoAtmosphere) :: candidate
    real(dp), allocatable :: pressure_model(:), density(:), mubar(:), usol_start(:,:)
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

    call create_atmosphere_candidate(self, candidate, err)
    if (allocated(err)) return

    call reset_press_temp_edd_profile(candidate%var)

    allocate(pressure_model(candidate%var%nz), density(candidate%var%nz), &
             mubar(candidate%var%nz), usol_start(candidate%dat%nq,candidate%var%nz))
    call map_atmosphere_p_to_grid(candidate%dat, candidate%var, pressure, &
                                  temperature, edd, mix, particle_radius, &
                                  trop_p=trop_p, pressure=pressure_model, &
                                  density=density, mubar=mubar, usol=usol_start, err=err)
    if (allocated(err)) return

    call finalize_atmosphere_initialization(candidate%dat, candidate%var, err)
    if (allocated(err)) return

    call prepare_atmosphere_candidate(candidate, usol_start, err)
    if (allocated(err)) return

    if (persistent_) then
      ! Optional dummy arguments may be forwarded directly to matching
      ! optional dummies. If either input was omitted, it remains absent in
      ! set_press_temp_edd_profile and its default behavior applies.
      call candidate%set_press_temp_edd_profile(pressure, temperature, edd, &
                                                trop_p=trop_p, &
                                                hydro_pressure=.true., &
                                                maintain_toa_pressure=maintain_toa_pressure_, &
                                                target_pressure=target_pressure, err=err)
      if (allocated(err)) return
    endif

    call commit_atmosphere_candidate(self, candidate, err)

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
