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
    use photochem_input, only: setup_atmosphere_from_file

    class(EvoAtmosphere), intent(inout) :: self
    character(len=*), intent(in) :: atmosphere_txt
    character(:), allocatable, intent(out) :: err

    type(EvoAtmosphere) :: candidate

    call create_atmosphere_candidate(self, candidate, err)
    if (allocated(err)) return

    call reset_press_temp_edd_profile(candidate%var)

    call setup_atmosphere_from_file(atmosphere_txt, candidate%dat, candidate%var, err)
    if (allocated(err)) return

    call prepare_atmosphere_candidate(candidate, err)
    if (allocated(err)) return

    call commit_atmosphere_candidate(self, candidate, err)

  end subroutine

  module subroutine initialize_atmosphere_z(self, z, temperature, edd, &
                                            surface_pressure, mix, &
                                            particle_radius, err)
    use photochem_input, only: map_atmosphere_z_to_grid

    class(EvoAtmosphere), intent(inout) :: self
    real(dp), intent(in) :: z(:), temperature(:), edd(:)
    real(dp), intent(in) :: surface_pressure
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    character(:), allocatable, intent(out) :: err

    type(EvoAtmosphere) :: candidate
    real(dp), allocatable :: pressure(:), density(:), mubar(:)

    call create_atmosphere_candidate(self, candidate, err)
    if (allocated(err)) return

    call reset_press_temp_edd_profile(candidate%var)

    allocate(pressure(candidate%var%nz), density(candidate%var%nz), &
             mubar(candidate%var%nz))
    call map_atmosphere_z_to_grid(candidate%dat, candidate%var, z, temperature, &
                                  edd, surface_pressure, mix, particle_radius, &
                                  pressure, density, mubar, err)
    if (allocated(err)) return

    candidate%var%top_atmos_from_file = .false.
    call prepare_atmosphere_candidate(candidate, err)
    if (allocated(err)) return

    call commit_atmosphere_candidate(self, candidate, err)

  end subroutine

  module subroutine initialize_atmosphere_p(self, pressure, temperature, edd, &
                                            mix, particle_radius, persistent, &
                                            trop_p, err)
    use photochem_input, only: map_atmosphere_p_to_grid

    class(EvoAtmosphere), intent(inout) :: self
    real(dp), intent(in) :: pressure(:), temperature(:), edd(:)
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    logical, optional, intent(in) :: persistent
    real(dp), optional, intent(in) :: trop_p
    character(:), allocatable, intent(out) :: err

    type(EvoAtmosphere) :: candidate
    real(dp), allocatable :: pressure_model(:), density(:), mubar(:)
    logical :: persistent_

    persistent_ = .false.
    if (present(persistent)) persistent_ = persistent
    if (.not. persistent_ .and. present(trop_p)) then
      err = '"trop_p" can only be specified when "persistent" is true.'
      return
    endif

    call create_atmosphere_candidate(self, candidate, err)
    if (allocated(err)) return

    call reset_press_temp_edd_profile(candidate%var)

    allocate(pressure_model(candidate%var%nz), density(candidate%var%nz), &
             mubar(candidate%var%nz))
    call map_atmosphere_p_to_grid(candidate%dat, candidate%var, pressure, &
                                  temperature, edd, mix, particle_radius, &
                                  pressure_model, density, mubar, err)
    if (allocated(err)) return

    candidate%var%top_atmos_from_file = .false.
    call prepare_atmosphere_candidate(candidate, err)
    if (allocated(err)) return

    if (persistent_) then
      if (present(trop_p)) then
        call candidate%set_press_temp_edd_profile(pressure, temperature, edd, &
                                                  trop_p=trop_p, &
                                                  hydro_pressure=.true., err=err)
      else
        call candidate%set_press_temp_edd_profile(pressure, temperature, edd, &
                                                  hydro_pressure=.true., err=err)
      endif
      if (allocated(err)) return
      candidate%var%usol_init = candidate%wrk%usol
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
    call clear_legacy_atmosphere_data(candidate%dat)
    call candidate%wrk%init(candidate%dat%nsp, candidate%dat%np, candidate%dat%nq, &
                            candidate%var%nz, candidate%dat%nrT, candidate%dat%kj, &
                            candidate%dat%nw)

  end subroutine

  subroutine prepare_atmosphere_candidate(candidate, err)
    type(EvoAtmosphere), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    call candidate%prep_atmosphere_unchecked(candidate%var%usol_init, err)
    if (allocated(err)) return

    ! Boundary conditions are imposed while preparing the atmosphere. Retain
    ! that prepared state as the canonical initial condition.
    candidate%var%usol_init = candidate%wrk%usol
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

  subroutine clear_legacy_atmosphere_data(dat)
    use photochem_types, only: PhotochemData

    type(PhotochemData), intent(inout) :: dat

    if (allocated(dat%z_file)) deallocate(dat%z_file)
    if (allocated(dat%T_file)) deallocate(dat%T_file)
    if (allocated(dat%edd_file)) deallocate(dat%edd_file)
    if (allocated(dat%den_file)) deallocate(dat%den_file)
    if (allocated(dat%mix_file)) deallocate(dat%mix_file)
    if (allocated(dat%particle_radius_file)) deallocate(dat%particle_radius_file)

  end subroutine

  module subroutine reset_press_temp_edd_profile(var)
    type(PhotochemVars), intent(inout) :: var

    var%press_temp_edd_profile%enabled = .false.
    var%press_temp_edd_profile%hydro_pressure = .true.
    var%press_temp_edd_profile%has_trop_p = .false.
    var%press_temp_edd_profile%trop_p = 0.0_dp
    if (allocated(var%press_temp_edd_profile%pressure)) &
      deallocate(var%press_temp_edd_profile%pressure)
    if (allocated(var%press_temp_edd_profile%temperature)) &
      deallocate(var%press_temp_edd_profile%temperature)
    if (allocated(var%press_temp_edd_profile%edd)) &
      deallocate(var%press_temp_edd_profile%edd)

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
