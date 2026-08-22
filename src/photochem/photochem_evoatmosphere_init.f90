submodule(photochem_evoatmosphere) photochem_evoatmosphere_init
  implicit none

  ! Aggregate construction and atmosphere-initialization implementation.
  
  ! Raw atmosphere-file data used only during legacy-file initialization.
  type :: AtmosphereFileProfile
    integer :: nlayer = 0
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: edd(:)
    real(dp), allocatable :: density(:)
    real(dp), allocatable :: mix(:,:)
    real(dp), allocatable :: particle_radius(:,:)
  end type

contains

  module subroutine AtmosphereState_allocate(self, dat, nz)
    class(AtmosphereState), intent(inout) :: self
    type(PhotochemData), intent(in) :: dat
    integer, intent(in) :: nz
    integer :: i

    allocate(self%z(nz), self%dz(nz), self%temperature(nz))
    allocate(self%edd(nz), self%particle_radius(dat%npq,nz), &
             self%usol(dat%nq,nz))
    allocate(self%grav(nz), self%xs_x_qy(nz,dat%kj,dat%nw))
    allocate(self%particle_xs(dat%np))
    if (dat%reverse) allocate(self%gibbs_energy(nz,dat%ng))

    do i = 1,dat%np
      if (dat%part_xs_file(i)%ThereIsData) then
        self%particle_xs(i)%ThereIsData = .true.
        allocate(self%particle_xs(i)%w0(nz,dat%nw))
        allocate(self%particle_xs(i)%qext(nz,dat%nw))
        allocate(self%particle_xs(i)%gt(nz,dat%nw))
      else
        self%particle_xs(i)%ThereIsData = .false.
      endif
    enddo

  end subroutine

  
  module function create_EvoAtmosphere(mechanism_file, settings_file, flux_file, atmosphere_txt, data_dir, err) result(self)
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphere) :: self

    self = create_EvoAtmosphere_static(mechanism_file, settings_file, &
                                       flux_file, data_dir, err)
    if (allocated(err)) return

    call self%initialize_from_atmosphere_file(atmosphere_txt, err)

  end function

  module function create_EvoAtmosphere_static(mechanism_file, settings_file, &
                                              flux_file, data_dir, err) result(self)
    use photochem_settings, only: PhotoSettings
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphere) :: self

    type(PhotoSettings) :: settings

    settings = PhotoSettings(settings_file, err)
    if (allocated(err)) return

    self%dat = PhotochemData(mechanism_file, settings, data_dir, err)
    if (allocated(err)) return

    self%var = PhotochemVars(self%dat, settings, flux_file, err)
    if (allocated(err)) return

    self%wrk = PhotochemWrk(self%dat%nsp, self%dat%np, self%dat%nq, &
                            self%var%nz, self%dat%nrT, self%dat%kj, &
                            self%dat%nw)
    self%atmosphere_initialized = .false.

  end function

  module subroutine initialize_from_atmosphere_file(self, atmosphere_txt, err)
    class(EvoAtmosphere), intent(inout) :: self
    character(len=*), intent(in) :: atmosphere_txt
    character(:), allocatable, intent(out) :: err

    type(AtmosphereState) :: state, previous_state
    logical :: was_initialized

    was_initialized = self%atmosphere_initialized
    if (was_initialized) then
      call previous_state%allocate(self%dat, self%var%nz)
      call copy_model_to_state(self, previous_state)
    endif

    call state%allocate(self%dat, self%var%nz)
    call copy_model_to_state(self, state)

    call map_atmosphere_file_to_grid(self%dat, atmosphere_txt, self%var%trop_alt, &
                                     state, err)
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

    ! Atmosphere-file initialization does not retain a persistent P-T-Kzz profile.
    call self%clear_press_temp_edd_profile(err)
    if (allocated(err)) return

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

  module subroutine initialize_atmosphere_z(self, z, temperature, edd, &
                                            surface_pressure, mix, &
                                            particle_radius, err)
    class(EvoAtmosphere), intent(inout) :: self
    real(dp), intent(in) :: z(:), temperature(:), edd(:)
    real(dp), intent(in) :: surface_pressure
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    character(:), allocatable, intent(out) :: err

    type(AtmosphereState) :: state, previous_state
    logical :: was_initialized

    was_initialized = self%atmosphere_initialized
    ! If initialized, then save the previous state
    if (was_initialized) then
      call previous_state%allocate(self%dat, self%var%nz)
      call copy_model_to_state(self, previous_state)
    endif

    call state%allocate(self%dat, self%var%nz)
    call copy_model_to_state(self, state)

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

  module subroutine copy_model_to_state(self, state)
    class(EvoAtmosphere), intent(in) :: self
    type(AtmosphereState), intent(inout) :: state

    ! Always preserve valid persistent settings
    state%press_temp_edd_profile = self%var%press_temp_edd_profile
    state%toa_pressure_maintenance = self%var%toa_pressure_maintenance

    if (.not. self%atmosphere_initialized) return

    ! Copy committed atmospheric state only after initialization
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

  end subroutine

  module subroutine copy_state_to_model(self, state)
    class(EvoAtmosphere), intent(inout) :: self
    type(AtmosphereState), intent(in) :: state

    self%var%press_temp_edd_profile = state%press_temp_edd_profile
    self%var%toa_pressure_maintenance = state%toa_pressure_maintenance

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

  end subroutine

  module subroutine initialize_atmosphere_p(self, pressure, temperature, edd, &
                                            mix, particle_radius, persistent, &
                                            trop_p, maintain_toa_pressure, &
                                            target_pressure, err)

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
    if (was_initialized) then
      call previous_state%allocate(self%dat, self%var%nz)
      call copy_model_to_state(self, previous_state)
    endif

    call state%allocate(self%dat, self%var%nz)
    call copy_model_to_state(self, state)

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

  module subroutine require_atmosphere_initialized(self, operation, err)
    class(EvoAtmosphere), intent(in) :: self
    character(len=*), intent(in) :: operation
    character(:), allocatable, intent(out) :: err

    if (.not. self%atmosphere_initialized) then
      err = 'EvoAtmosphere atmosphere is not initialized. Initialize the atmosphere before calling "'// &
            trim(operation)//'".'
    endif

  end subroutine
  

  subroutine read_atmosphere_file(atmosphere_txt, dat, profile, err)
    use futils, only: FileCloser
    use photochem_const, only: s_str_len
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(in) :: dat
    type(AtmosphereFileProfile), intent(out) :: profile
    character(:), allocatable, intent(out) :: err

    character(len=10000) :: line
    character(len=s_str_len) :: arr1(1000)
    character(len=s_str_len) :: arr11(1000)
    character(len=s_str_len),allocatable, dimension(:) :: labels
    integer :: ind(1)
    real(dp), allocatable :: temp(:,:)
    integer :: io, i, n, nn, ii
    type(FileCloser) :: file

    open(4, file=trim(atmosphere_txt),status='old',iostat=io)
    file%unit = 4
    if (io /= 0) then
      err = 'Can not open file '//trim(atmosphere_txt)
      return
    endif
    read(4,'(A)') line

    profile%nlayer = -1
    io = 0
    do while (io == 0)
      read(4,*,iostat=io)
      profile%nlayer = profile%nlayer + 1
    enddo

    allocate(profile%z(profile%nlayer))
    allocate(profile%temperature(profile%nlayer))
    allocate(profile%edd(profile%nlayer))
    allocate(profile%density(profile%nlayer))
    allocate(profile%mix(dat%nq, profile%nlayer))
    profile%z = 0.0_dp
    profile%temperature = 0.0_dp
    profile%edd = 0.0_dp
    profile%mix = 1.0e-40_dp
    if (dat%there_are_particles) then
      allocate(profile%particle_radius(dat%npq, profile%nlayer))
    endif

    rewind(4)
    read(4,'(A)') line
    n = 0
    nn = 0
    do i=1,1000
      read(line,*,iostat=io) arr1(1:i)
      if (io==-1) exit
      n = n+1
    enddo
    read(4,'(A)') line
    do i=1,1000
      read(line,*,iostat=io) arr11(1:i)
      if (io==-1) exit
      nn = nn+1
    enddo
    if (n /= nn) then
      err = 'There is a missing column label in the file '//trim(atmosphere_txt)
      return
    endif

    ! allocate memory
    allocate(labels(n))
    allocate(temp(n,profile%nlayer))
    rewind(4)
    read(4,'(A)') line
    read(line,*) (labels(i),i=1,n)

    ! First read in all the data into big array
    do i = 1,profile%nlayer
      read(4,*,iostat=io) (temp(ii,i),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in initial atmosphere in '//trim(atmosphere_txt)
        return
      endif
    enddo

    ! reads in mixing ratios
    do i=1,dat%nq
      ind = findloc(labels,dat%species_names(i))
      if (ind(1) /= 0) then
        profile%mix(i,:) = temp(ind(1),:)
      endif
    enddo

    if (dat%there_are_particles) then
      do i=1,dat%npq
        ind = findloc(labels,trim(dat%species_names(i))//"_r")
        if (ind(1) /= 0) then
          profile%particle_radius(i,:) = temp(ind(1),:)
        else
          ! did not find the data
          ! will set to 0.1 micron
          profile%particle_radius(i,:) = 1.0e-5_dp
        endif
      enddo
    endif

    ! reads in temperature
    ind = findloc(labels,'temp')
    if (ind(1) /= 0) then
      profile%temperature(:) = temp(ind(1),:)
    else
      err = '"temp" was not found in input file '//trim(atmosphere_txt)
      return
    endif

    ! reads in alt
    ind = findloc(labels,'alt')
    if (ind(1) /= 0) then
      profile%z(:) = temp(ind(1),:)*1.e5_dp ! convert to cm
    else
      err = '"alt" was not found in input file '//trim(atmosphere_txt)
      return
    endif

    ! reads in eddy diffusion
    ind = findloc(labels,'eddy')
    if (ind(1) /= 0) then
      profile%edd(:) = temp(ind(1),:)
    else
      err = '"eddy" was not found in input file '//trim(atmosphere_txt)
      return
    endif

    ! reads in density.
    ind = findloc(labels,'den')
    if (ind(1) /= 0) then
      profile%density(:) = temp(ind(1),:)
    else
      err = '"den" was not found in input file '//trim(atmosphere_txt)
      return
    endif

  end subroutine

  subroutine map_atmosphere_z_to_grid(dat, nz, trop_alt, z, temperature, &
                                             edd, surface_pressure, mix, &
                                             particle_radius, state, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_eqns, only: gravity, press_and_den, vertical_grid

    type(PhotochemData), intent(in) :: dat
    integer, intent(in) :: nz
    real(dp), intent(in) :: trop_alt
    real(dp), intent(in) :: z(:), temperature(:), edd(:)
    real(dp), intent(in) :: surface_pressure
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    type(AtmosphereState), intent(inout) :: state
    character(:), allocatable, intent(out) :: err

    real(dp), parameter :: mixing_ratio_floor = 1.0e-40_dp
    real(dp), allocatable :: gas_mix_normalized(:,:), gas_mix_model(:,:), particle_mix_model(:,:)
    real(dp), allocatable :: particle_mix_profile(:,:)
    real(dp), allocatable :: pressure(:), density(:), mubar(:)
    real(dp) :: gas_total, z_tolerance
    integer :: i, j, nprofile, ngas

    nprofile = size(z)
    ngas = dat%nq - dat%ng_1 + 1

    if (nz < 2) then
      err = 'Altitude initialization requires at least two model layers.'
      return
    endif
    if (nprofile < 2) then
      err = 'Altitude initialization requires at least two profile points.'
      return
    endif
    if (size(temperature) /= nprofile .or. size(edd) /= nprofile) then
      err = 'Altitude, temperature, and eddy-diffusion profiles must have the same length.'
      return
    endif
    if (size(mix,1) /= dat%nq .or. size(mix,2) /= nprofile) then
      err = 'The mixing-ratio array has the wrong shape.'
      return
    endif
    if (size(particle_radius,1) /= dat%npq .or. size(particle_radius,2) /= nprofile) then
      err = 'The particle-radius array has the wrong shape.'
      return
    endif

    if (.not. all(ieee_is_finite(z))) then
      err = 'Altitude profile contains a nonfinite value.'
      return
    endif
    if (.not. all(ieee_is_finite(temperature)) .or. any(temperature <= 0.0_dp)) then
      err = 'Temperature profile must contain only finite, positive values.'
      return
    endif
    if (.not. all(ieee_is_finite(edd)) .or. any(edd <= 0.0_dp)) then
      err = 'Eddy-diffusion profile must contain only finite, positive values.'
      return
    endif
    if (.not. ieee_is_finite(surface_pressure) .or. surface_pressure <= 0.0_dp) then
      err = 'Surface pressure must be finite and positive.'
      return
    endif
    if (.not. all(ieee_is_finite(mix)) .or. any(mix < 0.0_dp)) then
      err = 'Mixing ratios must contain only finite, nonnegative values.'
      return
    endif
    if (.not. all(ieee_is_finite(particle_radius)) .or. any(particle_radius <= 0.0_dp)) then
      err = 'Particle radii must contain only finite, positive values.'
      return
    endif
    if (any(z(2:) <= z(:nprofile-1))) then
      err = 'Altitude profile must be strictly increasing.'
      return
    endif
    z_tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(z(nprofile)))
    if (abs(z(1)) > z_tolerance) then
      err = 'The first altitude profile point must be zero.'
      return
    endif
    if (dat%gas_rainout .and. trop_alt >= z(nprofile)) then
      err = 'Tropopause altitude must be below the model-top altitude.'
      return
    endif

    allocate(gas_mix_normalized(ngas,nprofile))
    do j = 1,nprofile
      gas_total = sum(mix(dat%ng_1:dat%nq,j))
      if (.not. ieee_is_finite(gas_total) .or. gas_total <= 0.0_dp) then
        err = 'At least one gas mixing ratio must be positive at every altitude.'
        return
      endif
      gas_mix_normalized(:,j) = max(mix(dat%ng_1:dat%nq,j)/gas_total, mixing_ratio_floor)
      gas_mix_normalized(:,j) = gas_mix_normalized(:,j)/sum(gas_mix_normalized(:,j))
    enddo

    state%bottom_atmos = 0.0_dp
    state%top_atmos = z(nprofile)
    state%trop_alt = trop_alt
    call vertical_grid(state%bottom_atmos, state%top_atmos, state%z, state%dz)
    call gravity(dat%planet_radius, dat%planet_mass, state%z, state%grav)
    call interpolate_profile_1d(state%z, z, temperature, state%temperature, &
                                .false., .false., 'temperature', err)
    if (allocated(err)) return

    call interpolate_profile_1d(state%z, z, edd, state%edd, &
                                .true., .false., 'eddy diffusion', err)
    if (allocated(err)) return

    allocate(gas_mix_model(ngas,nz), pressure(nz), density(nz), mubar(nz))
    call interpolate_profiles_2d(state%z, z, gas_mix_normalized, gas_mix_model, &
                                 .true., .false., 'gas mixing ratios', err)
    if (allocated(err)) return
    do j = 1,nz
      gas_mix_model(:,j) = gas_mix_model(:,j)/sum(gas_mix_model(:,j))
      mubar(j) = sum(gas_mix_model(:,j)*dat%species_mass(dat%ng_1:dat%nq))
    enddo
    call press_and_den(state%temperature, state%grav, surface_pressure, &
                       state%dz, mubar, pressure, density)
    if (.not. all(ieee_is_finite(pressure)) .or. any(pressure <= 0.0_dp) .or. &
        .not. all(ieee_is_finite(density)) .or. any(density <= 0.0_dp)) then
      err = 'Hydrostatic integration produced an invalid pressure or density.'
      return
    endif

    state%usol = 0.0_dp
    do i = 1,ngas
      state%usol(dat%ng_1+i-1,:) = gas_mix_model(i,:)*density
    enddo

    if (dat%npq > 0) then
      allocate(particle_mix_profile(dat%npq,nprofile), &
               particle_mix_model(dat%npq,nz))
      particle_mix_profile = max(mix(:dat%npq,:), mixing_ratio_floor)
      call interpolate_profiles_2d(state%z, z, particle_mix_profile, &
                                   particle_mix_model, .true., .false., &
                                   'particle mixing ratios', err)
      if (allocated(err)) return

      call interpolate_profiles_2d(state%z, z, particle_radius, &
                                   state%particle_radius, .true., .false., &
                                   'particle radii', err)
      if (allocated(err)) return

      do i = 1,dat%npq
        state%usol(i,:) = particle_mix_model(i,:)*density
      enddo
    endif

  end subroutine

  subroutine map_atmosphere_p_to_grid(dat, nz, trop_alt_default, &
                                             profile_pressure, temperature, &
                                             edd, mix, particle_radius, trop_p, &
                                             state, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_const, only: k_boltz, N_avo
    use photochem_eqns, only: gravity

    type(PhotochemData), intent(in) :: dat
    integer, intent(in) :: nz
    real(dp), intent(in) :: trop_alt_default
    real(dp), intent(in) :: profile_pressure(:), temperature(:), edd(:)
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    real(dp), optional, intent(in) :: trop_p
    type(AtmosphereState), intent(inout) :: state
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: z(:), profile_mubar(:)
    real(dp) :: gas_total, inverse_radius, inverse_radius_new
    real(dp) :: inverse_radius_factor, delta_log_pressure
    real(dp) :: surface_z(1), surface_gravity(1)
    real(dp) :: trop_alt_array(1)
    real(dp) :: trop_alt
    integer :: i, nprofile, ierr

    nprofile = size(profile_pressure)

    ! These checks are needed before pressure and composition can safely be
    ! converted to altitude. The altitude mapper performs the remaining common
    ! profile validation before changing state.
    if (nprofile < 2) then
      err = 'Pressure initialization requires at least two profile points.'
      return
    endif
    if (size(temperature) /= nprofile .or. size(edd) /= nprofile) then
      err = 'Pressure, temperature, and eddy-diffusion profiles must have the same length.'
      return
    endif
    if (size(mix,1) /= dat%nq .or. size(mix,2) /= nprofile) then
      err = 'The mixing-ratio array has the wrong shape.'
      return
    endif
    if (.not. all(ieee_is_finite(profile_pressure)) .or. &
        any(profile_pressure <= 0.0_dp)) then
      err = 'Pressure profile must contain only finite, positive values.'
      return
    endif
    if (any(profile_pressure(2:) >= profile_pressure(:nprofile-1))) then
      err = 'Pressure profile must be strictly decreasing.'
      return
    endif
    if (.not. all(ieee_is_finite(temperature)) .or. any(temperature <= 0.0_dp)) then
      err = 'Temperature profile must contain only finite, positive values.'
      return
    endif
    if (.not. all(ieee_is_finite(mix)) .or. any(mix < 0.0_dp)) then
      err = 'Mixing ratios must contain only finite, nonnegative values.'
      return
    endif

    allocate(z(nprofile), profile_mubar(nprofile))
    do i = 1,nprofile
      gas_total = sum(mix(dat%ng_1:dat%nq,i))
      if (.not. ieee_is_finite(gas_total) .or. gas_total <= 0.0_dp) then
        err = 'At least one gas mixing ratio must be positive at every pressure.'
        return
      endif
      profile_mubar(i) = sum(mix(dat%ng_1:dat%nq,i)* &
                             dat%species_mass(dat%ng_1:dat%nq))/gas_total
      if (.not. ieee_is_finite(profile_mubar(i)) .or. profile_mubar(i) <= 0.0_dp) then
        err = 'Gas composition produced a nonfinite or nonpositive mean molecular weight.'
        return
      endif
    enddo

    ! For spherical gravity, hydrostatic balance can be integrated in inverse
    ! radius without a nonlinear solve:
    !
    !   d(1/r)/d(log(P)) = N_A*k_B*T/(mubar*g_surface*R_surface^2).
    !
    ! Apply the trapezoid rule to T/mubar between pressure knots. The configured
    ! planet radius is the radius at the lower pressure boundary.
    surface_z = 0.0_dp
    call gravity(dat%planet_radius, dat%planet_mass, surface_z, surface_gravity)
    if (.not. ieee_is_finite(surface_gravity(1)) .or. surface_gravity(1) <= 0.0_dp) then
      err = 'Could not compute finite, positive gravity at the lower boundary.'
      return
    endif

    inverse_radius_factor = N_avo*k_boltz/ &
                            (surface_gravity(1)*dat%planet_radius**2)
    inverse_radius = 1.0_dp/dat%planet_radius
    z(1) = 0.0_dp
    do i = 2,nprofile
      delta_log_pressure = log(profile_pressure(i))-log(profile_pressure(i-1))
      inverse_radius_new = inverse_radius + inverse_radius_factor*0.5_dp* &
                           (temperature(i-1)/profile_mubar(i-1) + &
                            temperature(i)/profile_mubar(i))*delta_log_pressure
      if (.not. ieee_is_finite(inverse_radius_new) .or. inverse_radius_new <= 0.0_dp) then
        err = 'Pressure profile extends beyond a finite hydrostatic altitude.'
        return
      endif

      z(i) = 1.0_dp/inverse_radius_new-dat%planet_radius
      if (.not. ieee_is_finite(z(i)) .or. z(i) <= z(i-1)) then
        err = 'Hydrostatic pressure integration did not produce increasing altitude.'
        return
      endif
      inverse_radius = inverse_radius_new
    enddo

    if (present(trop_p)) then
      if (.not. dat%gas_rainout) then
        err = '"trop_p" can only be supplied when gas rainout is enabled.'
        return
      endif
      if (.not. ieee_is_finite(trop_p) .or. trop_p <= 0.0_dp) then
        err = '"trop_p" must be finite and positive.'
        return
      endif
      if (trop_p > profile_pressure(1) .or. trop_p < profile_pressure(nprofile)) then
        err = '"trop_p" must lie within the supplied pressure profile.'
        return
      endif

      ! The pressure profile is descending while the interpolation abscissa
      ! must be ascending. Install the resulting altitude before the common
      ! altitude mapper validates the candidate tropopause.
      call interp([log(trop_p)], log(profile_pressure(nprofile:1:-1)), &
                  z(nprofile:1:-1), trop_alt_array, ierr=ierr)
      if (ierr /= 0 .or. .not. ieee_is_finite(trop_alt_array(1))) then
        err = 'Unable to determine the tropopause altitude from "trop_p".'
        return
      endif
      trop_alt = trop_alt_array(1)
    else
      trop_alt = trop_alt_default
    endif

    call map_atmosphere_z_to_grid(dat, nz, trop_alt, z, temperature, &
                                  edd, profile_pressure(1), mix, particle_radius, &
                                  state, err)
    if (allocated(err)) return

  end subroutine

  subroutine map_atmosphere_file_to_grid(dat, atmosphere_txt, trop_alt_default, state, err)
    use photochem_eqns, only: gravity, vertical_grid
    type(PhotochemData), intent(in) :: dat
    character(len=*), intent(in) :: atmosphere_txt
    real(dp), intent(in) :: trop_alt_default
    type(AtmosphereState), intent(inout) :: state
    character(:), allocatable, intent(out) :: err
    type(AtmosphereFileProfile) :: profile
    integer :: i
    real(dp), allocatable :: density(:), mix(:,:)

    call read_atmosphere_file(atmosphere_txt, dat, profile, err)
    if (allocated(err)) return

    if (profile%nlayer < 2) then
      err = 'Atmosphere file must contain at least two data rows.'
      return
    endif

    state%bottom_atmos = 0.0_dp
    state%top_atmos = profile%z(profile%nlayer) + &
                    0.5_dp*(profile%z(profile%nlayer) - profile%z(profile%nlayer-1))
    state%trop_alt = trop_alt_default

    if (state%top_atmos < state%bottom_atmos) then
      err = 'The top of the atmosphere must be bigger than the bottom'
      return
    endif

    if (dat%gas_rainout .and. state%trop_alt > state%top_atmos) then
      err = 'IOError: tropopause-altitude must be between the top and bottom of the atmosphere'
      return
    endif

    call vertical_grid(state%bottom_atmos, state%top_atmos, state%z, state%dz)
    call gravity(dat%planet_radius, dat%planet_mass, state%z, state%grav)

    if (size(state%usol,1) /= size(profile%mix,1) .or. &
        size(state%usol,2) /= size(state%z)) then
      err = 'The initial atmospheric-state array has the wrong shape.'
      return
    endif

    call interpolate_profile_1d(state%z, profile%z, profile%temperature, &
                                state%temperature, .false., .false., &
                                'temperature', err)
    if (allocated(err)) return

    call interpolate_profile_1d(state%z, profile%z, profile%edd, state%edd, &
                                .true., .false., 'eddy diffusion', err)
    if (allocated(err)) return

    allocate(density(size(state%z)), mix(size(profile%mix,1),size(state%z)))

    call interpolate_profile_1d(state%z, profile%z, profile%density, density, &
                                .true., .true., 'atmospheric density', err)
    if (allocated(err)) return

    call interpolate_profiles_2d(state%z, profile%z, profile%mix, mix, &
                                 .true., .false., 'mixing ratios', err)
    if (allocated(err)) return

    do i = 1,size(state%z)
      state%usol(:,i) = mix(:,i)*density(i)
    enddo

    if (dat%there_are_particles) then
      call interpolate_profiles_2d(state%z, profile%z, profile%particle_radius, &
                                   state%particle_radius, .true., .false., &
                                   'particle radii', err)
    endif
    if (allocated(err)) return

  end subroutine

  ! Interpolate one profile, optionally transforming its values to/from log10
  ! space. The profile values must be positive when log interpolation is used.
  subroutine interpolate_profile_1d(z_grid, profile_z, profile_values, &
                                    values_grid, log_interpolation, &
                                    linear_extrapolation, quantity, err)
    use futils, only: interp
    real(dp), intent(in) :: z_grid(:), profile_z(:), profile_values(:)
    real(dp), intent(out) :: values_grid(:)
    logical, intent(in) :: log_interpolation, linear_extrapolation
    character(*), intent(in) :: quantity
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: profile_values_work(:)
    integer :: ierr

    if (size(profile_z) /= size(profile_values) .or. &
        size(values_grid) /= size(z_grid)) then
      err = 'The one-dimensional interpolation arrays have incompatible shapes.'
      return
    endif

    allocate(profile_values_work(size(profile_values)))
    if (log_interpolation) then
      ! Log-space interpolation requires positive profile values. Callers
      ! that permit zeros should apply their physical floor before calling.
      profile_values_work = log10(profile_values)
    else
      profile_values_work = profile_values
    endif

    call interp(z_grid, profile_z, profile_values_work, values_grid, &
                linear_extrap=linear_extrapolation, ierr=ierr)
    if (ierr /= 0) then
      err = 'Unable to interpolate '//trim(quantity)//' onto the model grid.'
      return
    endif

    if (log_interpolation) then
      values_grid = 10.0_dp**values_grid
    endif

  end subroutine

  ! Apply the one-dimensional interpolation above independently to each row
  ! of a set of profiles sharing the same altitude grid.
  subroutine interpolate_profiles_2d(z_grid, profile_z, profile_values, &
                                     values_grid, log_interpolation, &
                                     linear_extrapolation, quantity, err)
    real(dp), intent(in) :: z_grid(:), profile_z(:), profile_values(:,:)
    real(dp), intent(out) :: values_grid(:,:)
    logical, intent(in) :: log_interpolation, linear_extrapolation
    character(*), intent(in) :: quantity
    character(:), allocatable, intent(out) :: err

    integer :: i

    if (size(profile_values,1) /= size(values_grid,1) .or. &
        size(profile_values,2) /= size(profile_z) .or. &
        size(values_grid,2) /= size(z_grid)) then
      err = 'The two-dimensional interpolation arrays have incompatible shapes.'
      return
    endif

    do i = 1,size(profile_values,1)
      call interpolate_profile_1d(z_grid, profile_z, profile_values(i,:), &
                                  values_grid(i,:), log_interpolation, &
                                  linear_extrapolation, quantity, err)
      if (allocated(err)) return
    enddo

  end subroutine

  ! Finalize derived atmospheric state after grid construction.
  module subroutine finalize_atmosphere_state(dat, state, err)
    use photochem_vars, only: refresh_temperature_dependent_vars, &
                              interp2particlexsdata
    type(PhotochemData), intent(in) :: dat
    type(AtmosphereState), intent(inout) :: state
    character(:), allocatable, intent(out) :: err

    call interp2particlexsdata(dat, state%particle_radius, state%particle_xs, err)
    if (allocated(err)) return

    call refresh_temperature_dependent_vars( &
      dat, state%temperature, state%z, state%bottom_atmos, state%top_atmos, &
      xs_x_qy=state%xs_x_qy, gibbs_energy=state%gibbs_energy, &
      trop_alt=state%trop_alt, trop_ind=state%trop_ind, err=err &
    )
    if (allocated(err)) return

  end subroutine

end submodule
