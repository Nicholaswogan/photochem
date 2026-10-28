  function allocate_evoatmospherestate() result(ptr) bind(c)
    type(c_ptr) :: ptr
    type(EvoAtmosphereState), pointer :: state
    allocate(state)
    ptr = c_loc(state)
  end function

  subroutine deallocate_evoatmospherestate(ptr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    type(EvoAtmosphereState), pointer :: state
    if (.not.c_associated(ptr)) return
    call c_f_pointer(ptr, state)
    deallocate(state)
  end subroutine

  subroutine evoatmospherestate_initialize_from_h5_wrapper(ptr, filename, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(c_char), intent(in) :: filename(*)
    character(c_char), intent(out) :: err(err_len+1)
    type(EvoAtmosphereState), pointer :: state
    character(:), allocatable :: filename_f, err_f

    allocate(character(len=len_cstring(filename)) :: filename_f)
    call copy_string_ctof(filename, filename_f)
    call c_f_pointer(ptr, state)
    call state%initialize_from_h5(filename_f, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) call copy_string_ftoc(err_f, err)
  end subroutine

  subroutine evoatmospherestate_save_to_h5_wrapper(ptr, filename, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(c_char), intent(in) :: filename(*)
    character(c_char), intent(out) :: err(err_len+1)
    type(EvoAtmosphereState), pointer :: state
    character(:), allocatable :: filename_f, err_f

    allocate(character(len=len_cstring(filename)) :: filename_f)
    call copy_string_ctof(filename, filename_f)
    call c_f_pointer(ptr, state)
    call state%save_to_h5(filename_f, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) call copy_string_ftoc(err_f, err)
  end subroutine

  subroutine evoatmospherestate_metadata_get(ptr, integers, reals, lengths) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: integers(10), lengths(3)
    real(c_double), intent(out) :: reals(12)
    type(EvoAtmosphereState), pointer :: state

    call c_f_pointer(ptr, state)
    ! Integer layout: nspecies, nq, ng, nsl, nz, np, enabled, mode, hydro, toa_enabled.
    integers = [size(state%species_names), state%nq, state%ng, state%nsl, size(state%z), &
                size(state%particle_radius,1), merge(1,0,state%press_temp_edd_profile%enabled), &
                state%press_temp_edd_profile%mode, merge(1,0,state%press_temp_edd_profile%hydro_pressure), &
                merge(1,0,state%toa_pressure_maintenance%enabled)]
    ! Real layout follows the named scalar fields in EvoAtmosphereState.pyx.
    reals = [state%planet_mass, state%bottom_atmos, state%top_atmos, state%trop_alt, state%planet_radius, &
             state%press_temp_edd_profile%trop_p, state%press_temp_edd_profile%temperature_tol, &
             state%press_temp_edd_profile%edd_tol, state%press_temp_edd_profile%extreme_factor, &
             state%toa_pressure_maintenance%target_pressure, state%toa_pressure_maintenance%pressure_factor, &
             state%toa_pressure_maintenance%extreme_pressure_factor]
    lengths = -1
    if (allocated(state%press_temp_edd_profile%pressure)) lengths(1) = size(state%press_temp_edd_profile%pressure)
    if (allocated(state%press_temp_edd_profile%temperature)) lengths(2) = size(state%press_temp_edd_profile%temperature)
    if (allocated(state%press_temp_edd_profile%edd)) lengths(3) = size(state%press_temp_edd_profile%edd)
  end subroutine

  subroutine evoatmospherestate_species_names_get(ptr, n, names) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: n
    character(c_char), intent(out) :: names(s_str_len*n)
    type(EvoAtmosphereState), pointer :: state
    integer :: i, j

    call c_f_pointer(ptr, state)
    do i = 1,n
      do j = 1,s_str_len
        names((i-1)*s_str_len+j) = state%species_names(i)(j:j)
      enddo
    enddo
  end subroutine

  subroutine evoatmospherestate_array_get(ptr, field, n, values) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: field, n
    real(c_double), intent(out) :: values(n)
    type(EvoAtmosphereState), pointer :: state

    call c_f_pointer(ptr, state)
    select case(field)
    case(1)
      values = state%z
    case(2)
      values = state%dz
    case(3)
      values = state%temperature
    case(4)
      values = state%edd
    case(5)
      values = reshape(state%particle_radius, [n])
    case(6)
      values = reshape(state%usol, [n])
    case(7)
      values = state%press_temp_edd_profile%pressure
    case(8)
      values = state%press_temp_edd_profile%temperature
    case(9)
      values = state%press_temp_edd_profile%edd
    end select
  end subroutine

  subroutine evoatmospherestate_create_wrapper(ptr, integers, reals, lengths, names, z, dz, temperature, edd, &
      particle_radius, usol, profile_pressure, profile_temperature, profile_edd, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: integers(10), lengths(3)
    real(c_double), intent(in) :: reals(12)
    character(c_char), intent(in) :: names(s_str_len*integers(1))
    real(c_double), intent(in) :: z(integers(5)), dz(integers(5)), temperature(integers(5)), edd(integers(5))
    real(c_double), intent(in) :: particle_radius(integers(6),integers(5)), usol(integers(2),integers(5))
    real(c_double), intent(in) :: profile_pressure(max(0,lengths(1))), profile_temperature(max(0,lengths(2)))
    real(c_double), intent(in) :: profile_edd(max(0,lengths(3)))
    character(c_char), intent(out) :: err(err_len+1)
    type(EvoAtmosphereState), pointer :: state
    type(EvoAtmosphereState) :: candidate
    type(PressureTempEddProfile) :: profile
    type(TOAPressureMaintenance) :: toa
    character(s_str_len), allocatable :: species_names(:)
    character(:), allocatable :: err_f
    integer :: i, j

    err(1) = c_null_char
    if (integers(1) < 3 .or. integers(2) < 1 .or. integers(5) < 2 .or. integers(6) < 0 .or. &
        any(lengths < -1) .or. any(integers([7,9,10]) < 0) .or. any(integers([7,9,10]) > 1)) then
      call copy_string_ftoc('Invalid atmospheric state dimensions or flags.', err)
      return
    endif
    allocate(species_names(integers(1)))
    do i = 1,integers(1)
      do j = 1,s_str_len
        species_names(i)(j:j) = names((i-1)*s_str_len+j)
      enddo
    enddo
    profile%enabled = integers(7) == 1
    profile%mode = integers(8)
    profile%hydro_pressure = integers(9) == 1
    profile%trop_p = reals(6)
    profile%temperature_tol = reals(7)
    profile%edd_tol = reals(8)
    profile%extreme_factor = reals(9)
    if (lengths(1) >= 0) profile%pressure = profile_pressure
    if (lengths(2) >= 0) profile%temperature = profile_temperature
    if (lengths(3) >= 0) profile%edd = profile_edd
    toa%enabled = integers(10) == 1
    toa%target_pressure = reals(10)
    toa%pressure_factor = reals(11)
    toa%extreme_pressure_factor = reals(12)
    candidate = EvoAtmosphereState(species_names, integers(2), integers(3), integers(4), reals(1), &
      reals(2), reals(3), reals(4), z, dz, temperature, edd, particle_radius, usol, reals(5), profile, toa, err_f)
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
      return
    endif
    call c_f_pointer(ptr, state)
    state = candidate
  end subroutine

