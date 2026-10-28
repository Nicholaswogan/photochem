submodule(photochem_evoatmosphere) photochem_evoatmosphere_state
  use h5fortran, only: hdf5_file, is_hdf5
  use hdf5
  use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_char, c_int, c_null_char
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none

  integer, parameter :: state_schema_version = 1

  interface
    function rename_file(old, new) bind(C, name="rename") result(status)
      import c_char, c_int
      character(c_char), intent(in) :: old(*), new(*)
      integer(c_int) :: status
    end function
  end interface

contains

  module function load_EvoAtmosphereState(filename, err) result(state)
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphereState) :: state

    call state%initialize_from_h5(filename, err)
  end function


  module subroutine EvoAtmosphereState_initialize_from_h5(self, filename, err)
    class(EvoAtmosphereState), intent(inout) :: self
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphereState), target :: raw
    type(EvoAtmosphereState) :: candidate
    type(hdf5_file) :: h

    if (len_trim(filename) == 0 .or. index(filename, c_null_char) /= 0) then
      err = 'An atmospheric state filename must be nonempty and contain no NUL characters.'
      return
    endif
    if (.not.is_hdf5(filename)) then
      err = 'Not a readable HDF5 atmospheric state file: '//filename
      return
    endif
    call h%open(filename, 'r')
    call read_state(h, raw, err)
    call h%close()
    if (allocated(err)) return

    candidate = EvoAtmosphereState( &
      raw%species_names, raw%nq, raw%ng, raw%nsl, raw%planet_mass, &
      raw%bottom_atmos, raw%top_atmos, raw%trop_alt, raw%z, raw%dz, &
      raw%temperature, raw%edd, raw%particle_radius, raw%usol, raw%planet_radius, &
      raw%press_temp_edd_profile, raw%toa_pressure_maintenance, err)
    if (allocated(err)) return

    ! Commit only after every read, close, and validation has succeeded.
    self%species_names = candidate%species_names
    self%nq = candidate%nq
    self%ng = candidate%ng
    self%nsl = candidate%nsl
    self%planet_mass = candidate%planet_mass
    self%bottom_atmos = candidate%bottom_atmos
    self%top_atmos = candidate%top_atmos
    self%trop_alt = candidate%trop_alt
    self%z = candidate%z
    self%dz = candidate%dz
    self%temperature = candidate%temperature
    self%edd = candidate%edd
    self%particle_radius = candidate%particle_radius
    self%usol = candidate%usol
    self%planet_radius = candidate%planet_radius
    self%press_temp_edd_profile = candidate%press_temp_edd_profile
    self%toa_pressure_maintenance = candidate%toa_pressure_maintenance
  end subroutine

  module subroutine EvoAtmosphereState_save_to_h5(self, filename, err)
    class(EvoAtmosphereState), target, intent(in) :: self
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err
    type(hdf5_file) :: h
    integer :: status, unit, attempt
    integer(int64) :: ticks
    character(64) :: suffix
    character(:), allocatable :: temporary
    logical :: reserved

    call validate_atmosphere_state(self, err)
    if (allocated(err)) return
    if (len_trim(filename) == 0 .or. index(filename, c_null_char) /= 0) then
      err = 'An atmospheric state filename must be nonempty and contain no NUL characters.'
      return
    endif
    ! Reserve a unique sibling without overwriting another writer's temporary.
    call system_clock(count=ticks)
    reserved = .false.
    do attempt = 1,100
      write(suffix, '(i0,a,i0)') ticks, '-', attempt
      temporary = trim(filename)//'.tmp-'//trim(suffix)
      open(newunit=unit, file=temporary, status='new', access='stream', iostat=status)
      if (status == 0) then
        close(unit, iostat=status)
        if (status == 0) reserved = .true.
        exit
      endif
    enddo
    if (.not.reserved) then
      err = 'Could not reserve a temporary atmospheric state file beside: '//filename
      return
    endif

    call h%open(temporary, 'w', comp_lvl=0)
    call write_state(h, self, err)
    call h%close()
    if (.not.allocated(err)) then
      status = rename_file(temporary//c_null_char, trim(filename)//c_null_char)
      if (status /= 0) err = 'Could not rename atmospheric state file to: '//filename
    endif
    if (allocated(err)) then
      open(newunit=unit, file=temporary, status='old', iostat=status)
      if (status == 0) then
        close(unit, status='delete', iostat=status)
        if (status /= 0) err = err//' Could not delete temporary file: '//temporary
      else
        err = err//' Could not reopen temporary file for cleanup: '//temporary
      endif
    endif
  end subroutine

  subroutine write_state(h, state, err)
    type(hdf5_file), intent(inout) :: h
    class(EvoAtmosphereState), target, intent(in) :: state
    character(:), allocatable, intent(out) :: err

    call h%write('schema_version', state_schema_version)
    ! Numeric datasets create the compatibility group before writing strings.
    call h%write('compatibility/nq', state%nq)
    call h%write('compatibility/ng', state%ng)
    call h%write('compatibility/nsl', state%nsl)
    call h%write('compatibility/planet_mass', state%planet_mass)
    call h%write('atmosphere/bottom_atmos', state%bottom_atmos)
    call h%write('atmosphere/top_atmos', state%top_atmos)
    call h%write('atmosphere/trop_alt', state%trop_alt)
    call h%write('atmosphere/planet_radius', state%planet_radius)
    call h%write('atmosphere/z', state%z)
    call h%write('atmosphere/dz', state%dz)
    call h%write('atmosphere/temperature', state%temperature)
    call h%write('atmosphere/edd', state%edd)
    call h%write('atmosphere/particle_radius', state%particle_radius)
    call h%write('atmosphere/usol', state%usol)
    call h%write('press_temp_edd_profile/enabled', merge(1,0,state%press_temp_edd_profile%enabled))
    call h%write('press_temp_edd_profile/mode', state%press_temp_edd_profile%mode)
    call h%write('press_temp_edd_profile/hydro_pressure', merge(1,0,state%press_temp_edd_profile%hydro_pressure))
    call h%write('press_temp_edd_profile/trop_p', state%press_temp_edd_profile%trop_p)
    call h%write('press_temp_edd_profile/temperature_tol', state%press_temp_edd_profile%temperature_tol)
    call h%write('press_temp_edd_profile/edd_tol', state%press_temp_edd_profile%edd_tol)
    call h%write('press_temp_edd_profile/extreme_factor', state%press_temp_edd_profile%extreme_factor)
    call write_optional_real(h, 'press_temp_edd_profile/pressure', state%press_temp_edd_profile%pressure)
    call write_optional_real(h, 'press_temp_edd_profile/temperature', state%press_temp_edd_profile%temperature)
    call write_optional_real(h, 'press_temp_edd_profile/edd', state%press_temp_edd_profile%edd)
    call h%write('toa_pressure_maintenance/enabled', merge(1,0,state%toa_pressure_maintenance%enabled))
    call h%write('toa_pressure_maintenance/target_pressure', state%toa_pressure_maintenance%target_pressure)
    call h%write('toa_pressure_maintenance/pressure_factor', state%toa_pressure_maintenance%pressure_factor)
    call h%write('toa_pressure_maintenance/extreme_pressure_factor', &
      state%toa_pressure_maintenance%extreme_pressure_factor)
    call write_species_names(h, state%species_names, err)
  end subroutine


  subroutine read_state(h, state, err)
    type(hdf5_file), intent(in) :: h
    type(EvoAtmosphereState), target, intent(out) :: state
    character(:), allocatable, intent(out) :: err
    integer :: version

    call read_integer(h, 'schema_version', version, err)
    if (allocated(err)) return
    if (version /= state_schema_version) then
      err = 'Unsupported atmospheric state HDF5 schema version.'
      return
    endif
    call read_species_names(h, state%species_names, err)
    call read_integer(h, 'compatibility/nq', state%nq, err)
    call read_integer(h, 'compatibility/ng', state%ng, err)
    call read_integer(h, 'compatibility/nsl', state%nsl, err)
    call read_real_scalar(h, 'compatibility/planet_mass', state%planet_mass, err)
    call read_real_scalar(h, 'atmosphere/bottom_atmos', state%bottom_atmos, err)
    call read_real_scalar(h, 'atmosphere/top_atmos', state%top_atmos, err)
    call read_real_scalar(h, 'atmosphere/trop_alt', state%trop_alt, err)
    call read_real_scalar(h, 'atmosphere/planet_radius', state%planet_radius, err)
    call read_real_array(h, 'atmosphere/z', state%z, err)
    call read_real_array(h, 'atmosphere/dz', state%dz, err)
    call read_real_array(h, 'atmosphere/temperature', state%temperature, err)
    call read_real_array(h, 'atmosphere/edd', state%edd, err)
    call read_real_array(h, 'atmosphere/particle_radius', state%particle_radius, err)
    call read_real_array(h, 'atmosphere/usol', state%usol, err)
    call read_flag(h, 'press_temp_edd_profile/enabled', state%press_temp_edd_profile%enabled, err)
    call read_integer(h, 'press_temp_edd_profile/mode', state%press_temp_edd_profile%mode, err)
    call read_flag(h, 'press_temp_edd_profile/hydro_pressure', state%press_temp_edd_profile%hydro_pressure, err)
    call read_real_scalar(h, 'press_temp_edd_profile/trop_p', state%press_temp_edd_profile%trop_p, err)
    call read_real_scalar(h, 'press_temp_edd_profile/temperature_tol', &
      state%press_temp_edd_profile%temperature_tol, err)
    call read_real_scalar(h, 'press_temp_edd_profile/edd_tol', state%press_temp_edd_profile%edd_tol, err)
    call read_real_scalar(h, 'press_temp_edd_profile/extreme_factor', state%press_temp_edd_profile%extreme_factor, err)
    call read_optional_real(h, 'press_temp_edd_profile/pressure', state%press_temp_edd_profile%pressure, err)
    call read_optional_real(h, 'press_temp_edd_profile/temperature', state%press_temp_edd_profile%temperature, err)
    call read_optional_real(h, 'press_temp_edd_profile/edd', state%press_temp_edd_profile%edd, err)
    call read_flag(h, 'toa_pressure_maintenance/enabled', state%toa_pressure_maintenance%enabled, err)
    call read_real_scalar(h, 'toa_pressure_maintenance/target_pressure', &
      state%toa_pressure_maintenance%target_pressure, err)
    call read_real_scalar(h, 'toa_pressure_maintenance/pressure_factor', &
      state%toa_pressure_maintenance%pressure_factor, err)
    call read_real_scalar(h, 'toa_pressure_maintenance/extreme_pressure_factor', &
      state%toa_pressure_maintenance%extreme_pressure_factor, err)
  end subroutine


  subroutine check_dataset(h, path, ndims, dtype, err)
    type(hdf5_file), intent(in) :: h
    character(*), intent(in) :: path
    integer, intent(in) :: ndims, dtype
    character(:), allocatable, intent(inout) :: err

    if (allocated(err)) return
    if (.not.h%exists(path)) then
      err = 'Missing atmospheric state dataset: '//path
      return
    endif
    if (h%ndims(path) /= ndims) then
      err = 'Incorrect atmospheric state dataset rank: '//path
      return
    endif
    if (h%class(path) /= dtype) err = 'Incorrect atmospheric state dataset type: '//path
  end subroutine

  subroutine read_integer(h, path, value, err)
    type(hdf5_file), intent(in) :: h
    character(*), intent(in) :: path
    integer, intent(out) :: value
    character(:), allocatable, intent(inout) :: err

    call check_dataset(h, path, 0, H5T_INTEGER_F, err)
    if (.not.allocated(err)) call h%read(path, value)
  end subroutine

  subroutine read_flag(h, path, value, err)
    type(hdf5_file), intent(in) :: h
    character(*), intent(in) :: path
    logical, intent(out) :: value
    character(:), allocatable, intent(inout) :: err
    integer :: flag

    call read_integer(h, path, flag, err)
    if (allocated(err)) return
    if (flag /= 0 .and. flag /= 1) then
      err = 'Atmospheric state flag must be 0 or 1: '//path
      return
    endif
    value = flag == 1
  end subroutine

  subroutine read_real_scalar(h, path, value, err)
    type(hdf5_file), intent(in) :: h
    character(*), intent(in) :: path
    real(dp), intent(out) :: value
    character(:), allocatable, intent(inout) :: err

    call check_dataset(h, path, 0, H5T_FLOAT_F, err)
    if (.not.allocated(err)) call h%read(path, value)
  end subroutine

  subroutine read_real_array(h, path, value, err)
    type(hdf5_file), intent(in) :: h
    character(*), intent(in) :: path
    real(dp), allocatable, intent(out) :: value(..)
    character(:), allocatable, intent(inout) :: err
    integer(hsize_t), allocatable :: dims(:)
    integer :: status

    call check_dataset(h, path, rank(value), H5T_FLOAT_F, err)
    if (allocated(err)) return
    call h%shape(path, dims)
    if (any(dims > int(huge(0),hsize_t))) then
      err = 'Atmospheric state dataset dimensions exceed the supported integer range: '//path
      return
    endif
    select rank(value)
    rank(1)
      allocate(value(dims(1)), stat=status)
      if (status == 0) call h%read(path, value)
    rank(2)
      allocate(value(dims(1),dims(2)), stat=status)
      if (status == 0) call h%read(path, value)
    rank default
      status = -1
    end select
    if (status /= 0) err = 'Could not allocate atmospheric state array: '//path
  end subroutine

  subroutine write_optional_real(h, path, value)
    type(hdf5_file), intent(inout) :: h
    character(*), intent(in) :: path
    real(dp), allocatable, intent(in) :: value(:)

    call h%write(path//'_allocated', merge(1,0,allocated(value)))
    if (allocated(value)) call h%write(path, value)
  end subroutine

  subroutine read_optional_real(h, path, value, err)
    type(hdf5_file), intent(in) :: h
    character(*), intent(in) :: path
    real(dp), allocatable, intent(out) :: value(:)
    character(:), allocatable, intent(inout) :: err
    logical :: present_array

    call read_flag(h, path//'_allocated', present_array, err)
    if (allocated(err)) return
    if (present_array) then
      call read_real_array(h, path, value, err)
    elseif (h%exists(path)) then
      err = 'Unexpected dataset for unallocated profile array: '//path
    endif
  end subroutine

  ! h5fortran has scalar strings but no string-array read/write interface.
  subroutine write_species_names(h, names, err)
    type(hdf5_file), intent(in) :: h
    character(s_str_len), target, intent(in) :: names(:)
    character(:), allocatable, intent(out) :: err
    integer(hid_t) :: dtype
    integer :: status

    call make_string_type(dtype, err)
    if (allocated(err)) return
    call write_dataset(h%lid, 'compatibility/species_names', dtype, &
      [int(size(names),hsize_t)], c_loc(names(1)(1:1)), err)
    call h5tclose_f(dtype, status)
    if (status < 0 .and. .not.allocated(err)) err = 'Could not close species string datatype.'
  end subroutine


  subroutine write_dataset(file_id, path, dtype, dims, buffer, err)
    integer(hid_t), intent(in) :: file_id, dtype
    character(*), intent(in) :: path
    integer(hsize_t), intent(in) :: dims(:)
    type(c_ptr), intent(in) :: buffer
    character(:), allocatable, intent(inout) :: err
    integer(hid_t) :: space_id, dataset_id
    integer :: status

    if (allocated(err)) return
    if (size(dims) == 0) then
      call h5screate_f(H5S_SCALAR_F, space_id, status)
    else
      call h5screate_simple_f(size(dims), dims, space_id, status)
    endif
    if (status < 0) then
      err = 'Could not create dataspace: '//path
      return
    endif
    call h5dcreate_f(file_id, path, dtype, space_id, dataset_id, status)
    if (status < 0) then
      err = 'Could not create dataset: '//path
    else
      if (all(dims > 0_hsize_t)) then
        call h5dwrite_f(dataset_id, dtype, buffer, status)
        if (status < 0) err = 'Could not write dataset: '//path
      endif
      call h5dclose_f(dataset_id, status)
      if (status < 0 .and. .not.allocated(err)) err = 'Could not close dataset: '//path
    endif
    call h5sclose_f(space_id, status)
    if (status < 0 .and. .not.allocated(err)) err = 'Could not close dataspace: '//path
  end subroutine

  subroutine close_dataset(dataset_id, path, err)
    integer(hid_t), intent(inout) :: dataset_id
    character(*), intent(in) :: path
    character(:), allocatable, intent(inout) :: err
    integer :: status

    if (dataset_id < 0) return
    call h5dclose_f(dataset_id, status)
    if (status < 0 .and. .not.allocated(err)) err = 'Could not close dataset: '//path
    dataset_id = -1
  end subroutine

  subroutine make_string_type(dtype, err)
    integer(hid_t), intent(out) :: dtype
    character(:), allocatable, intent(inout) :: err
    integer :: status, cleanup_status

    dtype = -1
    if (allocated(err)) return
    call h5tcopy_f(H5T_FORTRAN_S1, dtype, status)
    if (status >= 0) call h5tset_size_f(dtype, int(s_str_len,size_t), status)
    if (status < 0) then
      err = 'Could not construct species string datatype.'
      if (dtype >= 0) call h5tclose_f(dtype, cleanup_status)
      dtype = -1
    endif
  end subroutine

  subroutine read_species_names(h, names, err)
    type(hdf5_file), intent(in) :: h
    character(s_str_len), allocatable, target, intent(out) :: names(:)
    character(:), allocatable, intent(inout) :: err
    integer(hid_t) :: dataset_id, string_type, disk_type
    integer(hsize_t), allocatable :: dims(:)
    integer(size_t) :: bytes
    integer :: status
    logical :: variable_string

    call check_dataset(h, 'compatibility/species_names', 1, H5T_STRING_F, err)
    if (allocated(err)) return
    call h%shape('compatibility/species_names', dims)
    if (dims(1) < 3 .or. dims(1) > int(huge(0),hsize_t)) then
      err = 'Atmospheric state species list has an invalid length.'
      return
    endif
    call h5dopen_f(h%lid, 'compatibility/species_names', dataset_id, status)
    if (status < 0) then
      err = 'Could not open atmospheric state species names.'
      return
    endif
    disk_type = -1
    call h5dget_type_f(dataset_id, disk_type, status)
    if (status >= 0) call h5tis_variable_str_f(disk_type, variable_string, status)
    if (status >= 0) call h5tget_size_f(disk_type, bytes, status)
    if (status < 0) then
      err = 'Could not inspect species string datatype.'
    elseif (variable_string .or. bytes /= int(s_str_len,size_t)) then
      err = 'Atmospheric state species names require fixed-width strings of the supported length.'
    endif
    if (disk_type >= 0) then
      call h5tclose_f(disk_type, status)
      if (status < 0 .and. .not.allocated(err)) err = 'Could not close species string datatype.'
    endif
    if (.not.allocated(err)) then
      allocate(names(dims(1)), stat=status)
      if (status /= 0) then
        err = 'Could not allocate atmospheric state species names.'
      else
        call make_string_type(string_type, err)
        if (.not.allocated(err)) then
          call h5dread_f(dataset_id, string_type, c_loc(names(1)(1:1)), status)
          if (status < 0) err = 'Could not read atmospheric state species names.'
          call h5tclose_f(string_type, status)
          if (status < 0 .and. .not.allocated(err)) err = 'Could not close species string datatype.'
        endif
      endif
    endif
    call close_dataset(dataset_id, 'compatibility/species_names', err)
  end subroutine


end submodule
