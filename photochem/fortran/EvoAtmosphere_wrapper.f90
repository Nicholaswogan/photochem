  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  function allocate_evoatmosphere() result(ptr) bind(c)
    type(c_ptr) :: ptr
    type(EvoAtmosphere), pointer :: pc
    allocate(pc)
    ptr = c_loc(pc)
  end function
  
  subroutine deallocate_evoatmosphere(ptr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    type(EvoAtmosphere), pointer :: pc
    character(:), allocatable :: err_f
    
    call c_f_pointer(ptr, pc)
    deallocate(pc)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! subroutine wrappers  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine evoatmosphere_create_wrapper(ptr, mechanism_file, &
                                        settings_file, flux_file, &
                                        atmosphere_txt, atmosphere_txt_present, &
                                        data_dir, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(in) :: mechanism_file(*)
    character(kind=c_char), intent(in) :: settings_file(*)
    character(kind=c_char), intent(in) :: flux_file(*)
    character(kind=c_char), intent(in) :: atmosphere_txt(*)
    logical(c_bool), intent(in) :: atmosphere_txt_present
    character(kind=c_char), intent(in) :: data_dir(*)
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: mechanism_file_f
    character(len=:), allocatable :: settings_file_f
    character(len=:), allocatable :: flux_file_f
    character(len=:), allocatable :: atmosphere_txt_f
    character(len=:), allocatable :: data_dir_f
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(mechanism_file))::mechanism_file_f)
    allocate(character(len=len_cstring(settings_file))::settings_file_f)
    allocate(character(len=len_cstring(flux_file))::flux_file_f)
    allocate(character(len=len_cstring(data_dir))::data_dir_f)
    
    call copy_string_ctof(mechanism_file, mechanism_file_f)
    call copy_string_ctof(settings_file, settings_file_f)
    call copy_string_ctof(flux_file, flux_file_f)
    call copy_string_ctof(data_dir, data_dir_f)

    if (atmosphere_txt_present) then
      allocate(character(len=len_cstring(atmosphere_txt))::atmosphere_txt_f)
      call copy_string_ctof(atmosphere_txt, atmosphere_txt_f)
      pc = EvoAtmosphere(mechanism_file_f, settings_file_f, flux_file_f, &
                         atmosphere_txt_f, data_dir_f, err_f)
    else
      pc = EvoAtmosphere(mechanism_file_f, settings_file_f, flux_file_f, &
                         data_dir_f, err_f)
    endif
    
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_atmosphere_initialized_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    logical(c_bool), intent(out) :: val
    type(EvoAtmosphere), pointer :: pc

    call c_f_pointer(ptr, pc)
    val = pc%atmosphere_initialized
  end subroutine

  subroutine evoatmosphere_initialize_atmosphere_z_wrapper(ptr, nprofile, z, &
                                                            temperature, edd, &
                                                            surface_pressure, &
                                                            nq, mix, np, &
                                                            particle_radius, &
                                                            err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: nprofile, nq, np
    real(c_double), intent(in) :: z(nprofile), temperature(nprofile), edd(nprofile)
    real(c_double), intent(in) :: surface_pressure
    real(c_double), intent(in) :: mix(nq,nprofile), particle_radius(np,nprofile)
    character(kind=c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc

    call c_f_pointer(ptr, pc)
    call pc%initialize_atmosphere_z(z, temperature, edd, surface_pressure, &
                                    mix, particle_radius, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_initialize_atmosphere_p_wrapper(ptr, nprofile, &
                                                            pressure, &
                                                            temperature, edd, &
                                                            nq, mix, np, &
                                                            particle_radius, &
                                                            persistent, trop_p, &
                                                            trop_p_present, &
                                                            maintain_toa_pressure, &
                                                            maintain_toa_pressure_present, &
                                                            target_pressure, &
                                                            target_pressure_present, &
                                                            err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: nprofile, nq, np
    real(c_double), intent(in) :: pressure(nprofile), temperature(nprofile)
    real(c_double), intent(in) :: edd(nprofile)
    real(c_double), intent(in) :: mix(nq,nprofile), particle_radius(np,nprofile)
    logical(c_bool), intent(in) :: persistent, trop_p_present
    logical(c_bool), intent(in) :: maintain_toa_pressure, maintain_toa_pressure_present
    logical(c_bool), intent(in) :: target_pressure_present
    real(c_double), intent(in) :: trop_p, target_pressure
    character(kind=c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    logical :: persistent_f, maintain_toa_pressure_f

    call c_f_pointer(ptr, pc)
    persistent_f = persistent
    maintain_toa_pressure_f = maintain_toa_pressure
    ! The C ABI carries explicit presence flags for optional Fortran
    ! arguments. Preserve those flags when forwarding to the public method so
    ! invalid maintenance options can be diagnosed rather than silently
    ! ignored.
    if (trop_p_present) then
      if (target_pressure_present) then
        call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                        particle_radius, persistent=persistent_f, &
                                        trop_p=trop_p, &
                                        maintain_toa_pressure=maintain_toa_pressure_f, &
                                        target_pressure=target_pressure, err=err_f)
      elseif (maintain_toa_pressure_present) then
        call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                        particle_radius, persistent=persistent_f, &
                                        trop_p=trop_p, &
                                        maintain_toa_pressure=maintain_toa_pressure_f, &
                                        err=err_f)
      else
        call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                        particle_radius, persistent=persistent_f, &
                                        trop_p=trop_p, err=err_f)
      endif
    elseif (target_pressure_present) then
      if (maintain_toa_pressure_present) then
        call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                        particle_radius, persistent=persistent_f, &
                                        maintain_toa_pressure=maintain_toa_pressure_f, &
                                        target_pressure=target_pressure, err=err_f)
      else
        call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                        particle_radius, persistent=persistent_f, &
                                        target_pressure=target_pressure, err=err_f)
      endif
    elseif (maintain_toa_pressure_present) then
      call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                      particle_radius, persistent=persistent_f, &
                                      maintain_toa_pressure=maintain_toa_pressure_f, &
                                      err=err_f)
    else
      call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                      particle_radius, persistent=persistent_f, &
                                      err=err_f)
    endif
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_dat_get(ptr, ptr1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    type(c_ptr), intent(out) :: ptr1
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    ptr1 = c_loc(pc%dat)
  end subroutine

  subroutine evoatmosphere_var_get(ptr, ptr1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    type(c_ptr), intent(out) :: ptr1
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    ptr1 = c_loc(pc%var)
  end subroutine

  subroutine evoatmosphere_wrk_get(ptr, ptr1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    type(c_ptr), intent(out) :: ptr1
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    ptr1 = c_loc(pc%wrk)
  end subroutine

  subroutine evoatmosphere_prep_atmosphere_wrapper(ptr, nq, nz, usol, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol(nq, nz)
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    call pc%prep_atmosphere(usol, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_out2atmosphere_txt_wrapper(ptr, filename, number_of_decimals, overwrite, clip, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(in) :: filename(*)
    integer(c_int), intent(in) :: number_of_decimals
    logical(c_bool), intent(in) :: overwrite, clip
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: filename_f
    logical :: overwrite_f, clip_f
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(filename))::filename_f)
    call copy_string_ctof(filename, filename_f)
    overwrite_f = overwrite
    clip_f = clip
    
    call pc%out2atmosphere_txt(filename_f, number_of_decimals, overwrite_f, clip_f, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
    
  end subroutine

  subroutine evoatmosphere_gas_fluxes_wrapper(ptr, nq, surf_fluxes, top_fluxes, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: nq
    real(c_double), intent(out) :: surf_fluxes(nq)
    real(c_double), intent(out) :: top_fluxes(nq)
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    call pc%gas_fluxes(surf_fluxes,top_fluxes,err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_set_lower_bc_wrapper(ptr, species, bc_type, vdep, den, press, flux, height, missing, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(in) :: species(*)
    character(kind=c_char), intent(in) :: bc_type(*)
    real(c_double), intent(in) :: vdep
    real(c_double), intent(in) :: den
    real(c_double), intent(in) :: press
    real(c_double), intent(in) :: flux
    real(c_double), intent(in) :: height
    logical(c_bool), intent(in) :: missing
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: species_f
    character(len=:), allocatable :: bc_type_f
    
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(species))::species_f)
    allocate(character(len=len_cstring(bc_type))::bc_type_f)
    
    call copy_string_ctof(species, species_f)
    call copy_string_ctof(bc_type, bc_type_f)

    if (missing) then
      call pc%set_lower_bc(species_f, bc_type_f, err=err_f)
    else
      call pc%set_lower_bc(species_f, bc_type_f, vdep=vdep, den=den, press=press, flux=flux, height=height, err=err_f)
    endif
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  subroutine evoatmosphere_set_upper_bc_wrapper(ptr, species, bc_type, veff, flux, missing, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(in) :: species(*)
    character(kind=c_char), intent(in) :: bc_type(*)
    real(c_double), intent(in) :: veff
    real(c_double), intent(in) :: flux
    logical(c_bool), intent(in) :: missing
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: species_f
    character(len=:), allocatable :: bc_type_f
    
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(species))::species_f)
    allocate(character(len=len_cstring(bc_type))::bc_type_f)
    
    call copy_string_ctof(species, species_f)
    call copy_string_ctof(bc_type, bc_type_f)

    if (missing) then
      call pc%set_upper_bc(species_f, bc_type_f, err=err_f)
    else
      call pc%set_upper_bc(species_f, bc_type_f, veff=veff, flux=flux, err=err_f)
    endif
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_set_rate_fcn_wrapper(ptr, species_c, fcn_c, err) bind(c)
    use photochem_types, only: time_dependent_rate_fcn
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(in) :: species_c(*)
    type(c_funptr), value, intent(in) :: fcn_c
    character(kind=c_char), intent(out) :: err(err_len+1)

    type(EvoAtmosphere), pointer :: pc
    procedure(time_dependent_rate_fcn), pointer :: fcn_f
    character(:), allocatable :: species_f
    character(:), allocatable :: err_f

    ! Cast the void pointer to wrapped object.
    call c_f_pointer(ptr, pc)

    ! Copy c string to f string.
    allocate(character(len=len_cstring(species_c))::species_f)
    call copy_string_ctof(species_c, species_f)

    ! Convert c function pointer to a fortran function pointer (unsafe).
    call c_f_procpointer(fcn_c, fcn_f)

    ! Call the function.
    call pc%set_rate_fcn(species_f, fcn_f, err_f)

    ! Set error, if there is one.
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif 
  end subroutine

  subroutine evoatmosphere_set_temperature_wrapper(ptr, nz, temperature, &
                                                trop_alt, trop_alt_present, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: nz
    real(c_double), intent(in) :: temperature(nz)
    real(c_double), intent(in) :: trop_alt
    logical(c_bool), intent(in) :: trop_alt_present
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    if (trop_alt_present) then
      call pc%set_temperature(temperature, trop_alt, err_f)
    else
      call pc%set_temperature(temperature, err=err_f)
    endif
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
    
  end subroutine

  subroutine evoatmosphere_set_press_temp_edd_wrapper(ptr, P_dim1, P, T_dim1, T, edd_dim1, edd, &
                                                      trop_p, trop_p_present, &
                                                      hydro_pressure, hydro_pressure_present, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: P_dim1
    real(c_double), intent(in) :: P(P_dim1)
    integer(c_int), intent(in) :: T_dim1
    real(c_double), intent(in) :: T(T_dim1)
    integer(c_int), intent(in) :: edd_dim1
    real(c_double), intent(in) :: edd(edd_dim1)
    real(c_double), intent(in) :: trop_p
    logical(c_bool), intent(in) :: trop_p_present
    logical(c_bool), intent(in) :: hydro_pressure
    logical(c_bool), intent(in) :: hydro_pressure_present
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    logical :: hydro_pressure_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)

    hydro_pressure_f = hydro_pressure
    
    if (trop_p_present .and. hydro_pressure_present) then
      call pc%set_press_temp_edd(P, T, edd, trop_p=trop_p, hydro_pressure=hydro_pressure_f, err=err_f)
    elseif (trop_p_present .and. .not.hydro_pressure_present) then
      call pc%set_press_temp_edd(P, T, edd, trop_p=trop_p, err=err_f)
    elseif (.not.trop_p_present .and. hydro_pressure_present) then
      call pc%set_press_temp_edd(P, T, edd, hydro_pressure=hydro_pressure_f, err=err_f)
    else
      call pc%set_press_temp_edd(P, T, edd, err=err_f)
    endif
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
    
  end subroutine

  subroutine evoatmosphere_set_press_temp_edd_profile_wrapper(ptr, P_dim1, P, T_dim1, T, &
                                                      edd_dim1, edd, trop_p, trop_p_present, &
                                                      hydro_pressure, hydro_pressure_present, &
                                                      maintain_toa_pressure, maintain_toa_pressure_present, &
                                                      target_pressure, target_pressure_present, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: P_dim1
    real(c_double), intent(in) :: P(P_dim1)
    integer(c_int), intent(in) :: T_dim1
    real(c_double), intent(in) :: T(T_dim1)
    integer(c_int), intent(in) :: edd_dim1
    real(c_double), intent(in) :: edd(edd_dim1)
    real(c_double), intent(in) :: trop_p
    logical(c_bool), intent(in) :: trop_p_present
    logical(c_bool), intent(in) :: hydro_pressure
    logical(c_bool), intent(in) :: hydro_pressure_present
    logical(c_bool), intent(in) :: maintain_toa_pressure
    logical(c_bool), intent(in) :: maintain_toa_pressure_present
    real(c_double), intent(in) :: target_pressure
    logical(c_bool), intent(in) :: target_pressure_present
    character(kind=c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    logical :: hydro_pressure_f, maintain_toa_pressure_f
    real(c_double) :: target_pressure_f
    type(EvoAtmosphere), pointer :: pc

    call c_f_pointer(ptr, pc)
    hydro_pressure_f = .true.
    if (hydro_pressure_present) hydro_pressure_f = hydro_pressure
    maintain_toa_pressure_f = .true.
    if (maintain_toa_pressure_present) maintain_toa_pressure_f = maintain_toa_pressure
    target_pressure_f = 0.1_c_double
    if (target_pressure_present) target_pressure_f = target_pressure

    ! The non-tropopause optionals have explicit defaults, so they can be
    ! passed uniformly while preserving the meaningful presence semantics of
    ! trop_p (which controls whether a tropopause is configured).
    if (trop_p_present) then
      call pc%set_press_temp_edd_profile(P, T, edd, trop_p=trop_p, &
                                         hydro_pressure=hydro_pressure_f, &
                                         maintain_toa_pressure=maintain_toa_pressure_f, &
                                         target_pressure=target_pressure_f, err=err_f)
    else
      call pc%set_press_temp_edd_profile(P, T, edd, &
                                         hydro_pressure=hydro_pressure_f, &
                                         maintain_toa_pressure=maintain_toa_pressure_f, &
                                         target_pressure=target_pressure_f, err=err_f)
    endif
    err(1) = c_null_char
    if (allocated(err_f)) call copy_string_ftoc(err_f, err)

  end subroutine

  subroutine evoatmosphere_clear_press_temp_edd_profile_wrapper(ptr, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc

    call c_f_pointer(ptr, pc)
    call pc%clear_press_temp_edd_profile(err_f)
    err(1) = c_null_char
    if (allocated(err_f)) call copy_string_ftoc(err_f, err)

  end subroutine
  
  subroutine evoatmosphere_update_vertical_grid_wrapper(ptr, TOA_alt, TOA_alt_present, &
                                                     TOA_pressure, TOA_pressure_present, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: TOA_alt
    logical(c_bool), intent(in) :: TOA_alt_present
    real(c_double), intent(in) :: TOA_pressure
    logical(c_bool), intent(in) :: TOA_pressure_present
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    if (TOA_alt_present .and. .not.TOA_pressure_present) then
      call pc%update_vertical_grid(TOA_alt=TOA_alt, err=err_f)
    elseif (.not.TOA_alt_present .and. TOA_pressure_present) then
      call pc%update_vertical_grid(TOA_pressure=TOA_pressure, err=err_f)
    elseif (TOA_alt_present .and. TOA_pressure_present) then
      call pc%update_vertical_grid(TOA_alt=TOA_alt, TOA_pressure=TOA_pressure, err=err_f)
    elseif (.not.TOA_alt_present .and. .not.TOA_pressure_present) then
      call pc%update_vertical_grid(err=err_f)
    endif

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
    
  end subroutine

  subroutine evoatmosphere_regrid_prep_atmosphere_wrapper(ptr, nq, nz, usol, top_atmos, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol(nq, nz)
    real(c_double), intent(in) :: top_atmos
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    call pc%regrid_prep_atmosphere(usol, top_atmos, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  subroutine evoatmosphere_evolve_wrapper(ptr, filename, tstart, nq, nz, usol, nt, t_eval, overwrite, restart_from_file, &
                                          success, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(in) :: filename(*)
    real(c_double), intent(inout) :: tstart
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(inout) :: usol(nq, nz)
    integer(c_int), intent(in) :: nt
    real(c_double), intent(in) :: t_eval(nt)
    logical(c_bool), intent(in) :: overwrite
    logical(c_bool), intent(in) :: restart_from_file
    logical(c_bool), intent(out) :: success
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    logical :: overwrite_f, success_f, restart_from_file_f
    character(len=:), allocatable :: filename_f
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(filename))::filename_f)
    call copy_string_ctof(filename, filename_f)
    overwrite_f = overwrite
    restart_from_file_f = restart_from_file
    
    success_f = pc%evolve(filename_f, tstart, usol, t_eval, overwrite=overwrite_f, &
                          restart_from_file=restart_from_file_f, err=err_f)
    success = success_f
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_check_for_convergence_wrapper(ptr, converged, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    logical(c_bool), intent(out) :: converged
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    logical :: converged_f
  
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
  
    converged_f = pc%check_for_convergence(err_f)
    converged = converged_f

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_initialize_stepper_wrapper(ptr, nq, nz, usol_start, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol_start(nq, nz)
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
  
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    
    call pc%initialize_stepper(usol_start, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  function evoatmosphere_step_wrapper(ptr, err) result(tn) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(len=c_char), intent(out) :: err(err_len+1)
    real(c_double) :: tn
  
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
  
    tn = pc%step(err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end function
  
  subroutine evoatmosphere_destroy_stepper_wrapper(ptr, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
  
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
  
    call pc%destroy_stepper(err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_initialize_robust_stepper_wrapper(ptr, nq, nz, usol_start, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol_start(nq, nz)
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
  
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    
    call pc%initialize_robust_stepper(usol_start, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine evoatmosphere_robust_step_wrapper(ptr, give_up, converged, err) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    logical(c_bool), intent(out) :: give_up, converged
    character(len=c_char), intent(out) :: err(err_len+1)
    
    logical :: give_up_f, converged_f
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc

    call c_f_pointer(ptr, pc)
  
    call pc%robust_step(give_up_f, converged_f, err_f)

    give_up = give_up_f
    converged = converged_f

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  subroutine evoatmosphere_production_and_loss_wrapper(ptr, species, nq, nz, usol, pl_ptr, err) bind(c)
    use photochem, only: ProductionLoss
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(in) :: species(*)
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol(nq, nz)
    type(c_ptr), intent(out) :: pl_ptr
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: species_f
    character(:), allocatable :: err_f
    type(EvoAtmosphere), pointer :: pc
    type(ProductionLoss), pointer :: pl
    
    call c_f_pointer(ptr, pc)
    allocate(pl)
    
    allocate(character(len=len_cstring(species))::species_f)
    call copy_string_ctof(species, species_f)
    
    call pc%production_and_loss(species_f, usol, pl, err_f)
    if (allocated(err_f)) then
      deallocate(pl)
    else
      pl_ptr = c_loc(pl)
    endif
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine evoatmosphere_p_top_min_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    val = pc%P_top_min
  end subroutine

  subroutine evoatmosphere_p_top_min_set(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    pc%P_top_min = val
  end subroutine

  subroutine evoatmosphere_p_top_max_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    val = pc%P_top_max
  end subroutine

  subroutine evoatmosphere_p_top_max_set(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    pc%P_top_max = val
  end subroutine

  subroutine evoatmosphere_top_atmos_adjust_frac_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    val = pc%top_atmos_adjust_frac
  end subroutine

  subroutine evoatmosphere_top_atmos_adjust_frac_set(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(EvoAtmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    pc%top_atmos_adjust_frac = val
  end subroutine
