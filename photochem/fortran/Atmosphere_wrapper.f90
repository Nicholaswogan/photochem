module atmosphere_wrapper
  use photochem, only: Atmosphere, err_len, str_len
  use wrapper_utils, only: copy_string_ftoc, copy_string_ctof, len_cstring
  use iso_c_binding
  implicit none
  
contains
  
  !!!!!!!!!!!!!!!
  !!! version !!!
  !!!!!!!!!!!!!!!
  
  subroutine photochem_version_get(version_c) bind(c)
    use photochem, only: version 
    character(kind=c_char), intent(out) :: version_c(100+1)
    call copy_string_ftoc(version, version_c)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine allocate_atmosphere(ptr) bind(c)
    type(c_ptr), intent(out) :: ptr
    type(Atmosphere), pointer :: pc
    allocate(pc)
    ptr = c_loc(pc)
  end subroutine
  
  subroutine deallocate_atmosphere(ptr) bind(c)
    type(c_ptr), intent(in) :: ptr
    type(Atmosphere), pointer :: pc
    character(:), allocatable :: err_f
    call c_f_pointer(ptr, pc)
    deallocate(pc)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! subroutine wrappers  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine atmosphere_init_wrapper(ptr, data_dir, mechanism_file, &
                                     settings_file, flux_file, &
                                     atmosphere_txt, dat_ptr, &
                                     var_ptr, wrk_ptr , err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(in) :: data_dir(*)
    character(kind=c_char), intent(in) :: mechanism_file(*)
    character(kind=c_char), intent(in) :: settings_file(*)
    character(kind=c_char), intent(in) :: flux_file(*)
    character(kind=c_char), intent(in) :: atmosphere_txt(*)
    type(c_ptr), intent(out) :: dat_ptr, var_ptr, wrk_ptr
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: data_dir_f
    character(len=:), allocatable :: mechanism_file_f
    character(len=:), allocatable :: settings_file_f
    character(len=:), allocatable :: flux_file_f
    character(len=:), allocatable :: atmosphere_txt_f
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(data_dir))::data_dir_f)
    allocate(character(len=len_cstring(mechanism_file))::mechanism_file_f)
    allocate(character(len=len_cstring(settings_file))::settings_file_f)
    allocate(character(len=len_cstring(flux_file))::flux_file_f)
    allocate(character(len=len_cstring(atmosphere_txt))::atmosphere_txt_f)
    
    call copy_string_ctof(data_dir, data_dir_f)
    call copy_string_ctof(mechanism_file, mechanism_file_f)
    call copy_string_ctof(settings_file, settings_file_f)
    call copy_string_ctof(flux_file, flux_file_f)
    call copy_string_ctof(atmosphere_txt, atmosphere_txt_f)
    
    call pc%init(data_dir_f, &
                 mechanism_file_f, &
                 settings_file_f, &
                 flux_file_f, &
                 atmosphere_txt_f, &
                 err_f)
    
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
    dat_ptr = c_loc(pc%dat)
    var_ptr = c_loc(pc%var)
    wrk_ptr = c_loc(pc%wrk) 
  end subroutine
  
  subroutine atmosphere_photochemical_equilibrium_wrapper(ptr, success, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    logical(c_bool), intent(out) :: success
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    logical :: success_f
  
    type(Atmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
  
    call pc%photochemical_equilibrium(success_f, err_f)
    success = success_f

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  subroutine atmosphere_out2atmosphere_txt_wrapper(ptr, filename, overwrite, clip, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(in) :: filename(*)
    logical(c_bool), intent(in) :: overwrite, clip
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: filename_f
    logical :: overwrite_f, clip_f
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(filename))::filename_f)
    call copy_string_ctof(filename, filename_f)
    overwrite_f = overwrite
    clip_f = clip
    
    call pc%out2atmosphere_txt(filename_f, overwrite_f, clip_f, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
    
  end subroutine
  
  subroutine atmosphere_out2in_wrapper(ptr, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    call pc%out2in(err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
    
  end subroutine
  
  subroutine atmosphere_gas_fluxes_wrapper(ptr, surf_fluxes, top_fluxes, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: surf_fluxes(*)
    real(c_double), intent(out) :: top_fluxes(*)
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    call pc%gas_fluxes(surf_fluxes(1:pc%dat%nq),top_fluxes(1:pc%dat%nq),err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  subroutine atmosphere_set_lower_bc_wrapper(ptr, species, bc_type, vdep, mix, flux, height, missing, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(in) :: species(*)
    character(kind=c_char), intent(in) :: bc_type(*)
    real(c_double), intent(in) :: vdep
    real(c_double), intent(in) :: mix
    real(c_double), intent(in) :: flux
    real(c_double), intent(in) :: height
    logical(c_bool), intent(in) :: missing
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: species_f
    character(len=:), allocatable :: bc_type_f
    
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(species))::species_f)
    allocate(character(len=len_cstring(bc_type))::bc_type_f)
    
    call copy_string_ctof(species, species_f)
    call copy_string_ctof(bc_type, bc_type_f)

    if (missing) then
      call pc%set_lower_bc(species_f, bc_type_f, err=err_f)
    else
      call pc%set_lower_bc(species_f, bc_type_f, vdep=vdep, mix=mix, flux=flux, height=height, err=err_f)
    endif
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  subroutine atmosphere_set_upper_bc_wrapper(ptr, species, bc_type, veff, flux, missing, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(in) :: species(*)
    character(kind=c_char), intent(in) :: bc_type(*)
    real(c_double), intent(in) :: veff
    real(c_double), intent(in) :: flux
    logical(c_bool), intent(in) :: missing
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: species_f
    character(len=:), allocatable :: bc_type_f
    
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
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
  
  subroutine atmosphere_initialize_stepper_wrapper(ptr, nq, nz, usol_start, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol_start(nq, nz)
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
  
    type(Atmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
    
    call pc%initialize_stepper(usol_start, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  function atmosphere_step_wrapper(ptr, err) result(tn) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(len=c_char), intent(out) :: err(err_len+1)
    real(c_double) :: tn
  
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
  
    tn = pc%step(err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end function
  
  subroutine atmosphere_destroy_stepper_wrapper(ptr, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
  
    type(Atmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
  
    call pc%destroy_stepper(err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  subroutine atmosphere_production_and_loss_wrapper(ptr, species, nq, nz, usol, pl_ptr, err) bind(c)
    use photochem, only: ProductionLoss
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(in) :: species(*)
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol(nq, nz)
    type(c_ptr), intent(out) :: pl_ptr
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: species_f
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
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
  
  subroutine atmosphere_prep_atmosphere_wrapper(ptr, nq, nz, usol, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol(nq, nz)
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    
    call pc%prep_atmosphere(usol, err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  subroutine atmosphere_redox_conservation_wrapper(ptr, redox_factor, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: redox_factor
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    type(Atmosphere), pointer :: pc
    character(:), allocatable :: err_f
    
    call c_f_pointer(ptr, pc)
    
    redox_factor = pc%redox_conservation(err_f)
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
    
  end subroutine
  
  subroutine atmosphere_atom_conservation_wrapper(ptr, atom, con_ptr, err) bind(c)
    use photochem_types, only: AtomConservation
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(in) :: atom(*)
    type(c_ptr), intent(out) :: con_ptr
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    type(Atmosphere), pointer :: pc
    type(AtomConservation), pointer :: con
    character(:), allocatable :: err_f
    character(len=:), allocatable :: atom_f
    
    call c_f_pointer(ptr, pc)
    allocate(con)
    
    allocate(character(len=len_cstring(atom))::atom_f)
    call copy_string_ctof(atom, atom_f)
    
    
    con = pc%atom_conservation(atom_f, err_f)
    if (allocated(err_f)) then
      deallocate(con)
    else
      con_ptr = c_loc(con)
    endif

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif   
  end subroutine
  
  subroutine atmosphere_evolve_wrapper(ptr, filename, tstart, nq, nz, usol, nt, t_eval, overwrite, success, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(in) :: filename(*)
    real(c_double), intent(in) :: tstart
    integer(c_int), intent(in) :: nq, nz
    real(c_double), intent(in) :: usol(nq, nz)
    integer(c_int), intent(in) :: nt
    real(c_double), intent(in) :: t_eval(nt)
    logical(c_bool), intent(in) :: overwrite
    logical(c_bool), intent(out) :: success
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    logical :: overwrite_f, success_f
    character(len=:), allocatable :: filename_f
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(filename))::filename_f)
    call copy_string_ctof(filename, filename_f)
    overwrite_f = overwrite
    
    
    success_f = pc%evolve(filename_f, tstart, usol, t_eval, overwrite=overwrite_f, err=err_f)
    success = success_f
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine
  
  subroutine atmosphere_set_temperature_wrapper(ptr, nz, temperature, &
                                                trop_alt, trop_alt_present, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: nz
    real(c_double), intent(in) :: temperature(nz)
    real(c_double), intent(in) :: trop_alt
    logical(c_bool), intent(in) :: trop_alt_present
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(:), allocatable :: err_f
    type(Atmosphere), pointer :: pc
    
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

  subroutine atmosphere_set_photon_flux_fcn_wrapper(ptr, photon_flux_fcn_c) bind(c)
    use photochem_types, only: time_dependent_flux_fcn
    type(c_ptr), intent(in) :: ptr
    type(c_funptr), intent(in) :: photon_flux_fcn_c

    procedure(time_dependent_flux_fcn), pointer :: photon_flux_fcn_f

    type(Atmosphere), pointer :: pc

    call c_f_pointer(ptr, pc)
    call c_f_procpointer(photon_flux_fcn_c, photon_flux_fcn_f)
    call pc%set_photon_flux_fcn(photon_flux_fcn_f)

  end subroutine

end module