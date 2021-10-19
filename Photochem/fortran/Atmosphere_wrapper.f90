module atmosphere_wrapper
  use photochem, only: Atmosphere, err_len, str_len
  use wrapper_utils, only: copy_string_ftoc, copy_string_ctof, len_cstring
  use iso_c_binding
  implicit none
  
contains
  
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
    character(len=err_len) :: err_f
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
    
    err_f = ""
    call pc%init(data_dir_f, &
                 mechanism_file_f, &
                 settings_file_f, &
                 flux_file_f, &
                 atmosphere_txt_f, &
                 err_f)   
    call copy_string_ftoc(err_f, err)
    dat_ptr = c_loc(pc%dat)
    var_ptr = c_loc(pc%var)
    wrk_ptr = c_loc(pc%wrk) 
  end subroutine
  
  subroutine atmosphere_photochemical_equilibrium_wrapper(ptr, success, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    logical(4), intent(out) :: success
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(len=err_len) :: err_f
    logical :: success_f
  
    type(Atmosphere), pointer :: pc
    call c_f_pointer(ptr, pc)
  
    call pc%photochemical_equilibrium(success_f, err_f)
    success = success_f
    call copy_string_ftoc(err_f,err)
  end subroutine
  
  subroutine atmosphere_out2atmosphere_txt_wrapper(ptr, filename, overwrite, clip, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(in) :: filename(*)
    logical(4), intent(in) :: overwrite, clip
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=:), allocatable :: filename_f
    logical :: overwrite_f, clip_f
    character(len=err_len) :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    allocate(character(len=len_cstring(filename))::filename_f)
    call copy_string_ctof(filename, filename_f)
    overwrite_f = overwrite
    clip_f = clip
    
    err_f = ''
    call pc%out2atmosphere_txt(filename_f, overwrite_f, clip_f, err_f)
    call copy_string_ftoc(err_f,err)
    
  end subroutine
  
  subroutine atmosphere_out2in_wrapper(ptr, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=err_len) :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    err_f = ''
    call pc%out2in(err_f)
    call copy_string_ftoc(err_f,err)
    
  end subroutine
  
  subroutine atmosphere_surface_fluxes_wrapper(ptr, fluxes, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: fluxes(*)
    character(kind=c_char), intent(out) :: err(err_len+1)
    
    character(len=err_len) :: err_f
    type(Atmosphere), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    err_f = ''
    call pc%surface_fluxes(fluxes(1:pc%dat%nq),err_f)
    call copy_string_ftoc(err_f,err)
  end subroutine
  
end module