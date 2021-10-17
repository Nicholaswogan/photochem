module atmos_wrapper
  use iso_c_binding
  implicit none
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine allocate_photochem(ptr) bind(c)
    use Atmos, only: Photochem
    type(c_ptr), intent(out) :: ptr
    type(Photochem), pointer :: pc
    allocate(pc)
    ptr = c_loc(pc)
  end subroutine
  
  subroutine deallocate_photochem(ptr) bind(c)
    use Atmos, only: Photochem
    type(c_ptr), intent(in) :: ptr
    type(Photochem), pointer :: pc
    call c_f_pointer(ptr, pc)
    deallocate(pc)
  end subroutine
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! subroutine wrappers  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine photochem_init_wrapper(ptr, data_dir, mechanism_file, &
                                    settings_file, flux_file, &
                                    atmosphere_txt, err) bind(c)
    use Atmos, only: Photochem, str_len, err_len
    type(c_ptr), intent(in) :: ptr
    character(len=c_char), intent(in) :: data_dir(str_len+1)
    character(len=c_char), intent(in) :: mechanism_file(str_len+1)
    character(len=c_char), intent(in) :: settings_file(str_len+1)
    character(len=c_char), intent(in) :: flux_file(str_len+1)
    character(len=c_char), intent(in) :: atmosphere_txt(str_len+1)
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(len=str_len) :: data_dir_f
    character(len=str_len) :: mechanism_file_f
    character(len=str_len) :: settings_file_f
    character(len=str_len) :: flux_file_f
    character(len=str_len) :: atmosphere_txt_f
    character(len=err_len) :: err_f
    type(Photochem), pointer :: pc
    
    call c_f_pointer(ptr, pc)
    
    call copy_string_ctof(data_dir, data_dir_f)
    call copy_string_ctof(mechanism_file, mechanism_file_f)
    call copy_string_ctof(settings_file, settings_file_f)
    call copy_string_ctof(flux_file, flux_file_f)
    call copy_string_ctof(atmosphere_txt, atmosphere_txt_f)
    
    err_f = ""
    call pc%init(trim(data_dir_f), &
                 trim(mechanism_file_f), &
                 trim(settings_file_f), &
                 trim(flux_file_f), &
                 trim(atmosphere_txt_f), &
                 err_f)
                 
    call copy_string_ftoc(err_f, err)
  end subroutine
  
  subroutine photochem_photochemical_equilibrium_wrapper(ptr, success, err) bind(c)
    use Atmos, only: Photochem, err_len
    type(c_ptr), intent(in) :: ptr
    logical(c_bool), intent(out) :: success
    character(len=c_char), intent(out) :: err(err_len+1)
    
    character(len=err_len) :: err_f
    logical :: success_f
  
    type(Photochem), pointer :: pc
    call c_f_pointer(ptr, pc)
  
    call pc%photochemical_equilibrium(success_f, err_f)
    success = success_f
    call copy_string_ftoc(err_f,err)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!
  !!! Utilities  !!!
  !!!!!!!!!!!!!!!!!!
  
  subroutine copy_string_ctof(stringc,stringf)
  ! utility function to convert c string to fortran string
  character(len=*), intent(out) :: stringf
  character(c_char), intent(in) :: stringc(:)
  integer j
  stringf = ''
  char_loop: do j=1,min(size(stringc),len(stringf))
     if (stringc(j)==c_null_char) exit char_loop
     stringf(j:j) = stringc(j)
  end do char_loop
end subroutine copy_string_ctof

subroutine copy_string_ftoc(stringf,stringc)
  ! utility function to convert c string to fortran string
  character(len=*), intent(in) :: stringf
  character(c_char), intent(out) :: stringc(:)
  integer j,n
  n = len_trim(stringf)   
  do j=1,n    
    stringc(j) = stringf(j:j)   
  end do
  stringc(n+1) = c_null_char
end subroutine copy_string_ftoc
  
end module