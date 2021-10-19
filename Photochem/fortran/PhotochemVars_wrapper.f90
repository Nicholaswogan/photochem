module PhotochemVars_wrapper
  use photochem_const, only: s_str_len
  use photochem_types, only: PhotochemVars
  use iso_c_binding
  implicit none
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine allocate_photochemvars(ptr) bind(c)
    type(c_ptr), intent(out) :: ptr
    type(PhotochemVars), pointer :: dat
    allocate(dat)
    ptr = c_loc(dat)
  end subroutine
  
  subroutine deallocate_photochemvars(ptr) bind(c)
    type(c_ptr), intent(in) :: ptr
    type(PhotochemVars), pointer :: dat
    call c_f_pointer(ptr, dat)
    deallocate(dat)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine photochemvars_nz_get(ptr, nz) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: nz
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    nz = var%nz
  end subroutine
  
  subroutine photochemvars_usol_init_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    dim1 = size(var%usol_init,1)
    dim2 = size(var%usol_init,2)
  end subroutine
  
  subroutine photochemvars_usol_init_get(ptr, dim1, dim2, usol_init) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: usol_init(dim1, dim2)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    usol_init = var%usol_init
  end subroutine
  
  
end module