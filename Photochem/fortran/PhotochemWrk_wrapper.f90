module PhotochemWrk_wrapper
  use photochem_const, only: s_str_len
  use photochem_types, only: PhotochemWrk
  use iso_c_binding
  implicit none
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine allocate_photochemwrk(ptr) bind(c)
    type(c_ptr), intent(out) :: ptr
    type(PhotochemWrk), pointer :: wrk
    allocate(wrk)
    ptr = c_loc(wrk)
  end subroutine
  
  subroutine deallocate_photochemwrk(ptr) bind(c)
    type(c_ptr), intent(in) :: ptr
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    deallocate(wrk)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine photochemwrk_usol_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%usol,1)
    dim2 = size(wrk%usol,2)
  end subroutine
  
  subroutine photochemwrk_usol_get(ptr, dim1, dim2, usol) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: usol(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    usol = wrk%usol
  end subroutine
  
end module