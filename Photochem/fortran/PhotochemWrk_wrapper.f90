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
  
  subroutine photochemwrk_usol_set(ptr, dim1, dim2, usol) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(in) :: usol(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    wrk%usol = usol
  end subroutine
  
  subroutine photochemwrk_pressure_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%pressure,1)
  end subroutine
  
  subroutine photochemwrk_pressure_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%pressure
  end subroutine
  
  subroutine photochemwrk_density_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%density,1)
  end subroutine
  
  subroutine photochemwrk_density_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%density
  end subroutine

  subroutine photochemwrk_mubar_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%mubar,1)
  end subroutine
  
  subroutine photochemwrk_mubar_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%mubar
  end subroutine
  
  subroutine photochemwrk_prates_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%prates,1)
    dim2 = size(wrk%prates,2)
  end subroutine
  
  subroutine photochemwrk_prates_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%prates
  end subroutine
  
  subroutine photochemwrk_amean_grd_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%amean_grd,1)
    dim2 = size(wrk%amean_grd,2)
  end subroutine
  
  subroutine photochemwrk_amean_grd_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%amean_grd
  end subroutine
  
  subroutine photochemwrk_surf_radiance_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%surf_radiance,1)
  end subroutine
  
  subroutine photochemwrk_surf_radiance_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%surf_radiance
  end subroutine
  
end module