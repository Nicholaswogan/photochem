module PhotochemData_wrapper
  use photochem_const, only: s_str_len
  use photochem_types, only: PhotochemData
  use iso_c_binding
  implicit none
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine allocate_photochemdata(ptr) bind(c)
    type(c_ptr), intent(out) :: ptr
    type(PhotochemData), pointer :: dat
    allocate(dat)
    ptr = c_loc(dat)
  end subroutine
  
  subroutine deallocate_photochemdata(ptr) bind(c)
    type(c_ptr), intent(in) :: ptr
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    deallocate(dat)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine photochemdata_nq_get(ptr, nq) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: nq
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    nq = dat%nq
  end subroutine
  
  subroutine photochemdata_species_names_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%species_names)
  end subroutine
  
  subroutine photochemdata_species_names_get(ptr, dim1, species_names) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: species_names(dim1*s_str_len+1)
    type(PhotochemData), pointer :: dat
    
    
    integer :: i, j, k
    
    call c_f_pointer(ptr, dat)
    do i = 1,dim1
      do j = 1,s_str_len
        k = j + (i - 1) * s_str_len
        species_names(k) = dat%species_names(i)(j:j)
      enddo
    enddo
    species_names(dim1*s_str_len+1) = c_null_char
    
  end subroutine
  
  subroutine photochemdata_atoms_names_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%atoms_names)
  end subroutine
  
  subroutine photochemdata_atoms_names_get(ptr, dim1, names) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: names(dim1*s_str_len+1)
    type(PhotochemData), pointer :: dat
    
    
    integer :: i, j, k
    
    call c_f_pointer(ptr, dat)
    do i = 1,dim1
      do j = 1,s_str_len
        k = j + (i - 1) * s_str_len
        names(k) = dat%atoms_names(i)(j:j)
      enddo
    enddo
    names(dim1*s_str_len+1) = c_null_char
    
  end subroutine
  
  
  ! subroutine copy_string_f2d_to_c()
  
  
  
end module