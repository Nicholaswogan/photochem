
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine photochemdata_nq_get(ptr, nq) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: nq
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    nq = dat%nq
  end subroutine
  
  subroutine photochemdata_np_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    val = dat%np
  end subroutine
  
  subroutine photochemdata_nsp_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    val = dat%nsp
  end subroutine
  
  subroutine photochemdata_ng_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    val = dat%ng
  end subroutine
  
  subroutine photochemdata_nsl_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    val = dat%nsl
  end subroutine
  
  subroutine photochemdata_nll_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    val = dat%nll
  end subroutine
  
  subroutine photochemdata_nw_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    val = dat%nw
  end subroutine

  subroutine photochemdata_planet_mass_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    val = dat%planet_mass
  end subroutine

  subroutine photochemdata_planet_radius_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    val = dat%planet_radius
  end subroutine

  subroutine photochemdata_planet_radius_set(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dat%planet_radius = val
  end subroutine
  
  subroutine photochemdata_species_names_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%species_names)
  end subroutine
  
  subroutine photochemdata_species_names_get(ptr, dim1, species_names) bind(c)
    type(c_ptr), value, intent(in) :: ptr
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
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%atoms_names)
  end subroutine
  
  subroutine photochemdata_atoms_names_get(ptr, dim1, names) bind(c)
    type(c_ptr), value, intent(in) :: ptr
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
  
  subroutine photochemdata_reaction_equations_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%reaction_equations)
  end subroutine
  
  subroutine photochemdata_reaction_equations_get(ptr, dim1, names) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: names(dim1*m_str_len+1)
    type(PhotochemData), pointer :: dat
      
    integer :: i, j, k
    
    call c_f_pointer(ptr, dat)
    do i = 1,dim1
      do j = 1,m_str_len
        k = j + (i - 1) * m_str_len
        names(k) = dat%reaction_equations(i)(j:j)
      enddo
    enddo
    names(dim1*m_str_len+1) = c_null_char
    
  end subroutine
  
  subroutine photochemdata_photonums_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%photonums,1)
  end subroutine
  
  subroutine photochemdata_photonums_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    integer(c_int), intent(out) :: arr(dim1)
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    arr = dat%photonums
  end subroutine
  
  subroutine photochemdata_wavl_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%wavl,1)
  end subroutine
  
  subroutine photochemdata_wavl_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    arr = dat%wavl
  end subroutine
  
  subroutine photochemdata_species_mass_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%species_mass,1)
  end subroutine
  
  subroutine photochemdata_species_mass_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    arr = dat%species_mass
  end subroutine
  
  subroutine photochemdata_species_redox_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = size(dat%species_redox,1)
  end subroutine
  
  subroutine photochemdata_species_redox_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    arr = dat%species_redox
  end subroutine
  
  subroutine photochemdata_particle_sat_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    dim1 = 0
    if (.not.allocated(dat%particle_sat)) return
    dim1 = size(dat%particle_sat,1)
  end subroutine

  subroutine photochemdata_particle_sat_get(ptr, dim1, ptr1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    type(c_ptr), intent(out) :: ptr1(dim1)
    integer :: i
    type(PhotochemData), pointer :: dat
    call c_f_pointer(ptr, dat)
    do i = 1,dim1
      ptr1(i) = c_loc(dat%particle_sat(i))
    enddo
  end subroutine