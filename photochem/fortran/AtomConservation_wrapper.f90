  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine allocate_atomconservation(ptr) bind(c)
    type(c_ptr), intent(out) :: ptr
    type(AtomConservation), pointer :: con
    allocate(con)
    ptr = c_loc(con)
  end subroutine
  
  subroutine deallocate_atomconservation(ptr) bind(c)
    type(c_ptr), intent(in) :: ptr
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    deallocate(con)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine atomconservation_in_surf_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%in_surf
  end subroutine
  
  subroutine atomconservation_in_top_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%in_top
  end subroutine
  
  subroutine atomconservation_in_dist_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%in_dist
  end subroutine

  subroutine atomconservation_in_other_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%in_other
  end subroutine

  subroutine atomconservation_out_surf_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%out_surf
  end subroutine
  
  subroutine atomconservation_out_top_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%out_top
  end subroutine
  
  subroutine atomconservation_out_rain_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%out_rain
  end subroutine
  
  subroutine atomconservation_out_other_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%out_other
  end subroutine
  
  subroutine atomconservation_net_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%net
  end subroutine
  
  subroutine atomconservation_factor_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(AtomConservation), pointer :: con
    call c_f_pointer(ptr, con)
    val = con%factor
  end subroutine
