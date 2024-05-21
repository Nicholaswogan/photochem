  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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

! CondensationParameters

!~~ Getters and setters ~~!
subroutine condensationparameters_k_cond_get(ptr, val) bind(c)
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%k_cond
end subroutine

subroutine condensationparameters_k_cond_set(ptr, val) bind(c)
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  c%k_cond = val
end subroutine

subroutine condensationparameters_k_evap_get(ptr, val) bind(c)
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%k_evap
end subroutine

subroutine condensationparameters_k_evap_set(ptr, val) bind(c)
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  c%k_evap = val
end subroutine

subroutine condensationparameters_rhc_get(ptr, val) bind(c)
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%RHc
end subroutine

subroutine condensationparameters_rhc_set(ptr, val) bind(c)
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  c%RHc = val
end subroutine

subroutine condensationparameters_smooth_factor_get(ptr, val) bind(c)
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%smooth_factor
end subroutine

subroutine condensationparameters_smooth_factor_set(ptr, val) bind(c)
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  c%smooth_factor = val
end subroutine
