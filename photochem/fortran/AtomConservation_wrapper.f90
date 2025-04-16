! CondensationParameters

!~~ Getters and setters ~~!
subroutine condensationparameters_k_cond_get(ptr, val) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%k_cond
end subroutine

subroutine condensationparameters_k_cond_set(ptr, val) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  c%k_cond = val
end subroutine

subroutine condensationparameters_k_evap_get(ptr, val) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%k_evap
end subroutine

subroutine condensationparameters_k_evap_set(ptr, val) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  c%k_evap = val
end subroutine

subroutine condensationparameters_rhc_get(ptr, val) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%RHc
end subroutine

subroutine condensationparameters_rhc_set(ptr, val) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  c%RHc = val
end subroutine

subroutine condensationparameters_smooth_factor_get(ptr, val) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%smooth_factor
end subroutine

subroutine condensationparameters_smooth_factor_set(ptr, val) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(CondensationParameters), pointer :: c
  call c_f_pointer(ptr, c)
  c%smooth_factor = val
end subroutine

! SaturationData

subroutine saturationdata_sat_pressure_wrapper(ptr, T, Psat) bind(c)
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: T
  real(c_double), intent(out) :: Psat
  type(SaturationData), pointer :: sat
  call c_f_pointer(ptr, sat)
  Psat = sat%sat_pressure(T)
end subroutine