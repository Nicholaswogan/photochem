! ClimaRadtranWrk

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! getters and setters !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine climaradtranwrk_fup_a_get_size(ptr, dim1, dim2) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  dim1 = size(w%fup_a,1)
  dim2 = size(w%fup_a,2)
end subroutine

subroutine climaradtranwrk_fup_a_get(ptr, dim1, dim2, arr) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  arr = w%fup_a
end subroutine

subroutine climaradtranwrk_fdn_a_get_size(ptr, dim1, dim2) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  dim1 = size(w%fdn_a,1)
  dim2 = size(w%fdn_a,2)
end subroutine

subroutine climaradtranwrk_fdn_a_get(ptr, dim1, dim2, arr) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  arr = w%fdn_a
end subroutine

subroutine climaradtranwrk_fup_n_get_size(ptr, dim1) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  dim1 = size(w%fup_n,1)
end subroutine

subroutine climaradtranwrk_fup_n_get(ptr, dim1, arr) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  arr = w%fup_n
end subroutine

subroutine climaradtranwrk_fdn_n_get_size(ptr, dim1) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  dim1 = size(w%fdn_n,1)
end subroutine

subroutine climaradtranwrk_fdn_n_get(ptr, dim1, arr) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  arr = w%fdn_n
end subroutine

subroutine climaradtranwrk_amean_get_size(ptr, dim1, dim2) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  dim1 = size(w%amean,1)
  dim2 = size(w%amean,2)
end subroutine

subroutine climaradtranwrk_amean_get(ptr, dim1, dim2, arr) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  arr = w%amean
end subroutine

subroutine climaradtranwrk_tau_band_get_size(ptr, dim1, dim2) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  dim1 = size(w%tau_band,1)
  dim2 = size(w%tau_band,2)
end subroutine

subroutine climaradtranwrk_tau_band_get(ptr, dim1, dim2, arr) bind(c)
  use clima_radtran, only: ClimaRadtranWrk
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(ClimaRadtranWrk), pointer :: w
  call c_f_pointer(ptr, w)
  arr = w%tau_band
end subroutine
  