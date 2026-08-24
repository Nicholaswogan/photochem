! RTChannel

subroutine rtchannel_wavl_get_size(ptr, dim1) bind(c)
  use clima_radtran_types, only: RTChannel
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(RTChannel), pointer :: rtc
  call c_f_pointer(ptr, rtc)
  dim1 = size(rtc%wavl,1)
end subroutine

subroutine rtchannel_wavl_get(ptr, dim1, arr) bind(c)
  use clima_radtran_types, only: RTChannel
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(RTChannel), pointer :: rtc
  call c_f_pointer(ptr, rtc)
  arr = rtc%wavl
end subroutine

subroutine rtchannel_freq_get_size(ptr, dim1) bind(c)
  use clima_radtran_types, only: RTChannel
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(RTChannel), pointer :: rtc
  call c_f_pointer(ptr, rtc)
  dim1 = size(rtc%freq,1)
end subroutine

subroutine rtchannel_freq_get(ptr, dim1, arr) bind(c)
  use clima_radtran_types, only: RTChannel
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(RTChannel), pointer :: rtc
  call c_f_pointer(ptr, rtc)
  arr = rtc%freq
end subroutine