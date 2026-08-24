! futils

subroutine futils_rebin_wrapper(old_bins_len, old_bins, old_vals_len, old_vals, &
                                new_bins_len, new_bins, new_vals_len, new_vals, ierr) bind(c)
  use futils, only: rebin
  integer(c_int), intent(in) :: old_bins_len
  real(c_double), intent(in) :: old_bins(old_bins_len)
  integer(c_int), intent(in) :: old_vals_len
  real(c_double), intent(in) :: old_vals(old_vals_len)
  integer(c_int), intent(in) :: new_bins_len
  real(c_double), intent(in) :: new_bins(new_bins_len)
  integer(c_int), intent(in) :: new_vals_len
  real(c_double), intent(out) :: new_vals(new_vals_len)
  integer(c_int), intent(out) :: ierr
  call rebin(old_bins, old_vals, new_bins, new_vals, ierr)
end subroutine

subroutine futils_rebin_with_errors_wrapper(old_bins_len, old_bins, old_vals_len, old_vals, old_errs_len, old_errs, &
                                            new_bins_len, new_bins, new_vals_len, new_vals, new_errs_len, new_errs, ierr) bind(c)
  use futils, only: rebin_with_errors
  integer(c_int), intent(in) :: old_bins_len
  real(c_double), intent(in) :: old_bins(old_bins_len)
  integer(c_int), intent(in) :: old_vals_len
  real(c_double), intent(in) :: old_vals(old_vals_len)
  integer(c_int), intent(in) :: old_errs_len
  real(c_double), intent(in) :: old_errs(old_errs_len)
  integer(c_int), intent(in) :: new_bins_len
  real(c_double), intent(in) :: new_bins(new_bins_len)
  integer(c_int), intent(in) :: new_vals_len
  real(c_double), intent(out) :: new_vals(new_vals_len)
  integer(c_int), intent(in) :: new_errs_len
  real(c_double), intent(out) :: new_errs(new_errs_len)
  integer(c_int), intent(out) :: ierr
  call rebin_with_errors(old_bins, old_vals, old_errs, new_bins, new_vals, new_errs, ierr)
end subroutine

subroutine futils_rebin_error_message_wrapper(ierr, err) bind(c)
  use futils, only: rebin_error_message
  integer(c_int), intent(in) :: ierr
  character(kind=c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f

  err_f = rebin_error_message(ierr)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine
