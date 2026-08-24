cdef extern void futils_rebin_wrapper(
  int *old_bins_len, double *old_bins, 
  int *old_vals_len, double *old_vals, 
  int *new_bins_len, double *new_bins, 
  int *new_vals_len, double *new_vals, 
  int *ierr
)

cdef extern void futils_rebin_with_errors_wrapper(
  int *old_bins_len, double *old_bins, 
  int *old_vals_len, double *old_vals, 
  int *old_errs_len, double *old_errs, 
  int *new_bins_len, double *new_bins, 
  int *new_vals_len, double *new_vals, 
  int *new_errs_len, double *new_errs, 
  int *ierr
)

cdef extern void futils_rebin_error_message_wrapper(
  int *ierr, char *err
)
