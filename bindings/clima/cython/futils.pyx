cimport futils_pxd as f_pxd

cdef void rebin_error_message(int ierr):

  cdef char err[ERR_LEN+1]

  if ierr != 0:
    f_pxd.futils_rebin_error_message_wrapper(
      &ierr, 
      err
    )
    if len(err.strip()) > 0:
      raise Exception(err.decode().strip())

cpdef rebin(ndarray[double, ndim=1] old_bins, ndarray[double, ndim=1]  old_vals, ndarray[double, ndim=1]  new_bins):
  """Rebins `old_vals` defined on `old_bins` to `new_bins`. An example is
  rebinning a high resolution spectra of infrared emission of Earth 
  to a lower resolution. I have optimized the routine for downbinning
  data. Upbinning seems like an unlikely application.

  Parameters
  ----------
  old_bins : ndarray[double,ndim=1]
      Edges of bins for which old_vals are defined
  old_vals : ndarray[double,ndim=1]
      Values defined on old_bins.
  new_bins : ndarray[double,ndim=1]
      Edges of target bin that you want to rebin to.

  Returns
  -------
  new_vals : ndarray[double,ndim=1]
      Values defined on new_bins.
  """

  cdef int old_bins_len = old_bins.shape[0]
  cdef int old_vals_len = old_vals.shape[0]
  cdef int new_bins_len = new_bins.shape[0]
  cdef int new_vals_len = new_bins_len - 1
  cdef ndarray new_vals = np.empty(new_vals_len, np.double)
  cdef int ierr

  f_pxd.futils_rebin_wrapper(
    &old_bins_len, <double *> old_bins.data, 
    &old_vals_len, <double *> old_vals.data, 
    &new_bins_len, <double *> new_bins.data, 
    &new_vals_len, <double *> new_vals.data, 
    &ierr
  )
  
  rebin_error_message(ierr)

  return new_vals

cpdef rebin_with_errors(ndarray[double, ndim=1] old_bins, ndarray[double, ndim=1]  old_vals, ndarray[double, ndim=1]  old_errs, ndarray[double, ndim=1]  new_bins):
  """Rebins `old_vals` and `old_errs` defined on `old_bins` to `new_bins`. This function has
  the same behavior as `rebin`, except it also rebins errors.

  Parameters
  ----------
  old_bins : ndarray[double,ndim=1]
      Edges of bins for which old_vals are defined
  old_vals : ndarray[double,ndim=1]
      Values defined on old_bins.
  old_errs : ndarray[double,ndim=1]
      Standard deviations defined on old_bins.
  new_bins : ndarray[double,ndim=1]
      Edges of target bin that you want to rebin to.

  Returns
  -------
  new_vals : ndarray[double,ndim=1]
      Values defined on new_bins.
  new_errs : ndarray[double,ndim=1]
      Standard deviations defined on new_bins.
  """

  cdef int old_bins_len = old_bins.shape[0]
  cdef int old_vals_len = old_vals.shape[0]
  cdef int old_errs_len = old_errs.shape[0]
  cdef int new_bins_len = new_bins.shape[0]
  cdef int new_vals_len = new_bins_len - 1
  cdef ndarray new_vals = np.empty(new_vals_len, np.double)
  cdef int new_errs_len = new_bins_len - 1
  cdef ndarray new_errs = np.empty(new_errs_len, np.double)
  cdef int ierr

  f_pxd.futils_rebin_with_errors_wrapper(
    &old_bins_len, <double *> old_bins.data, 
    &old_vals_len, <double *> old_vals.data, 
    &old_errs_len, <double *> old_errs.data, 
    &new_bins_len, <double *> new_bins.data, 
    &new_vals_len, <double *> new_vals.data, 
    &new_errs_len, <double *> new_errs.data, 
    &ierr
  )

  rebin_error_message(ierr)

  return new_vals, new_errs