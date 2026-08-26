cimport futils_pxd as f_pxd

cdef void rebin_error_message(int ierr):

  cdef char err[ERR_LEN+1]

  if ierr != 0:
    f_pxd.futils_rebin_error_message_wrapper(
      &ierr, 
      err
    )
    if len(err.strip()) > 0:
      raise ClimaException(err.decode().strip())

cpdef rebin(ndarray[double, ndim=1] old_bins, ndarray[double, ndim=1]  old_vals, ndarray[double, ndim=1]  new_bins):
  """Rebin piecewise-constant values by overlap-weighted averaging.

  This routine is optimized for down-binning, such as reducing the resolution
  of a spectrum. Both edge arrays must be strictly increasing. Portions of new
  bins outside the old-bin extent receive no contribution.

  Parameters
  ----------
  old_bins : ndarray, shape (n_old + 1,)
      Original bin edges.
  old_vals : ndarray, shape (n_old,)
      Piecewise-constant values in the original bins.
  new_bins : ndarray, shape (n_new + 1,)
      Target bin edges.

  Returns
  -------
  ndarray, shape (n_new,)
      Overlap-weighted values in the target bins.

  Raises
  ------
  ClimaException
      If array lengths are inconsistent or either edge array is not strictly
      increasing.
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
  """Rebin piecewise-constant values and independent standard deviations.

  Values use the same overlap-weighted averaging as ``rebin``. Variances are
  propagated in quadrature. Unlike ``rebin``, the target-bin extent must lie
  within the original-bin extent.

  Parameters
  ----------
  old_bins : ndarray, shape (n_old + 1,)
      Original bin edges.
  old_vals : ndarray, shape (n_old,)
      Piecewise-constant values in the original bins.
  old_errs : ndarray, shape (n_old,)
      Nonnegative standard deviations in the original bins.
  new_bins : ndarray, shape (n_new + 1,)
      Target bin edges.

  Returns
  -------
  tuple[ndarray, ndarray]
      Re-binned values and standard deviations, each with shape ``(n_new,)``.

  Raises
  ------
  ClimaException
      If shapes or bin extents are inconsistent, edges are not strictly
      increasing, or an input standard deviation is negative.
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
