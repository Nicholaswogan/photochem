cimport ClimaRadtranWrk_pxd as rwrk_pxd

cdef class ClimaRadtranWrk:
  """Non-owning view of radiative-transfer results for one spectral channel.

  Vertical index 0 is the lower boundary and index ``nz`` is the
  top-of-atmosphere boundary. Spectral indices follow the corresponding
  ``RTChannel`` bins. Obtain this object from ``Radtran.wrk_ir`` or
  ``Radtran.wrk_sol`` and keep the parent climate model alive while using it.
  Every property returns a copy.
  """
  cdef rwrk_pxd.ClimaRadtranWrk *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass
  
  property fup_a:
    """ndarray, shape (nz + 1, nw): Upward spectral flux at boundaries (mW/m^2/Hz)."""
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_fup_a_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_fup_a_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property fdn_a:
    """ndarray, shape (nz + 1, nw): Downward spectral flux at boundaries (mW/m^2/Hz)."""
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_fdn_a_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_fdn_a_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property fup_n:
    """ndarray, shape (nz + 1,): Spectrally integrated upward flux (mW/m^2)."""
    def __get__(self):
      cdef int dim1
      rwrk_pxd.climaradtranwrk_fup_n_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rwrk_pxd.climaradtranwrk_fup_n_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property fdn_n:
    """ndarray, shape (nz + 1,): Spectrally integrated downward flux (mW/m^2)."""
    def __get__(self):
      cdef int dim1
      rwrk_pxd.climaradtranwrk_fdn_n_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rwrk_pxd.climaradtranwrk_fdn_n_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property amean:
    """ndarray, shape (nz + 1, nw): Mean photon intensity (photons/cm^2/s).

    This quantity is populated only for shortwave radiative transfer.
    """
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_amean_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_amean_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property tau_band:
    """ndarray, shape (nz, nw): Dimensionless band optical thickness."""
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_tau_band_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_tau_band_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
