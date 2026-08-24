cimport ClimaRadtranWrk_pxd as rwrk_pxd

cdef class ClimaRadtranWrk:
  cdef rwrk_pxd.ClimaRadtranWrk *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass
  
  property fup_a:
    """ndarray[double,ndim=2], size (nz+1,nw). mW/m^2/Hz in each wavelength 
    bin at the edges of the vertical grid. Upward radiation.
    """
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_fup_a_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_fup_a_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property fdn_a:
    """ndarray[double,ndim=2], size (nz+1,nw). mW/m^2/Hz in each wavelength 
    bin at the edges of the vertical grid. Downward radiation.
    """
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_fdn_a_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_fdn_a_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property fup_n:
    """ndarray[double,ndim=2], size (nz+1). `fup_a`, but integrated over
    wavelength
    """
    def __get__(self):
      cdef int dim1
      rwrk_pxd.climaradtranwrk_fup_n_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rwrk_pxd.climaradtranwrk_fup_n_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property fdn_n:
    """ndarray[double,ndim=2], size (nz+1). `fdn_a`, but integrated over
    wavelength
    """
    def __get__(self):
      cdef int dim1
      rwrk_pxd.climaradtranwrk_fdn_n_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rwrk_pxd.climaradtranwrk_fdn_n_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property amean:
    """ndarray[double,ndim=2], size (nz+1,nw). Mean intensity in photons/cm^2/s.
    Only used in solar radiative transfer.
    """
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_amean_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_amean_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property tau_band:
    "ndarray[double,ndim=2], size (nz,nw). Band optical thickness."
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_tau_band_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_tau_band_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
