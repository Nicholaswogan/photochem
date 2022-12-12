cimport PhotochemWrk_pxd as wrk_pxd

cdef class PhotochemWrk:
  cdef void *_ptr
  cdef bint _destroy

  def __cinit__(self, bint alloc = True):
    if alloc:
      wrk_pxd.allocate_photochemwrk(&self._ptr)
      self._destroy = True
    else:
      self._destroy = False

  def __dealloc__(self):
    if self._destroy:
      wrk_pxd.deallocate_photochemwrk(&self._ptr)
      self._ptr = NULL

  property tn:
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_tn_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      wrk_pxd.photochemwrk_tn_set(&self._ptr, &val)

  property usol:
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_usol_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_usol_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=2] usol_new_):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_usol_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray usol_new = np.asfortranarray(usol_new_)
      if usol_new.shape[0] != dim1 or usol_new.shape[1] != dim2:
        raise PhotoException("Input usol is the wrong size.")
      wrk_pxd.photochemwrk_usol_set(&self._ptr, &dim1, &dim2, <double *>usol_new.data)  
  
  property pressure:
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_pressure_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_pressure_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property density:
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_density_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_density_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property densities:
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_densities_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_densities_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property mubar:
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_mubar_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_mubar_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property prates:
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_prates_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_prates_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property amean_grd:
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_amean_grd_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_amean_grd_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property optical_depth:
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_optical_depth_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_optical_depth_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property surf_radiance:
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_surf_radiance_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_surf_radiance_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
  
