
cdef extern void allocate_photochemwrk(void *ptr)
cdef extern void deallocate_photochemwrk(void *ptr)

cdef extern void photochemwrk_usol_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_usol_get(void *ptr, int *dim1, int *dim2, double *usol)
cdef extern void photochemwrk_usol_set(void *ptr, int *dim1, int *dim2, double *usol)

cdef extern void photochemwrk_pressure_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_pressure_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_density_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_density_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_prates_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_prates_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_surf_radiance_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_surf_radiance_get(void *ptr, int *dim1, double *arr)

cdef class PhotochemWrk:
  cdef void *_ptr
  cdef bint _destroy

  def __cinit__(self, bint alloc = True):
    if alloc:
      allocate_photochemwrk(&self._ptr)
      self._destroy = True
    else:
      self._destroy = False

  def __dealloc__(self):
    if self._destroy:
      deallocate_photochemwrk(&self._ptr)
      self._ptr = NULL
      
  property usol:
    def __get__(self):
      cdef int dim1, dim2
      photochemwrk_usol_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      photochemwrk_usol_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=2] usol_new):
      cdef int dim1, dim2
      photochemwrk_usol_get_size(&self._ptr, &dim1, &dim2)
      if not np.isfortran(usol_new):
        raise PhotoException("Input usol must have fortran ordering.")
      if usol_new.shape[0] != dim1 or usol_new.shape[1] != dim2:
        raise PhotoException("Input usol is the wrong size.")
      photochemwrk_usol_set(&self._ptr, &dim1, &dim2, <double *>usol_new.data)  
  
  property pressure:
    def __get__(self):
      cdef int dim1
      photochemwrk_pressure_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemwrk_pressure_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property density:
    def __get__(self):
      cdef int dim1
      photochemwrk_density_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemwrk_density_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property prates:
    def __get__(self):
      cdef int dim1, dim2
      photochemwrk_prates_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      photochemwrk_prates_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property surf_radiance:
    def __get__(self):
      cdef int dim1
      photochemwrk_surf_radiance_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemwrk_surf_radiance_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
  
