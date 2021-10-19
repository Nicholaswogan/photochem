
cdef extern void allocate_photochemvars(void *ptr)
cdef extern void deallocate_photochemvars(void *ptr)

cdef extern void photochemvars_nz_get(void *ptr, int *nz)
cdef extern void photochemvars_usol_init_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemvars_usol_init_get(void *ptr, int *dim1, int *dim2, double *usol_init)

cdef class PhotochemVars:
  cdef void *_ptr
  cdef bint _destroy

  def __cinit__(self, bint alloc = True):
    if alloc:
      allocate_photochemvars(&self._ptr)
      self._destroy = True
    else:
      self._destroy = False

  def __dealloc__(self):
    if self._destroy:
      deallocate_photochemvars(&self._ptr)
      self._ptr = NULL
  
  property nz:
    def __get__(self):
      cdef int nz
      photochemvars_nz_get(&self._ptr, &nz)
      return nz
      
  property usol_init:
    def __get__(self):
      cdef int dim1, dim2
      photochemvars_usol_init_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      photochemvars_usol_init_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
      
  


    
    
    
    