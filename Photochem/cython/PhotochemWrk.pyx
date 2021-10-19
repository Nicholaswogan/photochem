
cdef extern void allocate_photochemwrk(void *ptr)
cdef extern void deallocate_photochemwrk(void *ptr)

cdef extern void photochemwrk_usol_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_usol_get(void *ptr, int *dim1, int *dim2, double *usol)

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
  
      
  


    
    
    
    