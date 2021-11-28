
cdef extern void allocate_atomconservation(void *ptr)
cdef extern void deallocate_atomconservation(void *ptr)

cdef extern void atomconservation_in_surf_get(void *ptr, double *val)
cdef extern void atomconservation_in_top_get(void *ptr, double *val)
cdef extern void atomconservation_in_dist_get(void *ptr, double *val)
cdef extern void atomconservation_out_surf_get(void *ptr, double *val)
cdef extern void atomconservation_out_top_get(void *ptr, double *val)
cdef extern void atomconservation_out_rain_get(void *ptr, double *val)
cdef extern void atomconservation_net_get(void *ptr, double *val)
cdef extern void atomconservation_factor_get(void *ptr, double *val)

cdef class AtomConservation:
  cdef void *_ptr

  def __cinit__(self, bint alloc = False):
    pass

  def __dealloc__(self):
    deallocate_atomconservation(&self._ptr)
    self._ptr = NULL 
      
  property in_surf:
    def __get__(self):
      cdef double val
      atomconservation_in_surf_get(&self._ptr, &val)
      return val
      
  property in_top:
    def __get__(self):
      cdef double val
      atomconservation_in_top_get(&self._ptr, &val)
      return val
      
  property in_dist:
    def __get__(self):
      cdef double val
      atomconservation_in_dist_get(&self._ptr, &val)
      return val
      
  property out_surf:
    def __get__(self):
      cdef double val
      atomconservation_out_surf_get(&self._ptr, &val)
      return val
      
  property out_top:
    def __get__(self):
      cdef double val
      atomconservation_out_top_get(&self._ptr, &val)
      return val
      
  property out_rain:
    def __get__(self):
      cdef double val
      atomconservation_out_rain_get(&self._ptr, &val)
      return val
      
  property net:
    def __get__(self):
      cdef double val
      atomconservation_net_get(&self._ptr, &val)
      return val

  property factor:
    def __get__(self):
      cdef double val
      atomconservation_factor_get(&self._ptr, &val)
      return val
