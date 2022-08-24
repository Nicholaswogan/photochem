
cimport AtomConservation_pxd as atom_pxd

cdef class AtomConservation:
  cdef void *_ptr

  def __cinit__(self, bint alloc = False):
    pass

  def __dealloc__(self):
    atom_pxd.deallocate_atomconservation(&self._ptr)
    self._ptr = NULL 
      
  property in_surf:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_in_surf_get(&self._ptr, &val)
      return val
      
  property in_top:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_in_top_get(&self._ptr, &val)
      return val
      
  property in_dist:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_in_dist_get(&self._ptr, &val)
      return val

  property in_other:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_in_other_get(&self._ptr, &val)
      return val
      
  property out_surf:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_out_surf_get(&self._ptr, &val)
      return val
      
  property out_top:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_out_top_get(&self._ptr, &val)
      return val
      
  property out_rain:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_out_rain_get(&self._ptr, &val)
      return val
      
  property out_other:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_out_other_get(&self._ptr, &val)
      return val
      
  property net:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_net_get(&self._ptr, &val)
      return val

  property factor:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_factor_get(&self._ptr, &val)
      return val
