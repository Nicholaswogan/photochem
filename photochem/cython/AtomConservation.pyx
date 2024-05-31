
cimport AtomConservation_pxd as atom_pxd

cdef class AtomConservation:
  """A class holding data describing the conservation of an atom
  in an atmosphere. This information is useful for determining if the model
  is at a steady-state. This class is produced from the `atom_conservation`
  routine.
  """

  cdef void *_ptr

  def __cinit__(self):
    self._ptr = NULL 

  def __dealloc__(self):
    atom_pxd.deallocate_atomconservation(self._ptr)
    self._ptr = NULL 
      
  property in_surf:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_in_surf_get(self._ptr, &val)
      return val
      
  property in_top:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_in_top_get(self._ptr, &val)
      return val
      
  property in_dist:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_in_dist_get(self._ptr, &val)
      return val

  property in_other:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_in_other_get(self._ptr, &val)
      return val
      
  property out_surf:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_out_surf_get(self._ptr, &val)
      return val
      
  property out_top:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_out_top_get(self._ptr, &val)
      return val
      
  property out_rain:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_out_rain_get(self._ptr, &val)
      return val
      
  property out_other:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_out_other_get(self._ptr, &val)
      return val
      
  property net:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_net_get(self._ptr, &val)
      return val

  property factor:
    def __get__(self):
      cdef double val
      atom_pxd.atomconservation_factor_get(self._ptr, &val)
      return val

cdef class CondensationParameters:

  cdef void *_ptr

  def __cinit__(self):
    self._ptr = NULL

  property k_cond:
    def __get__(self):
      cdef double val
      atom_pxd.condensationparameters_k_cond_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      atom_pxd.condensationparameters_k_cond_set(self._ptr, &val)

  property k_evap:
    def __get__(self):
      cdef double val
      atom_pxd.condensationparameters_k_evap_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      atom_pxd.condensationparameters_k_evap_set(self._ptr, &val)

  property RHc:
    def __get__(self):
      cdef double val
      atom_pxd.condensationparameters_rhc_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      atom_pxd.condensationparameters_rhc_set(self._ptr, &val)

  property smooth_factor:
    def __get__(self):
      cdef double val
      atom_pxd.condensationparameters_smooth_factor_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      atom_pxd.condensationparameters_smooth_factor_set(self._ptr, &val)
    
cdef class SaturationData:

  cdef void *_ptr

  def __cinit__(self):
    self._ptr = NULL

  def sat_pressure(self, double T):
    cdef double Psat
    atom_pxd.saturationdata_sat_pressure_wrapper(self._ptr, &T, &Psat)
    return Psat