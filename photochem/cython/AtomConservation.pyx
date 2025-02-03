
cimport AtomConservation_pxd as atom_pxd

cdef class CondensationParameters:

  cdef atom_pxd.CondensationParameters *_ptr

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

  cdef atom_pxd.SaturationData *_ptr

  def __cinit__(self):
    self._ptr = NULL

  def sat_pressure(self, double T):
    cdef double Psat
    atom_pxd.saturationdata_sat_pressure_wrapper(self._ptr, &T, &Psat)
    return Psat