
cimport AtomConservation_pxd as atom_pxd

cdef class CondensationParameters:
  """Parameters controlling finite-rate condensation and evaporation.

  Instances are borrowed views into an [EvoAtmosphere][photochem.EvoAtmosphere] model. Assigning
  a property changes the corresponding model setting; these objects should not
  be constructed directly.
  """

  cdef atom_pxd.CondensationParameters *_ptr

  def __cinit__(self):
    self._ptr = NULL

  property k_cond:
    """float. Dimensionless multiplier applied to ``Kzz/H**2``."""
    def __get__(self):
      cdef double val
      atom_pxd.condensationparameters_k_cond_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      atom_pxd.condensationparameters_k_cond_set(self._ptr, &val)

  property k_evap:
    """float. Dimensionless multiplier applied to ``wfall/H``."""
    def __get__(self):
      cdef double val
      atom_pxd.condensationparameters_k_evap_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      atom_pxd.condensationparameters_k_evap_set(self._ptr, &val)

  property RHc:
    """float. Relative humidity where condensation begins."""
    def __get__(self):
      cdef double val
      atom_pxd.condensationparameters_rhc_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      atom_pxd.condensationparameters_rhc_set(self._ptr, &val)

  property smooth_factor:
    """float. Fractional smoothing width used near thresholds to reduce stiffness."""
    def __get__(self):
      cdef double val
      atom_pxd.condensationparameters_smooth_factor_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      atom_pxd.condensationparameters_smooth_factor_set(self._ptr, &val)
    
cdef class SaturationData:
  """Saturation-vapor-pressure model for one condensable species.

  Instances are read-only borrowed views returned by
  [PhotochemData.particle_sat][photochem._photochem.PhotochemData.particle_sat]; they should not be constructed directly.
  """

  cdef atom_pxd.SaturationData *_ptr

  def __cinit__(self):
    self._ptr = NULL

  def sat_pressure(self, double T):
    """Evaluate saturation vapor pressure at temperature ``T``.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Saturation vapor pressure in dyn/cm^2.
    """
    cdef double Psat
    atom_pxd.saturationdata_sat_pressure_wrapper(self._ptr, &T, &Psat)
    return Psat
