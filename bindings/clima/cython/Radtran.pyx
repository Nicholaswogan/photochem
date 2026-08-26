cimport Radtran_pxd as rad_pxd

cdef class Radtran:
  """Non-owning view of an ``AdiabatClimate`` radiative-transfer model.

  Obtain this object from ``AdiabatClimate.rad``; it cannot be initialized
  independently. Keep the parent climate model alive while using this view.
  Array-valued properties return copies, while ``ir``, ``sol``, ``wrk_ir``,
  and ``wrk_sol`` return nested non-owning views.
  """
  cdef rad_pxd.Radtran *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass

  def set_bolometric_flux(self, double flux):
    """Set bolometric stellar flux by adjusting ``photon_scale_factor``.

    Parameters
    ----------
    flux : float
        Bolometric stellar flux at the planet (W/m^2).
    """
    rad_pxd.radtran_set_bolometric_flux_wrapper(self._ptr, &flux)

  def bolometric_flux(self):
    """Return the bolometric stellar flux at the planet.

    Returns 
    -------
    float
        Bolometric flux (W/m^2)
    """
    cdef double flux
    rad_pxd.radtran_bolometric_flux_wrapper(self._ptr, &flux)
    return flux

  def skin_temperature(self, double bond_albedo):
    """Return the radiative skin temperature for a Bond albedo.

    Parameters
    ----------
    bond_albedo : float
        Dimensionless Bond albedo.

    Returns 
    -------
    float
        Skin temperature (K).
    """
    cdef double T_skin
    rad_pxd.radtran_skin_temperature_wrapper(self._ptr, &bond_albedo, &T_skin)
    return T_skin

  def equilibrium_temperature(self, double bond_albedo):
    """Return the planetary equilibrium temperature for a Bond albedo.

    Parameters
    ----------
    bond_albedo : float
        Dimensionless Bond albedo.

    Returns 
    -------
    float
        Equilibrium temperature (K).
    """
    cdef double T_eq
    rad_pxd.radtran_equilibrium_temperature_wrapper(self._ptr, &bond_albedo, &T_eq)
    return T_eq

  def opacities2yaml(self):
    """Return a YAML fragment describing all configured model opacities.

    Returns 
    -------
    str
        YAML opacity configuration.
    """
    cdef int out_len
    cdef void *out_cp
    rad_pxd.radtran_opacities2yaml_wrapper_1(self._ptr, &out_len, &out_cp)
    cdef ndarray out_c = np.empty(out_len + 1, 'S1')
    rad_pxd.radtran_opacities2yaml_wrapper_2(self._ptr, &out_cp, &out_len, <char *>out_c.data)
    return out_c[:-1].tobytes().decode()

  def set_custom_optical_properties(self, ndarray[double, ndim=1] wv, ndarray[double, ndim=1] P, ndarray[double, ndim=2] dtau_dz, ndarray[double, ndim=2] w0, ndarray[double, ndim=2] g0):
    """Set pressure-dependent custom optical properties.

    Pressure must be strictly decreasing and wavelength strictly increasing;
    both grids must be positive. The profiles are copied into interpolation
    data and included in subsequent opacity calculations.

    Parameters
    ----------
    wv : ndarray, shape (n_wavelength,)
        Wavelengths (nm).
    P : ndarray, shape (n_pressure,)
        Pressure grid (dyn/cm^2).
    dtau_dz : ndarray, shape (n_pressure, n_wavelength)
        Optical depth per altitude (1/cm).
    w0 : ndarray, shape (n_pressure, n_wavelength)
        Dimensionless single-scattering albedo.
    g0 : ndarray, shape (n_pressure, n_wavelength)
        Dimensionless asymmetry parameter.
    """
    cdef int dim_wv = wv.shape[0]
    cdef int dim_P = P.shape[0]
    cdef int dim1_dtau_dz = dtau_dz.shape[0]
    cdef int dim2_dtau_dz = dtau_dz.shape[1]
    dtau_dz = np.asfortranarray(dtau_dz)
    cdef int dim1_w0 = w0.shape[0]
    cdef int dim2_w0 = w0.shape[1]
    w0 = np.asfortranarray(w0)
    cdef int dim1_g0 = g0.shape[0]
    cdef int dim2_g0 = g0.shape[1]
    g0 = np.asfortranarray(g0)

    cdef char err[ERR_LEN+1]

    rad_pxd.radtran_set_custom_optical_properties(
      self._ptr, &dim_wv, <double *> wv.data, &dim_P, <double *> P.data,  
      &dim1_dtau_dz, &dim2_dtau_dz, <double *> dtau_dz.data, 
      &dim1_w0, &dim2_w0, <double *> w0.data, 
      &dim1_g0, &dim2_g0, <double *> g0.data, 
      err
    )

    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def unset_custom_optical_properties(self):
    """Remove custom optical properties installed by ``set_custom_optical_properties``."""
    rad_pxd.radtran_unset_custom_optical_properties(self._ptr)
    
  property zenith_u:
    """ndarray, shape (n_zenith,): Writable cosines of solar zenith angles.

    Getting this property returns a copy; assign the complete array to update
    the model.
    """
    def __get__(self):
      cdef int dim1
      rad_pxd.radtran_zenith_u_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rad_pxd.radtran_zenith_u_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] arr):
      cdef int dim1
      rad_pxd.radtran_zenith_u_get_size(self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise ClimaException('zenith_u must have shape (number_of_zeniths,).')
      rad_pxd.radtran_zenith_u_set(self._ptr, &dim1, <double *>arr.data)

  property surface_albedo:
    """ndarray, shape (sol.nw,): Writable surface albedo in each solar bin.

    Getting this property returns a copy; assign the complete array to update
    the model.
    """
    def __get__(self):
      cdef int dim1
      rad_pxd.radtran_surface_albedo_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rad_pxd.radtran_surface_albedo_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] arr):
      cdef int dim1
      rad_pxd.radtran_surface_albedo_get_size(self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise ClimaException('surface_albedo must match the solar wavelength grid.')
      rad_pxd.radtran_surface_albedo_set(self._ptr, &dim1, <double *>arr.data)

  property surface_emissivity:
    """ndarray, shape (ir.nw,): Writable surface emissivity in each longwave bin.

    Getting this property returns a copy; assign the complete array to update
    the model.
    """
    def __get__(self):
      cdef int dim1
      rad_pxd.radtran_surface_emissivity_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rad_pxd.radtran_surface_emissivity_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] arr):
      cdef int dim1
      rad_pxd.radtran_surface_emissivity_get_size(self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise ClimaException('surface_emissivity must match the IR wavelength grid.')
      rad_pxd.radtran_surface_emissivity_set(self._ptr, &dim1, <double *>arr.data)

  property has_hard_surface:
    """bool: Whether to use the hard-surface lower thermal boundary condition.

    False selects the gas-giant no-hard-surface diffusion boundary.
    """
    def __get__(self):
      cdef cbool val
      rad_pxd.radtran_has_hard_surface_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      rad_pxd.radtran_has_hard_surface_set(self._ptr, &val)

  property photon_scale_factor:
    """float: Writable dimensionless multiplier applied to the stellar spectrum."""
    def __get__(self):
      cdef double val
      rad_pxd.radtran_photon_scale_factor_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      rad_pxd.radtran_photon_scale_factor_set(self._ptr, &val)

  property ir_tau_min:
    """float: Writable dimensionless thin-layer guard for longwave two-stream terms."""
    def __get__(self):
      cdef double val
      rad_pxd.radtran_ir_tau_min_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      rad_pxd.radtran_ir_tau_min_set(self._ptr, &val)

  property ir:
    """RTChannel: Non-owning view of longwave spectral-grid metadata."""
    def __get__(self):
      var = RTChannel()
      rad_pxd.radtran_ir_get(self._ptr, &var._ptr)
      return var

  property sol:
    """RTChannel: Non-owning view of shortwave spectral-grid metadata."""
    def __get__(self):
      var = RTChannel()
      rad_pxd.radtran_sol_get(self._ptr, &var._ptr)
      return var

  property wrk_ir:
    """ClimaRadtranWrk: Non-owning view of current longwave results."""
    def __get__(self):
      var = ClimaRadtranWrk()
      rad_pxd.radtran_wrk_ir_get(self._ptr, &var._ptr)
      return var

  property wrk_sol:
    """ClimaRadtranWrk: Non-owning view of current shortwave results."""
    def __get__(self):
      var = ClimaRadtranWrk()
      rad_pxd.radtran_wrk_sol_get(self._ptr, &var._ptr)
      return var
  


  
