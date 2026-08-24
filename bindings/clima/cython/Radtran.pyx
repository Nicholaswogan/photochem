cimport Radtran_pxd as rad_pxd

cdef class Radtran:
  cdef rad_pxd.Radtran *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass

  def set_bolometric_flux(self, double flux):
    """Sets the bolometric stellar flux by adjusting the `photon_scale_factor`.

    Parameters
    ----------
    flux : float
        Bolometric flux (W/m^2)
    """
    rad_pxd.radtran_set_bolometric_flux_wrapper(self._ptr, &flux)

  def bolometric_flux(self):
    """The bolometric stellar flux at the planet in W/m^2

    Returns 
    -------
    flux : float
        Bolometric flux (W/m^2)
    """
    cdef double flux
    rad_pxd.radtran_bolometric_flux_wrapper(self._ptr, &flux)
    return flux

  def skin_temperature(self, double bond_albedo):
    """The skin temperature

    Parameters
    ----------
    bond_albedo : float
        The bond albedo of a planet

    Returns 
    -------
    float
        The skin temperature
    """
    cdef double T_skin
    rad_pxd.radtran_skin_temperature_wrapper(self._ptr, &bond_albedo, &T_skin)
    return T_skin

  def equilibrium_temperature(self, double bond_albedo):
    """The equilibrium temperature

    Parameters
    ----------
    bond_albedo : float
        The bond albedo of a planet

    Returns 
    -------
    float
        The equilibrium temperature
    """
    cdef double T_eq
    rad_pxd.radtran_equilibrium_temperature_wrapper(self._ptr, &bond_albedo, &T_eq)
    return T_eq

  def opacities2yaml(self):
    """Returns a yaml string representing all opacities in the model.

    Returns 
    -------
    str
        string representing all opacities.
    """
    cdef int out_len
    cdef void *out_cp
    rad_pxd.radtran_opacities2yaml_wrapper_1(self._ptr, &out_len, &out_cp)
    cdef ndarray out_c = np.empty(out_len + 1, 'S1')
    rad_pxd.radtran_opacities2yaml_wrapper_2(self._ptr, &out_cp, &out_len, <char *>out_c.data)
    return out_c[:-1].tobytes().decode()

  def set_custom_optical_properties(self, ndarray[double, ndim=1] wv, ndarray[double, ndim=1] P, ndarray[double, ndim=2] dtau_dz, ndarray[double, ndim=2] w0, ndarray[double, ndim=2] g0):
    """Sets custom optical properties

    Parameters
    ----------
    wv : ndarray[double,ndim=1]
        Array of of wavelengths in nm
    P : ndarray[double,ndim=1]
        Array of pressures in dynes/cm^2. Must be decreasing.
    dtau_dz : ndarray[double,ndim=2]
        Optical depth per altitude (1/cm). Shape (len(P), len(wv)).
    w0 : ndarray[double,ndim=2]
        Single scattering albedo. Shape (len(P), len(wv)).
    g0 : ndarray[double,ndim=2]
        Asymmetry parameter. Shape (len(P), len(wv)).
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
    "Unsets custom optical properties set with `set_custom_optical_properties`."
    rad_pxd.radtran_unset_custom_optical_properties(self._ptr)
    
  property zenith_u:
    "ndarray[double,ndim=1]. cosine of the zenith angle in radians."
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
        raise ClimaException('"zenith_u" is the wrong size')
      rad_pxd.radtran_zenith_u_set(self._ptr, &dim1, <double *>arr.data)

  property surface_albedo:
    "ndarray[double,ndim=1]. The surface albedo in each solar wavelength bin"
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
        raise ClimaException('"surface_albedo" is the wrong size')
      rad_pxd.radtran_surface_albedo_set(self._ptr, &dim1, <double *>arr.data)

  property surface_emissivity:
    "ndarray[double,ndim=1]. The surface emissivity in each IR wavelength bin"
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
        raise ClimaException('"surface_emissivity" is the wrong size')
      rad_pxd.radtran_surface_emissivity_set(self._ptr, &dim1, <double *>arr.data)

  property has_hard_surface:
    "bool. If true use hard-surface lower thermal BC; if false use gas-giant no-hard-surface BC."
    def __get__(self):
      cdef cbool val
      rad_pxd.radtran_has_hard_surface_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      rad_pxd.radtran_has_hard_surface_set(self._ptr, &val)

  property photon_scale_factor:
    """float. A scale factor that is applied to `photons_sol` so that
    bolometric luminosity can be easily changed."""
    def __get__(self):
      cdef double val
      rad_pxd.radtran_photon_scale_factor_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      rad_pxd.radtran_photon_scale_factor_set(self._ptr, &val)

  property ir_tau_min:
    "float. Thin-layer optical-depth guard used in IR two-stream source terms."
    def __get__(self):
      cdef double val
      rad_pxd.radtran_ir_tau_min_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      rad_pxd.radtran_ir_tau_min_set(self._ptr, &val)

  property ir:
    "The RTChannel for longwave radiative transfer"
    def __get__(self):
      var = RTChannel()
      rad_pxd.radtran_ir_get(self._ptr, &var._ptr)
      return var

  property sol:
    "The RTChannel for shortwave radiative transfer"
    def __get__(self):
      var = RTChannel()
      rad_pxd.radtran_sol_get(self._ptr, &var._ptr)
      return var

  property wrk_ir:
    "The ClimaRadtranWrk for longwave radiative transfer"
    def __get__(self):
      var = ClimaRadtranWrk()
      rad_pxd.radtran_wrk_ir_get(self._ptr, &var._ptr)
      return var

  property wrk_sol:
    "The ClimaRadtranWrk for shortwave radiative transfer"
    def __get__(self):
      var = ClimaRadtranWrk()
      rad_pxd.radtran_wrk_sol_get(self._ptr, &var._ptr)
      return var
  


  
