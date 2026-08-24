
cimport AdiabatClimate_pxd as wa_pxd
  
cdef class AdiabatClimate:
  """This class is a 1-D climate model. The code draws multispecies 
  pseudoadiabats upward (Eq. 1 in Graham et al. 2021, PSJ), connected
  to an assumed isothermal stratosphere. The model also can perform radiative
  transfer and do a non-linear solve for the incoming solar radiation
  that balances the outgoing longwave radiation.
  """

  cdef wa_pxd.AdiabatClimate *_ptr
  cdef cbool _init_called

  def __cinit__(self, *args, **kwargs):
    self._init_called = False
    self._ptr = wa_pxd.allocate_adiabatclimate()

  def __dealloc__(self):
    wa_pxd.deallocate_adiabatclimate(self._ptr)

  def __getattribute__(self, name):
    if not self._init_called:
      raise ClimaException('The "__init__" method of AdiabatClimate has not been called.')
    return super().__getattribute__(name)

  def __setattr__(self, name, value):
    if not self._init_called:
      raise ClimaException('The "__init__" method of AdiabatClimate has not been called.')
    PyObject_GenericSetAttr(self, name, value)

  def __init__(self, str species_file, str settings_file, 
                     str flux_file, data_dir = None, cbool double_radiative_grid = True):           
    """Initializes `AdiabatClimate`

    Parameters
    ----------
    species_file : str
        The input species .yaml file
    settings_file : str
        The input settings .yaml file
    flux_file : str
        The input stellar flux file
    data_dir : str, optional
        The directory where all data is stored.
    double_radiative_grid : bool, optional
        If True (default), radiative transfer uses a refined grid where each
        physical layer is split into two RT layers, plus two ghost RT layers
        above the model top to improve TOA numerical stability.
    """

    self._init_called = True

    if data_dir == None:
      data_dir_ = photochem_clima_data.DATA_DIR
    else:
      data_dir_ = data_dir
    
    # convert strings to char
    cdef bytes species_file_b = pystring2cstring(species_file)
    cdef char *species_file_c = species_file_b
    cdef bytes settings_file_b = pystring2cstring(settings_file)
    cdef char *settings_file_c = settings_file_b
    cdef bytes flux_file_b = pystring2cstring(flux_file)
    cdef char *flux_file_c = flux_file_b
    cdef bytes data_dir_b = pystring2cstring(data_dir_)
    cdef char *data_dir_c = data_dir_b
    cdef char err[ERR_LEN+1]
    
    # Initialize
    wa_pxd.adiabatclimate_create_wrapper(self._ptr, species_file_c,
                                         settings_file_c, flux_file_c, data_dir_c, &double_radiative_grid,
                                         err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def make_profile(self, double T_surf, ndarray[double, ndim=1] P_i_surf):
    """Constructs an atmosphere using a multispecies pseudoadiabat (Eq. 1 in Graham et al. 2021, PSJ)
    troposphere connected to an isothermal stratosphere.

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    P_i_surf : ndarray[double,ndim=1]
        Array of surface pressures of each species (dynes/cm^2)
    """
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_profile_wrapper(self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def make_profile_dry(self, ndarray[double, ndim=1] P, ndarray[double, ndim=1] T, ndarray[double, ndim=2] f_i):
    """Given a P, T and mixing ratios, this function will update all atmospheric variables
    (except self%P_trop and self%convecting_with_below) to reflect these inputs. 
    The atmosphere is assumed to be dry (no condensibles). Any gas exceeding saturation 
    will not be altered.

    Parameters
    ----------
    P : ndarray[double,ndim=1]
        Pressure (dynes/cm^2). The first element is the surface.
    T : ndarray[double,ndim=1]
        Temperature (K) defined on `P`.
    f_i : ndarray[double,ndim=2]
        Mixing ratios defined on `P` of shape (size(P),ng).
    """
    cdef int dim_P = P.shape[0]
    cdef int dim_T = T.shape[0]
    cdef ndarray f_i_ = np.asfortranarray(f_i)
    cdef int dim1_f_i = f_i_.shape[0]
    cdef int dim2_f_i = f_i_.shape[1]
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_profile_dry_wrapper(
      self._ptr, 
      &dim_P, <double *>P.data,
      &dim_T, <double *>T.data,
      &dim1_f_i, &dim2_f_i, <double *>f_i_.data,
      err
    )
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    
  def make_column(self, double T_surf, ndarray[double, ndim=1] N_i_surf):
    """Similar to `make_profile`, but instead the input is column reservoirs 
    of each gas. 

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    N_i_surf : ndarray[double,ndim=1]
        Array of columns of each species (mol/cm^2)
    """
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_column_wrapper(self._ptr, &T_surf,
    &ng, <double *>N_i_surf.data, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
  
  def make_profile_bg_gas(self, double T_surf, ndarray[double, ndim=1] P_i_surf, double P_surf, str bg_gas):
    """Similar to `make_profile`, but instead imposes a background gas and fixed
    surface pressure. We do a non-linear solve for the background gas pressure
    so that the desired total surface pressure is satisfied.

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    P_i_surf : ndarray[double,ndim=1]
        Array of surface pressures of each species (dynes/cm^2)
    P_surf : float
        The surface pressure (dynes/cm^2)
    bg_gas : str
        The name of the background gas
    """
    cdef int ng = P_i_surf.shape[0]
    cdef bytes bg_gas_b = pystring2cstring(bg_gas)
    cdef char *bg_gas_c = bg_gas_b
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_profile_bg_gas_wrapper(self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, &P_surf, bg_gas_c, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def TOA_fluxes(self, double T_surf, ndarray[double, ndim=1] P_i_surf):
    """Calls `make_profile`, then does radiative transfer on the constructed atmosphere

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    P_i_surf : ndarray[double,ndim=1]
        Array of surface pressures of each species (dynes/cm^2)

    Returns
    -------
    tuple
        The tuple has two elements. The first is the incoming solar radiation at the top
        of the atmosphere, and the second is the outgoing longwave radiation at the top
        of the atmosphere
    """
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.adiabatclimate_toa_fluxes_wrapper(self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, &ISR, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR

  def TOA_fluxes_column(self, double T_surf, ndarray[double, ndim=1] N_i_surf):
    """Calls `make_column`, then does radiative transfer on the constructed atmosphere

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    N_i_surf : ndarray[double,ndim=1]
        Array of columns of each species (mol/cm^2)

    Returns
    -------
    tuple
        The tuple has two elements. The first is the incoming solar radiation at the top
        of the atmosphere, and the second is the outgoing longwave radiation at the top
        of the atmosphere
    """
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.adiabatclimate_toa_fluxes_column_wrapper(self._ptr, &T_surf,
    &ng, <double *>N_i_surf.data, &ISR, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR

  def TOA_fluxes_bg_gas(self, double T_surf, ndarray[double, ndim=1] P_i_surf, double P_surf, str bg_gas):
    """Calls `make_profile_bg_gas`, then does radiative transfer on the constructed atmosphere

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    P_i_surf : ndarray[double,ndim=1]
        Array of surface pressures of each species (dynes/cm^2)
    P_surf : float
        The surface pressure (dynes/cm^2)
    bg_gas : str
        The name of the background gas

    Returns
    -------
    tuple
        The tuple has two elements. The first is the incoming solar radiation at the top
        of the atmosphere, and the second is the outgoing longwave radiation at the top
        of the atmosphere
    """
    cdef int ng = P_i_surf.shape[0]
    cdef bytes bg_gas_b = pystring2cstring(bg_gas)
    cdef char *bg_gas_c = bg_gas_b
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.adiabatclimate_toa_fluxes_bg_gas_wrapper(self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, &P_surf, bg_gas_c, &ISR, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR

  def TOA_fluxes_dry(self, ndarray[double, ndim=1] P, ndarray[double, ndim=1] T, ndarray[double, ndim=2] f_i):
    """Calls `make_profile_dry`, then does radiative transfer on the constructed atmosphere

    Parameters
    ----------
    P : ndarray[double,ndim=1]
        Pressure (dynes/cm^2). The first element is the surface.
    T : ndarray[double,ndim=1]
        Temperature (K) defined on `P`.
    f_i : ndarray[double,ndim=2]
        Mixing ratios defined on `P` of shape (size(P),ng).
    """
    cdef int dim_P = P.shape[0]
    cdef int dim_T = T.shape[0]
    cdef ndarray f_i_ = np.asfortranarray(f_i)
    cdef int dim1_f_i = f_i_.shape[0]
    cdef int dim2_f_i = f_i_.shape[1]
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.adiabatclimate_toa_fluxes_dry_wrapper(
      self._ptr, 
      &dim_P, <double *>P.data,
      &dim_T, <double *>T.data,
      &dim1_f_i, &dim2_f_i, <double *>f_i_.data,
      &ISR, &OLR,
      err
    )
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR
  
  def surface_temperature(self, ndarray[double, ndim=1] P_i_surf, double T_guess = 280):
    """Does a non-linear solve for the surface temperature that balances incoming solar
    and outgoing longwave radiation. Uses `make_profile`.

    Parameters
    ----------
    P_i_surf : ndarray[double,ndim=1]
        Array of surface pressures of each species (dynes/cm^2)
    T_guess : float, optional
        Guess for the surface temperature (K)

    Returns
    -------
    float
        The surface temperature at an equilibrium climate state
    """
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.adiabatclimate_surface_temperature_wrapper(self._ptr, 
    &ng, <double *>P_i_surf.data, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  def surface_temperature_column(self, ndarray[double, ndim=1] N_i_surf, double T_guess = 280):
    """Does a non-linear solve for the surface temperature that balances incoming solar
    and outgoing longwave radiation. Uses `make_column`.

    Parameters
    ----------
    N_i_surf : ndarray[double,ndim=1]
        Array of columns of each species (mol/cm^2)
    T_guess : float, optional
        Guess for the surface temperature (K)

    Returns
    -------
    float
        The surface temperature at an equilibrium climate state
    """
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.adiabatclimate_surface_temperature_column_wrapper(self._ptr, 
    &ng, <double *>N_i_surf.data, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  def set_ocean_solubility_fcn(self, str species, object fcn):
    """Sets a function for describing how gases dissolve in a liquid ocean.

    Parameters
    ----------
    species : str
        Species that the ocean is made of
    fcn : function
        A Numba cfunc that describes the solubility of other gases in the ocean
    """
    cdef bytes species_b = pystring2cstring(species)
    cdef char *species_c = species_b
    cdef char err[ERR_LEN+1]
    cdef uintptr_t fcn_l
    cdef wa_pxd.ocean_solubility_fcn fcn_c

    if fcn is None:
      fcn_l = 0
      fcn_c = NULL
    else:
      argtypes = (ct.c_double, ct.c_int32, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_void_p)
      restype = None
      fcn_argtypes = fcn.ctypes.argtypes
      if not (len(fcn_argtypes) == 5 and all([fcn_argtypes[i] == argtypes[i] for i in range(4)])):
        raise ClimaException("The callback function has the wrong argument types.")
      if not fcn.ctypes.restype == restype:
        raise ClimaException("The callback function has the wrong return type.")

      fcn_l = fcn.address
      fcn_c = <wa_pxd.ocean_solubility_fcn> fcn_l

    wa_pxd.adiabatclimate_set_ocean_solubility_fcn_wrapper(self._ptr, species_c, fcn_c, err)
    if len(err.strip()) > 0:
       raise ClimaException(err.decode("utf-8").strip())

  def surface_temperature_bg_gas(self, ndarray[double, ndim=1] P_i_surf, double P_surf, str bg_gas, double T_guess = 280):
    """Similar to surface_temperature. The difference is that this function imposes
    a background gas and fixed surface pressure.

    Parameters
    ----------
    P_i_surf : ndarray[double,ndim=1]
        Array of surface pressures of each species (dynes/cm^2)
    P_surf : float
        The surface pressure (dynes/cm^2)
    bg_gas : str
        The name of the background gas
    T_guess : float, optional
        Guess for the surface temperature (K)

    Returns
    -------
    float
        The surface temperature at an equilibrium climate state
    """
    cdef int ng = P_i_surf.shape[0]
    cdef bytes bg_gas_b = pystring2cstring(bg_gas)
    cdef char *bg_gas_c = bg_gas_b
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.adiabatclimate_surface_temperature_bg_gas_wrapper(self._ptr, 
    &ng, <double *>P_i_surf.data, &P_surf, bg_gas_c, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  def set_particle_density_and_radii(self, ndarray[double, ndim=1] P, ndarray[double, ndim=2] pdensities, ndarray[double, ndim=2] pradii):
    """Sets particle densities and radii.

    Parameters
    ----------
    P_i_surf : ndarray[double,ndim=1]
        Array of pressures in dynes/cm^2
    pdensities : ndarray[double,ndim=2]
        Particle densities in particles/cm^3 at each pressure and for each 
        particle in the model. Shape (nz, np).
    pradii : ndarray[double,ndim=2]
        Particle radii in cm at each pressure and for each particle
        in the model. Shape (nz, np).
    """
    cdef int dim_P = P.shape[0]
    cdef int dim1_pdensities = pdensities.shape[0]
    cdef int dim2_pdensities = pdensities.shape[1]
    pdensities = np.asfortranarray(pdensities)
    cdef int dim1_pradii = pradii.shape[0]
    cdef int dim2_pradii = pradii.shape[1]
    pradii = np.asfortranarray(pradii)
    cdef char err[ERR_LEN+1]

    wa_pxd.adiabatclimate_set_particle_density_and_radii(
      self._ptr, &dim_P, <double *> P.data,  
      &dim1_pdensities, &dim2_pdensities, <double *> pdensities.data, 
      &dim1_pradii, &dim2_pradii, <double *> pradii.data, 
      err
    )

    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def RCE(self, ndarray[double, ndim=1] P_i_surf, double T_surf_guess, ndarray[double, ndim=1] T_guess,
          convecting_with_below = None, custom_dry_mix = None):
    """Compute full radiative-convective equilibrium.

    Parameters
    ----------
    P_i_surf : ndarray[double,ndim=1]
        Array of surface pressures of each species (dynes/cm^2)
    T_surf_guess : float
        A guess for the surface temperature (K)
    T_guess : ndarray[double,ndim=1]
        A guess for the temperature in each atmospheric layer (K)
    convecting_with_below : ndarray[bool,ndim=1], optional
        An array describing a guess for the radiative vs. convective 
        regions of the atmosphere
    custom_dry_mix : dict of {str: ndarray[double,ndim=1]}, optional
        Allows the user to specify vertically inhomogenous mixing ratios for dry species. 
        The dictionary must contain a "pressure" key specifying a pressure grid in dynes/cm^2.
        The other keys are mixing ratio profiles for dry species at each pressure point. 
        The total custom surface pressure is determined by:

        ```python
        P_dry_surf = 0.0
        for key in custom_dry_mix:
            if key != 'pressure':
                ind = self.species_names.index(key)
                P_dry_surf += P_i_surf[ind]
        ```

        The underlying code interpolates `custom_dry_mix` to the appropriate pressure levels.
        Note, these mixing ratios are not absolute volume mixing ratios. They are instead
        gas concentrations relative to all gases specified in `custom_dry_mix`.

    Returns
    -------
    bool
        Whether the routine converged or not.
    """
    cdef int ng = P_i_surf.shape[0]
    cdef int dim_T_guess = T_guess.shape[0]

    cdef ndarray[cbool, ndim=1] convecting_with_below_ = np.array([False],dtype=np.bool_)
    cdef cbool convecting_with_below_present
    if convecting_with_below is None:
      convecting_with_below_present = False
    else:
      convecting_with_below_present = True
      convecting_with_below_ = convecting_with_below
    cdef int dim_convecting_with_below = convecting_with_below_.shape[0]

    # Custom mixing ratio
    cdef cbool custom_present = False
    # initialize arrays
    cdef int dim_sp_custom = 1
    cdef ndarray sp_custom_c = np.zeros(dim_sp_custom*S_STR_LEN + 1, 'S1')
    cdef int dim_P_custom = 1
    cdef ndarray[double, ndim=1] P_custom = np.empty(dim_P_custom, np.double)
    cdef int dim1_mix_custom = 1
    cdef int dim2_mix_custom = 1
    cdef ndarray[double, ndim=2, mode='fortran'] mix_custom = np.empty((dim1_mix_custom,dim2_mix_custom), np.double)
    cdef list species

    if custom_dry_mix is not None:
      custom_present = True

      species = list(custom_dry_mix.keys())
      if 'pressure' not in species:
        raise ClimaException('`pressure` must be a key in `custom_dry_mix`')
      species.remove('pressure')
      dim_sp_custom = len(species)
      sp_custom_c = list2cstring(species, S_STR_LEN)

      P_custom = custom_dry_mix['pressure']
      dim_P_custom = len(P_custom)

      dim1_mix_custom = dim_P_custom
      dim2_mix_custom = dim_sp_custom
      mix_custom = np.empty((dim_P_custom,dim_sp_custom), np.double, order='F')
      for i,sp in enumerate(species):
        if sp != 'pressure':
          mix_custom[:,i] = custom_dry_mix[sp]

    cdef cbool converged
    cdef char err[ERR_LEN+1]

    wa_pxd.adiabatclimate_rce_wrapper(
      self._ptr, &ng, <double *>P_i_surf.data, &T_surf_guess, &dim_T_guess, <double *>T_guess.data, 
      &convecting_with_below_present, &dim_convecting_with_below, <cbool *>convecting_with_below_.data,
      &custom_present, &dim_sp_custom, <char*> sp_custom_c.data, &dim_P_custom, <double *> P_custom.data, 
      &dim1_mix_custom, &dim2_mix_custom, <double *> mix_custom.data, 
      &converged, err
    )
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return converged

  def to_regular_grid(self):
    "Re-grids atmosphere so that each grid cell is equally spaced in altitude."
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_to_regular_grid_wrapper(self._ptr, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def out2atmosphere_txt(self, str filename, ndarray[double, ndim=1] eddy, int number_of_decimals=5, cbool overwrite = False, cbool clip = True):
    """Saves state of the atmosphere to a file.

    Parameters
    ----------
    filename : str
        Output filename
    eddy : ndarray[double,ndim=1]
        Array of eddy diffusions (cm^2/s) to write to the output file. This is
        useful for coupling to the photochemical model.
    number_of_decimals : int
        Number of decimals
    overwrite : bool, optional
        If true, then output file can be overwritten, by default False
    clip : bool, optional
        If true, then mixing ratios are clipped at a very small 
        positive number, by default False
    """  
    cdef int nz = eddy.shape[0]
    cdef bytes filename_b = pystring2cstring(filename)
    cdef char *filename_c = filename_b
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_out2atmosphere_txt_wrapper(self._ptr, filename_c, &nz, <double *>eddy.data, &number_of_decimals, &overwrite, &clip, err)  
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def heat_redistribution_parameters(self):
    """For considering a tidally locked planet. This function computes key parameters for
    Equation (10) in Koll (2022), ApJ. The function must be called after calling a function
    `self.TOA_fluxes`, because it uses atmosphere properties, and radiative properties.

    Returns
    -------
    tuple
        The tuple has three elements: `(tau_LW, k_term, f_term)`. 
        `tau_LW` is the long wavelength optical depth.
        `k_term` is the k term in Equation (10).
        `f_term` is the heat redistribution parameter, f in Equation (10)
    """
    cdef char err[ERR_LEN+1];
    cdef double tau_LW, k_term, f_term;
    wa_pxd.adiabatclimate_heat_redistribution_parameters_wrapper(self._ptr, &tau_LW, &k_term, &f_term, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return tau_LW, k_term, f_term

  property P_top:
    "float. Pressure of the top of the atmosphere (dynes/cm^2)"
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_top_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_p_top_set(self._ptr, &val)

  property reference_pressure:
    """float. Reference pressure for which `planet_radius` is defined (dynes/cm^2).
    If <= 0, `planet_radius` is interpreted at `P_surf`.
    """
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_reference_pressure_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_reference_pressure_set(self._ptr, &val)

  property T_trop:
    "float. Tropopause temperature (K)"
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_t_trop_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_t_trop_set(self._ptr, &val)

  property use_make_column_P_guess:
    """bool. If True, then any function that calls `make_column` will
    use the initial guess in `self.make_column_P_guess`
    """
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_use_make_column_p_guess_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_use_make_column_p_guess_set(self._ptr, &val)

  property make_column_P_guess:
    """ndarray[double,ndim=1], shape (ng). Initial guess for surface pressure 
    of all gases for `make_column`.
    """
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_make_column_p_guess_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_make_column_p_guess_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] arr):
      cdef int dim1
      wa_pxd.adiabatclimate_make_column_p_guess_get_size(self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise ClimaException('"make_column_P_guess" is the wrong size')
      wa_pxd.adiabatclimate_make_column_p_guess_set(self._ptr, &dim1, <double *>arr.data)

  property solve_for_T_trop:
    """bool. If True, then Tropopause temperature is non-linearly solved for such that
    it matches the skin temperature. The initial guess will always be self.T_trop.
    """
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_solve_for_t_trop_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_solve_for_t_trop_set(self._ptr, &val)

  property albedo_fcn:
    """Callback that sets the surface albedo based on the surface temperature.
    This can be used to parameterize the ice-albedo feedback.
    """
    def __set__(self, object fcn):
      cdef uintptr_t fcn_l
      cdef wa_pxd.temp_dependent_albedo_fcn fcn_c
      if fcn is None:
        fcn_l = 0
        fcn_c = NULL
      else:
        argtypes = (ct.c_double,)
        restype = ct.c_double
        if not fcn.ctypes.argtypes == argtypes:
          raise ClimaException("The callback function has the wrong argument types.")
        if not fcn.ctypes.restype == restype:
          raise ClimaException("The callback function has the wrong return type.")
        fcn_l = fcn.address
        fcn_c = <wa_pxd.temp_dependent_albedo_fcn> fcn_l

      wa_pxd.adiabatclimate_albedo_fcn_set(self._ptr, fcn_c)

  property RH:
    "ndarray[double,ndim=1], shape (ng). Relative humidity of each species"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_rh_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_rh_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] arr):
      cdef int dim1
      wa_pxd.adiabatclimate_rh_get_size(self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise ClimaException('"RH" is the wrong size')
      wa_pxd.adiabatclimate_rh_set(self._ptr, &dim1, <double *>arr.data)

  property ocean_args_p:
    "int or NoneType. Pointer to data that is passed to ocean solubility functions."
    def __set__(self, object p_int):
      cdef uintptr_t p1
      cdef void * p
      if p_int is None:
        p = NULL
      else:
        p1 = p_int
        p = <void *>p1
      wa_pxd.adiabatclimate_ocean_args_p_set(self._ptr, p)

  property tidally_locked_dayside:
    """bool. If True, then will attempt to compute the climate corresponding to the
    observed dayside temperature of a tidally locked planet.
    """
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_tidally_locked_dayside_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_tidally_locked_dayside_set(self._ptr, &val)

  property L:
    "float. Circulation's horizontal scale (cm)."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_l_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_l_set(self._ptr, &val)

  property chi:
    "float. Heat engine efficiency term (no units)."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_chi_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_chi_set(self._ptr, &val)

  property n_LW:
    "float. = 1 or 2 (no units)."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_n_lw_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_n_lw_set(self._ptr, &val)

  property Cd:
    "float. Drag coefficient (no units)."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_cd_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_cd_set(self._ptr, &val)

  property surface_heat_flow:
    "float. surface heat flow (mW/m^2)."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_surface_heat_flow_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_surface_heat_flow_set(self._ptr, &val)

  property species_names:
    "List, shape (ng). The name of each species in the model"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_species_names_get_size(self._ptr, &dim1)
      cdef ndarray species_names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      wa_pxd.adiabatclimate_species_names_get(self._ptr, &dim1, <char *>species_names_c.data)
      return c2stringarr(species_names_c, S_STR_LEN, dim1)

  property particle_names:
    "List, shape (np). The name of each particle in the model"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_particle_names_get_size(self._ptr, &dim1)
      cdef ndarray particle_names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      wa_pxd.adiabatclimate_particle_names_get(self._ptr, &dim1, <char *>particle_names_c.data)
      return c2stringarr(particle_names_c, S_STR_LEN, dim1)

  property rad:
    "Radtran object that does radiative transfer"
    def __get__(self):
      var = Radtran()
      wa_pxd.adiabatclimate_rad_get(self._ptr, &var._ptr)
      return var
  
  property convecting_with_below:
    """ndarray[bool,ndim=1], shape (nz). If True, then the layer below 
    is convecting with the current layer. Index 1 determines if the 
    first atomspheric layer is convecting with the ground.
    """
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_convecting_with_below_get_size(self._ptr, &dim1)
      cdef ndarray[cbool, ndim=1] arr = np.empty(dim1, bool)
      wa_pxd.adiabatclimate_convecting_with_below_get(self._ptr, &dim1, <cbool *>arr.data)
      return arr

  property lapse_rate:
    "ndarray[double,ndim=1], shape (nz). The true lapse rate (dlnT/dlnP)."
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_lapse_rate_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_lapse_rate_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property lapse_rate_intended:
    "ndarray[double,ndim=1], shape (nz). The computed lapse rate (dlnT/dlnP)."
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_lapse_rate_intended_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_lapse_rate_intended_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property convective_newton_step_size:
    "float. Fraction of the Newton step used in convective classification (0..1)."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_convective_newton_step_size_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_convective_newton_step_size_set(self._ptr, &val)
  
  property convective_max_boundary_shift:
    "int. Boundary-motion limiter for convective mask updates. If < 0, no limiter is applied."
    def __get__(self):
      cdef int val
      wa_pxd.adiabatclimate_convective_max_boundary_shift_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      wa_pxd.adiabatclimate_convective_max_boundary_shift_set(self._ptr, &val)

  property prevent_overconvection:
    "bool. If true, shrink convective tops when a strong inversion exists above."
    def __get__(self):
      cdef bint val
      wa_pxd.adiabatclimate_prevent_overconvection_get(self._ptr, &val)
      return bool(val)
    def __set__(self, bint val):
      wa_pxd.adiabatclimate_prevent_overconvection_set(self._ptr, &val)

  property require_mode2:
    "bool. If true, require passing through mode 2 when mode 1 converges before optional mode 3 polishing."
    def __get__(self):
      cdef bint val
      wa_pxd.adiabatclimate_require_mode2_get(self._ptr, &val)
      return bool(val)
    def __set__(self, bint val):
      wa_pxd.adiabatclimate_require_mode2_set(self._ptr, &val)

  property rtol:
    "float. Relative tolerance of integration."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_rtol_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_rtol_set(self._ptr, &val)
  
  property atol:
    "float. Absolute tolerance of integration."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_atol_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_atol_set(self._ptr, &val)

  property tol_make_column:
    "float. Tolerance for nonlinear solve in make_column."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_tol_make_column_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_tol_make_column_set(self._ptr, &val)
    
  property epsj:
    "float. Perturbation for the jacobian"
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_epsj_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_epsj_set(self._ptr, &val)

  property xtol_rc:
    "float. xtol for RC equilibrium"
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_xtol_rc_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_xtol_rc_set(self._ptr, &val)

  property dt_increment:
    "float. Multiplicative growth factor for PTC timestep updates."
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_dt_increment_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_dt_increment_set(self._ptr, &val)

  property rce_solve_strategy:
    "int. RCE nonlinear solve strategy selector."
    def __get__(self):
      cdef int val
      wa_pxd.adiabatclimate_rce_solve_strategy_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      wa_pxd.adiabatclimate_rce_solve_strategy_set(self._ptr, &val)

  property max_rc_iters:
    "int. Max number of iterations in the RCE routine"
    def __get__(self):
      cdef int val
      wa_pxd.adiabatclimate_max_rc_iters_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      wa_pxd.adiabatclimate_max_rc_iters_set(self._ptr, &val)

  property max_rc_iters_convection:
    """int. Max number of iterations for which convective layers can
    be converged to radiative layers in the RCE routine
    """
    def __get__(self):
      cdef int val
      wa_pxd.adiabatclimate_max_rc_iters_convection_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      wa_pxd.adiabatclimate_max_rc_iters_convection_set(self._ptr, &val)

  property verbose:
    "bool. verbosity"
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_verbose_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_verbose_set(self._ptr, &val)

  property f_i_surf:
    "ndarray[double,ndim=1], shape (ng). Surface mixing ratios."
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_f_i_surf_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_f_i_surf_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property P_surf:
    "float. Surface pressure (dynes/cm^2)"
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_surf_get(self._ptr, &val)
      return val

  property P_trop:
    "float. Tropopause pressure (dynes/cm^2)"
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_trop_get(self._ptr, &val)
      return val

  property P:
    "ndarray[double,ndim=1], shape (nz). Pressure in each grid cell (dynes/cm^2)"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_p_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_p_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property T_surf:
    "float. Surface temperature (K)"
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_t_surf_get(self._ptr, &val)
      return val

  property T:
    "ndarray[double,ndim=1], shape (nz). Temperature in each grid cell (K)"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_t_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_t_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property f_i:
    """ndarray[double,ndim=2], shape (nz,ng). Mixing ratios of each species at 
    each atmospheric layer
    """
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_f_i_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_f_i_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property z:
    "ndarray[double,ndim=1], shape (nz). Altitude at the center of the grid cell (cm)"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_z_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_z_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property dz:
    "ndarray[double,ndim=1], shape (nz). Thickness of each grid cell (cm)"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_dz_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_dz_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property gravity_surf:
    "double. Surface gravity (cm/s^2)"
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_gravity_surf_get(self._ptr, &val)
      return val

  property gravity:
    "ndarray[double,ndim=1], shape (nz). Gravity at grid-cell centers (cm/s^2)"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_gravity_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_gravity_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property densities:
    "ndarray[double,ndim=2], shape (nz,ng). Densities in each grid cell (molecules/cm^3)"
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_densities_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_densities_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property pdensities:
    "ndarray[double,ndim=2], shape (nz,np). Particle densities in particles/cm^3."
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_pdensities_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_pdensities_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property pradii:
    "ndarray[double,ndim=2], shape (nz,np). Particle radii in cm."
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_pradii_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_pradii_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property N_atmos:
    "ndarray[double,ndim=1], shape (ng). Reservoir of gas in the atmosphere (mol/cm^2)"
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_n_atmos_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_n_atmos_get(self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property N_surface:
    """ndarray[double,ndim=1], shape (ng). Reservoir of gas on surface (mol/cm^2).
    Strictly consistent when `reference_pressure <= 0` (i.e., `planet_radius` at `P_surf`).
    """
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_n_surface_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_n_surface_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property N_ocean:
    """ndarray[double,ndim=2], shape (ng,ng). Reservoir of gas dissolved in oceans (mol/cm^2).
    There can be multiple oceans. The gases dissolved in ocean made of species 0 is given by `N_ocean[:,0]`.
    Strictly consistent when `reference_pressure <= 0` (i.e., `planet_radius` at `P_surf`).
    """
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_n_ocean_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1,dim2), np.double, order="f")
      wa_pxd.adiabatclimate_n_ocean_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
 
