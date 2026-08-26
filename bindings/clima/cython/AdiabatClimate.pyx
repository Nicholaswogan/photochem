
cimport AdiabatClimate_pxd as wa_pxd
  
cdef class AdiabatClimate:
  """One-dimensional multispecies adiabatic and radiative-convective climate model.

  The model constructs a pseudoadiabatic troposphere connected to an
  isothermal stratosphere and can compute top-of-atmosphere radiative fluxes
  or solve for radiative-convective equilibrium. Arrays indexed by gas use
  ``species_names`` order; arrays indexed by particle use ``particle_names``.
  Atmospheric profiles run from the surface upward.

  Array-valued properties return copies. The ``rad`` property and its nested
  channel and work objects are non-owning views; keep this model alive while
  using them.
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
    """Initialize an ``AdiabatClimate`` model from its input files.

    Construction allocates the model and configurable defaults but does not
    construct a physical atmospheric profile. Profile-building methods update
    the atmospheric state.

    Parameters
    ----------
    species_file : str
        Species YAML file.
    settings_file : str
        Climate-settings YAML file.
    flux_file : str
        Stellar-flux text file.
    data_dir : str, optional
        Directory containing radiative-transfer data. By default, use the
        installed ``photochem_clima_data`` package.
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
    """Construct a multispecies pseudoadiabat connected to an isothermal stratosphere.

    The pseudoadiabat follows Equation 1 of Graham et al. (2021, PSJ).
    Atmospheric state fields, particle profiles, reservoirs, lapse rates, and
    the convective mask are updated on success.

    Parameters
    ----------
    T_surf : float
        Surface temperature (K).
    P_i_surf : ndarray, shape (ng,)
        Surface partial pressures in ``species_names`` order (dyn/cm^2).
    """
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_profile_wrapper(self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def make_profile_dry(self, ndarray[double, ndim=1] P, ndarray[double, ndim=1] T, ndarray[double, ndim=2] f_i):
    """Construct a dry atmosphere from prescribed pressure, temperature, and composition.

    Condensation is not applied, so supersaturated gas is retained. All
    atmospheric state fields except ``P_trop`` and ``convecting_with_below``
    are updated on success.

    Parameters
    ----------
    P : ndarray, shape (n,)
        Pressure ordered from the surface upward (dyn/cm^2).
    T : ndarray, shape (n,)
        Temperature defined on ``P`` (K).
    f_i : ndarray, shape (n, ng)
        Gas volume mixing ratios in ``species_names`` order.
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
    """Construct a profile whose total gas reservoirs match prescribed columns.

    This method repeatedly calls ``make_profile`` while solving for surface
    partial pressures. On success, ``make_column_P_guess`` is updated to the
    solved partial pressures.

    Parameters
    ----------
    T_surf : float
        Surface temperature (K).
    N_i_surf : ndarray, shape (ng,)
        Total gas reservoirs in ``species_names`` order (mol/cm^2), including
        atmospheric, surface, and ocean reservoirs.
    """
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_column_wrapper(self._ptr, &T_surf,
    &ng, <double *>N_i_surf.data, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
  
  def make_profile_bg_gas(self, double T_surf, ndarray[double, ndim=1] P_i_surf, double P_surf, str bg_gas):
    """Construct a profile with a named background gas and fixed surface pressure.

    The background gas partial pressure is adjusted until the profile's total
    surface pressure equals ``P_surf``. Other entries of ``P_i_surf`` are retained.

    Parameters
    ----------
    T_surf : float
        Surface temperature (K).
    P_i_surf : ndarray, shape (ng,)
        Initial surface partial pressures in ``species_names`` order (dyn/cm^2).
    P_surf : float
        Target total surface pressure (dyn/cm^2).
    bg_gas : str
        Background-gas name. It must occur in ``species_names``.
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
    """Construct a pseudoadiabatic profile and return its top-of-atmosphere fluxes.

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    P_i_surf : ndarray, shape (ng,)
        Surface partial pressures in ``species_names`` order (dyn/cm^2).

    Returns
    -------
    tuple[float, float]
        Incoming stellar radiation and outgoing longwave radiation at the top
        of the atmosphere, respectively (mW/m^2).
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
    """Construct a column-reservoir profile and return top-of-atmosphere fluxes.

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    N_i_surf : ndarray, shape (ng,)
        Total gas reservoirs in ``species_names`` order (mol/cm^2).

    Returns
    -------
    tuple[float, float]
        Incoming stellar radiation and outgoing longwave radiation at the top
        of the atmosphere, respectively (mW/m^2).
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
    """Construct a fixed-pressure background-gas profile and return TOA fluxes.

    Parameters
    ----------
    T_surf : float
        The surface temperature (K)
    P_i_surf : ndarray, shape (ng,)
        Initial surface partial pressures in ``species_names`` order (dyn/cm^2).
    P_surf : float
        Target total surface pressure (dyn/cm^2).
    bg_gas : str
        The name of the background gas

    Returns
    -------
    tuple[float, float]
        Incoming stellar radiation and outgoing longwave radiation at the top
        of the atmosphere, respectively (mW/m^2).
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
    """Construct a dry prescribed profile and return top-of-atmosphere fluxes.

    Parameters
    ----------
    P : ndarray, shape (n,)
        Pressure ordered from the surface upward (dyn/cm^2).
    T : ndarray, shape (n,)
        Temperature defined on ``P`` (K).
    f_i : ndarray, shape (n, ng)
        Gas volume mixing ratios in ``species_names`` order.

    Returns
    -------
    tuple[float, float]
        Incoming stellar radiation and outgoing longwave radiation at the top
        of the atmosphere, respectively (mW/m^2).
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
    """Solve for radiative-equilibrium surface temperature using ``make_profile``.

    If ``solve_for_T_trop`` is true, the tropopause temperature is solved
    simultaneously. If ``tidally_locked_dayside`` is true, the dayside
    redistribution correction is included. The converged atmospheric and
    radiative state remains stored in the model.

    Parameters
    ----------
    P_i_surf : ndarray, shape (ng,)
        Surface partial pressures in ``species_names`` order (dyn/cm^2).
    T_guess : float, optional
        Initial surface-temperature guess (K), by default 280.

    Returns
    -------
    float
        Equilibrium surface temperature (K).
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
    """Solve for radiative-equilibrium surface temperature using ``make_column``.

    Other solver behavior matches ``surface_temperature``.

    Parameters
    ----------
    N_i_surf : ndarray, shape (ng,)
        Total gas reservoirs in ``species_names`` order (mol/cm^2).
    T_guess : float, optional
        Initial surface-temperature guess (K), by default 280.

    Returns
    -------
    float
        Equilibrium surface temperature (K).
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
    """Set the solubility callback for an ocean made from a modeled species.

    The callback receives surface temperature, the number of gases, gas partial
    pressures in bar, an output array for solubilities in mol/kg, and
    ``ocean_args_p``. Keep the compiled callback alive while the model uses it.

    Parameters
    ----------
    species : str
        Name of the species that makes up the ocean.
    fcn : numba.cfunc or None
        Compiled callback with C signature
        ``void(double, int32, double*, double*, void*)``. Pass ``None`` to
        clear the callback for this ocean.
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
    """Solve for surface temperature with a background gas and fixed pressure.

    Profile construction uses ``make_profile_bg_gas``; other solver behavior
    matches ``surface_temperature``.

    Parameters
    ----------
    P_i_surf : ndarray, shape (ng,)
        Initial surface partial pressures in ``species_names`` order (dyn/cm^2).
    P_surf : float
        Target total surface pressure (dyn/cm^2).
    bg_gas : str
        The name of the background gas
    T_guess : float, optional
        Initial surface-temperature guess (K), by default 280.

    Returns
    -------
    float
        Equilibrium surface temperature (K).
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
    """Set pressure-dependent particle number-density and radius profiles.

    Inputs are retained as log-pressure interpolation functions and sampled
    onto the model grid whenever a profile is constructed. ``P`` must be
    strictly decreasing. Values outside its range are held at the nearest
    endpoint. Zero values are represented internally by a small positive floor.

    Parameters
    ----------
    P : ndarray, shape (n,)
        Pressure grid (dyn/cm^2).
    pdensities : ndarray, shape (n, np)
        Particle number densities (particles/cm^3).
    pradii : ndarray, shape (n, np)
        Particle radii (cm).
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
    """Solve for full radiative-convective equilibrium.

    The solution and convective classification are stored in the atmospheric
    state properties. The return value reports convergence; invalid inputs and
    model failures raise ``ClimaException``.

    Parameters
    ----------
    P_i_surf : ndarray, shape (ng,)
        Surface partial pressures in ``species_names`` order (dyn/cm^2).
    T_surf_guess : float
        Initial surface-temperature guess (K).
    T_guess : ndarray, shape (nz,)
        Initial layer-temperature guess (K).
    convecting_with_below : ndarray of bool, shape (nz,), optional
        Initial convective connectivity. Element 0 indicates whether the lowest
        atmospheric layer convects with the surface.
    custom_dry_mix : dict[str, ndarray], optional
        Vertically inhomogeneous mixing ratios for dry gases. The dictionary
        must contain a ``"pressure"`` key specifying a strictly decreasing
        pressure grid in dyn/cm^2.
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
        True if the RCE solve converged.
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
    """Regrid the current atmosphere to equal-width altitude layers.

    Gas number densities are conservatively rebinned, temperature is
    interpolated, and pressure and mixing ratios are recomputed. This method
    mutates the stored grid and atmospheric state.
    """
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_to_regular_grid_wrapper(self._ptr, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def out2atmosphere_txt(self, str filename, ndarray[double, ndim=1] eddy, int number_of_decimals=5, cbool overwrite = False, cbool clip = True):
    """Write the current atmosphere in Photochem's legacy text format.

    This method first calls ``to_regular_grid`` and therefore mutates the stored
    atmosphere, including when a later file validation or write error occurs.
    Output altitude is in km, pressure in bar, total gas density in
    molecules/cm^3, temperature in K, and eddy diffusion in cm^2/s.

    Parameters
    ----------
    filename : str
        Output filename
    eddy : ndarray, shape (nz,)
        Eddy diffusion to write (cm^2/s).
    number_of_decimals : int, optional
        Number of decimal places, from 2 through 17; by default 5.
    overwrite : bool, optional
        Allow replacement of an existing file, by default False.
    clip : bool, optional
        Clip output mixing ratios below ``1.0e-40``, by default True.
    """  
    cdef int nz = eddy.shape[0]
    cdef bytes filename_b = pystring2cstring(filename)
    cdef char *filename_c = filename_b
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_out2atmosphere_txt_wrapper(self._ptr, filename_c, &nz, <double *>eddy.data, &number_of_decimals, &overwrite, &clip, err)  
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def heat_redistribution_parameters(self):
    """Compute tidally locked heat-redistribution parameters.

    This implements Equation 10 of Koll (2022, ApJ). Call a top-of-atmosphere
    flux method first because this calculation uses the current atmospheric
    state and radiative-transfer results.

    Returns
    -------
    tuple[float, float, float]
        ``(tau_LW, k_term, f_term)``: the dimensionless Planck-weighted
        longwave optical depth, Equation 10's dimensionless ``k`` term, and
        dimensionless heat-redistribution factor ``f``.
    """
    cdef char err[ERR_LEN+1];
    cdef double tau_LW, k_term, f_term;
    wa_pxd.adiabatclimate_heat_redistribution_parameters_wrapper(self._ptr, &tau_LW, &k_term, &f_term, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return tau_LW, k_term, f_term

  property P_top:
    """float: Writable target pressure at the model top (dyn/cm^2)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_top_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_p_top_set(self._ptr, &val)

  property reference_pressure:
    """float: Writable pressure at which the configured planet radius is defined.

    Units are dyn/cm^2. If nonpositive, the radius is defined at ``P_surf``.
    """
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_reference_pressure_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_reference_pressure_set(self._ptr, &val)

  property T_trop:
    """float: Writable isothermal-stratosphere temperature (K)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_t_trop_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_t_trop_set(self._ptr, &val)

  property use_make_column_P_guess:
    """bool: Whether ``make_column`` first uses ``make_column_P_guess``.
    """
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_use_make_column_p_guess_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_use_make_column_p_guess_set(self._ptr, &val)

  property make_column_P_guess:
    """ndarray, shape (ng,): Writable surface-partial-pressure guess.

    Values are in dyn/cm^2 and follow ``species_names`` order. Getting this
    property returns a copy; assign the complete array to update the model.
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
        raise ClimaException('make_column_P_guess must have shape (number_of_gases,).')
      wa_pxd.adiabatclimate_make_column_p_guess_set(self._ptr, &dim1, <double *>arr.data)

  property solve_for_T_trop:
    """bool: Whether equilibrium solvers also solve the tropopause temperature.

    When enabled, the target is the radiative skin temperature and ``T_trop``
    supplies the initial guess.
    """
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_solve_for_t_trop_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_solve_for_t_trop_set(self._ptr, &val)

  property albedo_fcn:
    """numba.cfunc or None: Write-only temperature-dependent albedo callback.

    The callback can parameterize ice-albedo feedback and must have C signature
    ``double(double)`` with temperature in K and dimensionless albedo returned.
    Keep the compiled callback alive while the model uses it. Assign ``None``
    to clear it.
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
    """ndarray, shape (ng,): Writable relative humidity of each gas.

    Values follow ``species_names`` order. Getting this property returns a
    copy; assign the complete array to update the model.
    """
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
        raise ClimaException('RH must have shape (number_of_gases,).')
      wa_pxd.adiabatclimate_rh_set(self._ptr, &dim1, <double *>arr.data)

  property ocean_args_p:
    """int or None: Write-only address passed to ocean-solubility callbacks.

    The pointed-to data and its lifetime are managed by the caller. Assign
    ``None`` to pass a null pointer.
    """
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
    """bool: Whether equilibrium solvers apply tidally locked dayside redistribution.
    """
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_tidally_locked_dayside_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_tidally_locked_dayside_set(self._ptr, &val)

  property L:
    """float: Writable circulation horizontal scale (cm)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_l_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_l_set(self._ptr, &val)

  property chi:
    """float: Writable dimensionless heat-engine efficiency."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_chi_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_chi_set(self._ptr, &val)

  property n_LW:
    """float: Writable dimensionless longwave-opacity exponent, normally 1 or 2."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_n_lw_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_n_lw_set(self._ptr, &val)

  property Cd:
    """float: Writable dimensionless drag coefficient."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_cd_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_cd_set(self._ptr, &val)

  property surface_heat_flow:
    """float: Writable heat flow from the surface into the atmosphere (mW/m^2)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_surface_heat_flow_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_surface_heat_flow_set(self._ptr, &val)

  property species_names:
    """list[str]: Gas names in the ordering used by gas-indexed arrays."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_species_names_get_size(self._ptr, &dim1)
      cdef ndarray species_names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      wa_pxd.adiabatclimate_species_names_get(self._ptr, &dim1, <char *>species_names_c.data)
      return c2stringarr(species_names_c, S_STR_LEN, dim1)

  property particle_names:
    """list[str]: Particle names in the ordering used by particle-indexed arrays."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_particle_names_get_size(self._ptr, &dim1)
      cdef ndarray particle_names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      wa_pxd.adiabatclimate_particle_names_get(self._ptr, &dim1, <char *>particle_names_c.data)
      return c2stringarr(particle_names_c, S_STR_LEN, dim1)

  property rad:
    """Radtran: Non-owning view of the model's radiative-transfer object.

    Keep this ``AdiabatClimate`` alive while using the returned object.
    """
    def __get__(self):
      var = Radtran()
      wa_pxd.adiabatclimate_rad_get(self._ptr, &var._ptr)
      return var
  
  property convecting_with_below:
    """ndarray, shape (nz,): Read-only copy of convective connectivity.

    A true element means that the corresponding layer convects with the layer
    below. Element 0 refers to the lowest layer and the surface.
    """
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_convecting_with_below_get_size(self._ptr, &dim1)
      cdef ndarray[cbool, ndim=1] arr = np.empty(dim1, bool)
      wa_pxd.adiabatclimate_convecting_with_below_get(self._ptr, &dim1, <cbool *>arr.data)
      return arr

  property lapse_rate:
    """ndarray, shape (nz,): Read-only copy of actual ``dln(T)/dln(P)``."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_lapse_rate_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_lapse_rate_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property lapse_rate_intended:
    """ndarray, shape (nz,): Read-only copy of adiabatic target ``dln(T)/dln(P)``."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_lapse_rate_intended_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_lapse_rate_intended_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property convective_newton_step_size:
    """float: Writable fraction of the Newton step used for convective classification."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_convective_newton_step_size_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_convective_newton_step_size_set(self._ptr, &val)

  property convective_hysteresis_frac_on:
    """float: Writable fractional threshold for making RCE layers convective.

    The applied threshold is the maximum of ``convective_hysteresis_min`` and
    this value times the magnitude of ``lapse_rate_intended``.
    """
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_convective_hysteresis_frac_on_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_convective_hysteresis_frac_on_set(self._ptr, &val)

  property convective_hysteresis_frac_off:
    """float: Writable fractional threshold for making RCE layers radiative.

    The applied threshold is the maximum of ``convective_hysteresis_min`` and
    this value times the magnitude of ``lapse_rate_intended``.
    """
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_convective_hysteresis_frac_off_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_convective_hysteresis_frac_off_set(self._ptr, &val)

  property convective_hysteresis_min:
    """float: Writable absolute ``dln(T)/dln(P)`` hysteresis threshold."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_convective_hysteresis_min_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_convective_hysteresis_min_set(self._ptr, &val)
  
  property convective_max_boundary_shift:
    """int: Writable convective-boundary motion limit in layers.

    Negative values disable the limiter.
    """
    def __get__(self):
      cdef int val
      wa_pxd.adiabatclimate_convective_max_boundary_shift_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      wa_pxd.adiabatclimate_convective_max_boundary_shift_set(self._ptr, &val)

  property prevent_overconvection:
    """bool: Whether strong inversions shrink the tops of convective zones."""
    def __get__(self):
      cdef bint val
      wa_pxd.adiabatclimate_prevent_overconvection_get(self._ptr, &val)
      return bool(val)
    def __set__(self, bint val):
      wa_pxd.adiabatclimate_prevent_overconvection_set(self._ptr, &val)

  property require_mode2:
    """bool: Whether RCE must pass through mode 2 before optional mode 3 polishing."""
    def __get__(self):
      cdef bint val
      wa_pxd.adiabatclimate_require_mode2_get(self._ptr, &val)
      return bool(val)
    def __set__(self, bint val):
      wa_pxd.adiabatclimate_require_mode2_set(self._ptr, &val)

  property rtol:
    """float: Writable relative tolerance for adiabatic-profile integrations."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_rtol_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_rtol_set(self._ptr, &val)
  
  property atol:
    """float: Writable absolute tolerance for adiabatic-profile integrations."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_atol_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_atol_set(self._ptr, &val)

  property tol_make_column:
    """float: Writable nonlinear-solver tolerance used by ``make_column``."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_tol_make_column_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_tol_make_column_set(self._ptr, &val)
    
  property epsj:
    """float: Writable temperature perturbation for the RCE Jacobian (K)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_epsj_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_epsj_set(self._ptr, &val)

  property xtol_rc:
    """float: Writable nonlinear-solver step tolerance for RCE."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_xtol_rc_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_xtol_rc_set(self._ptr, &val)

  property dt_increment:
    """float: Writable multiplicative growth factor for PTC timestep updates."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_dt_increment_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_dt_increment_set(self._ptr, &val)

  property rce_solve_strategy:
    """int: Writable RCE nonlinear-solver strategy.

    Use ``RCE_SOLVE_HYBRJ_ONLY``, ``RCE_SOLVE_PTC_THEN_HYBRJ``, or
    ``RCE_SOLVE_HYBRJ_THEN_PTC_THEN_HYBRJ`` from ``photochem.clima``.
    """
    def __get__(self):
      cdef int val
      wa_pxd.adiabatclimate_rce_solve_strategy_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      wa_pxd.adiabatclimate_rce_solve_strategy_set(self._ptr, &val)

  property max_rc_iters:
    """int: Writable maximum number of outer RCE iterations."""
    def __get__(self):
      cdef int val
      wa_pxd.adiabatclimate_max_rc_iters_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      wa_pxd.adiabatclimate_max_rc_iters_set(self._ptr, &val)

  property max_rc_iters_convection:
    """int: Writable RCE iterations allowed to convert convective layers to radiative."""
    def __get__(self):
      cdef int val
      wa_pxd.adiabatclimate_max_rc_iters_convection_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      wa_pxd.adiabatclimate_max_rc_iters_convection_set(self._ptr, &val)

  property compute_solar_in_jac:
    """bool: Whether the RCE Jacobian recomputes shortwave radiative transfer.

    When false, the Jacobian reuses the base-state shortwave calculation for
    each temperature perturbation.
    """
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_compute_solar_in_jac_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_compute_solar_in_jac_set(self._ptr, &val)

  property verbose:
    """bool: Whether RCE writes progress information."""
    def __get__(self):
      cdef cbool val
      wa_pxd.adiabatclimate_verbose_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      wa_pxd.adiabatclimate_verbose_set(self._ptr, &val)

  property f_i_surf:
    """ndarray, shape (ng,): Read-only copy of surface gas volume mixing ratios."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_f_i_surf_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_f_i_surf_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property P_surf:
    """float: Surface pressure (dyn/cm^2)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_surf_get(self._ptr, &val)
      return val

  property P_trop:
    """float: Tropopause pressure (dyn/cm^2)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_trop_get(self._ptr, &val)
      return val

  property P:
    """ndarray, shape (nz,): Read-only copy of layer-center pressure (dyn/cm^2)."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_p_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_p_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property T_surf:
    """float: Surface temperature (K)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_t_surf_get(self._ptr, &val)
      return val

  property T:
    """ndarray, shape (nz,): Read-only copy of layer-center temperature (K)."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_t_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_t_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property f_i:
    """ndarray, shape (nz, ng): Read-only copy of gas volume mixing ratios."""
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_f_i_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_f_i_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property z:
    """ndarray, shape (nz,): Read-only copy of layer-center altitude (cm)."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_z_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_z_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property dz:
    """ndarray, shape (nz,): Read-only copy of layer thickness (cm)."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_dz_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_dz_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property gravity_surf:
    """float: Surface gravitational acceleration (cm/s^2)."""
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_gravity_surf_get(self._ptr, &val)
      return val

  property gravity:
    """ndarray, shape (nz,): Read-only copy of layer-center gravity (cm/s^2)."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_gravity_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_gravity_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property densities:
    """ndarray, shape (nz, ng): Read-only copy of gas number densities (molecules/cm^3)."""
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_densities_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_densities_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property pdensities:
    """ndarray, shape (nz, np): Read-only copy of particle number densities (particles/cm^3)."""
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_pdensities_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_pdensities_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property pradii:
    """ndarray, shape (nz, np): Read-only copy of particle radii (cm)."""
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_pradii_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_pradii_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property N_atmos:
    """ndarray, shape (ng,): Read-only copy of atmospheric gas reservoirs (mol/cm^2)."""
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_n_atmos_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_n_atmos_get(self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property N_surface:
    """ndarray, shape (ng,): Read-only copy of surface gas reservoirs (mol/cm^2).

    Strictly consistent when ``reference_pressure <= 0``, meaning that the
    configured planet radius is defined at ``P_surf``.
    """
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_n_surface_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_n_surface_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property N_ocean:
    """ndarray, shape (ng, ng): Read-only copy of ocean gas reservoirs (mol/cm^2).

    ``N_ocean[:, j]`` contains gases dissolved in the ocean made from species
    ``j``. Values are strictly consistent when ``reference_pressure <= 0``.
    """
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_n_ocean_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1,dim2), np.double, order="f")
      wa_pxd.adiabatclimate_n_ocean_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
 
