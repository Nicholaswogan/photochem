
cimport EvoAtmosphere_pxd as ea_pxd

cdef class EvoAtmosphere:
  """A photochemical model which assumes no background gas. Once initialized,
  this class can integrate an atmosphere forward in time to a steady
  state.
  """

  cdef ea_pxd.EvoAtmosphere *_ptr
  cdef bool _init_called

  def __cinit__(self, *args, **kwargs):
    self._init_called = False
    self._ptr = ea_pxd.allocate_evoatmosphere()

  def __dealloc__(self):
    ea_pxd.deallocate_evoatmosphere(self._ptr)

  def __getattribute__(self, name):
    if not self._init_called:
      raise PhotoException('The "__init__" method of EvoAtmosphere has not been called.')
    return super().__getattribute__(name)

  def __setattr__(self, name, value):
    if not self._init_called:
      raise PhotoException('The "__init__" method of EvoAtmosphere has not been called.')
    PyObject_GenericSetAttr(self, name, value)

  def __init__(self, str mechanism_file, str settings_file,
               str flux_file, str atmosphere_txt=None, data_dir=None):
    """Configure a photochemical model and optionally initialize its atmosphere.

    When ``atmosphere_txt`` is omitted, static model setup is completed but
    atmospheric operations remain unavailable until one of the explicit
    atmosphere initialization methods succeeds.

    Parameters
    ----------
    mechanism_file : str
        Path to the reaction mechanism file (yaml format).
    settings_file : str
        Path to the settings file (yaml format).
    flux_file : str
        Path to the file describing the stellar flux.
    atmosphere_txt : str, optional
        Path to the file containing altitude, total number density, temperature, 
        eddy diffusion, initial concentrations of each gas (mixing ratios), 
        and particle radii. If omitted, the model is left uninitialized.
    data_dir : str, optional
        Path to the data directory containing photolysis cross sections and other data
        needed to run the model
    """
    self._init_called = True

    if data_dir == None:
      data_dir_ = photochem_clima_data.DATA_DIR
    else:
      data_dir_ = data_dir
    
    # convert strings to char
    cdef bytes mechanism_file_b = pystring2cstring(mechanism_file)
    cdef char *mechanism_file_c = mechanism_file_b
    cdef bytes settings_file_b = pystring2cstring(settings_file)
    cdef char *settings_file_c = settings_file_b
    cdef bytes flux_file_b = pystring2cstring(flux_file)
    cdef char *flux_file_c = flux_file_b
    cdef bytes atmosphere_txt_b = b''
    cdef char *atmosphere_txt_c = atmosphere_txt_b
    cdef bool atmosphere_txt_present = atmosphere_txt is not None
    cdef bytes data_dir_b = pystring2cstring(data_dir_)
    cdef char *data_dir_c = data_dir_b
    cdef char err[ERR_LEN+1]
    
    # Configure static model state and optionally initialize the atmosphere.
    if atmosphere_txt_present:
      atmosphere_txt_b = pystring2cstring(atmosphere_txt)
      atmosphere_txt_c = atmosphere_txt_b
    ea_pxd.evoatmosphere_create_wrapper(
      self._ptr, mechanism_file_c, settings_file_c, flux_file_c,
      atmosphere_txt_c, &atmosphere_txt_present, data_dir_c, err
    )
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  property atmosphere_initialized:
    """Whether atmosphere-dependent model state has been initialized."""
    def __get__(self):
      cdef bool val
      ea_pxd.evoatmosphere_atmosphere_initialized_get(self._ptr, &val)
      return val

  cdef tuple _prepare_atmosphere_composition(self, int nprofile, dict mix,
                                              particle_radius,
                                              double default_mix,
                                              double default_particle_radius):
    """Convert species-keyed profiles to arrays in Fortran species order."""
    cdef int nq = self.dat.nq
    cdef int nparticles = self.dat.np
    cdef ndarray mix_ = np.full((nq,nprofile), default_mix, dtype=np.double,
                                order='F')
    cdef ndarray particle_radius_ = np.full(
        (nparticles,nprofile), default_particle_radius, dtype=np.double,
        order='F'
    )
    cdef ndarray explicit_gas_total = np.zeros(nprofile, dtype=np.double)
    cdef ndarray profile
    cdef int i

    if not np.isfinite(default_mix) or default_mix < 0.0:
      raise PhotoException('default_mix must be finite and nonnegative.')
    if nparticles > 0 and (not np.isfinite(default_particle_radius) or
                           default_particle_radius <= 0.0):
      raise PhotoException('default_particle_radius must be finite and positive.')

    species_names = self.dat.species_names[:nq]
    particle_names = species_names[:nparticles]
    gas_names = species_names[nparticles:]
    unknown_mix = set(mix) - set(species_names)
    if unknown_mix:
      raise PhotoException('Unknown species in mix: '+', '.join(sorted(unknown_mix)))
    if not (set(mix) & set(gas_names)):
      raise PhotoException('mix must explicitly contain at least one gas species.')

    for i, name in enumerate(species_names):
      if name in mix:
        profile = np.asarray(mix[name], dtype=np.double)
        if profile.ndim != 1 or profile.size != nprofile:
          raise PhotoException("Mixing-ratio profile for '"+name+"' has the wrong shape.")
        if not np.all(np.isfinite(profile)) or np.any(profile < 0.0):
          raise PhotoException("Mixing-ratio profile for '"+name+"' must be finite and nonnegative.")
        mix_[i,:] = profile
        if i >= nparticles:
          explicit_gas_total += profile
    if np.any(explicit_gas_total <= 0.0):
      raise PhotoException('Explicitly supplied gases must have a positive total mixing ratio at every profile point.')

    if particle_radius is not None:
      unknown_radii = set(particle_radius) - set(particle_names)
      if unknown_radii:
        raise PhotoException('Unknown particles in particle_radius: '+', '.join(sorted(unknown_radii)))
      for i, name in enumerate(particle_names):
        if name in particle_radius:
          profile = np.asarray(particle_radius[name], dtype=np.double)
          if profile.ndim != 1 or profile.size != nprofile:
            raise PhotoException("Particle-radius profile for '"+name+"' has the wrong shape.")
          if not np.all(np.isfinite(profile)) or np.any(profile <= 0.0):
            raise PhotoException("Particle-radius profile for '"+name+"' must be finite and positive.")
          particle_radius_[i,:] = profile

    return mix_, particle_radius_

  def initialize_atmosphere_z(self, ndarray[double, ndim=1] z,
                              ndarray[double, ndim=1] temperature,
                              ndarray[double, ndim=1] edd,
                              double surface_pressure, dict mix,
                              particle_radius=None,
                              double default_mix=1.0e-40,
                              double default_particle_radius=1.0e-5):
    """Initialize the atmosphere from altitude-based profiles.

    Temperature is interpolated linearly in altitude. Eddy diffusion,
    mixing ratios, and particle radii are interpolated in log10 space. Gas
    mixing ratios are normalized internally, and total number density is
    derived by hydrostatic integration upward from ``surface_pressure``.

    The first altitude must be zero, and the final altitude defines the top of
    the model domain. ``mix`` may contain any subset of evolved species;
    omitted species use ``default_mix``. At least one gas must be supplied
    explicitly and the explicitly supplied gases must have a positive total
    mixing ratio at every altitude. Likewise, ``particle_radius`` may contain
    any subset of particles, with omitted radii using
    ``default_particle_radius``.

    Fixed-density and fixed-partial-pressure lower boundary conditions
    override the corresponding bottom-layer mixing ratios. Successful
    initialization destroys any active integrator and disables a persistent
    pressure-based temperature--eddy profile. If initialization fails, the
    existing atmospheric state and integrator are retained.

    Parameters
    ----------
    z : ndarray, shape (nprofile,)
        Altitude profile knots in cm, including both domain edges.
    temperature : ndarray, shape (nprofile,)
        Temperature at each altitude knot in K.
    edd : ndarray, shape (nprofile,)
        Eddy diffusion at each altitude knot in cm^2/s.
    surface_pressure : float
        Total pressure at the lower domain edge in dyn/cm^2.
    mix : dict[str, ndarray]
        Mixing-ratio profiles keyed by species name. Unspecified species use
        ``default_mix``.
    particle_radius : dict[str, ndarray], optional
        Particle-radius profiles in cm, keyed by particle name. Unspecified
        particles use ``default_particle_radius``.
    default_mix : float, optional
        Mixing ratio assigned to species omitted from ``mix``. The default is
        1.0e-40.
    default_particle_radius : float, optional
        Radius in cm assigned to particles omitted from ``particle_radius``.
        The default is 1.0e-5 cm, matching legacy atmosphere-file behavior.
    """
    cdef ndarray z_ = np.ascontiguousarray(z, dtype=np.double)
    cdef ndarray temperature_ = np.ascontiguousarray(temperature, dtype=np.double)
    cdef ndarray edd_ = np.ascontiguousarray(edd, dtype=np.double)
    cdef int nprofile = z_.size
    cdef int nq
    cdef int nparticles
    cdef ndarray mix_
    cdef ndarray particle_radius_
    cdef char err[ERR_LEN+1]

    if temperature_.size != nprofile or edd_.size != nprofile:
      raise PhotoException('z, temperature, and edd must have the same length.')

    mix_, particle_radius_ = self._prepare_atmosphere_composition(
      nprofile, mix, particle_radius, default_mix, default_particle_radius
    )
    nq = mix_.shape[0]
    nparticles = particle_radius_.shape[0]

    ea_pxd.evoatmosphere_initialize_atmosphere_z_wrapper(
      self._ptr, &nprofile, <double *>z_.data,
      <double *>temperature_.data, <double *>edd_.data, &surface_pressure,
      &nq, <double *>mix_.data, &nparticles,
      <double *>particle_radius_.data, err
    )
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def initialize_atmosphere_p(self, ndarray[double, ndim=1] pressure,
                              ndarray[double, ndim=1] temperature,
                              ndarray[double, ndim=1] edd, dict mix,
                              particle_radius=None, bint persistent=False,
                              tropopause_pressure=None,
                              double default_mix=1.0e-40,
                              double default_particle_radius=1.0e-5,
                              maintain_toa_pressure=None, target_pressure=None):
    """Initialize the atmosphere from pressure-based profiles.

    The pressure profile must be strictly decreasing. Its first and last
    points define the lower and upper boundaries of the model domain, and the
    configured planet radius is interpreted at the lower boundary. Altitude
    is constructed hydrostatically using gas composition only.

    ``mix`` and ``particle_radius`` use the same partial-dictionary convention
    as :meth:`initialize_atmosphere_z`. Omitted species and particle radii use
    ``default_mix`` and ``default_particle_radius``, respectively.

    By default, pressure is used only to construct the initial atmosphere. If
    ``persistent`` is true, temperature and eddy diffusion are retained as
    functions of hydrostatic pressure and reapplied during every atmospheric
    preparation. Composition remains part of the evolving ODE state. When gas
    rainout and persistence are both enabled, ``tropopause_pressure`` is
    required.

    Persistent initialization enables approximate TOA-pressure maintenance by
    default. Use ``maintain_toa_pressure=False`` to disable it, or supply
    ``target_pressure`` in dyn/cm^2 to select its target. These options are
    only valid with ``persistent=True``.

    The requested initial pressure grid is retained even when its TOA pressure
    lies outside the maintenance band. If a robust stepper is initialized,
    its startup preflight performs any required regrid before CVODE starts.

    Successful initialization destroys any active integrator. If
    initialization fails, the existing atmosphere, persistent profile, and
    integrator are retained.

    Parameters
    ----------
    pressure : ndarray, shape (nprofile,)
        Strictly decreasing pressure knots in dyn/cm^2, including both domain
        boundaries.
    temperature : ndarray, shape (nprofile,)
        Temperature at each pressure knot in K.
    edd : ndarray, shape (nprofile,)
        Eddy diffusion at each pressure knot in cm^2/s.
    mix : dict[str, ndarray]
        Mixing-ratio profiles keyed by species name. Unspecified species use
        ``default_mix``.
    particle_radius : dict[str, ndarray], optional
        Particle-radius profiles in cm, keyed by particle name. Unspecified
        particles use ``default_particle_radius``.
    persistent : bool, optional
        Retain temperature and eddy diffusion as functions of hydrostatic
        pressure. The default is false.
    tropopause_pressure : float, optional
        Tropopause pressure in dyn/cm^2. May be supplied only when
        ``persistent`` is true and is required when gas rainout is enabled.
    default_mix : float, optional
        Mixing ratio assigned to species omitted from ``mix``. The default is
        1.0e-40.
    default_particle_radius : float, optional
        Radius in cm assigned to particles omitted from ``particle_radius``.
        The default is 1.0e-5 cm, matching legacy atmosphere-file behavior.
    maintain_toa_pressure : bool, optional
        Enable approximate TOA-pressure maintenance when ``persistent`` is
        true. Defaults to true for persistent initialization.
    target_pressure : float, optional
        Target TOA pressure for approximate maintenance in dyn/cm^2. The
        default is 0.1 dyn/cm^2.
    """
    cdef ndarray pressure_ = np.ascontiguousarray(pressure, dtype=np.double)
    cdef ndarray temperature_ = np.ascontiguousarray(temperature, dtype=np.double)
    cdef ndarray edd_ = np.ascontiguousarray(edd, dtype=np.double)
    cdef int nprofile = pressure_.size
    cdef int nq
    cdef int nparticles
    cdef ndarray mix_
    cdef ndarray particle_radius_
    cdef bool persistent_ = persistent
    cdef double tropopause_pressure_ = 0.0
    cdef bool tropopause_pressure_present = False
    cdef bool maintain_toa_pressure_ = True
    cdef bool maintain_toa_pressure_present = False
    cdef double target_pressure_ = 0.1
    cdef bool target_pressure_present = False
    cdef char err[ERR_LEN+1]

    if temperature_.size != nprofile or edd_.size != nprofile:
      raise PhotoException('pressure, temperature, and edd must have the same length.')
    if tropopause_pressure is not None:
      tropopause_pressure_present = True
      tropopause_pressure_ = tropopause_pressure
    if maintain_toa_pressure is not None:
      maintain_toa_pressure_present = True
      maintain_toa_pressure_ = maintain_toa_pressure
    if target_pressure is not None:
      target_pressure_present = True
      target_pressure_ = target_pressure

    mix_, particle_radius_ = self._prepare_atmosphere_composition(
      nprofile, mix, particle_radius, default_mix, default_particle_radius
    )
    nq = mix_.shape[0]
    nparticles = particle_radius_.shape[0]

    ea_pxd.evoatmosphere_initialize_atmosphere_p_wrapper(
      self._ptr, &nprofile, <double *>pressure_.data,
      <double *>temperature_.data, <double *>edd_.data,
      &nq, <double *>mix_.data, &nparticles,
      <double *>particle_radius_.data, &persistent_,
      &tropopause_pressure_, &tropopause_pressure_present,
      &maintain_toa_pressure_, &maintain_toa_pressure_present,
      &target_pressure_, &target_pressure_present, err
    )
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    
  property dat:
    """The PhotochemData class. Data in this class almost never changes after the
    class is initialized.
    """
    def __get__(self):
      dat = PhotochemData()
      ea_pxd.evoatmosphere_dat_get(self._ptr, &dat._ptr)
      return dat
      
  property var:
    """The PhotochemVars class. Data in this class can change between photochemical 
    integrations.
    """
    def __get__(self):
      var = PhotochemVars()
      ea_pxd.evoatmosphere_var_get(self._ptr, &var._ptr)
      return var
      
  property wrk:
    """The PhotochemWrk class. Data in this class changes during each step of 
    integration.
    """
    def __get__(self):
      wrk = PhotochemWrkEvo()
      ea_pxd.evoatmosphere_wrk_get(self._ptr, &wrk._ptr)
      return wrk

  def prep_atmosphere(self, ndarray[double, ndim=2] usol):
    """Given `usol`, the densities of each species in the atmosphere,
    this subroutine calculates reaction rates, photolysis rates, etc.
    and puts this information into self.wrk. self.wrk contains all the
    information needed for `dochem` to compute chemistry.

    Parameters
    ----------
    usol : ndarray[double,ndim=2]
        Number densities (molecules/cm^3)
    """
    cdef char err[ERR_LEN+1]
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef ndarray usol_ = np.asfortranarray(usol)
    if usol_.shape[0] != nq or usol_.shape[1] != nz:
      raise PhotoException("Input usol is the wrong size.")
      
    ea_pxd.evoatmosphere_prep_atmosphere_wrapper(self._ptr, &nq, &nz, <double *>usol_.data, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def out2atmosphere_txt(self,str filename, int number_of_decimals=5, bool overwrite = False, bool clip = True):
    """Saves state of the atmosphere using the concentrations in self.wrk.usol.

    Parameters
    ----------
    filename : str
        Output filename
    number_of_decimals : int, optional
        Number of decimals
    overwrite : bool, optional
        If true, then output file can be overwritten, by default False
    clip : bool, optional
        If true, then mixing ratios are clipped at a very small 
        positive number, by default False
    """   
    cdef bytes filename_b = pystring2cstring(filename)
    cdef char *filename_c = filename_b
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_out2atmosphere_txt_wrapper(self._ptr, filename_c, &number_of_decimals, &overwrite, &clip, err)  
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def gas_fluxes(self):
    """Computes gas fluxes at model boundaries in order to maintain
    current atmospheric concentrations. Uses the densities stored in
    self.wrk.usol.

    Returns
    -------
    tuple
        First element are the surface fluxes, and the second are top-of-atmosphere
        fluxes. Units molecules/cm^2/s.
    """
    cdef ndarray surf_fluxes = np.empty(self.dat.nq, np.double)
    cdef ndarray top_fluxes = np.empty(self.dat.nq, np.double)
    cdef int nq = self.dat.nq
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_gas_fluxes_wrapper(self._ptr, &nq, <double *>surf_fluxes.data, <double *>top_fluxes.data, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    surface = {}
    names = self.dat.species_names
    for i in range(self.dat.nq):
      surface[names[i]] = surf_fluxes[i]
    top = {}
    names = self.dat.species_names
    for i in range(self.dat.nq):
      top[names[i]] = top_fluxes[i]
    return surface, top

  def set_lower_bc(self, str species, str bc_type, vdep = None, den = None, press = None,
                            flux = None, height = None):
    """Sets a lower boundary condition.

    Parameters
    ----------
    species : str
        Species to set boundary condition
    bc_type : str
        Boundary condition type
    vdep : float, optional
        Deposition velocity (cm/s)
    den : float, optional
        Density (molecules/cm^3)
    press : float, optional
        Pressure (dynes/cm^2)
    flux : float, optional
        Flux (molecules/cm^2/s)
    height : float, optional
        Height in atmosphere (km)
    """
    cdef bytes species_b = pystring2cstring(species)
    cdef char *species_c = species_b
    cdef bytes bc_type_b = pystring2cstring(bc_type)
    cdef char *bc_type_c = bc_type_b
    cdef double vdep_c = 0
    cdef double den_c = 0
    cdef double press_c = 0
    cdef double flux_c = 0
    cdef double height_c = 0
    cdef bool missing = False
    if bc_type == 'vdep':
      if vdep == None:
        missing = True
      else:
        vdep_c = vdep
    elif bc_type == 'den':
      if den == None:
        missing = True
      else:
        den_c = den
    elif bc_type == 'press':
      if press == None:
        missing = True
      else:
        press_c = press
    elif bc_type == 'flux':
      if flux == None:
        missing = True
      else:
        flux_c = flux
    elif bc_type == 'vdep + dist flux':
      if vdep == None or flux == None or height == None:
        missing = True
      else:
        vdep_c = vdep
        flux_c = flux
        height_c = height
    elif bc_type == 'Moses':
      pass
      
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_set_lower_bc_wrapper(self._ptr, species_c, bc_type_c, 
                                      &vdep_c, &den_c, &press_c, &flux_c, &height_c, &missing, err);
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
  
  def set_upper_bc(self, str species, str bc_type, veff = None,flux = None):
    """Sets upper boundary condition.

    Parameters
    ----------
    species : str
        Species to set boundary condition
    bc_type : str
        Boundary condition type
    veff : float, optional
        effusion velocity (cm/s)
    flux : float, optional
        Flux (molecules/cm^2/s)
    """
    cdef bytes species_b = pystring2cstring(species)
    cdef char *species_c = species_b
    cdef bytes bc_type_b = pystring2cstring(bc_type)
    cdef char *bc_type_c = bc_type_b
    cdef double veff_c = 0
    cdef double flux_c = 0
    cdef bool missing = False
    if bc_type == 'veff':
      if veff == None:
        missing = True
      else:
        veff_c = veff
    elif bc_type == 'flux':
      if flux == None:
        missing = True
      else:
        flux_c = flux
      
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_set_upper_bc_wrapper(self._ptr, species_c, bc_type_c, 
                                      &veff_c, &flux_c, &missing, err);
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def set_rate_fcn(self, str species, object fcn):
    """Sets a function describing a custom rate for a species.
    This could be useful for modeling external processes not in the
    model.

    Parameters
    ----------
    species : str
        Species name
    fcn : function
        A Numba cfunc that describes the time-dependent rate
    """
    cdef bytes species_b = pystring2cstring(species)
    cdef char *species_c = species_b
    cdef char err[ERR_LEN+1]
    cdef uintptr_t fcn_l
    cdef ea_pxd.time_dependent_rate_fcn fcn_c

    if fcn is None:
      fcn_l = 0
      fcn_c = NULL
    else:
      argtypes = (ct.c_double, ct.c_int32, ct.POINTER(ct.c_double))
      restype = None
      if not fcn.ctypes.argtypes == argtypes:
        raise PhotoException("The callback function has the wrong argument types.")
      if not fcn.ctypes.restype == restype:
        raise PhotoException("The callback function has the wrong return type.")

      fcn_l = fcn.address
      fcn_c = <ea_pxd.time_dependent_rate_fcn> fcn_l
      
    ea_pxd.evoatmosphere_set_rate_fcn_wrapper(self._ptr, species_c, fcn_c, err)
    if len(err.strip()) > 0:
       raise PhotoException(err.decode("utf-8").strip())

  def set_temperature(self, ndarray[double, ndim=1] temperature, trop_alt = None):
    """Changes the temperature profile.

    This method cannot be used while a CVODE stepper is initialized or while
    a persistent pressure-temperature-eddy profile is enabled. Call
    :meth:`destroy_stepper` or :meth:`clear_press_temp_edd_profile` first,
    respectively.

    Parameters
    ----------
    temperature : ndarray[double,ndim=1]
        new temperature at each atomspheric layer
    trop_alt : float, optional
        Tropopause altitude (cm). Only necessary if gas rainout is enabled.
    """
    cdef char err[ERR_LEN+1]
    cdef int nz = temperature.size
    
    cdef double trop_alt_ = 0.0
    cdef bool trop_alt_present = False
    if trop_alt != None:
      trop_alt_present = True
      trop_alt_ = trop_alt
      
    ea_pxd.evoatmosphere_set_temperature_wrapper(self._ptr, &nz, <double *>temperature.data, 
                                       &trop_alt_, &trop_alt_present, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
  
  def set_press_temp_edd(self, ndarray[double, ndim=1] P, ndarray[double, ndim=1] T, ndarray[double, ndim=1] edd, trop_p = None, hydro_pressure = None):
    """Map pressure-temperature and pressure-eddy-diffusion profiles onto
    the model's current altitude grid.

    Temperature is interpolated linearly in log10 pressure. Eddy diffusion
    is interpolated linearly in log10 pressure-log10 eddy space. If the input
    profile does not reach the model surface, its deepest two points are used
    to extrapolate it to the surface. Values above the input profile are held
    constant.

    By default, the mapping uses hydrostatic pressure and solves for it from
    the bottom of the atmosphere upward. The alternative uses the actual gas
    pressure, ``density * k_boltz * T``.

    On success, this method updates the model temperature, eddy diffusion,
    and derived atmospheric working state. It does not evolve atmospheric
    species number densities. This method cannot be used while a CVODE
    stepper is initialized or while a persistent pressure-temperature-eddy
    profile is enabled; call :meth:`destroy_stepper` or
    :meth:`clear_press_temp_edd_profile` first, respectively.
    
    Parameters
    ----------
    P : ndarray of float64, shape (n,)
        Strictly decreasing pressure profile in dynes/cm^2. Must contain at
        least two finite, positive values. ``T`` and ``edd`` must have the
        same shape.
    T : ndarray of float64, shape (n,)
        Temperature corresponding to ``P``, in K. Must be finite and
        positive.
    edd : ndarray of float64, shape (n,)
        Eddy diffusion corresponding to ``P``, in cm^2/s. Must be finite and
        positive.
    trop_p : float, optional
        Tropopause pressure in dynes/cm^2. Required when gas rainout is
        enabled. When supplied, the mapped pressure must decrease strictly
        with altitude so the tropopause altitude is unambiguous.
    hydro_pressure : bool, default=True
        Use hydrostatic pressure if True. If False, use the actual gas
        pressure computed from number density and temperature.

    Raises
    ------
    PhotoException
        If the profiles are invalid or cannot be mapped onto the altitude
        grid.
    """
    cdef char err[ERR_LEN+1]
    cdef ndarray P_ = np.ascontiguousarray(P)
    cdef ndarray T_ = np.ascontiguousarray(T)
    cdef ndarray edd_ = np.ascontiguousarray(edd)
    cdef int P_dim1 = P_.size
    cdef int T_dim1 = T_.size
    cdef int edd_dim1 = edd_.size
    
    cdef double trop_p_ = 0.0
    cdef bool trop_p_present = False
    if trop_p != None:
      trop_p_present = True
      trop_p_ = trop_p

    cdef bool hydro_pressure_ = True
    cdef bool hydro_pressure_present = False
    if hydro_pressure != None:
      hydro_pressure_present = True
      hydro_pressure_ = hydro_pressure
      
    ea_pxd.evoatmosphere_set_press_temp_edd_wrapper(
      self._ptr, &P_dim1, <double *>P_.data,
      &T_dim1, <double *>T_.data, &edd_dim1, <double *>edd_.data,
      &trop_p_, &trop_p_present, &hydro_pressure_,
      &hydro_pressure_present, err
    )
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def set_press_temp_edd_profile(self, ndarray[double, ndim=1] P, ndarray[double, ndim=1] T, ndarray[double, ndim=1] edd,
                                 trop_p = None, hydro_pressure = None,
                                 maintain_toa_pressure = None, target_pressure = None):
    """Prescribe persistent pressure-temperature and pressure-eddy profiles.

    The profiles are mapped onto the current altitude grid immediately and
    during every subsequent atmospheric preparation, including every ODE
    right-hand-side evaluation. Temperature and eddy diffusion therefore
    remain functions of atmospheric pressure as composition evolves.

    Temperature is interpolated linearly in log10 pressure. Eddy diffusion
    is interpolated linearly in log10 pressure-log10 eddy space. If the input
    profile does not reach the model surface, its deepest two points are used
    to extrapolate it to the surface. Values above the input profile are held
    constant. Vertical-grid updates preserve and remap the profile.

    This method cannot be called while a stepper is initialized. Call
    :meth:`destroy_stepper` first. While persistent mode is enabled,
    :meth:`set_temperature` and :meth:`set_press_temp_edd` cannot be used;
    call :meth:`clear_press_temp_edd_profile` first. Persistent profiles enable
    approximate TOA-pressure maintenance by default; use
    ``maintain_toa_pressure=False`` to disable it.

    Parameters
    ----------
    P : ndarray of float64, shape (n,)
        Strictly decreasing pressure profile in dynes/cm^2. Must contain at
        least two finite, positive values. ``T`` and ``edd`` must have the
        same shape.
    T : ndarray of float64, shape (n,)
        Temperature corresponding to ``P``, in K. Must be finite and
        positive.
    edd : ndarray of float64, shape (n,)
        Eddy diffusion corresponding to ``P``, in cm^2/s. Must be finite and
        positive.
    trop_p : float, optional
        Tropopause pressure in dynes/cm^2. Required when gas rainout is
        enabled. When supplied, the mapped pressure must decrease strictly
        with altitude so the tropopause altitude is unambiguous.
    hydro_pressure : bool, default=True
        Use hydrostatic pressure if True. If False, use actual gas pressure,
        ``density * k_boltz * T``.
    maintain_toa_pressure : bool, default=True
        Enable approximate TOA-pressure maintenance during robust stepping.
        This option is only meaningful while the persistent profile is active.
    target_pressure : float, optional, default=0.1
        Target TOA pressure for approximate maintenance, in dynes/cm^2. Must
        be finite and positive when supplied.
    """
    cdef char err[ERR_LEN+1]
    cdef ndarray P_ = np.ascontiguousarray(P)
    cdef ndarray T_ = np.ascontiguousarray(T)
    cdef ndarray edd_ = np.ascontiguousarray(edd)
    cdef int P_dim1 = P_.size
    cdef int T_dim1 = T_.size
    cdef int edd_dim1 = edd_.size

    cdef double trop_p_ = 0.0
    cdef bool trop_p_present = False
    if trop_p != None:
      trop_p_present = True
      trop_p_ = trop_p

    cdef bool hydro_pressure_ = True
    cdef bool hydro_pressure_present = False
    if hydro_pressure != None:
      hydro_pressure_present = True
      hydro_pressure_ = hydro_pressure

    cdef bool maintain_toa_pressure_ = True
    cdef bool maintain_toa_pressure_present = False
    if maintain_toa_pressure != None:
      maintain_toa_pressure_present = True
      maintain_toa_pressure_ = maintain_toa_pressure

    cdef double target_pressure_ = 0.1
    cdef bool target_pressure_present = False
    if target_pressure != None:
      target_pressure_present = True
      target_pressure_ = target_pressure

    ea_pxd.evoatmosphere_set_press_temp_edd_profile_wrapper(
      self._ptr, &P_dim1, <double *>P_.data,
      &T_dim1, <double *>T_.data, &edd_dim1, <double *>edd_.data,
      &trop_p_, &trop_p_present, &hydro_pressure_,
      &hydro_pressure_present, &maintain_toa_pressure_,
      &maintain_toa_pressure_present, &target_pressure_,
      &target_pressure_present, err
    )
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def clear_press_temp_edd_profile(self):
    """Disable the persistent pressure-temperature-eddy profile.

    The most recently mapped altitude-based temperature and eddy-diffusion
    profiles remain in place. This also disables approximate TOA-pressure
    maintenance. This method cannot be called while a stepper is initialized;
    call :meth:`destroy_stepper` first.
    """
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_clear_press_temp_edd_profile_wrapper(self._ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    
  def update_vertical_grid(self, TOA_alt = None, TOA_pressure = None):
    """Re-does the vertical grid so that the pressure at the top of the
    atmosphere is at `TOA_alt` or `TOA_pressure`. If the TOA needs to be raised above the current
    TOA, then the function constantly extrapolates mixing ratios, temperature,
    eddy diffusion, and particle radii.

    Parameters
    ----------
    TOA_alt : float
        New top of atmosphere altitude (cm)
    TOA_pressure : float
        New top of atmosphere pressure (dynes/cm^2)
    """
    cdef char err[ERR_LEN+1]

    cdef double TOA_alt_ = 0.0
    cdef bool TOA_alt_present = False
    if TOA_alt != None:
      TOA_alt_present = True
      TOA_alt_ = TOA_alt

    cdef double TOA_pressure_ = 0.0
    cdef bool TOA_pressure_present = False
    if TOA_pressure != None:
      TOA_pressure_present = True
      TOA_pressure_ = TOA_pressure

    ea_pxd.evoatmosphere_update_vertical_grid_wrapper(self._ptr, &TOA_alt_, &TOA_alt_present, 
                                                  &TOA_pressure_, &TOA_pressure_present, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def evolve(self, str filename, double tstart, ndarray[double, ndim=2] usol, ndarray[double, ndim=1] t_eval, bool overwrite = False, bool restart_from_file = False):
    """Evolve atmosphere through time on the current fixed vertical grid,
    and save output in a binary Fortran file. Grid changes must be requested
    explicitly with :meth:`update_vertical_grid` or handled by the robust
    stepper's TOA-maintenance mode.

    Parameters
    ----------
    filename : str
        Filename to save results.
    tstart : float
        start time in seconds
    usol : ndarray[double,ndim=2]
        Initial number densities (molecules/cm^3)
    t_eval : ndarray[double,ndim=1]
        times to evaluate the solution
    overwrite : bool
        If true, then overwrites pre-existing files with `filename`
    restart_from_file : bool
        If true, then the integration restarts from the input file.

    Returns
    -------
    bool
        If True, then integration was successful.
    """
    cdef bytes filename_b = pystring2cstring(filename)
    cdef char *filename_c = filename_b
    cdef char err[ERR_LEN+1]
    cdef bool success
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef int nt = t_eval.size
    cdef ndarray usol_ = np.asfortranarray(usol)
    if usol_.shape[0] != nq or usol_.shape[1] != nz:
      raise PhotoException("Input usol is the wrong size.")
      
    ea_pxd.evoatmosphere_evolve_wrapper(self._ptr, filename_c, &tstart, &nq, &nz, <double *>usol_.data, &nt, <double *>t_eval.data, &overwrite, &restart_from_file, &success, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return success

  def check_for_convergence(self):
    "Determines if integration has converged to photochemical steady-state."  
    cdef bool converged
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_check_for_convergence_wrapper(self._ptr, &converged, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return converged

  def initialize_stepper(self, ndarray[double, ndim=2] usol_start):
    """Initializes an integration starting at `usol_start`.

    Parameters
    ----------
    usol_start : ndarray[double,ndim=2]
        Initial number densities (molecules/cm^3)
    """   
    cdef char err[ERR_LEN+1]
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef ndarray usol_start_ = np.asfortranarray(usol_start)
    if usol_start_.shape[0] != nq or usol_start_.shape[1] != nz:
      raise PhotoException("Input usol_start is the wrong size.")
    ea_pxd.evoatmosphere_initialize_stepper_wrapper(self._ptr, &nq, &nz,  <double *>usol_start_.data, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    
  def step(self):
    """Takes one internal integration step. Function `initialize_stepper`.
    must have been called before this

    Returns
    -------
    float
        Current time in the integration.
    """
    cdef char err[ERR_LEN+1]
    cdef double tn = ea_pxd.evoatmosphere_step_wrapper(self._ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return tn
    
  def destroy_stepper(self):
    "Deallocates memory created during `initialize_stepper`"
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_destroy_stepper_wrapper(self._ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def initialize_robust_stepper(self, ndarray[double, ndim=2] usol_start):
    """Initializes a robust integration starting at `usol_start`.

    The integration begins at time zero and resets the total accepted-step and
    failed-step counters. Optional approximate TOA-pressure maintenance is
    configured through ``self.var.toa_pressure_maintenance`` and requires a
    persistent pressure-based temperature/eddy-diffusion profile. When it is
    enabled, the initial composition is prepared and the model top is brought
    inside the configured pressure band before CVODE starts. This preflight is
    performed here so pressure-based initialization can retain its requested
    domain endpoints.

    Parameters
    ----------
    usol_start : ndarray[double,ndim=2]
        Initial number densities (molecules/cm^3)
    """   
    cdef char err[ERR_LEN+1]
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef ndarray usol_start_ = np.asfortranarray(usol_start)
    if usol_start_.shape[0] != nq or usol_start_.shape[1] != nz:
      raise PhotoException("Input usol_start is the wrong size.")
    ea_pxd.evoatmosphere_initialize_robust_stepper_wrapper(self._ptr, &nq, &nz,  <double *>usol_start_.data, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    
  def robust_step(self):
    """Takes one internal robust integration step. Function `initialize_robust_stepper`
    must have been called before this.

    A failed step is recovered from the last committed state without advancing
    logical time. Scheduled CVODE restarts preserve logical time and total
    counters, but discard segment-local convergence history. When
    ``self.var.toa_pressure_maintenance.enabled`` is true, accepted steps may
    trigger a pressure-targeted vertical-grid update after chemistry has
    converged; successful updates restart CVODE while preserving the total
    integration counters.

    Returns
    -------
    give_up : bool
         If True, then the algorithm thinks it is time to give up.
    converged : bool
        If True, then the integration has converged to a steady state.
    """
    cdef bool give_up, converged
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_robust_step_wrapper(self._ptr, &give_up, &converged, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return give_up, converged

  def find_steady_state(self):
    """Integrates using a robust stepper until a steady state has been achieved.

    Returns
    -------
    converged : bool
        If True, then the integration has converged to a steady state.
    """

    cdef bool converged, give_up
    converged = False

    self.initialize_robust_stepper(self.wrk.usol)

    while True:
      give_up, converged = self.robust_step()

      if give_up:
        converged = False
        return converged
      
      if converged:
        return converged

      PyErr_CheckSignals()

  def production_and_loss(self, str species, ndarray[double, ndim=2] usol):
    """Computes the production and loss of input `species`.
    See ProductionLoss object in photochem_types.f90.

    Parameters
    ----------
    species : str
        name of species
    usol : ndarray[double,ndim=2]
        Number densities (molecules/cm^3)

    Returns
    -------
    object
        Type describing production and loss of species
    """
    cdef bytes species_b = pystring2cstring(species)
    cdef char *species_c = species_b
    cdef char err[ERR_LEN+1]
    
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef ndarray usol_ = np.asfortranarray(usol)
    if usol_.shape[0] != nq or usol_.shape[1] != nz:
      raise PhotoException("Input usol is the wrong size.")
      
    cdef pl_pxd.ProductionLoss *pl_ptr
    ea_pxd.evoatmosphere_production_and_loss_wrapper(self._ptr, species_c, &nq, &nz, <double *>usol_.data, &pl_ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    pl = ProductionLoss()
    pl._ptr = pl_ptr
    return pl

  def mole_fraction_dict(self):
    """Makes a dictionary describing the atmospheric composition and structure
    using the densities in `self.wrk.usol`

    Returns
    -------
    dict
        Atmospheric composition and structure.
        - Key "alt" is altitude in cm
        - Key "temp" is temperature in K
        - Key "pressure" is pressure in dynes/cm^2
        - Key "density" is number density in molecules/cm^3
        - There is a key for each species giving its volume mixing ratio
    """   
    out = {}
    out['alt'] = self.var.z
    out['temp'] = self.var.temperature
    out['pressure'] = self.wrk.pressure
    out['density'] = self.wrk.density
    names = self.dat.species_names[:-2] # all but hv and M
    cdef ndarray usol = self.wrk.usol
    cdef ndarray densities = self.wrk.densities
    cdef ndarray density = self.wrk.density
    cdef int nq = self.dat.nq
    cdef int i
    # Evolving species
    for i in range(nq):
      out[names[i]] = usol[i,:]/density
    # Short lived species
    for i in range(self.dat.nsl):
      out[names[nq+i]] = densities[nq+i,:]/density
    return out
