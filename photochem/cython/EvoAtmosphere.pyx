
cimport EvoAtmosphere_pxd as ea_pxd

cdef class EvoAtmosphere:
  """A photochemical model which assumes no background gas. Once initialized,
  this class can integrate an atmosphere forward in time to a steady
  state. The model can also, optionally, self-consistently evolve climate for
  terrestrial planets.
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
               str flux_file, str atmosphere_txt, data_dir = None):           
    """Initializes the photochemical model.

    Parameters
    ----------
    mechanism_file : str
        Path to the reaction mechanism file (yaml format).
    settings_file : str
        Path to the settings file (yaml format).
    flux_file : str
        Path to the file describing the stellar flux.
    atmosphere_txt : str
        Path to the file containing altitude, total number density, temperature, 
        eddy diffusion, initial concentrations of each gas (mixing ratios), 
        and particle radii.
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
    cdef bytes atmosphere_txt_b = pystring2cstring(atmosphere_txt)
    cdef char *atmosphere_txt_c = atmosphere_txt_b
    cdef bytes data_dir_b = pystring2cstring(data_dir_)
    cdef char *data_dir_c = data_dir_b
    cdef char err[ERR_LEN+1]
    
    # Initialize
    ea_pxd.evoatmosphere_create_wrapper(self._ptr, mechanism_file_c,
                                       settings_file_c, flux_file_c,
                                       atmosphere_txt_c, data_dir_c, err)
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

    Parameters
    ----------
    temperature : ndarray[double,ndim=1]
        new temperature at each atomspheric layer
    trop_alt : float, optional
        Tropopause altitude (cm). Only necessary if rainout == True, 
        or fix_water_in_trop == True.
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
    """Given an input P, T, and edd, the code will find the temperature and eddy diffusion profile
    on the current altitude-grid that matches the inputs.
    
    Parameters
    ----------
    P : ndarray[double,ndim=1]
        Pressure (dynes/cm^2)
    T : ndarray[double,ndim=1]
        Temperature (K)
    edd : ndarray[double,ndim=1]
        Eddy diffusion (cm^2/s)
    trop_p : float, optional
        Tropopause pressure (dynes/cm^2). Only necessary if rainout == True, 
        or fix_water_in_trop == True.
    hydro_pressure : bool, optional
        If True, then use hydrostatic pressure. If False then use the
        actual pressure in the atmosphere. Default is True.
    """
    cdef char err[ERR_LEN+1]
    cdef int P_dim1 = P.size
    cdef int T_dim1 = T.size
    cdef int edd_dim1 = edd.size
    
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
      
    ea_pxd.evoatmosphere_set_press_temp_edd_wrapper(self._ptr, &P_dim1, <double *>P.data, &T_dim1, <double *>T.data, &edd_dim1, <double *>edd.data,
                                               &trop_p_, &trop_p_present, 
                                               &hydro_pressure_, &hydro_pressure_present, err)
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

  def regrid_prep_atmosphere(self, ndarray[double, ndim=2] usol, double top_atmos):
    """This subroutine calculates re-grids the model so that the top of the model domain 
    is at `top_atmos` the computes reaction rates, photolysis rates, etc.
    and puts this information into self.wrk. self.wrk contains all the
    information needed for `dochem` to compute chemistry.

    Parameters
    ----------
    usol : ndarray[double,ndim=2]
        The number densities (molecules/cm^3)
    top_atmos : float
        The top of the model domain (cm)
    """
    cdef char err[ERR_LEN+1]
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef ndarray usol_ = np.asfortranarray(usol)
    if usol_.shape[0] != nq or usol_.shape[1] != nz:
      raise PhotoException("Input usol is the wrong size.")
      
    ea_pxd.evoatmosphere_regrid_prep_atmosphere_wrapper(self._ptr, &nq, &nz, <double *>usol_.data, &top_atmos, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    
  def evolve(self, str filename, double tstart, ndarray[double, ndim=2] usol, ndarray[double, ndim=1] t_eval, bool overwrite = False, bool restart_from_file = False):
    """Evolve atmosphere through time, and saves output in a 
    binary Fortran file.

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

  property T_surf:
    """double. The surface temperature (K). Only relevent when doing time-dependent
    photochemical-climate simulation.
    """
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_t_surf_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_t_surf_set(self._ptr, &val)

  property T_trop:
    "double. Assumed tropopause temperature for climate calculations (K)."
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_t_trop_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_t_trop_set(self._ptr, &val)

  property albedo_fcn:
    """A function describing a temperature dependent surface albedo.
    This is useful for modeling the ice-albedo feedback.
    """
    def __set__(self, object fcn):
      cdef uintptr_t fcn_l
      cdef ea_pxd.temp_dependent_albedo_fcn fcn_c

      if fcn is None:
        fcn_l = 0
        fcn_c = NULL
      else:
        argtypes = (ct.c_double,)
        restype = ct.c_double
        if not fcn.ctypes.argtypes == argtypes:
          raise PhotoException("The callback function has the wrong argument types.")
        if not fcn.ctypes.restype == restype:
          raise PhotoException("The callback function has the wrong return type.")

        fcn_l = fcn.address
        fcn_c = <ea_pxd.temp_dependent_albedo_fcn> fcn_l
      
      ea_pxd.evoatmosphere_albedo_fcn_set(self._ptr, fcn_c)

  property P_top_min:
    """double. When running the `evolve` routine, this determines
    the minimum pressure of the top of the model domain (bars). 
    If the pressure gets smaller than this value, then the integration will stop
    and re-grid the model domain before continuing integration, so that the 
    top of the atmosphere has a bigger pressure than `P_top_min`.
    """
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_p_top_min_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_p_top_min_set(self._ptr, &val)

  property P_top_max:
    """double. When running the `evolve` routine, this determines
    the maximum pressure of the top of the model domain (bars). 
    If the pressure gets larger than this value, then the integration will stop
    and re-grid the model domain before continuing integration, so that the 
    top of the atmosphere has a smaller pressure than `P_top_max`.
    """
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_p_top_max_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_p_top_max_set(self._ptr, &val)

  property top_atmos_adjust_frac:
    """Sets the fractional amount that the top of the model domain changes
    when integration is haulted by `P_top_min` or `P_top_max`
    """
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_top_atmos_adjust_frac_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_top_atmos_adjust_frac_set(self._ptr, &val)
    
