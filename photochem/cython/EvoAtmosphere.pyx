
cimport EvoAtmosphere_pxd as ea_pxd

cdef class EvoAtmosphere:
  """A photochemical model which assumes no background gas. Once initialized,
  this class can integrate an atmosphere forward in time, and can investigate
  the reactions producing and destroying each molecule. The model can also
  optionally self-consistently evolve climate.
  """

  cdef void *_ptr
  cdef void *_dat_ptr
  cdef void *_var_ptr
  cdef void *_wrk_ptr

  def __init__(self, mechanism_file = None, settings_file = None, 
                     flux_file = None, atmosphere_txt = None, data_dir = None):           
    # Allocate memory
    ea_pxd.allocate_evoatmosphere(&self._ptr)

    if data_dir == None:
      data_dir_ = os.path.dirname(os.path.realpath(__file__))+'/data'
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
    ea_pxd.evoatmosphere_create_wrapper(&self._ptr, mechanism_file_c,
                                       settings_file_c, flux_file_c,
                                       atmosphere_txt_c, data_dir_c, 
                                       &self._dat_ptr, &self._var_ptr, &self._wrk_ptr,
                                       err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def __dealloc__(self):
    ea_pxd.deallocate_evoatmosphere(&self._ptr);
    
  property dat:
    """The PhotochemData class. Data in this class almost never changes after the
    `Atmosphere` class is initialized.
    """
    def __get__(self):
      dat = PhotochemData(alloc = False)
      dat._ptr = self._dat_ptr
      return dat
      
  property var:
    """The PhotochemVars class. Data in this class can change between photochemical 
    integrations.
    """
    def __get__(self):
      var = PhotochemVars(alloc = False)
      var._ptr = self._var_ptr
      return var
      
  property wrk:
    """The PhotochemWrk class. Data in this class changes during each step of 
    integration.
    """
    def __get__(self):
      wrk = PhotochemWrk(alloc = False)
      wrk._ptr = self._wrk_ptr
      return wrk

  def out2atmosphere_txt(self,filename = None, bool overwrite = False, bool clip = True):
    """Saves state of the atmosphere using the mixing ratios in self.wrk.usol.

    Parameters
    ----------
    filename : str
        Output filename
    overwrite : bool, optional
        If true, then output file can be overwritten, by default False
    clip : bool, optional
        If true, then mixing ratios are clipped at a very small 
        positive number, by default False
    """   
    cdef bytes filename_b = pystring2cstring(filename)
    cdef char *filename_c = filename_b
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_out2atmosphere_txt_wrapper(&self._ptr, filename_c, &overwrite, &clip, err)  
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
        fluxes.
    """
    cdef ndarray surf_fluxes = np.empty(self.dat.nq, np.double)
    cdef ndarray top_fluxes = np.empty(self.dat.nq, np.double)
    cdef int nq = self.dat.nq
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_gas_fluxes_wrapper(&self._ptr, &nq, <double *>surf_fluxes.data, <double *>top_fluxes.data, err)
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
      
    ea_pxd.evoatmosphere_regrid_prep_atmosphere_wrapper(&self._ptr, &nq, &nz, <double *>usol_.data, &top_atmos, err)
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
      
    ea_pxd.evoatmosphere_evolve_wrapper(&self._ptr, filename_c, &tstart, &nq, &nz, <double *>usol_.data, &nt, <double *>t_eval.data, &overwrite, &restart_from_file, &success, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return success

  def initialize_stepper(self, ndarray[double, ndim=2] usol_start):
    """Initializes an integration starting at `usol_start`.

    Parameters
    ----------
    usol_start : ndarray[double,ndim=2]
        Initial mixing ratios
    """   
    cdef char err[ERR_LEN+1]
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef ndarray usol_start_ = np.asfortranarray(usol_start)
    if usol_start_.shape[0] != nq or usol_start_.shape[1] != nz:
      raise PhotoException("Input usol_start is the wrong size.")
    ea_pxd.evoatmosphere_initialize_stepper_wrapper(&self._ptr, &nq, &nz,  <double *>usol_start_.data, err)
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
    cdef double tn = ea_pxd.evoatmosphere_step_wrapper(&self._ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return tn
    
  def destroy_stepper(self):
    "Deallocates memory created during `initialize_stepper`"
    cdef char err[ERR_LEN+1]
    ea_pxd.evoatmosphere_destroy_stepper_wrapper(&self._ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def production_and_loss(self, str species, ndarray[double, ndim=2] usol, double top_atmos):
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
      
    cdef void *pl_ptr
    ea_pxd.evoatmosphere_production_and_loss_wrapper(&self._ptr, species_c, &nq, &nz, <double *>usol_.data, &top_atmos, &pl_ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    pl = ProductionLoss()
    pl._ptr = pl_ptr
    return pl

  property T_surf:
    "double. The surface temperature (K)"
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_t_surf_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_t_surf_set(&self._ptr, &val)

  property T_trop:
    "double. The tropopause temperature."
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_t_trop_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_t_trop_set(&self._ptr, &val)

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
      
      ea_pxd.evoatmosphere_albedo_fcn_set(&self._ptr, fcn_c)

  property P_top_min:
    """double. When running the `evolve` routine, this determines
    the minimum pressure of the top of the model domain (bars). 
    If the pressure gets smaller than this value, then the integration will stop
    and re-grid the model domain before continuing integration, so that the 
    top of the atmosphere has a bigger pressure than `P_top_min`.
    """
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_p_top_min_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_p_top_min_set(&self._ptr, &val)

  property P_top_max:
    """double. When running the `evolve` routine, this determines
    the maximum pressure of the top of the model domain (bars). 
    If the pressure gets larger than this value, then the integration will stop
    and re-grid the model domain before continuing integration, so that the 
    top of the atmosphere has a smaller pressure than `P_top_max`.
    """
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_p_top_max_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_p_top_max_set(&self._ptr, &val)

  property top_atmos_adjust_frac:
    """Sets the fractional amount that the top of the model domain changes
    when integration is haulted by `P_top_min` or `P_top_max`
    """
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_top_atmos_adjust_frac_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_top_atmos_adjust_frac_set(&self._ptr, &val)
    
