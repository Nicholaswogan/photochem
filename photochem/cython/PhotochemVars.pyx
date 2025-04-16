cimport PhotochemVars_pxd as var_pxd

cdef class PhotochemVars:
  """This class contains data that can change between independent
  model integrations.
  """

  cdef var_pxd.PhotochemVars *_ptr

  def __cinit__(self):
    self._ptr = NULL
  
  property nz:
    "The number of vertical atmospheric layers"
    def __get__(self):
      cdef int nz
      var_pxd.photochemvars_nz_get(self._ptr, &nz)
      return nz

  property top_atmos:
    "The top of the model domain (cm)"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_top_atmos_get(self._ptr, &val)
      return val

  property bottom_atmos:
    "The bottom of the model domain (cm)"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_bottom_atmos_get(self._ptr, &val)
      return val
      
  property usol_init:
    """ndarray[double,dim=2], shape (nq,nz). Contains the initial concentration
    of atmospheric species (molecules/cm^3).
    """
    def __get__(self):
      cdef int dim1, dim2
      var_pxd.photochemvars_usol_init_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      var_pxd.photochemvars_usol_init_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property particle_radius:
    "ndarray[double,dim=2], shape (npq,nz). cm"
    def __get__(self):
      cdef int dim1, dim2
      var_pxd.photochemvars_particle_radius_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      var_pxd.photochemvars_particle_radius_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=2] arr):
      cdef int dim1, dim2
      var_pxd.photochemvars_particle_radius_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr_ = np.asfortranarray(arr)
      if arr_.shape[0] != dim1 or arr_.shape[1] != dim2:
        raise PhotoException("Input array is the wrong size.")
      var_pxd.photochemvars_particle_radius_set(self._ptr, &dim1, &dim2, <double *>arr_.data)

  property diurnal_fac:
    "double. Default is 0.5, to account for half planet facing the sun."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_diurnal_fac_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_diurnal_fac_set(self._ptr, &val)

  property trop_alt:
    "double. Tropopause altitude."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_trop_alt_get(self._ptr, &val)
      return val

  property trop_ind:
    "int. Tropopause index."
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_trop_ind_get(self._ptr, &val)
      return val

  property relative_humidity:
    "double. Relative humidity of H2O in the troposphere."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_relative_humidity_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_relative_humidity_set(self._ptr, &val)

  property H2O_cond_params:
    "CondensationParameters. H2O condensation rate parameters."
    def __get__(self):
      cdef atom_pxd.CondensationParameters *ptr1
      var_pxd.photochemvars_h2o_cond_params_get(self._ptr, &ptr1)
      val = CondensationParameters()
      val._ptr = ptr1
      return val
  
  property z:
    "ndarray[double,dim=1], shape (nz). The altitude of the center of each atmopsheric layer (cm)"
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_z_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_z_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property photon_flux_fcn:
    "A function for altering the photon flux over time"
    def __set__(self, object fcn):
      cdef uintptr_t fcn_l
      cdef var_pxd.time_dependent_flux_fcn fcn_c
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
        fcn_c = <var_pxd.time_dependent_flux_fcn> fcn_l

      var_pxd.photochemvars_photon_flux_fcn_set(self._ptr, fcn_c)

  property cond_params:
    """list, shape (np). Parameters describing condensation and evaporation rates and
    the RH needed for condensation.
    """
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_cond_params_get_size(self._ptr, &dim1)
      cdef atom_pxd.CondensationParameters **arrp = <atom_pxd.CondensationParameters **> malloc(dim1 * sizeof(atom_pxd.CondensationParameters *))
      var_pxd.photochemvars_cond_params_get(self._ptr, &dim1, arrp)
      arr1 = []
      for i in range(dim1):
        tmp = CondensationParameters()
        tmp._ptr = arrp[i]
        arr1.append(tmp)
      free(arrp)
      return arr1

  property temperature:
    "ndarray[double,dim=1], shape (nz). The temperature of each atmospheric layer (K)"
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_temperature_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_temperature_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property edd:
    "ndarray[double,dim=1], shape (nz). The eddy diffusion of each atmospheric layer (cm^2/s)"
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_edd_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_edd_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] edd_new_):
      cdef int dim1
      var_pxd.photochemvars_edd_get_size(self._ptr, &dim1)
      cdef ndarray edd_new = np.asfortranarray(edd_new_)
      if edd_new.shape[0] != dim1:
        raise PhotoException("Input edd is the wrong size.")
      var_pxd.photochemvars_edd_set(self._ptr, &dim1, <double *>edd_new.data)

  property custom_binary_diffusion_fcn:
    "A function for specifying a custom binary diffusion parameter (b_ij)"
    def __set__(self, object fcn):
      cdef uintptr_t fcn_l
      cdef var_pxd.binary_diffusion_fcn fcn_c
      if fcn is None:
        fcn_l = 0
        fcn_c = NULL
      else:
        argtypes = (ct.c_double, ct.c_double, ct.c_double)
        restype = ct.c_double
        if not fcn.ctypes.argtypes == argtypes:
          raise PhotoException("The callback function has the wrong argument types.")
        if not fcn.ctypes.restype == restype:
          raise PhotoException("The callback function has the wrong return type.")
        fcn_l = fcn.address
        fcn_c = <var_pxd.binary_diffusion_fcn> fcn_l

      var_pxd.photochemvars_custom_binary_diffusion_fcn_set(self._ptr, fcn_c)
      
  property photon_flux:
    "ndarray[double,dim=1], shape (nw). photon/cm^2/s in each wavelength bin hitting planet."
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_photon_flux_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_photon_flux_get(self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property grav:
    "ndarray[double,dim=1], shape (nz). The gravitational acceleration at the center of each grid cell (cm/s^2)."
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_grav_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_grav_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property at_photo_equilibrium:
    "bool. If True, then the model is at photochemical equilibrium."
    def __get__(self):
      cdef bool at_photo_equilibrium
      var_pxd.photochemvars_at_photo_equilibrium_get(self._ptr, &at_photo_equilibrium)
      return at_photo_equilibrium
      
  property surface_pressure:
    "double. The surface pressure in bars."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_surface_pressure_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_surface_pressure_set(self._ptr, &val)

  property tauc:
    "ndarray[double,dim=2], shape (nz,nw). Custom optical depth in each layer."
    def __get__(self):
      cdef int dim1, dim2
      var_pxd.photochemvars_tauc_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      var_pxd.photochemvars_tauc_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=2] arr_):
      cdef int dim1, dim2
      var_pxd.photochemvars_tauc_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.asfortranarray(arr_)
      if arr.shape[0] != dim1 or arr.shape[1] != dim2:
        raise PhotoException("Input is the wrong size.")
      var_pxd.photochemvars_tauc_set(self._ptr, &dim1, &dim2, <double *>arr.data)  

  property w0c:
    "ndarray[double,dim=2], shape (nz,nw). Custom single scattering albedo."
    def __get__(self):
      cdef int dim1, dim2
      var_pxd.photochemvars_w0c_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      var_pxd.photochemvars_w0c_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=2] arr_):
      cdef int dim1, dim2
      var_pxd.photochemvars_w0c_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.asfortranarray(arr_)
      if arr.shape[0] != dim1 or arr.shape[1] != dim2:
        raise PhotoException("Input is the wrong size.")
      var_pxd.photochemvars_w0c_set(self._ptr, &dim1, &dim2, <double *>arr.data)  

  property g0c:
    "ndarray[double,dim=2], shape (nz,nw). Custom asymmetry parameter."
    def __get__(self):
      cdef int dim1, dim2
      var_pxd.photochemvars_g0c_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      var_pxd.photochemvars_g0c_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=2] arr_):
      cdef int dim1, dim2
      var_pxd.photochemvars_g0c_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.asfortranarray(arr_)
      if arr.shape[0] != dim1 or arr.shape[1] != dim2:
        raise PhotoException("Input is the wrong size.")
      var_pxd.photochemvars_g0c_set(self._ptr, &dim1, &dim2, <double *>arr.data)  

  property max_error_reinit_attempts:
    """int. number of times to reinitialize CVODE when it returns
    a potentially recoverable error.
    """
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_max_error_reinit_attempts_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_max_error_reinit_attempts_set(self._ptr, &val)
  
  property rtol:
    "double. Integration relative tolerance."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_rtol_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_rtol_set(self._ptr, &val)
      
  property atol:
    """double. Integration absolute tolerance. If autodiff == .true., then the model
    works better when atol is smaller (e.g., atol = ~1.0e-18).
    """
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_atol_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_atol_set(self._ptr, &val)

  property mxsteps:
    "int. Max number of steps before integrator will give up."
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_mxsteps_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_mxsteps_set(self._ptr, &val)
      
  property equilibrium_time:
    "double. Atomsphere considered in equilibrium if integrations reaches this time (seconds)"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_equilibrium_time_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_equilibrium_time_set(self._ptr, &val)

  property conv_hist_factor:
    """double. For convergence checking. Considers mixing ratio change between t_now and time 
    t = t_now*conv_hist_factor to see if atmosphere is changing.
    """
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_conv_hist_factor_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_conv_hist_factor_set(self._ptr, &val)

  property conv_min_mix:
    "double. Minimum mixing ratio considered in convergence checking."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_conv_min_mix_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_conv_min_mix_set(self._ptr, &val)

  property conv_longdy:
    """double. Threshold normalized change in mixing ratios for converchecking check.
    A reasonable value is ~1.0e-2.
    """
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_conv_longdy_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_conv_longdy_set(self._ptr, &val)

  property conv_longdydt:
    """double. Threshold normalized change in mixing ratios per time change for 
    convergence checking.
    """
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_conv_longdydt_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_conv_longdydt_set(self._ptr, &val)

  property max_dt:
    """double. Maximum time step size (seconds).
    """
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_max_dt_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_max_dt_set(self._ptr, &val)

  property autodiff:
    """bool. If True, then the chemistry terms of the Jacobian are computed uses 
    foward mode automatic differentiation.
    """
    def __get__(self):
      cdef bool val
      var_pxd.photochemvars_autodiff_get(self._ptr, &val)
      return val
    def __set__(self, bool val):
      var_pxd.photochemvars_autodiff_set(self._ptr, &val)

  property epsj:
    "double. Perturbation for finite difference Jacobian calculation, when autodiff == .false."
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_epsj_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_epsj_set(self._ptr, &val)
  
  property verbose:
    "int. 0 == no printing. 1 == some printing. 2 == bunch of printing."
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_verbose_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_verbose_set(self._ptr, &val)

  property fast_arbitrary_rate:
    "double. arbitrary rate that is fast (1/s). Used for keeping H2O at saturation in troposphere"
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_fast_arbitrary_rate_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_fast_arbitrary_rate_set(self._ptr, &val)

  property upwind_molec_diff:
    """bool. If True, then the code uses a 1st order upwind method for the advective molecular
    diffusion terms instead of a centered scheme. This permits stability (at the cost 
    of accuracy) for atmospheres with strong molcular advection in the upper atmosphere.
    """
    def __get__(self):
      cdef bool val
      var_pxd.photochemvars_upwind_molec_diff_get(self._ptr, &val)
      return val
    def __set__(self, bool val):
      var_pxd.photochemvars_upwind_molec_diff_set(self._ptr, &val)

  property nerrors_before_giveup:
    "int. Number of integration errors before giving up completely"
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_nerrors_before_giveup_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_nerrors_before_giveup_set(self._ptr, &val)

  property nsteps_before_conv_check:
    "int. Number of steps to take before checking for convergence"
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_nsteps_before_conv_check_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_nsteps_before_conv_check_set(self._ptr, &val)

  property nsteps_before_reinit:
    "int. Number of steps before reinitializing the integration"
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_nsteps_before_reinit_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_nsteps_before_reinit_set(self._ptr, &val)

  property nsteps_before_giveup:
    "int. Number of total steps to take before giving up."
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_nsteps_before_giveup_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_nsteps_before_giveup_set(self._ptr, &val)