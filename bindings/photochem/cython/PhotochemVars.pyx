cimport PhotochemVars_pxd as var_pxd

cdef class TOAPressureMaintenance:
  """Live configuration for optional robust-stepper TOA pressure maintenance.

  Instances are borrowed views owned by an [EvoAtmosphere][photochem.EvoAtmosphere] model.
  Assigning a property immediately changes that model's configuration.
  """

  cdef var_pxd.TOAPressureMaintenance *_ptr

  def __cinit__(self):
    self._ptr = NULL

  property enabled:
    """bool. Enable automatic TOA pressure maintenance."""
    def __get__(self):
      cdef bool val
      var_pxd.taopressuremaintenance_enabled_get(self._ptr, &val)
      return val
    def __set__(self, bool val):
      var_pxd.taopressuremaintenance_enabled_set(self._ptr, &val)

  property target_pressure:
    """float. Target TOA pressure in dyn/cm^2."""
    def __get__(self):
      cdef double val
      var_pxd.taopressuremaintenance_target_pressure_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.taopressuremaintenance_target_pressure_set(self._ptr, &val)

  property pressure_factor:
    """double. Acceptable multiplicative pressure band around the target."""
    def __get__(self):
      cdef double val
      var_pxd.taopressuremaintenance_pressure_factor_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.taopressuremaintenance_pressure_factor_set(self._ptr, &val)

  property nsteps_between_updates:
    """int. Minimum accepted robust steps between maintenance attempts."""
    def __get__(self):
      cdef int val
      var_pxd.taopressuremaintenance_nsteps_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.taopressuremaintenance_nsteps_set(self._ptr, &val)

  property max_failures:
    """int. Failed maintenance attempts allowed before integration stops."""
    def __get__(self):
      cdef int val
      var_pxd.taopressuremaintenance_max_failures_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.taopressuremaintenance_max_failures_set(self._ptr, &val)

cdef class PhotochemVars:
  """Prepared atmospheric state and mutable model configuration.

  This object is a live borrowed view available as [EvoAtmosphere.var][photochem.EvoAtmosphere.var];
  it should not be constructed directly. Array-valued getters return copies,
  while assignments to settable properties update the model. Atmospheric
  profile values become valid only after ``EvoAtmosphere.atmosphere_initialized``
  is true.
  """

  cdef var_pxd.PhotochemVars *_ptr

  def __cinit__(self):
    self._ptr = NULL
  
  property nz:
    """int. Number of vertical atmospheric layers."""
    def __get__(self):
      cdef int nz
      var_pxd.photochemvars_nz_get(self._ptr, &nz)
      return nz

  property top_atmos:
    """float. Altitude of the top of the model domain in cm."""
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_top_atmos_get(self._ptr, &val)
      return val

  property bottom_atmos:
    """float. Altitude of the bottom of the model domain in cm."""
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_bottom_atmos_get(self._ptr, &val)
      return val
      
  property particle_radius:
    """ndarray[float], shape (np, nz). Particle radii in cm."""
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
        raise PhotoException("particle_radius must have shape (np, nz).")
      var_pxd.photochemvars_particle_radius_set(self._ptr, &dim1, &dim2, <double *>arr_.data)

  property diurnal_fac:
    """float. Factor multiplying photon flux in photolysis-rate calculations.

    The default 0.5 accounts for half the planet facing the star.
    """
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_diurnal_fac_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_diurnal_fac_set(self._ptr, &val)

  property trop_alt:
    """float. Tropopause altitude in cm."""
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_trop_alt_get(self._ptr, &val)
      return val

  property trop_ind:
    """int. Internal tropopause index using the model's Fortran convention."""
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_trop_ind_get(self._ptr, &val)
      return val

  property rainfall_rate:
    """float. Rainfall relative to Earth's average of 1.1e17 molecules/cm^2/s."""
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_rainfall_rate_get(self._ptr, &val)
      return val

  property z:
    """ndarray[float], shape (nz,). Layer-center altitudes in cm."""
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_z_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_z_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property dz:
    """ndarray[float], shape (nz,). Atmospheric layer thicknesses in cm."""
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_dz_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_dz_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property photon_flux_fcn:
    """Numba callback that supplies stellar photon flux with time.

    Assign ``None`` to disable the callback. Otherwise supply a compiled C
    callback with signature ``void(double time, int32 nw, double* flux)``;
    the callback must fill ``flux`` with ``nw`` values in photons/cm^2/s.
    """
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
    """list[CondensationParameters], length ``np``. Condensation and evaporation parameters.

    Entries are live borrowed views; assigning their properties changes the
    model configuration.
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
    """ndarray[float], shape (nz,). Layer temperatures in K."""
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_temperature_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_temperature_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property edd:
    """ndarray[float], shape (nz,). Eddy diffusion coefficients in cm^2/s."""
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
        raise PhotoException("edd must have shape (nz,).")
      var_pxd.photochemvars_edd_set(self._ptr, &dim1, <double *>edd_new.data)

  property custom_binary_diffusion_fcn:
    """Numba callback for a custom binary diffusion parameter.

    Assign ``None`` to restore the built-in calculation. Otherwise supply a
    compiled C callback with signature ``double(double mu_i, double mubar,
    double T)``. The molar masses are in g/mol, temperature is in K, and the
    callback returns the binary diffusion parameter in molecules cm^-1 s^-1.
    """
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
    """ndarray[float], shape (nw,). Incident photons/cm^2/s per wavelength bin."""
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_photon_flux_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_photon_flux_get(self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property grav:
    """ndarray[float], shape (nz,). Layer-center gravity in cm/s^2."""
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_grav_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_grav_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property tauc:
    """ndarray[float], shape (nz, nw). Custom layer optical depth."""
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
        raise PhotoException("tauc must have shape (nz, nw).")
      var_pxd.photochemvars_tauc_set(self._ptr, &dim1, &dim2, <double *>arr.data)  

  property w0c:
    """ndarray[float], shape (nz, nw). Custom single-scattering albedo."""
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
        raise PhotoException("w0c must have shape (nz, nw).")
      var_pxd.photochemvars_w0c_set(self._ptr, &dim1, &dim2, <double *>arr.data)  

  property g0c:
    """ndarray[float], shape (nz, nw). Custom scattering asymmetry parameter."""
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
        raise PhotoException("g0c must have shape (nz, nw).")
      var_pxd.photochemvars_g0c_set(self._ptr, &dim1, &dim2, <double *>arr.data)  

  property max_error_reinit_attempts:
    """int. Number of times to reinitialize CVODE when it returns
    a potentially recoverable error during [EvoAtmosphere.evolve][photochem.EvoAtmosphere.evolve].
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
    """float. Dimensionless absolute tolerance in mixing ratio.

    The model multiplies this value by each layer's hydrostatic number density
    to construct the CVODE absolute-tolerance vector.
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
    """float. Time in s at which an integration is considered converged."""
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
    """double. Threshold normalized change in mixing ratios for convergence.
    A reasonable value is ~1.0e-2.
    """
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_conv_longdy_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_conv_longdy_set(self._ptr, &val)

  property conv_longdydt:
    """double. Threshold mixing-ratio change rate for convergence, in s^-1."""
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

  property initial_dt:
    """double. Initial CVODE time step (seconds)."""
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_initial_dt_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_initial_dt_set(self._ptr, &val)

  property max_err_test_failures:
    """int. Maximum CVODE error-test failures allowed per attempted step."""
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_max_err_test_failures_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_max_err_test_failures_set(self._ptr, &val)

  property max_order:
    """int. Maximum order used by CVODE's BDF method."""
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_max_order_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_max_order_set(self._ptr, &val)

  property jacobian_method:
    """int. Method used to compute the chemistry Jacobian.

    ``1`` selects automatic differentiation, ``2`` selects finite differences,
    and ``3`` selects the default analytical implementation.
    """
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_jacobian_method_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_jacobian_method_set(self._ptr, &val)

  property epsj:
    "double. Relative perturbation used by the finite-difference Jacobian method."
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

  property upwind_molec_diff:
    """bool. If True, then the code uses a 1st order upwind method for the advective molecular
    diffusion terms instead of a centered scheme. This permits stability (at the cost 
    of accuracy) for atmospheres with strong molecular advection in the upper atmosphere.
    """
    def __get__(self):
      cdef bool val
      var_pxd.photochemvars_upwind_molec_diff_get(self._ptr, &val)
      return val
    def __set__(self, bool val):
      var_pxd.photochemvars_upwind_molec_diff_set(self._ptr, &val)

  property nerrors_before_giveup:
    """int. Number of failed-step recovery restarts allowed. The next
    integration error makes ``robust_step`` give up.
    """
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_nerrors_before_giveup_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_nerrors_before_giveup_set(self._ptr, &val)

  property nsteps_before_conv_check:
    """int. Accepted steps after initialization or restart to take before
    checking the non-time convergence criteria.
    """
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_nsteps_before_conv_check_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_nsteps_before_conv_check_set(self._ptr, &val)

  property nsteps_before_reinit:
    """int. Accepted steps per integration segment. At this count CVODE is
    restarted and segment-local convergence history is discarded.
    """
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_nsteps_before_reinit_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_nsteps_before_reinit_set(self._ptr, &val)

  property nsteps_before_giveup:
    "int. Maximum total accepted steps before ``robust_step`` gives up."
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_nsteps_before_giveup_get(self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_nsteps_before_giveup_set(self._ptr, &val)

  property reinit_min_density:
    """double. Minimum number density retained when robust integration
    restarts after a recoverable solver failure (molecules/cm^3).
    """
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_reinit_min_density_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_reinit_min_density_set(self._ptr, &val)

  property toa_pressure_maintenance:
    """TOAPressureMaintenance. Live automatic TOA-pressure settings.

    Automatic maintenance requires a persistent pressure-based temperature
    and eddy-diffusion profile.
    """
    def __get__(self):
      cdef TOAPressureMaintenance maintenance = TOAPressureMaintenance()
      var_pxd.photochemvars_toa_pressure_maintenance_get(self._ptr, &maintenance._ptr)
      return maintenance
