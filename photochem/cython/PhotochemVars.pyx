cimport PhotochemVars_pxd as var_pxd

cdef class PhotochemVars:
  cdef void *_ptr
  cdef bint _destroy

  def __cinit__(self, bint alloc = True):
    if alloc:
      var_pxd.allocate_photochemvars(&self._ptr)
      self._destroy = True
    else:
      self._destroy = False

  def __dealloc__(self):
    if self._destroy:
      var_pxd.deallocate_photochemvars(&self._ptr)
      self._ptr = NULL
  
  property nz:
    def __get__(self):
      cdef int nz
      var_pxd.photochemvars_nz_get(&self._ptr, &nz)
      return nz

  property top_atmos:
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_top_atmos_get(&self._ptr, &val)
      return val

  property bottom_atmos:
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_bottom_atmos_get(&self._ptr, &val)
      return val
      
  property usol_init:
    def __get__(self):
      cdef int dim1, dim2
      var_pxd.photochemvars_usol_init_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      var_pxd.photochemvars_usol_init_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property z:
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_z_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_z_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property temperature:
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_temperature_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_temperature_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property edd:
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_edd_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_edd_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property photon_flux:
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_photon_flux_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_photon_flux_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property grav:
    def __get__(self):
      cdef int dim1
      var_pxd.photochemvars_grav_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      var_pxd.photochemvars_grav_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property at_photo_equilibrium:
    def __get__(self):
      cdef bool at_photo_equilibrium
      var_pxd.photochemvars_at_photo_equilibrium_get(&self._ptr, &at_photo_equilibrium)
      return at_photo_equilibrium
      
  property surface_pressure:
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_surface_pressure_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_surface_pressure_set(&self._ptr, &val)

  property max_error_reinit_attempts:
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_max_error_reinit_attempts_get(&self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_max_error_reinit_attempts_set(&self._ptr, &val)
  
  property rtol:
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_rtol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_rtol_set(&self._ptr, &val)
      
  property atol:
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_atol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_atol_set(&self._ptr, &val)

  property mxsteps:
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_mxsteps_get(&self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_mxsteps_set(&self._ptr, &val)
      
  property equilibrium_time:
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_equilibrium_time_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_equilibrium_time_set(&self._ptr, &val)
  
  property verbose:
    def __get__(self):
      cdef int val
      var_pxd.photochemvars_verbose_get(&self._ptr, &val)
      return val
    def __set__(self, int val):
      var_pxd.photochemvars_verbose_set(&self._ptr, &val)

  property fast_arbitrary_rate:
    def __get__(self):
      cdef double val
      var_pxd.photochemvars_fast_arbitrary_rate_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      var_pxd.photochemvars_fast_arbitrary_rate_set(&self._ptr, &val)

    
  