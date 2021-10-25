
cdef extern void allocate_photochemvars(void *ptr)
cdef extern void deallocate_photochemvars(void *ptr)

cdef extern void photochemvars_nz_get(void *ptr, int *nz)

cdef extern void photochemvars_at_photo_equilibrium_get(void *ptr, bint *at_photo_equilibrium)

cdef extern void photochemvars_usol_init_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemvars_usol_init_get(void *ptr, int *dim1, int *dim2, double *usol_init)

cdef extern void photochemvars_temperature_get_size(void *ptr, int *dim1)
cdef extern void photochemvars_temperature_get(void *ptr, int *dim1, double *temperature)

cdef extern void photochemvars_grav_get_size(void *ptr, int *dim1)
cdef extern void photochemvars_grav_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemvars_z_get_size(void *ptr, int *dim1)
cdef extern void photochemvars_z_get(void *ptr, int *dim1, double *z)

cdef extern void photochemvars_rtol_get(void *ptr, double *val)
cdef extern void photochemvars_rtol_set(void *ptr, double *val)

cdef extern void photochemvars_atol_get(void *ptr, double *val)
cdef extern void photochemvars_atol_set(void *ptr, double *val)

cdef extern void photochemvars_mxsteps_get(void *ptr, int *val)
cdef extern void photochemvars_mxsteps_set(void *ptr, int *val)

cdef extern void photochemvars_equilibrium_time_get(void *ptr, double *val)
cdef extern void photochemvars_equilibrium_time_set(void *ptr, double *val)

cdef extern void photochemvars_verbose_get(void *ptr, int *val)
cdef extern void photochemvars_verbose_set(void *ptr, int *val)

cdef class PhotochemVars:
  cdef void *_ptr
  cdef bint _destroy

  def __cinit__(self, bint alloc = True):
    if alloc:
      allocate_photochemvars(&self._ptr)
      self._destroy = True
    else:
      self._destroy = False

  def __dealloc__(self):
    if self._destroy:
      deallocate_photochemvars(&self._ptr)
      self._ptr = NULL
  
  property nz:
    def __get__(self):
      cdef int nz
      photochemvars_nz_get(&self._ptr, &nz)
      return nz
      
  property usol_init:
    def __get__(self):
      cdef int dim1, dim2
      photochemvars_usol_init_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      photochemvars_usol_init_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property z:
    def __get__(self):
      cdef int dim1
      photochemvars_z_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemvars_z_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property temperature:
    def __get__(self):
      cdef int dim1
      photochemvars_temperature_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemvars_temperature_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property grav:
    def __get__(self):
      cdef int dim1
      photochemvars_grav_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemvars_grav_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property at_photo_equilibrium:
    def __get__(self):
      cdef bint at_photo_equilibrium
      photochemvars_at_photo_equilibrium_get(&self._ptr, &at_photo_equilibrium)
      return at_photo_equilibrium
      
  property rtol:
    def __get__(self):
      cdef double val
      photochemvars_rtol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      photochemvars_rtol_set(&self._ptr, &val)
      
  property atol:
    def __get__(self):
      cdef double val
      photochemvars_atol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      photochemvars_atol_set(&self._ptr, &val)

  property mxsteps:
    def __get__(self):
      cdef int val
      photochemvars_mxsteps_get(&self._ptr, &val)
      return val
    def __set__(self, int val):
      photochemvars_mxsteps_set(&self._ptr, &val)
      
  property equilibrium_time:
    def __get__(self):
      cdef double val
      photochemvars_equilibrium_time_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      photochemvars_equilibrium_time_set(&self._ptr, &val)
  
  property verbose:
    def __get__(self):
      cdef int val
      photochemvars_verbose_get(&self._ptr, &val)
      return val
    def __set__(self, int val):
      photochemvars_verbose_set(&self._ptr, &val)

    
    
    
    