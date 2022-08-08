
cimport EvoAtmosphere_pxd as ea_pxd

cdef class EvoAtmosphere:
  cdef void *_ptr
  cdef void *_dat_ptr
  cdef void *_var_ptr
  cdef void *_wrk_ptr

  def __init__(self, mechanism_file = None, settings_file = None, 
                     flux_file = None, atmosphere_txt = None):           
    # Allocate memory
    ea_pxd.allocate_evoatmosphere(&self._ptr)
    
    # convert strings to char
    cdef bytes data_dir_b = pystring2cstring(os.path.dirname(os.path.realpath(__file__))+'/data')
    cdef char *data_dir_c = data_dir_b
    cdef bytes mechanism_file_b = pystring2cstring(mechanism_file)
    cdef char *mechanism_file_c = mechanism_file_b
    cdef bytes settings_file_b = pystring2cstring(settings_file)
    cdef char *settings_file_c = settings_file_b
    cdef bytes flux_file_b = pystring2cstring(flux_file)
    cdef char *flux_file_c = flux_file_b
    cdef bytes atmosphere_txt_b = pystring2cstring(atmosphere_txt)
    cdef char *atmosphere_txt_c = atmosphere_txt_b
    cdef char err[ERR_LEN+1]
    
    # Initialize
    ea_pxd.evoatmosphere_init_wrapper(&self._ptr, data_dir_c, mechanism_file_c,
                                       settings_file_c, flux_file_c,
                                       atmosphere_txt_c, 
                                       &self._dat_ptr, &self._var_ptr, &self._wrk_ptr,
                                       err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def __dealloc__(self):
    ea_pxd.deallocate_evoatmosphere(&self._ptr);
    
  property dat:
    def __get__(self):
      dat = PhotochemData(alloc = False)
      dat._ptr = self._dat_ptr
      return dat
      
  property var:
    def __get__(self):
      var = PhotochemVars(alloc = False)
      var._ptr = self._var_ptr
      return var
      
  property wrk:
    def __get__(self):
      wrk = PhotochemWrk(alloc = False)
      wrk._ptr = self._wrk_ptr
      return wrk
    
  def regrid_prep_atmosphere(self, ndarray[double, ndim=2] usol_, double top_atmos):
    cdef char err[ERR_LEN+1]
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef ndarray usol = np.asfortranarray(usol_)
    if usol.shape[0] != nq or usol.shape[1] != nz:
      raise PhotoException("Input usol is the wrong size.")
      
    ea_pxd.evoatmosphere_regrid_prep_atmosphere_wrapper(&self._ptr, &nq, &nz, <double *>usol.data, &top_atmos, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    
  def evolve(self, str filename, double tstart, ndarray[double, ndim=2] usol_, ndarray[double, ndim=1] t_eval, bool overwrite = False):
    cdef bytes filename_b = pystring2cstring(filename)
    cdef char *filename_c = filename_b
    cdef char err[ERR_LEN+1]
    cdef bool success
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    cdef int nt = t_eval.size
    cdef ndarray usol = np.asfortranarray(usol_)
    if usol.shape[0] != nq or usol.shape[1] != nz:
      raise PhotoException("Input usol is the wrong size.")
      
    ea_pxd.evoatmosphere_evolve_wrapper(&self._ptr, filename_c, &tstart, &nq, &nz, <double *>usol.data, &nt, <double *>t_eval.data, &overwrite, &success, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return success

  property T_surf:
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_t_surf_get(&self._ptr, &val)
      return val

  property T_trop:
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_t_trop_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_t_trop_set(&self._ptr, &val)

  property P_top_min:
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_p_top_min_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_p_top_min_set(&self._ptr, &val)

  property P_top_max:
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_p_top_max_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_p_top_max_set(&self._ptr, &val)

  property top_atmos_adjust_frac:
    def __get__(self):
      cdef double val
      ea_pxd.evoatmosphere_top_atmos_adjust_frac_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      ea_pxd.evoatmosphere_top_atmos_adjust_frac_set(&self._ptr, &val)
    
