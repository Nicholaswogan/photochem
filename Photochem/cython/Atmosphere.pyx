
# allocate and destroy
cdef extern void allocate_atmosphere(void *ptr);
cdef extern void deallocate_atmosphere(void *ptr);

# subroutines
cdef extern void atmosphere_init_wrapper(void *ptr, char *data_dir, char *mechanism_file,
                                        char *settings_file, char *flux_file,
                                        char *atmosphere_txt, void *dat_ptr, void *var_ptr,
                                        void *wrk_ptr, char *err);
cdef extern void atmosphere_photochemical_equilibrium_wrapper(void *ptr, bint *success, char* err)   
cdef extern void atmosphere_out2atmosphere_txt_wrapper(void *ptr, char *filename, bint *overwrite, bint *clip, char *err)                         
cdef extern void atmosphere_out2in_wrapper(void *ptr, char *err)
cdef extern void atmosphere_surface_fluxes_wrapper(void *ptr, double *fluxes, char *err)

class PhotoException(Exception):
    pass
    
cdef class Atmosphere:
  cdef void *_ptr
  cdef void *_dat_ptr
  cdef void *_var_ptr
  cdef void *_wrk_ptr
  
  def __init__(self, mechanism_file = None, settings_file = None, 
                     flux_file = None, atmosphere_txt = None):           
    # Allocate memory
    allocate_atmosphere(&self._ptr)
    
    # convert strings to char
    cdef bytes data_dir_b = pystring2cstring(os.path.dirname(os.path.realpath(__file__))+'/../data')
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
    atmosphere_init_wrapper(&self._ptr, data_dir_c, mechanism_file_c,
                                       settings_file_c, flux_file_c,
                                       atmosphere_txt_c, 
                                       &self._dat_ptr, &self._var_ptr, &self._wrk_ptr,
                                       err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())

  def __dealloc__(self):
    deallocate_atmosphere(&self._ptr);
    
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
      
  def photochemical_equilibrium(self):
    cdef bint success
    cdef char err[ERR_LEN+1]
    atmosphere_photochemical_equilibrium_wrapper(&self._ptr, &success, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return success
    
  def out2atmosphere_txt(self,filename = None, bint overwrite = False, bint clip = True):
    cdef bytes filename_b = pystring2cstring(filename)
    cdef char *filename_c = filename_b
    cdef char err[ERR_LEN+1]
    atmosphere_out2atmosphere_txt_wrapper(&self._ptr, filename_c, &overwrite, &clip, err)  
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
      
  def out2in(self):
    cdef char err[ERR_LEN+1]
    atmosphere_out2in_wrapper(&self._ptr, err)  
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
      
  def out_dict(self):
    if not self.var.at_photo_equilibrium:
      raise PhotoException("Must integrate to photochemical equilibrium before making output dictionary")
    out = {}
    out['alt'] = self.var.z/1e5
    out['temp'] = self.var.temperature
    out['pressure'] = self.wrk.pressure
    out['density'] = self.wrk.density
    names = self.dat.species_names
    cdef ndarray usol = self.wrk.usol
    cdef int i
    for i in range(self.dat.nq):
      out[names[i]] = usol[i,:]
    return out
    
  def surface_fluxes(self):
    cdef ndarray fluxes = np.empty(self.dat.nq, np.double)
    cdef char err[ERR_LEN+1]
    atmosphere_surface_fluxes_wrapper(&self._ptr, <double *>fluxes.data, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    out = {}
    names = self.dat.species_names
    for i in range(self.dat.nq):
      out[names[i]] = fluxes[i]
    return out
      
cdef pystring2cstring(str pystring):
  # add a null c char, and convert to byes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring
  
    
    
    
    
    