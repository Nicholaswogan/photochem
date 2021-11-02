
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
cdef extern void atmosphere_set_lower_bc_wrapper(void *ptr, char *species, char *bc_type, 
                                                    double *vdep, double *mix, double *flux, double *height, bint *missing, char *err)
cdef extern void atmosphere_set_upper_bc_wrapper(void *ptr, char *species, 
                                                    char *bc_type, double *veff, double *flux, bint *missing, char *err)

cdef extern void atmosphere_initialize_stepper_wrapper(void *ptr, int *nq, int *nz, double *usol_start, char *err)
cdef extern double atmosphere_step_wrapper(void *ptr, char *err)
cdef extern void atmosphere_destroy_stepper_wrapper(void *ptr, char *err)

cdef extern void atmosphere_production_and_loss_wrapper(void *ptr, char *species, int *nq, 
                                                        int *nz, double *usol, void *pl_ptr, char *err)

cdef pystring2cstring(str pystring):
  # add a null c char, and convert to byes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring

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
    
  def initialize_stepper(self, ndarray[double, ndim=2] usol_start):
    cdef char err[ERR_LEN+1]
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    if usol_start.shape[0] != nq or usol_start.shape[1] != nz:
      raise PhotoException("Input usol_start is the wrong size.")
    atmosphere_initialize_stepper_wrapper(&self._ptr, &nq, &nz,  <double *>usol_start.data, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    
  def step(self):
    cdef char err[ERR_LEN+1]
    cdef double tn = atmosphere_step_wrapper(&self._ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return tn
    
  def destroy_stepper(self):
    cdef char err[ERR_LEN+1]
    atmosphere_destroy_stepper_wrapper(&self._ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
      
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
      
  def mole_fraction_dict(self):
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
    
  def set_lower_bc(self, species = None, bc_type = None, vdep = None, mix = None,
                            flux = None, height = None):
    
    cdef bytes species_b = pystring2cstring(species)
    cdef char *species_c = species_b
    cdef bytes bc_type_b = pystring2cstring(bc_type)
    cdef char *bc_type_c = bc_type_b
    cdef double vdep_c = 0
    cdef double mix_c = 0
    cdef double flux_c = 0
    cdef double height_c = 0
    cdef bint missing = False
    if bc_type == 'vdep':
      if vdep == None:
        missing = True
      else:
        vdep_c = vdep
    elif bc_type == 'mix':
      if mix == None:
        missing = True
      else:
        mix_c = mix
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
    atmosphere_set_lower_bc_wrapper(&self._ptr, species_c, bc_type_c, 
                                      &vdep_c, &mix_c, &flux_c, &height_c, &missing, err);
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
  
  def set_upper_bc(self, species = None, bc_type = None, veff = None,flux = None):
    
    cdef bytes species_b = pystring2cstring(species)
    cdef char *species_c = species_b
    cdef bytes bc_type_b = pystring2cstring(bc_type)
    cdef char *bc_type_c = bc_type_b
    cdef double veff_c = 0
    cdef double flux_c = 0
    cdef bint missing = False
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
    atmosphere_set_upper_bc_wrapper(&self._ptr, species_c, bc_type_c, 
                                      &veff_c, &flux_c, &missing, err);
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
      
  def production_and_loss(self, str species, ndarray[double, ndim=2] usol):
    cdef bytes species_b = pystring2cstring(species)
    cdef char *species_c = species_b
    cdef char err[ERR_LEN+1]
    
    cdef int nq = self.dat.nq
    cdef int nz = self.var.nz
    if not np.isfortran(usol):
      raise PhotoException("Input usol must have fortran ordering.")
    if usol.shape[0] != nq or usol.shape[1] != nz:
      raise PhotoException("Input usol is the wrong size.")
      
    cdef void *pl_ptr
    atmosphere_production_and_loss_wrapper(&self._ptr, species_c, &nq, &nz, <double *>usol.data, &pl_ptr, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    pl = ProductionLoss()
    pl._ptr = pl_ptr
    return pl


