import numpy as __np
from numpy cimport ndarray 
import os

cdef extern void allocate_atmosphere(void *ptr);
cdef extern void deallocate_atmosphere(void *ptr);

cdef extern void atmosphere_init_wrapper(void *ptr, char *data_dir, char *mechanism_file,
                                        char *settings_file, char *flux_file,
                                        char *atmosphere_txt, char *err);
cdef extern void atmosphere_photochemical_equilibrium_wrapper(void *ptr, bint *success, char* err)                               

class PhotoException(Exception):
    pass
    
cdef class Atmosphere:
  cdef void *_ptr
  
  def __init__(self, mechanism_file = None, settings_file = None, 
                     flux_file = None, atmosphere_txt = None):           
    allocate_atmosphere(&self._ptr)
    data_dir_b = (os.path.dirname(os.path.realpath(__file__))+'/../data').encode('utf-8')
    cdef char *data_dir_c = data_dir_b
    mechanism_file_b = mechanism_file.encode('utf-8')
    cdef char *mechanism_file_c = mechanism_file_b
    settings_file_b = settings_file.encode('utf-8')
    cdef char *settings_file_c = settings_file_b
    flux_file_b = flux_file.encode('utf-8')
    cdef char *flux_file_c = flux_file_b
    atmosphere_txt_b = atmosphere_txt.encode('utf-8')
    cdef char *atmosphere_txt_c = atmosphere_txt_b
    cdef char err[1024+1]
    atmosphere_init_wrapper(&self._ptr, data_dir_c, mechanism_file_c,
                                       settings_file_c, flux_file_c,
                                       atmosphere_txt_c, err)
    if len(err.strip()) > 0:
      deallocate_atmosphere(&self._ptr)
      raise PhotoException(err.decode("utf-8").strip())

  def __dealloc__(self):
    deallocate_atmosphere(&self._ptr);
    
  def photochemical_equilibrium(self):
    cdef bint success
    cdef char err[1024+1]
    atmosphere_photochemical_equilibrium_wrapper(&self._ptr, &success, err)
    if len(err.strip()) > 0:
      raise PhotoException(err.decode("utf-8").strip())
    return success
    
    
    
    
    
    
    
    