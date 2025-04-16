from numpy cimport ndarray
from libcpp cimport bool
from libc.stdlib cimport malloc, free
from libc.stdint cimport uintptr_t
from cpython.object cimport PyObject_GenericSetAttr
from cpython.exc cimport PyErr_CheckSignals
import numpy as np
import ctypes as ct
import os
import photochem_clima_data

DEF S_STR_LEN = 20;
DEF M_STR_LEN = 100;
DEF STR_LEN = 1024;
DEF ERR_LEN = 1024;

include "EvoAtmosphere.pyx"
include "PhotochemData.pyx"
include "PhotochemVars.pyx"
include "PhotochemWrk.pyx"
include "ProductionLoss.pyx"
include "AtomConservation.pyx"

# version
cdef extern void photochem_version_get(char *version_c)
  
def _photochem_version():
  cdef char version_c[100+1]
  photochem_version_get(version_c)
  return version_c.decode("utf-8").strip()
  
__version__ = _photochem_version()

# utils
cdef pystring2cstring(str pystring):
  # add a null c char, and convert to byes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring

cdef c2stringarr(ndarray c_str_arr, int str_len, int arr_len):  
  bs = c_str_arr[:-1].tobytes()
  return [bs[i:i+str_len].decode().strip() for i in range(0, str_len*arr_len, str_len)]

class PhotoException(Exception):
    pass
