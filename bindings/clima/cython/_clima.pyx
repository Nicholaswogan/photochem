from numpy cimport ndarray, uint8_t
from libc.stdint cimport uintptr_t
from libcpp cimport bool as cbool
from cpython.object cimport PyObject_GenericSetAttr
import numpy as np
import ctypes as ct
import os
import photochem_clima_data

DEF S_STR_LEN = 20;
DEF ERR_LEN = 1024;

include "futils.pyx"
include "RTChannel.pyx"
include "ClimaRadtranWrk.pyx"
include "Radtran.pyx"
include "AdiabatClimate.pyx"

# version
cdef extern void clima_version_get(char *version_c)

def _clima_version():
  cdef char version_c[100+1]
  clima_version_get(version_c)
  return version_c.decode("utf-8").strip()
  
__version__ = _clima_version()

# utils
cdef pystring2cstring(str pystring):
  # add a null c char, and convert to byes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring

cdef c2stringarr(ndarray c_str_arr, int str_len, int arr_len):  
  bs = c_str_arr[:-1].tobytes()
  return [bs[i:i+str_len].decode().strip() for i in range(0, str_len*arr_len, str_len)]

cdef list2cstring(list arr, int str_len):
  arr_c = np.zeros(len(arr)*str_len + 1, 'S1')
  for i in range(len(arr)):
    if len(arr[i]) > str_len:
          raise Exception('Failed to convert Python list to a C string')
    arr_c[i*str_len:(i+1)*str_len] = b' '
    arr_c[i*str_len:i*str_len+len(arr[i])] = np.array([elem.encode('utf-8') for elem in arr[i]])
  return arr_c

class ClimaException(Exception):
    pass
