from numpy cimport ndarray
from libcpp cimport bool
import numpy as np
import ctypes as ct
import os

DEF S_STR_LEN = 20;
DEF M_STR_LEN = 100;
DEF STR_LEN = 1024;
DEF ERR_LEN = 1024;

include "Atmosphere.pyx"
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
