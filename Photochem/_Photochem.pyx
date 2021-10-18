from numpy cimport ndarray
import numpy as np
import os

cdef extern from "Photochem.h":
  cdef int S_STR_LEN "S_STR_LEN"

include "Atmosphere.pyx"
include "PhotochemData.pyx"
include "PhotochemVars.pyx"

    
    
    
    
    