from numpy cimport ndarray
import numpy as np
import os

cdef int S_STR_LEN = 20;
cdef int STR_LEN = 1024;
cdef int ERR_LEN = 1024;

include "Atmosphere.pyx"
include "PhotochemData.pyx"
include "PhotochemVars.pyx"

    
    
    
    
    