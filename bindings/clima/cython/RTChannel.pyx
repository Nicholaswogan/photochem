cimport RTChannel_pxd as rtc_pxd

cdef class RTChannel:
  cdef rtc_pxd.RTChannel *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass

  property wavl:
    "ndarray[double,ndim=1], shape (nw+1). The edges of the radiative transfer bins in wavelength (nm)"
    def __get__(self):
      cdef int dim1
      rtc_pxd.rtchannel_wavl_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rtc_pxd.rtchannel_wavl_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property freq:
    "ndarray[double,ndim=1], shape (nw+1). The edges of the radiative transfer bins in frequency (1/s)"
    def __get__(self):
      cdef int dim1
      rtc_pxd.rtchannel_freq_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rtc_pxd.rtchannel_freq_get(self._ptr, &dim1, <double *>arr.data)
      return arr