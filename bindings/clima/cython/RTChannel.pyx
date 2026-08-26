cimport RTChannel_pxd as rtc_pxd

cdef class RTChannel:
  """Non-owning view of spectral-grid metadata for one radiative channel.

  Obtain this object from ``Radtran.ir`` or ``Radtran.sol`` and keep the parent
  climate model alive while using it. Array-valued properties return copies.
  """
  cdef rtc_pxd.RTChannel *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass

  property wavl:
    """ndarray, shape (nw + 1,): Read-only copy of spectral-bin edges (nm)."""
    def __get__(self):
      cdef int dim1
      rtc_pxd.rtchannel_wavl_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rtc_pxd.rtchannel_wavl_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property freq:
    """ndarray, shape (nw + 1,): Read-only copy of spectral-bin edges (Hz)."""
    def __get__(self):
      cdef int dim1
      rtc_pxd.rtchannel_freq_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rtc_pxd.rtchannel_freq_get(self._ptr, &dim1, <double *>arr.data)
      return arr
