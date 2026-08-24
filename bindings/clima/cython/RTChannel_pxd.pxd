
cdef extern from *:
  struct RTChannel:
    pass

# getters and setters

cdef extern void rtchannel_wavl_get_size(RTChannel *ptr, int *dim1)
cdef extern void rtchannel_wavl_get(RTChannel *ptr, int *dim1, double *arr)

cdef extern void rtchannel_freq_get_size(RTChannel *ptr, int *dim1)
cdef extern void rtchannel_freq_get(RTChannel *ptr, int *dim1, double *arr)