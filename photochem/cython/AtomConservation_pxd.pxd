from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern void deallocate_atomconservation(void *ptr)

cdef extern void atomconservation_in_surf_get(void *ptr, double *val)
cdef extern void atomconservation_in_top_get(void *ptr, double *val)
cdef extern void atomconservation_in_dist_get(void *ptr, double *val)
cdef extern void atomconservation_in_other_get(void *ptr, double *val)
cdef extern void atomconservation_out_surf_get(void *ptr, double *val)
cdef extern void atomconservation_out_top_get(void *ptr, double *val)
cdef extern void atomconservation_out_rain_get(void *ptr, double *val)
cdef extern void atomconservation_out_other_get(void *ptr, double *val)
cdef extern void atomconservation_net_get(void *ptr, double *val)
cdef extern void atomconservation_factor_get(void *ptr, double *val)

# CondensationParameters
cdef extern void condensationparameters_k_cond_get(void *ptr, double *val)
cdef extern void condensationparameters_k_cond_set(void *ptr, double *val)

cdef extern void condensationparameters_k_evap_get(void *ptr, double *val)
cdef extern void condensationparameters_k_evap_set(void *ptr, double *val)

cdef extern void condensationparameters_rhc_get(void *ptr, double *val)
cdef extern void condensationparameters_rhc_set(void *ptr, double *val)

cdef extern void condensationparameters_smooth_factor_get(void *ptr, double *val)
cdef extern void condensationparameters_smooth_factor_set(void *ptr, double *val)