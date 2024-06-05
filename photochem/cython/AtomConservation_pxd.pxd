from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct AtomConservation:
    pass
  struct CondensationParameters:
    pass
  struct SaturationData:
    pass

cdef extern void deallocate_atomconservation(AtomConservation *ptr)

cdef extern void atomconservation_in_surf_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_in_top_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_in_dist_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_in_other_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_out_surf_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_out_top_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_out_rain_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_out_other_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_net_get(AtomConservation *ptr, double *val)
cdef extern void atomconservation_factor_get(AtomConservation *ptr, double *val)

# CondensationParameters
cdef extern void condensationparameters_k_cond_get(CondensationParameters *ptr, double *val)
cdef extern void condensationparameters_k_cond_set(CondensationParameters *ptr, double *val)

cdef extern void condensationparameters_k_evap_get(CondensationParameters *ptr, double *val)
cdef extern void condensationparameters_k_evap_set(CondensationParameters *ptr, double *val)

cdef extern void condensationparameters_rhc_get(CondensationParameters *ptr, double *val)
cdef extern void condensationparameters_rhc_set(CondensationParameters *ptr, double *val)

cdef extern void condensationparameters_smooth_factor_get(CondensationParameters *ptr, double *val)
cdef extern void condensationparameters_smooth_factor_set(CondensationParameters *ptr, double *val)

# SaturationData
cdef extern void saturationdata_sat_pressure_wrapper(SaturationData *ptr, double *T, double *Psat)