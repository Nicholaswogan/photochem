from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct CondensationParameters:
    pass
  struct SaturationData:
    pass

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