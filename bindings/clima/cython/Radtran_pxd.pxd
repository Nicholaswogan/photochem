cimport RTChannel_pxd as rtc_pxd
cimport ClimaRadtranWrk_pxd as rwrk_pxd
from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass
cdef extern from *:
  struct Radtran:
    pass

cdef extern void radtran_set_bolometric_flux_wrapper(Radtran *ptr, double *flux)

cdef extern void radtran_bolometric_flux_wrapper(Radtran *ptr, double *flux)

cdef extern void radtran_skin_temperature_wrapper(Radtran *ptr, double *bond_albedo, double *T_skin)

cdef extern void radtran_equilibrium_temperature_wrapper(Radtran *ptr, double *bond_albedo, double *T_eq)

cdef extern void radtran_opacities2yaml_wrapper_1(Radtran *ptr, int *out_len, void **out_cp)
cdef extern void radtran_opacities2yaml_wrapper_2(Radtran *ptr, void **out_cp, int *out_len, char* out_c)

cdef extern void radtran_set_custom_optical_properties(
  Radtran *ptr, int *dim_wv, double *wv, int *dim_P, double *P,  
  int *dim1_dtau_dz, int *dim2_dtau_dz, double *dtau_dz, 
  int *dim1_w0, int *dim2_w0, double *w0, 
  int *dim1_g0, int *dim2_g0, double *g0, 
  char *err
)
cdef extern void radtran_unset_custom_optical_properties(
  Radtran *ptr
)

# getters and setters
cdef extern void radtran_zenith_u_get_size(Radtran *ptr, int *dim1)
cdef extern void radtran_zenith_u_get(Radtran *ptr, int *dim1, double *arr)
cdef extern void radtran_zenith_u_set(Radtran *ptr, int *dim1, double *arr)

cdef extern void radtran_surface_albedo_get_size(Radtran *ptr, int *dim1)
cdef extern void radtran_surface_albedo_get(Radtran *ptr, int *dim1, double *arr)
cdef extern void radtran_surface_albedo_set(Radtran *ptr, int *dim1, double *arr)

cdef extern void radtran_surface_emissivity_get_size(Radtran *ptr, int *dim1)
cdef extern void radtran_surface_emissivity_get(Radtran *ptr, int *dim1, double *arr)
cdef extern void radtran_surface_emissivity_set(Radtran *ptr, int *dim1, double *arr)

cdef extern void radtran_has_hard_surface_get(Radtran *ptr, bool *val)
cdef extern void radtran_has_hard_surface_set(Radtran *ptr, bool *val)

cdef extern void radtran_photon_scale_factor_get(Radtran *ptr, double *val)
cdef extern void radtran_photon_scale_factor_set(Radtran *ptr, double *val)

cdef extern void radtran_ir_tau_min_get(Radtran *ptr, double *val)
cdef extern void radtran_ir_tau_min_set(Radtran *ptr, double *val)

cdef extern void radtran_ir_get(Radtran *ptr, rtc_pxd.RTChannel **ptr1);

cdef extern void radtran_sol_get(Radtran *ptr, rtc_pxd.RTChannel **ptr1);

cdef extern void radtran_wrk_ir_get(Radtran *ptr, rwrk_pxd.ClimaRadtranWrk **ptr1);

cdef extern void radtran_wrk_sol_get(Radtran *ptr, rwrk_pxd.ClimaRadtranWrk **ptr1);
