from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern void allocate_photochemvars(void *ptr)
cdef extern void deallocate_photochemvars(void *ptr)

cdef extern void photochemvars_nz_get(void *ptr, int *nz)

cdef extern void photochemvars_top_atmos_get(void *ptr, double *val)

cdef extern void photochemvars_bottom_atmos_get(void *ptr, double *val)

cdef extern void photochemvars_at_photo_equilibrium_get(void *ptr, bool *at_photo_equilibrium)

cdef extern void photochemvars_usol_init_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemvars_usol_init_get(void *ptr, int *dim1, int *dim2, double *usol_init)

cdef extern void photochemvars_temperature_get_size(void *ptr, int *dim1)
cdef extern void photochemvars_temperature_get(void *ptr, int *dim1, double *temperature)

cdef extern void photochemvars_edd_get_size(void *ptr, int *dim1)
cdef extern void photochemvars_edd_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemvars_photon_flux_get_size(void *ptr, int *dim1)
cdef extern void photochemvars_photon_flux_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemvars_grav_get_size(void *ptr, int *dim1)
cdef extern void photochemvars_grav_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemvars_z_get_size(void *ptr, int *dim1)
cdef extern void photochemvars_z_get(void *ptr, int *dim1, double *z)

cdef extern void photochemvars_surface_pressure_get(void *ptr, double *val)
cdef extern void photochemvars_surface_pressure_set(void *ptr, double *val)

cdef extern void photochemvars_max_error_reinit_attempts_get(void *ptr, int *val)
cdef extern void photochemvars_max_error_reinit_attempts_set(void *ptr, int *val)

cdef extern void photochemvars_rtol_get(void *ptr, double *val)
cdef extern void photochemvars_rtol_set(void *ptr, double *val)

cdef extern void photochemvars_atol_get(void *ptr, double *val)
cdef extern void photochemvars_atol_set(void *ptr, double *val)

cdef extern void photochemvars_mxsteps_get(void *ptr, int *val)
cdef extern void photochemvars_mxsteps_set(void *ptr, int *val)

cdef extern void photochemvars_equilibrium_time_get(void *ptr, double *val)
cdef extern void photochemvars_equilibrium_time_set(void *ptr, double *val)

cdef extern void photochemvars_verbose_get(void *ptr, int *val)
cdef extern void photochemvars_verbose_set(void *ptr, int *val)

cdef extern void photochemvars_fast_arbitrary_rate_get(void *ptr, double *val)
cdef extern void photochemvars_fast_arbitrary_rate_set(void *ptr, double *val)