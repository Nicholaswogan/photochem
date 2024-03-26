from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern void allocate_photochemwrk(void *ptr)
cdef extern void deallocate_photochemwrk(void *ptr)

cdef extern void photochemwrk_nsteps_get(void *ptr, int *val)

cdef extern void photochemwrk_t_history_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_t_history_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_mix_history_get_size(void *ptr, int *dim1, int *dim2, int *dim3)
cdef extern void photochemwrk_mix_history_get(void *ptr, int *dim1, int *dim2, int *dim3, double *arr)

cdef extern void photochemwrk_longdy_get(void *ptr, double *val)

cdef extern void photochemwrk_longdydt_get(void *ptr, double *val)

cdef extern void photochemwrk_tn_get(void *ptr, double *val)
cdef extern void photochemwrk_tn_set(void *ptr, double *val)

cdef extern void photochemwrk_usol_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_usol_get(void *ptr, int *dim1, int *dim2, double *usol)
cdef extern void photochemwrk_usol_set(void *ptr, int *dim1, int *dim2, double *usol)

cdef extern void photochemwrk_pressure_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_pressure_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_density_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_density_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_densities_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_densities_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_mubar_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_mubar_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_prates_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_prates_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_amean_grd_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_amean_grd_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_optical_depth_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_optical_depth_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_surf_radiance_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_surf_radiance_get(void *ptr, int *dim1, double *arr)

# PhotochemWrkEvo
cdef extern void photochemwrkevo_pressure_hydro_get_size(void *ptr, int *dim1)
cdef extern void photochemwrkevo_pressure_hydro_get(void *ptr, int *dim1, double *arr)