from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from "PhotochemWrk.h":
  struct PhotochemWrk:
    pass
  struct PhotochemWrkEvo:
    pass

cdef extern void photochemwrk_nsteps_total_get(PhotochemWrk *ptr, int *val)

cdef extern void photochemwrk_nsteps_get(PhotochemWrk *ptr, int *val)

cdef extern void photochemwrk_t_history_get_size(PhotochemWrk *ptr, int *dim1)
cdef extern void photochemwrk_t_history_get(PhotochemWrk *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_mix_history_get_size(PhotochemWrk *ptr, int *dim1, int *dim2, int *dim3)
cdef extern void photochemwrk_mix_history_get(PhotochemWrk *ptr, int *dim1, int *dim2, int *dim3, double *arr)

cdef extern void photochemwrk_longdy_get(PhotochemWrk *ptr, double *val)

cdef extern void photochemwrk_longdydt_get(PhotochemWrk *ptr, double *val)

cdef extern void photochemwrk_tn_get(PhotochemWrk *ptr, double *val)
cdef extern void photochemwrk_tn_set(PhotochemWrk *ptr, double *val)

cdef extern void photochemwrk_usol_get_size(PhotochemWrk *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_usol_get(PhotochemWrk *ptr, int *dim1, int *dim2, double *usol)
cdef extern void photochemwrk_usol_set(PhotochemWrk *ptr, int *dim1, int *dim2, double *usol)

cdef extern void photochemwrk_pressure_get_size(PhotochemWrk *ptr, int *dim1)
cdef extern void photochemwrk_pressure_get(PhotochemWrk *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_density_get_size(PhotochemWrk *ptr, int *dim1)
cdef extern void photochemwrk_density_get(PhotochemWrk *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_densities_get_size(PhotochemWrk *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_densities_get(PhotochemWrk *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_rx_rates_get_size(PhotochemWrk *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_rx_rates_get(PhotochemWrk *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_mubar_get_size(PhotochemWrk *ptr, int *dim1)
cdef extern void photochemwrk_mubar_get(PhotochemWrk *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_prates_get_size(PhotochemWrk *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_prates_get(PhotochemWrk *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_amean_grd_get_size(PhotochemWrk *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_amean_grd_get(PhotochemWrk *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_optical_depth_get_size(PhotochemWrk *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_optical_depth_get(PhotochemWrk *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemwrk_surf_radiance_get_size(PhotochemWrk *ptr, int *dim1)
cdef extern void photochemwrk_surf_radiance_get(PhotochemWrk *ptr, int *dim1, double *arr)

# PhotochemWrkEvo
cdef extern void photochemwrkevo_pressure_hydro_get_size(PhotochemWrkEvo *ptr, int *dim1)
cdef extern void photochemwrkevo_pressure_hydro_get(PhotochemWrkEvo *ptr, int *dim1, double *arr)