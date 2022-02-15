from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern void allocate_photochemwrk(void *ptr)
cdef extern void deallocate_photochemwrk(void *ptr)

cdef extern void photochemwrk_usol_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void photochemwrk_usol_get(void *ptr, int *dim1, int *dim2, double *usol)
cdef extern void photochemwrk_usol_set(void *ptr, int *dim1, int *dim2, double *usol)

cdef extern void photochemwrk_pressure_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_pressure_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemwrk_density_get_size(void *ptr, int *dim1)
cdef extern void photochemwrk_density_get(void *ptr, int *dim1, double *arr)

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