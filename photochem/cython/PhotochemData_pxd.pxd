from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern void allocate_photochemdata(void *ptr)
cdef extern void deallocate_photochemdata(void *ptr)

cdef extern void photochemdata_nq_get(void *ptr, int *nq)
cdef extern void photochemdata_np_get(void *ptr, int *nq)
cdef extern void photochemdata_nsp_get(void *ptr, int *nq)
cdef extern void photochemdata_ng_get(void *ptr, int *nq)
cdef extern void photochemdata_nsl_get(void *ptr, int *nq)
cdef extern void photochemdata_nll_get(void *ptr, int *nq)
cdef extern void photochemdata_nw_get(void *ptr, int *nq)

cdef extern void photochemdata_planet_mass_get(void *ptr, double *val)

cdef extern void photochemdata_planet_radius_get(void *ptr, double *val)

cdef extern void photochemdata_species_names_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_species_names_get(void *ptr, int *dim1, char* species_names)

cdef extern void photochemdata_atoms_names_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_atoms_names_get(void *ptr, int *dim1, char* names)

cdef extern void photochemdata_reaction_equations_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_reaction_equations_get(void *ptr, int *dim1, char* names)

cdef extern void photochemdata_photonums_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_photonums_get(void *ptr, int *dim1, int *arr)

cdef extern void photochemdata_wavl_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_wavl_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemdata_species_mass_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_species_mass_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemdata_species_redox_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_species_redox_get(void *ptr, int *dim1, double *arr)

