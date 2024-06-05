cimport AtomConservation_pxd as atom_pxd
from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct PhotochemData:
    pass

cdef extern void photochemdata_nq_get(PhotochemData *ptr, int *nq)
cdef extern void photochemdata_np_get(PhotochemData *ptr, int *nq)
cdef extern void photochemdata_nsp_get(PhotochemData *ptr, int *nq)
cdef extern void photochemdata_ng_get(PhotochemData *ptr, int *nq)
cdef extern void photochemdata_nsl_get(PhotochemData *ptr, int *nq)
cdef extern void photochemdata_nll_get(PhotochemData *ptr, int *nq)
cdef extern void photochemdata_nw_get(PhotochemData *ptr, int *nq)

cdef extern void photochemdata_planet_mass_get(PhotochemData *ptr, double *val)

cdef extern void photochemdata_planet_radius_get(PhotochemData *ptr, double *val)
cdef extern void photochemdata_planet_radius_set(PhotochemData *ptr, double *val)

cdef extern void photochemdata_species_names_get_size(PhotochemData *ptr, int *dim1)
cdef extern void photochemdata_species_names_get(PhotochemData *ptr, int *dim1, char* species_names)

cdef extern void photochemdata_atoms_names_get_size(PhotochemData *ptr, int *dim1)
cdef extern void photochemdata_atoms_names_get(PhotochemData *ptr, int *dim1, char* names)

cdef extern void photochemdata_reaction_equations_get_size(PhotochemData *ptr, int *dim1)
cdef extern void photochemdata_reaction_equations_get(PhotochemData *ptr, int *dim1, char* names)

cdef extern void photochemdata_photonums_get_size(PhotochemData *ptr, int *dim1)
cdef extern void photochemdata_photonums_get(PhotochemData *ptr, int *dim1, int *arr)

cdef extern void photochemdata_wavl_get_size(PhotochemData *ptr, int *dim1)
cdef extern void photochemdata_wavl_get(PhotochemData *ptr, int *dim1, double *arr)

cdef extern void photochemdata_species_mass_get_size(PhotochemData *ptr, int *dim1)
cdef extern void photochemdata_species_mass_get(PhotochemData *ptr, int *dim1, double *arr)

cdef extern void photochemdata_species_redox_get_size(PhotochemData *ptr, int *dim1)
cdef extern void photochemdata_species_redox_get(PhotochemData *ptr, int *dim1, double *arr)

cdef extern void photochemdata_particle_sat_get_size(PhotochemData *ptr, int *dim1)
cdef extern void photochemdata_particle_sat_get(PhotochemData *ptr, int *dim1, atom_pxd.SaturationData **ptr1)

