cdef extern from *:
  struct EvoAtmosphereState:
    pass

cdef extern EvoAtmosphereState *allocate_evoatmospherestate()
cdef extern void deallocate_evoatmospherestate(EvoAtmosphereState *ptr)
cdef extern void evoatmospherestate_initialize_from_h5_wrapper(EvoAtmosphereState *ptr, char *filename, char *err)
cdef extern void evoatmospherestate_save_to_h5_wrapper(EvoAtmosphereState *ptr, char *filename, char *err)
cdef extern void evoatmospherestate_metadata_get(EvoAtmosphereState *ptr, int *integers, double *reals, int *lengths)
cdef extern void evoatmospherestate_species_names_get(EvoAtmosphereState *ptr, int *n, char *names)
cdef extern void evoatmospherestate_array_get(EvoAtmosphereState *ptr, int *field, int *n, double *values)
cdef extern void evoatmospherestate_create_wrapper(EvoAtmosphereState *ptr, int *integers, double *reals,
    int *lengths, char *names, double *z, double *dz, double *temperature, double *edd,
    double *particle_radius, double *usol, double *profile_pressure, double *profile_temperature,
    double *profile_edd, char *err)

