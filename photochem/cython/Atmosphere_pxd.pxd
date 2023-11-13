from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

# callback signatures
ctypedef void (*time_dependent_rate_fcn)(double tn, int nz, double *rate)

# allocate and destroy
cdef extern void allocate_atmosphere(void *ptr);
cdef extern void deallocate_atmosphere(void *ptr);

# subroutines
cdef extern void atmosphere_create_wrapper(void *ptr, char *mechanism_file,
                                        char *settings_file, char *flux_file,
                                        char *atmosphere_txt, char *data_dir, void *dat_ptr, void *var_ptr,
                                        void *wrk_ptr, char *err);

cdef extern void atmosphere_photochemical_equilibrium_wrapper(void *ptr, bool *success, char* err)
  
cdef extern void atmosphere_out2atmosphere_txt_wrapper(void *ptr, char *filename, bool *overwrite, bool *clip, char *err)                         
cdef extern void atmosphere_out2in_wrapper(void *ptr, char *err)
cdef extern void atmosphere_gas_fluxes_wrapper(void *ptr, double *surf_fluxes, double *top_fluxes, char *err)
cdef extern void atmosphere_set_lower_bc_wrapper(void *ptr, char *species, char *bc_type, 
                                                    double *vdep, double *mix, double *flux, double *height, bool *missing, char *err)
cdef extern void atmosphere_set_upper_bc_wrapper(void *ptr, char *species, 
                                                    char *bc_type, double *veff, double *flux, bool *missing, char *err)

cdef extern void atmosphere_initialize_stepper_wrapper(void *ptr, int *nq, int *nz, double *usol_start, char *err)
cdef extern double atmosphere_step_wrapper(void *ptr, char *err)
cdef extern void atmosphere_destroy_stepper_wrapper(void *ptr, char *err)

cdef extern void atmosphere_production_and_loss_wrapper(void *ptr, char *species, int *nq, 
                                                        int *nz, double *usol, void *pl_ptr, char *err)

cdef extern void atmosphere_prep_atmosphere_wrapper(void *ptr, int *nq, int *nz, double *usol, char *err)

cdef extern void atmosphere_redox_conservation_wrapper(void *ptr, double *redox_factor, char *err)
cdef extern void atmosphere_atom_conservation_wrapper(void *ptr, char *atom, void *con_ptr, char *err)

cdef extern void atmosphere_evolve_wrapper(void *ptr, char *filename, 
                double *tstart, int *nq, int *nz, double *usol, 
                int *nt, double *t_eval, bool *overwrite, bool *success, char *err)

cdef extern void atmosphere_set_press_temp_edd_wrapper(void *ptr, int *P_dim1, double *P, int *T_dim1, double *T, int *edd_dim1, double *edd,
                                                      double *trop_p, bool *trop_p_present, char *err)

cdef extern void atmosphere_set_temperature_wrapper(void *ptr, int *nz, double *temperature, 
                                                    double *trop_alt, bool *trop_alt_present, char *err)

cdef extern void atmosphere_set_rate_fcn_wrapper(void *ptr, char *species_c, time_dependent_rate_fcn fcn, char *err)

cdef extern void atmosphere_update_vertical_grid_wrapper(void *ptr, double *toa_alt, bool *toa_alt_present,
                                                         double *toa_pressure, bool *toa_pressure_present, char *err)