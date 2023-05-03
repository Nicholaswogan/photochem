from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

# callback signatures
ctypedef double (*temp_dependent_albedo_fcn)(double T_surf)

# allocate and destroy
cdef extern void allocate_evoatmosphere(void *ptr);
cdef extern void deallocate_evoatmosphere(void *ptr);

# subroutines
cdef extern void evoatmosphere_create_wrapper(void *ptr, char *data_dir, char *mechanism_file,
                                            char *settings_file, char *flux_file,
                                            char *atmosphere_txt, void *dat_ptr, void *var_ptr,
                                            void *wrk_ptr, char *err);

cdef extern void evoatmosphere_out2atmosphere_txt_wrapper(void *ptr, char *filename, bool *overwrite, bool *clip, char *err)   

cdef extern void evoatmosphere_regrid_prep_atmosphere_wrapper(void *ptr, int *nq, int *nz, double *usol, double *top_atmos, char *err)

cdef extern void evoatmosphere_evolve_wrapper(void *ptr, char *filename, 
                double *tstart, int *nq, int *nz, double *usol, 
                int *nt, double *t_eval, bool *overwrite, bool *restart_from_file, bool *success, char *err)

cdef extern void evoatmosphere_production_and_loss_wrapper(void *ptr, char *species, int *nq, 
                                                        int *nz, double *usol, double *top_atmos, void *pl_ptr, char *err)

cdef extern void evoatmosphere_set_albedo_fcn_wrapper(void *ptr, temp_dependent_albedo_fcn fcn)

# getters and setters
cdef extern void evoatmosphere_t_surf_get(void *ptr, double *val)
cdef extern void evoatmosphere_t_surf_set(void *ptr, double *val)

cdef extern void evoatmosphere_t_trop_get(void *ptr, double *val)
cdef extern void evoatmosphere_t_trop_set(void *ptr, double *val)

cdef extern void evoatmosphere_p_top_min_get(void *ptr, double *val)
cdef extern void evoatmosphere_p_top_min_set(void *ptr, double *val)

cdef extern void evoatmosphere_p_top_max_get(void *ptr, double *val)
cdef extern void evoatmosphere_p_top_max_set(void *ptr, double *val)

cdef extern void evoatmosphere_top_atmos_adjust_frac_get(void *ptr, double *val)
cdef extern void evoatmosphere_top_atmos_adjust_frac_set(void *ptr, double *val)
