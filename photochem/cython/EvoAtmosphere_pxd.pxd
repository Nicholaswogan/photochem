cimport ProductionLoss_pxd as pl_pxd
cimport PhotochemData_pxd as dat_pxd
cimport PhotochemVars_pxd as var_pxd
cimport PhotochemWrk_pxd as wrk_pxd
from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct EvoAtmosphere:
    pass

# callback signatures
ctypedef double (*temp_dependent_albedo_fcn)(double T_surf)
ctypedef void (*time_dependent_rate_fcn)(double tn, int nz, double *rate)

# allocate and destroy
cdef extern EvoAtmosphere *allocate_evoatmosphere();
cdef extern void deallocate_evoatmosphere(EvoAtmosphere *ptr);

# subroutines
cdef extern void evoatmosphere_create_wrapper(EvoAtmosphere *ptr, char *mechanism_file,
                                            char *settings_file, char *flux_file,
                                            char *atmosphere_txt, char *data_dir, char *err);

cdef extern void evoatmosphere_dat_get(EvoAtmosphere *ptr, dat_pxd.PhotochemData **ptr1)
cdef extern void evoatmosphere_var_get(EvoAtmosphere *ptr, var_pxd.PhotochemVars **ptr1)
cdef extern void evoatmosphere_wrk_get(EvoAtmosphere *ptr, wrk_pxd.PhotochemWrk **ptr1)

cdef extern void evoatmosphere_prep_atmosphere_wrapper(EvoAtmosphere *ptr, int *nq, int *nz, double *usol, char *err)

cdef extern void evoatmosphere_out2atmosphere_txt_wrapper(EvoAtmosphere *ptr, char *filename, int *number_of_decimals, bool *overwrite, bool *clip, char *err)   
cdef extern void evoatmosphere_gas_fluxes_wrapper(EvoAtmosphere *ptr, int *nq, double *surf_fluxes, double *top_fluxes, char *err)
cdef extern void evoatmosphere_set_lower_bc_wrapper(EvoAtmosphere *ptr, char *species, char *bc_type, 
                                                    double *vdep, double *den, double *press, double *flux, double *height, bool *missing, char *err)
cdef extern void evoatmosphere_set_upper_bc_wrapper(EvoAtmosphere *ptr, char *species, 
                                                    char *bc_type, double *veff, double *flux, bool *missing, char *err)
cdef extern void evoatmosphere_set_rate_fcn_wrapper(EvoAtmosphere *ptr, char *species_c, time_dependent_rate_fcn fcn, char *err)
cdef extern void evoatmosphere_set_temperature_wrapper(EvoAtmosphere *ptr, int *nz, double *temperature, 
                                                    double *trop_alt, bool *trop_alt_present, char *err)
cdef extern void evoatmosphere_set_press_temp_edd_wrapper(EvoAtmosphere *ptr, int *P_dim1, double *P, int *T_dim1, double *T, int *edd_dim1, double *edd,
                                                      double *trop_p, bool *trop_p_present, 
                                                      bool *hydro_pressure, bool *hydro_pressure_present, char *err)
cdef extern void evoatmosphere_update_vertical_grid_wrapper(EvoAtmosphere *ptr, double *toa_alt, bool *toa_alt_present,
                                                         double *toa_pressure, bool *toa_pressure_present, char *err)
cdef extern void evoatmosphere_regrid_prep_atmosphere_wrapper(EvoAtmosphere *ptr, int *nq, int *nz, double *usol, double *top_atmos, char *err)

cdef extern void evoatmosphere_evolve_wrapper(EvoAtmosphere *ptr, char *filename, 
                double *tstart, int *nq, int *nz, double *usol, 
                int *nt, double *t_eval, bool *overwrite, bool *restart_from_file, bool *success, char *err)

cdef extern void evoatmosphere_check_for_convergence_wrapper(EvoAtmosphere *ptr, bool *converged, char* err)

cdef extern void evoatmosphere_initialize_stepper_wrapper(EvoAtmosphere *ptr, int *nq, int *nz, double *usol_start, char *err)

cdef extern double evoatmosphere_step_wrapper(EvoAtmosphere *ptr, char *err)

cdef extern void evoatmosphere_destroy_stepper_wrapper(EvoAtmosphere *ptr, char *err)

cdef extern void evoatmosphere_initialize_robust_stepper_wrapper(EvoAtmosphere *ptr, int *nq, int *nz, double *usol_start, char *err)

cdef extern void evoatmosphere_robust_step_wrapper(EvoAtmosphere *ptr, bool *give_up, bool *converged, char *err)

cdef extern void evoatmosphere_production_and_loss_wrapper(EvoAtmosphere *ptr, char *species, int *nq, 
                                                        int *nz, double *usol, pl_pxd.ProductionLoss **pl_ptr, char *err)

# getters and setters
cdef extern void evoatmosphere_t_surf_get(EvoAtmosphere *ptr, double *val)
cdef extern void evoatmosphere_t_surf_set(EvoAtmosphere *ptr, double *val)

cdef extern void evoatmosphere_t_trop_get(EvoAtmosphere *ptr, double *val)
cdef extern void evoatmosphere_t_trop_set(EvoAtmosphere *ptr, double *val)

cdef extern void evoatmosphere_albedo_fcn_set(EvoAtmosphere *ptr, temp_dependent_albedo_fcn fcn)

cdef extern void evoatmosphere_p_top_min_get(EvoAtmosphere *ptr, double *val)
cdef extern void evoatmosphere_p_top_min_set(EvoAtmosphere *ptr, double *val)

cdef extern void evoatmosphere_p_top_max_get(EvoAtmosphere *ptr, double *val)
cdef extern void evoatmosphere_p_top_max_set(EvoAtmosphere *ptr, double *val)

cdef extern void evoatmosphere_top_atmos_adjust_frac_get(EvoAtmosphere *ptr, double *val)
cdef extern void evoatmosphere_top_atmos_adjust_frac_set(EvoAtmosphere *ptr, double *val)
