cimport AtomConservation_pxd as atom_pxd
cimport ProductionLoss_pxd as pl_pxd
cimport PhotochemData_pxd as dat_pxd
cimport PhotochemVars_pxd as var_pxd
cimport PhotochemWrk_pxd as wrk_pxd
from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct Atmosphere:
    pass

# callback signatures
ctypedef void (*time_dependent_rate_fcn)(double tn, int nz, double *rate)

# allocate and destroy
cdef extern Atmosphere *allocate_atmosphere();
cdef extern void deallocate_atmosphere(Atmosphere *ptr);

# subroutines
cdef extern void atmosphere_create_wrapper(Atmosphere *ptr, char *mechanism_file,
                                        char *settings_file, char *flux_file,
                                        char *atmosphere_txt, char *data_dir, char *err);

cdef extern void atmosphere_dat_get(Atmosphere *ptr, dat_pxd.PhotochemData **ptr1)
cdef extern void atmosphere_var_get(Atmosphere *ptr, var_pxd.PhotochemVars **ptr1)
cdef extern void atmosphere_wrk_get(Atmosphere *ptr, wrk_pxd.PhotochemWrk **ptr1)

cdef extern void atmosphere_check_for_convergence_wrapper(Atmosphere *ptr, bool *converged, char* err)
cdef extern void atmosphere_photochemical_equilibrium_wrapper(Atmosphere *ptr, bool *success, char* err)
  
cdef extern void atmosphere_out2atmosphere_txt_wrapper(Atmosphere *ptr, char *filename, int *number_of_decimals, bool *overwrite, bool *clip, char *err)                         
cdef extern void atmosphere_out2in_wrapper(Atmosphere *ptr, char *err)
cdef extern void atmosphere_gas_fluxes_wrapper(Atmosphere *ptr, int *nq, double *surf_fluxes, double *top_fluxes, char *err)
cdef extern void atmosphere_set_lower_bc_wrapper(Atmosphere *ptr, char *species, char *bc_type, 
                                                    double *vdep, double *mix, double *flux, double *height, bool *missing, char *err)
cdef extern void atmosphere_set_upper_bc_wrapper(Atmosphere *ptr, char *species, 
                                                    char *bc_type, double *veff, double *flux, bool *missing, char *err)

cdef extern void atmosphere_initialize_stepper_wrapper(Atmosphere *ptr, int *nq, int *nz, double *usol_start, char *err)
cdef extern double atmosphere_step_wrapper(Atmosphere *ptr, char *err)
cdef extern void atmosphere_destroy_stepper_wrapper(Atmosphere *ptr, char *err)

cdef extern void atmosphere_production_and_loss_wrapper(Atmosphere *ptr, char *species, int *nq, 
                                                        int *nz, double *usol, pl_pxd.ProductionLoss **pl_ptr, char *err)

cdef extern void atmosphere_prep_atmosphere_wrapper(Atmosphere *ptr, int *nq, int *nz, double *usol, char *err)

cdef extern void atmosphere_redox_conservation_wrapper(Atmosphere *ptr, double *redox_factor, char *err)
cdef extern void atmosphere_atom_conservation_wrapper(Atmosphere *ptr, char *atom, atom_pxd.AtomConservation **con_ptr, char *err)

cdef extern void atmosphere_evolve_wrapper(Atmosphere *ptr, char *filename, 
                double *tstart, int *nq, int *nz, double *usol, 
                int *nt, double *t_eval, bool *overwrite, bool *success, char *err)

cdef extern void atmosphere_set_press_temp_edd_wrapper(Atmosphere *ptr, int *P_dim1, double *P, int *T_dim1, double *T, int *edd_dim1, double *edd,
                                                      double *trop_p, bool *trop_p_present, char *err)

cdef extern void atmosphere_set_temperature_wrapper(Atmosphere *ptr, int *nz, double *temperature, 
                                                    double *trop_alt, bool *trop_alt_present, char *err)

cdef extern void atmosphere_set_rate_fcn_wrapper(Atmosphere *ptr, char *species_c, time_dependent_rate_fcn fcn, char *err)

cdef extern void atmosphere_update_vertical_grid_wrapper(Atmosphere *ptr, double *toa_alt, bool *toa_alt_present,
                                                         double *toa_pressure, bool *toa_pressure_present, char *err)