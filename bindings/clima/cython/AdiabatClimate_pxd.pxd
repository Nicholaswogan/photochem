cimport Radtran_pxd as rad_pxd
from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct AdiabatClimate:
    pass

# callback signatures
ctypedef double (*temp_dependent_albedo_fcn)(double T_surf)
ctypedef void (*ocean_solubility_fcn)(double T_surf, int ng, double *P_i, double *m_i, void *args_p)

# allocate and destroy
cdef extern AdiabatClimate *allocate_adiabatclimate();
cdef extern void deallocate_adiabatclimate(AdiabatClimate *ptr);

# wrappers for functions
cdef extern void adiabatclimate_create_wrapper(AdiabatClimate *ptr, char *species_file,
                                                    char *settings_file, char *flux_file, char *data_dir, bool *double_radiative_grid, char *err);

cdef extern void adiabatclimate_make_profile_wrapper(AdiabatClimate *ptr, double *T_surf, int *ng, 
                                      double *P_i_surf, char *err)

cdef extern void adiabatclimate_make_column_wrapper(AdiabatClimate *ptr, double *T_surf, int *ng, 
                                      double *N_i_surf, char *err)

cdef extern void adiabatclimate_make_profile_bg_gas_wrapper(AdiabatClimate *ptr, double *T_surf, int *ng, 
                                      double *P_i_surf, double *P_surf, char *bg_gas, char *err)

cdef extern void adiabatclimate_make_profile_dry_wrapper(
  AdiabatClimate *ptr, 
  int *dim_P, double *P,
  int *dim_T, double *T,
  int *dim1_f_i, int *dim2_f_i, double *f_i,
  char *err
)

cdef extern void adiabatclimate_toa_fluxes_wrapper(AdiabatClimate *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *ISR, double *OLR, char *err)

cdef extern void adiabatclimate_toa_fluxes_column_wrapper(AdiabatClimate *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *ISR, double *OLR, char *err)

cdef extern void adiabatclimate_toa_fluxes_bg_gas_wrapper(AdiabatClimate *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *P_surf, char *bg_gas, double *ISR, 
                                     double *OLR, char *err)

cdef extern void adiabatclimate_toa_fluxes_dry_wrapper(
  AdiabatClimate *ptr, 
  int *dim_P, double *P,
  int *dim_T, double *T,
  int *dim1_f_i, int *dim2_f_i, double *f_i,
  double *ISR, double *OLR,
  char *err
)

cdef extern void adiabatclimate_surface_temperature_wrapper(AdiabatClimate *ptr, int *ng, 
                                      double *P_i_surf, double *T_guess, double *T_surf, char *err)

cdef extern void adiabatclimate_surface_temperature_column_wrapper(AdiabatClimate *ptr, int *ng, 
                                      double *N_i_surf, double *T_guess, double *T_surf, char *err)

cdef extern void adiabatclimate_surface_temperature_bg_gas_wrapper(AdiabatClimate *ptr, int *ng, 
                                      double *P_i_surf, double *P_surf, char *bg_gas, 
                                      double *T_guess, double *T_surf, char *err)

cdef extern void adiabatclimate_set_particle_density_and_radii(
  AdiabatClimate *ptr, int *dim_P, double *P,  
  int *dim1_pdensities, int *dim2_pdensities, double *pdensities, 
  int *dim1_pradii, int *dim2_pradii, double *pradii, 
  char *err
)

cdef extern void adiabatclimate_rce_wrapper(
  AdiabatClimate *ptr, int *ng, double *P_i_surf, double *T_surf_guess, int *dim_T_guess, double *T_guess, 
  bool *convecting_with_below_present, int *dim_convecting_with_below, bool *convecting_with_below,
  bool *custom_present, int *dim_sp_custom, char* sp_custom, int *dim_P_custom, double *P_custom, 
  int *dim1_mix_custom, int *dim2_mix_custom, double *mix_custom,
  bool *converged, char *err
)

cdef extern void adiabatclimate_set_ocean_solubility_fcn_wrapper(AdiabatClimate *ptr, char *species_c, ocean_solubility_fcn fcn, char *err)

cdef extern void adiabatclimate_to_regular_grid_wrapper(AdiabatClimate *ptr, char *err)

cdef extern void adiabatclimate_out2atmosphere_txt_wrapper(AdiabatClimate *ptr, char *filename, int *nz, double *eddy, int *number_of_decimals, bool *overwrite, bool *clip, char *err)

cdef extern void adiabatclimate_heat_redistribution_parameters_wrapper(AdiabatClimate *ptr, double *tau_LW, double *k_term, double *f_term, char *err)

# getters and setters
cdef extern void adiabatclimate_p_top_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_p_top_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_reference_pressure_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_reference_pressure_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_t_trop_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_t_trop_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_use_make_column_p_guess_get(AdiabatClimate *ptr, bool *val)
cdef extern void adiabatclimate_use_make_column_p_guess_set(AdiabatClimate *ptr, bool *val)

cdef extern void adiabatclimate_make_column_p_guess_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_make_column_p_guess_get(AdiabatClimate *ptr, int *dim1, double *val)
cdef extern void adiabatclimate_make_column_p_guess_set(AdiabatClimate *ptr, int *dim1, double *val)

cdef extern void adiabatclimate_solve_for_t_trop_get(AdiabatClimate *ptr, bool *val)
cdef extern void adiabatclimate_solve_for_t_trop_set(AdiabatClimate *ptr, bool *val)

cdef extern void adiabatclimate_albedo_fcn_set(AdiabatClimate *ptr, temp_dependent_albedo_fcn fcn)

cdef extern void adiabatclimate_rh_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_rh_get(AdiabatClimate *ptr, int *dim1, double *val)
cdef extern void adiabatclimate_rh_set(AdiabatClimate *ptr, int *dim1, double *val)

cdef extern void adiabatclimate_ocean_args_p_set(AdiabatClimate *ptr, void *p)

cdef extern void adiabatclimate_tidally_locked_dayside_get(AdiabatClimate *ptr, bool *val)
cdef extern void adiabatclimate_tidally_locked_dayside_set(AdiabatClimate *ptr, bool *val)

cdef extern void adiabatclimate_l_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_l_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_chi_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_chi_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_n_lw_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_n_lw_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_cd_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_cd_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_surface_heat_flow_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_surface_heat_flow_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_species_names_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_species_names_get(AdiabatClimate *ptr, int *dim1, char* species_names)

cdef extern void adiabatclimate_particle_names_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_particle_names_get(AdiabatClimate *ptr, int *dim1, char* particle_names)

cdef extern void adiabatclimate_rad_get(AdiabatClimate *ptr, rad_pxd.Radtran **ptr1)

cdef extern void adiabatclimate_convecting_with_below_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_convecting_with_below_get(AdiabatClimate *ptr, int *dim1, bool *arr)

cdef extern void adiabatclimate_lapse_rate_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_lapse_rate_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_lapse_rate_intended_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_lapse_rate_intended_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_convective_newton_step_size_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_convective_newton_step_size_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_convective_max_boundary_shift_get(AdiabatClimate *ptr, int *val)
cdef extern void adiabatclimate_convective_max_boundary_shift_set(AdiabatClimate *ptr, int *val)

cdef extern void adiabatclimate_prevent_overconvection_get(AdiabatClimate *ptr, bint *val)
cdef extern void adiabatclimate_prevent_overconvection_set(AdiabatClimate *ptr, bint *val)
cdef extern void adiabatclimate_require_mode2_get(AdiabatClimate *ptr, bint *val)
cdef extern void adiabatclimate_require_mode2_set(AdiabatClimate *ptr, bint *val)

cdef extern void adiabatclimate_rtol_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_rtol_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_atol_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_atol_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_tol_make_column_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_tol_make_column_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_epsj_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_epsj_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_xtol_rc_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_xtol_rc_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_dt_increment_get(AdiabatClimate *ptr, double *val)
cdef extern void adiabatclimate_dt_increment_set(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_rce_solve_strategy_get(AdiabatClimate *ptr, int *val)
cdef extern void adiabatclimate_rce_solve_strategy_set(AdiabatClimate *ptr, int *val)

cdef extern void adiabatclimate_max_rc_iters_get(AdiabatClimate *ptr, int *val)
cdef extern void adiabatclimate_max_rc_iters_set(AdiabatClimate *ptr, int *val)

cdef extern void adiabatclimate_max_rc_iters_convection_get(AdiabatClimate *ptr, int *val)
cdef extern void adiabatclimate_max_rc_iters_convection_set(AdiabatClimate *ptr, int *val)

cdef extern void adiabatclimate_verbose_get(AdiabatClimate *ptr, bool *val)
cdef extern void adiabatclimate_verbose_set(AdiabatClimate *ptr, bool *val)

cdef extern void adiabatclimate_f_i_surf_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_f_i_surf_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_p_surf_get(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_p_trop_get(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_p_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_p_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_t_surf_get(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_t_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_t_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_f_i_get_size(AdiabatClimate *ptr, int *dim1, int *dim2)
cdef extern void adiabatclimate_f_i_get(AdiabatClimate *ptr, int *dim1, int *dim2, double *arr)

cdef extern void adiabatclimate_z_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_z_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_dz_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_dz_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_gravity_surf_get(AdiabatClimate *ptr, double *val)

cdef extern void adiabatclimate_gravity_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_gravity_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_densities_get_size(AdiabatClimate *ptr, int *dim1, int *dim2)
cdef extern void adiabatclimate_densities_get(AdiabatClimate *ptr, int *dim1, int *dim2, double *arr)

cdef extern void adiabatclimate_pdensities_get_size(AdiabatClimate *ptr, int *dim1, int *dim2)
cdef extern void adiabatclimate_pdensities_get(AdiabatClimate *ptr, int *dim1, int *dim2, double *arr)

cdef extern void adiabatclimate_pradii_get_size(AdiabatClimate *ptr, int *dim1, int *dim2)
cdef extern void adiabatclimate_pradii_get(AdiabatClimate *ptr, int *dim1, int *dim2, double *arr)

cdef extern void adiabatclimate_n_atmos_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_n_atmos_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_n_surface_get_size(AdiabatClimate *ptr, int *dim1)
cdef extern void adiabatclimate_n_surface_get(AdiabatClimate *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_n_ocean_get_size(AdiabatClimate *ptr, int *dim1, int *dim2)
cdef extern void adiabatclimate_n_ocean_get(AdiabatClimate *ptr, int *dim1, int *dim2, double *arr)
