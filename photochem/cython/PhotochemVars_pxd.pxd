cimport AtomConservation_pxd as atom_pxd
from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct PhotochemVars:
    pass

# callback signatures
ctypedef void (*time_dependent_flux_fcn)(double tn, int nw, double *photon_flux)
ctypedef double (*binary_diffusion_fcn)(double mu_i, double mubar, double T)

cdef extern void photochemvars_nz_get(PhotochemVars *ptr, int *nz)

cdef extern void photochemvars_top_atmos_get(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_bottom_atmos_get(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_at_photo_equilibrium_get(PhotochemVars *ptr, bool *at_photo_equilibrium)

cdef extern void photochemvars_usol_init_get_size(PhotochemVars *ptr, int *dim1, int *dim2)
cdef extern void photochemvars_usol_init_get(PhotochemVars *ptr, int *dim1, int *dim2, double *usol_init)

cdef extern void photochemvars_particle_radius_get_size(PhotochemVars *ptr, int *dim1, int *dim2)
cdef extern void photochemvars_particle_radius_get(PhotochemVars *ptr, int *dim1, int *dim2, double *arr)
cdef extern void photochemvars_particle_radius_set(PhotochemVars *ptr, int *dim1, int *dim2, double *arr)

cdef extern void photochemvars_diurnal_fac_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_diurnal_fac_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_trop_alt_get(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_trop_ind_get(PhotochemVars *ptr, int *val)

cdef extern void photochemvars_relative_humidity_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_relative_humidity_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_h2o_cond_params_get(PhotochemVars *ptr, atom_pxd.CondensationParameters **ptr1)

cdef extern void photochemvars_photon_flux_fcn_set(PhotochemVars *ptr, time_dependent_flux_fcn fcn)

cdef extern void photochemvars_cond_params_get_size(PhotochemVars *ptr, int *dim1)
cdef extern void photochemvars_cond_params_get(PhotochemVars *ptr, int *dim1, atom_pxd.CondensationParameters **ptr1)

cdef extern void photochemvars_temperature_get_size(PhotochemVars *ptr, int *dim1)
cdef extern void photochemvars_temperature_get(PhotochemVars *ptr, int *dim1, double *temperature)

cdef extern void photochemvars_edd_get_size(PhotochemVars *ptr, int *dim1)
cdef extern void photochemvars_edd_get(PhotochemVars *ptr, int *dim1, double *arr)
cdef extern void photochemvars_edd_set(PhotochemVars *ptr, int *dim1, double *arr)

cdef extern void photochemvars_custom_binary_diffusion_fcn_set(PhotochemVars *ptr, binary_diffusion_fcn fcn)

cdef extern void photochemvars_photon_flux_get_size(PhotochemVars *ptr, int *dim1)
cdef extern void photochemvars_photon_flux_get(PhotochemVars *ptr, int *dim1, double *arr)

cdef extern void photochemvars_grav_get_size(PhotochemVars *ptr, int *dim1)
cdef extern void photochemvars_grav_get(PhotochemVars *ptr, int *dim1, double *arr)

cdef extern void photochemvars_z_get_size(PhotochemVars *ptr, int *dim1)
cdef extern void photochemvars_z_get(PhotochemVars *ptr, int *dim1, double *z)

cdef extern void photochemvars_surface_pressure_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_surface_pressure_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_tauc_get_size(PhotochemVars *ptr, int *dim1, int *dim2)
cdef extern void photochemvars_tauc_get(PhotochemVars *ptr, int *dim1, int *dim2, double *val)
cdef extern void photochemvars_tauc_set(PhotochemVars *ptr, int *dim1, int *dim2, double *val)

cdef extern void photochemvars_w0c_get_size(PhotochemVars *ptr, int *dim1, int *dim2)
cdef extern void photochemvars_w0c_get(PhotochemVars *ptr, int *dim1, int *dim2, double *val)
cdef extern void photochemvars_w0c_set(PhotochemVars *ptr, int *dim1, int *dim2, double *val)

cdef extern void photochemvars_g0c_get_size(PhotochemVars *ptr, int *dim1, int *dim2)
cdef extern void photochemvars_g0c_get(PhotochemVars *ptr, int *dim1, int *dim2, double *val)
cdef extern void photochemvars_g0c_set(PhotochemVars *ptr, int *dim1, int *dim2, double *val)

cdef extern void photochemvars_max_error_reinit_attempts_get(PhotochemVars *ptr, int *val)
cdef extern void photochemvars_max_error_reinit_attempts_set(PhotochemVars *ptr, int *val)

cdef extern void photochemvars_rtol_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_rtol_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_atol_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_atol_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_mxsteps_get(PhotochemVars *ptr, int *val)
cdef extern void photochemvars_mxsteps_set(PhotochemVars *ptr, int *val)

cdef extern void photochemvars_equilibrium_time_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_equilibrium_time_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_conv_hist_factor_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_conv_hist_factor_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_conv_min_mix_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_conv_min_mix_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_conv_longdy_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_conv_longdy_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_conv_longdydt_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_conv_longdydt_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_max_dt_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_max_dt_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_autodiff_get(PhotochemVars *ptr, bool *val)
cdef extern void photochemvars_autodiff_set(PhotochemVars *ptr, bool *val)

cdef extern void photochemvars_epsj_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_epsj_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_verbose_get(PhotochemVars *ptr, int *val)
cdef extern void photochemvars_verbose_set(PhotochemVars *ptr, int *val)

cdef extern void photochemvars_fast_arbitrary_rate_get(PhotochemVars *ptr, double *val)
cdef extern void photochemvars_fast_arbitrary_rate_set(PhotochemVars *ptr, double *val)

cdef extern void photochemvars_upwind_molec_diff_get(PhotochemVars *ptr, bool *val)
cdef extern void photochemvars_upwind_molec_diff_set(PhotochemVars *ptr, bool *val)

cdef extern void photochemvars_nerrors_before_giveup_get(PhotochemVars *ptr, int *val)
cdef extern void photochemvars_nerrors_before_giveup_set(PhotochemVars *ptr, int *val)

cdef extern void photochemvars_nsteps_before_conv_check_get(PhotochemVars *ptr, int *val)
cdef extern void photochemvars_nsteps_before_conv_check_set(PhotochemVars *ptr, int *val)

cdef extern void photochemvars_nsteps_before_reinit_get(PhotochemVars *ptr, int *val)
cdef extern void photochemvars_nsteps_before_reinit_set(PhotochemVars *ptr, int *val)

cdef extern void photochemvars_nsteps_before_giveup_get(PhotochemVars *ptr, int *val)
cdef extern void photochemvars_nsteps_before_giveup_set(PhotochemVars *ptr, int *val)