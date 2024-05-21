
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine photochemvars_nz_get(ptr, nz) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: nz
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    nz = var%nz
  end subroutine
  
  subroutine photochemvars_top_atmos_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%top_atmos
  end subroutine

  subroutine photochemvars_bottom_atmos_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%bottom_atmos
  end subroutine
  
  subroutine photochemvars_at_photo_equilibrium_get(ptr, at_photo_equilibrium) bind(c)
    type(c_ptr), intent(in) :: ptr
    logical(c_bool), intent(out) :: at_photo_equilibrium
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    at_photo_equilibrium = var%at_photo_equilibrium
  end subroutine
  
  subroutine photochemvars_usol_init_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    dim1 = size(var%usol_init,1)
    dim2 = size(var%usol_init,2)
  end subroutine
  
  subroutine photochemvars_usol_init_get(ptr, dim1, dim2, usol_init) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: usol_init(dim1, dim2)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    usol_init = var%usol_init
  end subroutine

  subroutine photochemvars_trop_alt_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%trop_alt
  end subroutine

  subroutine photochemvars_trop_ind_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%trop_ind
  end subroutine

  subroutine photochemvars_relative_humidity_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%relative_humidity
  end subroutine
  
  subroutine photochemvars_relative_humidity_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%relative_humidity = val
  end subroutine

  subroutine photochemvars_h2o_cond_params_get(ptr, ptr1) bind(c)
    type(c_ptr), intent(in) :: ptr
    type(c_ptr), intent(out) :: ptr1
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    ptr1 = c_loc(var%H2O_cond_params)
  end subroutine

  subroutine photochemvars_photon_flux_fcn_set(ptr, photon_flux_fcn_c) bind(c)
    use photochem_types, only: time_dependent_flux_fcn
    type(c_ptr), intent(in) :: ptr
    type(c_funptr), value, intent(in) :: photon_flux_fcn_c
  
    procedure(time_dependent_flux_fcn), pointer :: photon_flux_fcn_f
    type(PhotochemVars), pointer :: var
  
    call c_f_pointer(ptr, var)
    call c_f_procpointer(photon_flux_fcn_c, photon_flux_fcn_f)
    var%photon_flux_fcn => photon_flux_fcn_f
  
  end subroutine

  subroutine photochemvars_cond_params_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    dim1 = size(var%cond_params,1)
  end subroutine

  subroutine photochemvars_cond_params_get(ptr, dim1, ptr1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    type(c_ptr), intent(out) :: ptr1(dim1)
    integer :: i
    type(CondensationParameters), pointer :: t1_p(:)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    t1_p => var%cond_params
    do i = 1,dim1
      ptr1(i) = c_loc(t1_p(i))
    enddo
  end subroutine
  
  subroutine photochemvars_temperature_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    dim1 = size(var%temperature,1)
  end subroutine
  
  subroutine photochemvars_temperature_get(ptr, dim1, temperature) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: temperature(dim1)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    temperature = var%temperature
  end subroutine
  
  subroutine photochemvars_edd_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    dim1 = size(var%edd,1)
  end subroutine
  
  subroutine photochemvars_edd_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    arr = var%edd
  end subroutine

  subroutine photochemvars_edd_set(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(in) :: arr(dim1)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%edd = arr
  end subroutine

  subroutine photochemvars_custom_binary_diffusion_fcn_set(ptr, fcn_c) bind(c)
    use photochem_types, only: binary_diffusion_fcn
    type(c_ptr), intent(in) :: ptr
    type(c_funptr), value, intent(in) :: fcn_c
  
    procedure(binary_diffusion_fcn), pointer :: fcn_f
    type(PhotochemVars), pointer :: var
  
    call c_f_pointer(ptr, var)
    call c_f_procpointer(fcn_c, fcn_f)
    var%custom_binary_diffusion_fcn => fcn_f
  
  end subroutine
  
  subroutine photochemvars_photon_flux_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    dim1 = size(var%photon_flux,1)
  end subroutine
  
  subroutine photochemvars_photon_flux_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    arr = var%photon_flux
  end subroutine
  
  subroutine photochemvars_grav_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    dim1 = size(var%grav,1)
  end subroutine
  
  subroutine photochemvars_grav_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    arr = var%grav
  end subroutine
  
  subroutine photochemvars_z_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    dim1 = size(var%z,1)
  end subroutine
  
  subroutine photochemvars_z_get(ptr, dim1, z) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: z(dim1)
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    z = var%z
  end subroutine
  
  subroutine photochemvars_surface_pressure_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%surface_pressure
  end subroutine
  
  subroutine photochemvars_surface_pressure_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%surface_pressure = val
  end subroutine

  subroutine photochemvars_max_error_reinit_attempts_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%max_error_reinit_attempts
  end subroutine
  
  subroutine photochemvars_max_error_reinit_attempts_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%max_error_reinit_attempts = val
  end subroutine
  
  subroutine photochemvars_rtol_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%rtol
  end subroutine
  
  subroutine photochemvars_rtol_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%rtol = val
  end subroutine
  
  subroutine photochemvars_atol_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%atol
  end subroutine
  
  subroutine photochemvars_atol_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%atol = val
  end subroutine
  
  subroutine photochemvars_mxsteps_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%mxsteps
  end subroutine
  
  subroutine photochemvars_mxsteps_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%mxsteps = val
  end subroutine
  
  subroutine photochemvars_equilibrium_time_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%equilibrium_time
  end subroutine
  
  subroutine photochemvars_equilibrium_time_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%equilibrium_time = val
  end subroutine
  
  subroutine photochemvars_conv_hist_factor_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%conv_hist_factor
  end subroutine
  
  subroutine photochemvars_conv_hist_factor_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%conv_hist_factor = val
  end subroutine

  subroutine photochemvars_conv_min_mix_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%conv_min_mix
  end subroutine
  
  subroutine photochemvars_conv_min_mix_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%conv_min_mix = val
  end subroutine

  subroutine photochemvars_conv_longdy_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%conv_longdy
  end subroutine
  
  subroutine photochemvars_conv_longdy_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%conv_longdy = val
  end subroutine

  subroutine photochemvars_conv_longdydt_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%conv_longdydt
  end subroutine
  
  subroutine photochemvars_conv_longdydt_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%conv_longdydt = val
  end subroutine

  subroutine photochemvars_autodiff_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    logical(c_bool), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%autodiff
  end subroutine

  subroutine photochemvars_autodiff_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    logical(c_bool), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%autodiff = val
  end subroutine

  subroutine photochemvars_epsj_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%epsj
  end subroutine
  
  subroutine photochemvars_epsj_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%epsj = val
  end subroutine

  subroutine photochemvars_verbose_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%verbose
  end subroutine
  
  subroutine photochemvars_verbose_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%verbose = val
  end subroutine

  subroutine photochemvars_fast_arbitrary_rate_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    val = var%fast_arbitrary_rate
  end subroutine
  
  subroutine photochemvars_fast_arbitrary_rate_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemVars), pointer :: var
    call c_f_pointer(ptr, var)
    var%fast_arbitrary_rate = val
  end subroutine
  