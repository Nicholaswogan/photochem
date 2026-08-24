! AdiabatClimate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! allocator and destroyer !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function allocate_adiabatclimate() result(ptr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr) :: ptr
  type(AdiabatClimate), pointer :: pc
  allocate(pc)
  ptr = c_loc(pc)
end function

subroutine deallocate_adiabatclimate(ptr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  type(AdiabatClimate), pointer :: pc
  character(:), allocatable :: err_f
  call c_f_pointer(ptr, pc)
  deallocate(pc)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine wrappers  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adiabatclimate_create_wrapper(ptr, species_file, &
                                  settings_file, flux_file, data_dir, double_radiative_grid, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  character(kind=c_char), intent(in) :: species_file(*)
  character(kind=c_char), intent(in) :: settings_file(*)
  character(kind=c_char), intent(in) :: flux_file(*)
  character(kind=c_char), intent(in) :: data_dir(*)
  logical(c_bool), optional, intent(in) :: double_radiative_grid
  character(kind=c_char), intent(out) :: err(err_len+1)
  
  character(len=:), allocatable :: species_file_f
  character(len=:), allocatable :: settings_file_f
  character(len=:), allocatable :: flux_file_f
  character(len=:), allocatable :: data_dir_f
  logical :: double_radiative_grid_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c
  
  call c_f_pointer(ptr, c)
  
  allocate(character(len=len_cstring(species_file))::species_file_f)
  allocate(character(len=len_cstring(settings_file))::settings_file_f)
  allocate(character(len=len_cstring(flux_file))::flux_file_f)
  allocate(character(len=len_cstring(data_dir))::data_dir_f)
  double_radiative_grid_f = double_radiative_grid
  
  call copy_string_ctof(species_file, species_file_f)
  call copy_string_ctof(settings_file, settings_file_f)
  call copy_string_ctof(flux_file, flux_file_f)
  call copy_string_ctof(data_dir, data_dir_f)
  
  c = AdiabatClimate(species_file_f, &
                     settings_file_f, &
                     flux_file_f, &
                     data_dir_f, &
                     double_radiative_grid_f, &
                     err_f)
  
  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif
end subroutine

subroutine adiabatclimate_make_profile_wrapper(ptr, T_surf, ng, P_i_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%make_profile(T_surf, P_i_surf, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_make_column_wrapper(ptr, T_surf, ng, N_i_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: N_i_surf(ng)
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%make_column(T_surf, N_i_surf, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_make_profile_bg_gas_wrapper(ptr, T_surf, ng, P_i_surf, P_surf, bg_gas, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: P_surf
  character(kind=c_char), intent(in) :: bg_gas(*)
  character(c_char), intent(out) :: err(err_len+1)

  character(len=:), allocatable :: bg_gas_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  allocate(character(len=len_cstring(bg_gas))::bg_gas_f)
  call copy_string_ctof(bg_gas, bg_gas_f)

  call c%make_profile_bg_gas(T_surf, P_i_surf, P_surf, bg_gas_f, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_make_profile_dry_wrapper(ptr, dim_P, P, dim_T, T, dim1_f_i, dim2_f_i, f_i, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim_P
  real(c_double), intent(in) :: P(dim_P)
  integer(c_int), intent(in) :: dim_T
  real(c_double), intent(in) :: T(dim_T)
  integer(c_int), intent(in) :: dim1_f_i, dim2_f_i
  real(c_double), intent(in) :: f_i(dim1_f_i,dim2_f_i)
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)
  call c%make_profile_dry(P, T, f_i, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_toa_fluxes_wrapper(ptr, T_surf, ng, P_i_surf, ISR, OLR, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(out) :: ISR, OLR
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%TOA_fluxes(T_surf, P_i_surf, ISR, OLR, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_toa_fluxes_column_wrapper(ptr, T_surf, ng, N_i_surf, ISR, OLR, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: N_i_surf(ng)
  real(c_double), intent(out) :: ISR, OLR
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%TOA_fluxes_column(T_surf, N_i_surf, ISR, OLR, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_toa_fluxes_bg_gas_wrapper(ptr, T_surf, ng, P_i_surf, P_surf, bg_gas, ISR, OLR, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: P_surf
  character(kind=c_char), intent(in) :: bg_gas(*)
  real(c_double), intent(out) :: ISR, OLR
  character(c_char), intent(out) :: err(err_len+1)

  character(len=:), allocatable :: bg_gas_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  allocate(character(len=len_cstring(bg_gas))::bg_gas_f)
  call copy_string_ctof(bg_gas, bg_gas_f)

  call c%TOA_fluxes_bg_gas(T_surf, P_i_surf, P_surf, bg_gas_f, ISR, OLR, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_toa_fluxes_dry_wrapper(ptr, dim_P, P, dim_T, T, dim1_f_i, dim2_f_i, f_i, ISR, OLR, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim_P
  real(c_double), intent(in) :: P(dim_P)
  integer(c_int), intent(in) :: dim_T
  real(c_double), intent(in) :: T(dim_T)
  integer(c_int), intent(in) :: dim1_f_i, dim2_f_i
  real(c_double), intent(in) :: f_i(dim1_f_i,dim2_f_i)
  real(c_double), intent(out) :: ISR, OLR
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)
  call c%TOA_fluxes_dry(P, T, f_i, ISR, OLR, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_surface_temperature_wrapper(ptr, ng, P_i_surf, T_guess, T_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: T_guess
  real(c_double), intent(out) :: T_surf
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  T_surf = c%surface_temperature(P_i_surf, T_guess, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_surface_temperature_column_wrapper(ptr, ng, N_i_surf, T_guess, T_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: N_i_surf(ng)
  real(c_double), intent(in) :: T_guess
  real(c_double), intent(out) :: T_surf
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  T_surf = c%surface_temperature_column(N_i_surf, T_guess, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_surface_temperature_bg_gas_wrapper(ptr, ng, P_i_surf, P_surf, bg_gas, T_guess, T_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: P_surf
  character(kind=c_char), intent(in) :: bg_gas(*)
  real(c_double), intent(in) :: T_guess
  real(c_double), intent(out) :: T_surf
  character(c_char), intent(out) :: err(err_len+1)

  character(len=:), allocatable :: bg_gas_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  allocate(character(len=len_cstring(bg_gas))::bg_gas_f)
  call copy_string_ctof(bg_gas, bg_gas_f)

  T_surf = c%surface_temperature_bg_gas(P_i_surf, P_surf, bg_gas_f, T_guess, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_set_particle_density_and_radii(ptr, dim_P, P,  dim1_pdensities, dim2_pdensities, pdensities, &
                                                         dim1_pradii, dim2_pradii, pradii, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim_P
  real(c_double), intent(in) :: P(dim_P)
  integer(c_int), intent(in) :: dim1_pdensities, dim2_pdensities
  real(c_double), intent(in) :: pdensities(dim1_pdensities, dim2_pdensities)
  integer(c_int), intent(in) :: dim1_pradii, dim2_pradii
  real(c_double), intent(in) :: pradii(dim1_pradii, dim2_pradii)
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)
  call c%set_particle_density_and_radii(P, pdensities, pradii, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_rce_wrapper(ptr, ng, P_i_surf, T_surf_guess, dim_T_guess, T_guess, &
                                      convecting_with_below_present, dim_convecting_with_below, convecting_with_below, &
                                      custom_present, dim_sp_custom, sp_custom, dim_P_custom, P_custom, &
                                      dim1_mix_custom, dim2_mix_custom, mix_custom, &
                                      converged, err) bind(c)
  use clima, only: AdiabatClimate, s_str_len
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: T_surf_guess
  integer(c_int), intent(in) :: dim_T_guess
  real(c_double), intent(in) :: T_guess(dim_T_guess)
  logical(c_bool), intent(in) :: convecting_with_below_present
  integer(c_int), intent(in) :: dim_convecting_with_below
  logical(c_bool), intent(in) :: convecting_with_below(dim_convecting_with_below)
  logical(c_bool), intent(in) :: custom_present
  integer(c_int), intent(in) :: dim_sp_custom
  character(c_char), intent(out) :: sp_custom(dim_sp_custom*s_str_len+1)
  integer(c_int), intent(in) :: dim_P_custom
  real(c_double), intent(in) :: P_custom(dim_P_custom)
  integer(c_int), intent(in) :: dim1_mix_custom, dim2_mix_custom
  real(c_double), intent(in) :: mix_custom(dim1_mix_custom,dim2_mix_custom)
  logical(c_bool), intent(out) :: converged
  character(c_char), intent(out) :: err(err_len+1)

  character(s_str_len), allocatable :: sp_custom_f(:)
  integer :: i, j, k
  logical, allocatable :: convecting_with_below_f(:)
  logical :: converged_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  allocate(convecting_with_below_f(dim_convecting_with_below))
  convecting_with_below_f = convecting_with_below

  allocate(sp_custom_f(dim_sp_custom))
  do i = 1,dim_sp_custom
    do j = 1,s_str_len
      k = j + (i - 1) * s_str_len
      sp_custom_f(i)(j:j) = sp_custom(k)
    enddo
  enddo
  
  converged_f = .false.
  if (convecting_with_below_present .and. custom_present) then
    converged_f = c%RCE(P_i_surf, T_surf_guess, T_guess, convecting_with_below_f, sp_custom_f, P_custom, mix_custom, err_f)
  elseif (convecting_with_below_present .and. .not.custom_present) then
    converged_f = c%RCE(P_i_surf, T_surf_guess, T_guess, convecting_with_below_f, err=err_f)
  elseif (.not.convecting_with_below_present .and. custom_present) then
    converged_f = c%RCE(P_i_surf, T_surf_guess, T_guess, sp_custom=sp_custom_f, P_custom=P_custom, mix_custom=mix_custom, err=err_f)
  else
    converged_f = c%RCE(P_i_surf, T_surf_guess, T_guess, err=err_f)
  endif
  converged = converged_f
  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_set_ocean_solubility_fcn_wrapper(ptr, species_c, fcn_c, err) bind(c)
  use clima, only: AdiabatClimate
  use clima, only: ocean_solubility_fcn
  type(c_ptr), value, intent(in) :: ptr
  character(kind=c_char), intent(in) :: species_c(*)
  type(c_funptr), value, intent(in) :: fcn_c
  character(kind=c_char), intent(out) :: err(err_len+1)

  type(AdiabatClimate), pointer :: c
  procedure(ocean_solubility_fcn), pointer :: fcn_f
  character(:), allocatable :: species_f
  character(:), allocatable :: err_f

  ! Cast the void pointer to wrapped object.
  call c_f_pointer(ptr, c)

  ! Copy c string to f string.
  allocate(character(len=len_cstring(species_c))::species_f)
  call copy_string_ctof(species_c, species_f)

  ! Convert c function pointer to a fortran function pointer.
  call c_f_procpointer(fcn_c, fcn_f)

  ! Call the function.
  call c%set_ocean_solubility_fcn(species_f, fcn_f, err_f)

  ! Set error, if there is one.
  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif 
end subroutine

subroutine adiabatclimate_to_regular_grid_wrapper(ptr, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%to_regular_grid(err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_out2atmosphere_txt_wrapper(ptr, filename, nz, eddy, number_of_decimals, overwrite, clip, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  character(kind=c_char), intent(in) :: filename(*)
  integer(c_int), intent(in) :: nz
  real(c_double), intent(in) :: eddy(nz)
  integer(c_int), intent(in) :: number_of_decimals
  logical(c_bool), intent(in) :: overwrite, clip
  character(kind=c_char), intent(out) :: err(err_len+1)
  
  character(len=:), allocatable :: filename_f
  logical :: overwrite_f, clip_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c
  
  call c_f_pointer(ptr, c)
  
  allocate(character(len=len_cstring(filename))::filename_f)
  call copy_string_ctof(filename, filename_f)
  overwrite_f = overwrite
  clip_f = clip
  
  call c%out2atmosphere_txt(filename_f, eddy, number_of_decimals, overwrite_f, clip_f, err_f)
  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif
  
end subroutine

subroutine adiabatclimate_heat_redistribution_parameters_wrapper(ptr, tau_LW, k_term, f_term, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: tau_LW
  real(c_double), intent(out) :: k_term
  real(c_double), intent(out) :: f_term
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%heat_redistribution_parameters(tau_LW, k_term, f_term, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! getters and setters !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adiabatclimate_p_top_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%P_top
end subroutine

subroutine adiabatclimate_p_top_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%P_top = val
end subroutine

subroutine adiabatclimate_reference_pressure_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%reference_pressure
end subroutine

subroutine adiabatclimate_reference_pressure_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%reference_pressure = val
end subroutine

subroutine adiabatclimate_t_trop_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%T_trop
end subroutine

subroutine adiabatclimate_t_trop_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%T_trop = val
end subroutine

subroutine adiabatclimate_use_make_column_p_guess_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%use_make_column_P_guess
end subroutine

subroutine adiabatclimate_use_make_column_p_guess_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%use_make_column_P_guess = val
end subroutine

subroutine adiabatclimate_make_column_p_guess_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%make_column_P_guess,1)
end subroutine

subroutine adiabatclimate_make_column_p_guess_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%make_column_P_guess
end subroutine

subroutine adiabatclimate_make_column_p_guess_set(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(in) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%make_column_P_guess = arr
end subroutine

subroutine adiabatclimate_solve_for_t_trop_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%solve_for_T_trop
end subroutine

subroutine adiabatclimate_solve_for_t_trop_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%solve_for_T_trop = val
end subroutine

subroutine adiabatclimate_albedo_fcn_set(ptr, albedo_fcn_c) bind(c)
  use clima, only: AdiabatClimate
  use clima, only: temp_dependent_albedo_fcn
  type(c_ptr), value, intent(in) :: ptr
  type(c_funptr), value, intent(in) :: albedo_fcn_c

  procedure(temp_dependent_albedo_fcn), pointer :: albedo_fcn_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)
  call c_f_procpointer(albedo_fcn_c, albedo_fcn_f)
  c%albedo_fcn => albedo_fcn_f

end subroutine

subroutine adiabatclimate_rh_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%RH,1)
end subroutine

subroutine adiabatclimate_rh_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%RH
end subroutine

subroutine adiabatclimate_rh_set(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(in) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%RH = arr
end subroutine

subroutine adiabatclimate_ocean_args_p_set(ptr, p) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  type(c_ptr), value, intent(in) :: p
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%ocean_args_p = p
end subroutine

subroutine adiabatclimate_tidally_locked_dayside_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%tidally_locked_dayside
end subroutine

subroutine adiabatclimate_tidally_locked_dayside_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%tidally_locked_dayside = val
end subroutine

subroutine adiabatclimate_l_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%L
end subroutine

subroutine adiabatclimate_l_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%L = val
end subroutine

subroutine adiabatclimate_chi_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%chi
end subroutine

subroutine adiabatclimate_chi_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%chi = val
end subroutine

subroutine adiabatclimate_n_lw_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%n_LW
end subroutine

subroutine adiabatclimate_n_lw_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%n_LW = val
end subroutine

subroutine adiabatclimate_cd_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%Cd
end subroutine

subroutine adiabatclimate_cd_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%Cd = val
end subroutine

subroutine adiabatclimate_surface_heat_flow_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%surface_heat_flow
end subroutine

subroutine adiabatclimate_surface_heat_flow_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%surface_heat_flow = val
end subroutine

subroutine adiabatclimate_species_names_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%species_names)
end subroutine

subroutine adiabatclimate_species_names_get(ptr, dim1, species_names) bind(c)
  use clima, only: AdiabatClimate, s_str_len
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  character(kind=c_char), intent(out) :: species_names(dim1*s_str_len+1)
  type(AdiabatClimate), pointer :: c
  
  integer :: i, j, k
  
  call c_f_pointer(ptr, c)
  do i = 1,dim1
    do j = 1,s_str_len
      k = j + (i - 1) * s_str_len
      species_names(k) = c%species_names(i)(j:j)
    enddo
  enddo
  species_names(dim1*s_str_len+1) = c_null_char
  
end subroutine

subroutine adiabatclimate_particle_names_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%particle_names)
end subroutine

subroutine adiabatclimate_particle_names_get(ptr, dim1, particle_names) bind(c)
  use clima, only: AdiabatClimate, s_str_len
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  character(kind=c_char), intent(out) :: particle_names(dim1*s_str_len+1)
  type(AdiabatClimate), pointer :: c
  
  integer :: i, j, k
  
  call c_f_pointer(ptr, c)
  do i = 1,dim1
    do j = 1,s_str_len
      k = j + (i - 1) * s_str_len
      particle_names(k) = c%particle_names(i)(j:j)
    enddo
  enddo
  particle_names(dim1*s_str_len+1) = c_null_char
  
end subroutine

subroutine adiabatclimate_rad_get(ptr, ptr1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  ptr1 = c_loc(c%rad)
end subroutine

subroutine adiabatclimate_convecting_with_below_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%convecting_with_below,1)
end subroutine

subroutine adiabatclimate_convecting_with_below_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  logical(c_bool), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%convecting_with_below
end subroutine

subroutine adiabatclimate_lapse_rate_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%lapse_rate,1)
end subroutine

subroutine adiabatclimate_lapse_rate_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%lapse_rate
end subroutine

subroutine adiabatclimate_lapse_rate_intended_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%lapse_rate_intended,1)
end subroutine

subroutine adiabatclimate_lapse_rate_intended_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%lapse_rate_intended
end subroutine

subroutine adiabatclimate_convective_newton_step_size_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%convective_newton_step_size
end subroutine

subroutine adiabatclimate_convective_newton_step_size_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%convective_newton_step_size = val
end subroutine

subroutine adiabatclimate_convective_max_boundary_shift_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%convective_max_boundary_shift
end subroutine

subroutine adiabatclimate_convective_max_boundary_shift_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%convective_max_boundary_shift = val
end subroutine

subroutine adiabatclimate_prevent_overconvection_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%prevent_overconvection
end subroutine

subroutine adiabatclimate_prevent_overconvection_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%prevent_overconvection = val
end subroutine

subroutine adiabatclimate_require_mode2_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%require_mode2
end subroutine

subroutine adiabatclimate_require_mode2_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%require_mode2 = val
end subroutine

subroutine adiabatclimate_rtol_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%rtol
end subroutine

subroutine adiabatclimate_rtol_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%rtol = val
end subroutine

subroutine adiabatclimate_atol_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%atol
end subroutine

subroutine adiabatclimate_atol_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%atol = val
end subroutine

subroutine adiabatclimate_tol_make_column_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%tol_make_column
end subroutine

subroutine adiabatclimate_tol_make_column_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%tol_make_column = val
end subroutine

subroutine adiabatclimate_epsj_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%epsj
end subroutine

subroutine adiabatclimate_epsj_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%epsj = val
end subroutine

subroutine adiabatclimate_xtol_rc_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%xtol_rc
end subroutine

subroutine adiabatclimate_xtol_rc_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%xtol_rc = val
end subroutine

subroutine adiabatclimate_dt_increment_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%dt_increment
end subroutine

subroutine adiabatclimate_dt_increment_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%dt_increment = val
end subroutine

subroutine adiabatclimate_rce_solve_strategy_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%rce_solve_strategy
end subroutine

subroutine adiabatclimate_rce_solve_strategy_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%rce_solve_strategy = val
end subroutine

subroutine adiabatclimate_max_rc_iters_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%max_rc_iters
end subroutine

subroutine adiabatclimate_max_rc_iters_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%max_rc_iters = val
end subroutine

subroutine adiabatclimate_max_rc_iters_convection_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%max_rc_iters_convection
end subroutine

subroutine adiabatclimate_max_rc_iters_convection_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%max_rc_iters_convection = val
end subroutine

subroutine adiabatclimate_verbose_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%verbose
end subroutine

subroutine adiabatclimate_verbose_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  logical(c_bool), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%verbose = val
end subroutine

subroutine adiabatclimate_f_i_surf_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%f_i_surf,1)
end subroutine

subroutine adiabatclimate_f_i_surf_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%f_i_surf
end subroutine

subroutine adiabatclimate_p_surf_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%P_surf
end subroutine

subroutine adiabatclimate_p_trop_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%P_trop
end subroutine

subroutine adiabatclimate_p_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%P,1)
end subroutine

subroutine adiabatclimate_p_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%P
end subroutine

subroutine adiabatclimate_t_surf_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%T_surf
end subroutine

subroutine adiabatclimate_t_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%T,1)
end subroutine

subroutine adiabatclimate_t_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%T
end subroutine

subroutine adiabatclimate_f_i_get_size(ptr, dim1, dim2) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%f_i,1)
  dim2 = size(c%f_i,2)
end subroutine

subroutine adiabatclimate_f_i_get(ptr, dim1, dim2, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%f_i
end subroutine

subroutine adiabatclimate_z_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%z,1)
end subroutine

subroutine adiabatclimate_z_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%z
end subroutine

subroutine adiabatclimate_dz_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%dz,1)
end subroutine

subroutine adiabatclimate_dz_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%dz
end subroutine

subroutine adiabatclimate_gravity_surf_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%gravity_surf
end subroutine

subroutine adiabatclimate_gravity_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%gravity,1)
end subroutine

subroutine adiabatclimate_gravity_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%gravity
end subroutine

subroutine adiabatclimate_densities_get_size(ptr, dim1, dim2) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%densities,1)
  dim2 = size(c%densities,2)
end subroutine

subroutine adiabatclimate_densities_get(ptr, dim1, dim2, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%densities
end subroutine

subroutine adiabatclimate_pdensities_get_size(ptr, dim1, dim2) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%pdensities,1)
  dim2 = size(c%pdensities,2)
end subroutine

subroutine adiabatclimate_pdensities_get(ptr, dim1, dim2, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%pdensities
end subroutine

subroutine adiabatclimate_pradii_get_size(ptr, dim1, dim2) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%pradii,1)
  dim2 = size(c%pradii,2)
end subroutine

subroutine adiabatclimate_pradii_get(ptr, dim1, dim2, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%pradii
end subroutine

subroutine adiabatclimate_n_atmos_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%N_atmos,1)
end subroutine

subroutine adiabatclimate_n_atmos_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%N_atmos
end subroutine

subroutine adiabatclimate_n_surface_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%N_surface,1)
end subroutine

subroutine adiabatclimate_n_surface_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%N_surface
end subroutine

subroutine adiabatclimate_n_ocean_get_size(ptr, dim1, dim2) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%N_ocean,1)
  dim2 = size(c%N_ocean,2)
end subroutine

subroutine adiabatclimate_n_ocean_get(ptr, dim1, dim2, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%N_ocean
end subroutine
