
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine photochemwrk_nsteps_total_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    val = wrk%nsteps_total
  end subroutine

  subroutine photochemwrk_nsteps_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: val
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    val = wrk%nsteps
  end subroutine

  subroutine photochemwrk_t_history_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%t_history,1)
  end subroutine
  
  subroutine photochemwrk_t_history_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%t_history
  end subroutine

  subroutine photochemwrk_mix_history_get_size(ptr, dim1, dim2, dim3) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2, dim3
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%mix_history,1)
    dim2 = size(wrk%mix_history,2)
    dim3 = size(wrk%mix_history,3)
  end subroutine
  
  subroutine photochemwrk_mix_history_get(ptr, dim1, dim2, dim3, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2, dim3
    real(c_double), intent(out) :: arr(dim1, dim2, dim3)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%mix_history
  end subroutine

  subroutine photochemwrk_longdy_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    val = wrk%longdy
  end subroutine

  subroutine photochemwrk_longdydt_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    val = wrk%longdydt
  end subroutine

  subroutine photochemwrk_tn_get(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    val = wrk%tn
  end subroutine
  
  subroutine photochemwrk_tn_set(ptr, val) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    wrk%tn = val
  end subroutine
  
  subroutine photochemwrk_usol_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%usol,1)
    dim2 = size(wrk%usol,2)
  end subroutine
  
  subroutine photochemwrk_usol_get(ptr, dim1, dim2, usol) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: usol(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    usol = wrk%usol
  end subroutine
  
  subroutine photochemwrk_usol_set(ptr, dim1, dim2, usol) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(in) :: usol(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    wrk%usol = usol
  end subroutine
  
  subroutine photochemwrk_pressure_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%pressure,1)
  end subroutine
  
  subroutine photochemwrk_pressure_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%pressure
  end subroutine
  
  subroutine photochemwrk_density_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%density,1)
  end subroutine
  
  subroutine photochemwrk_density_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%density
  end subroutine
  
  subroutine photochemwrk_densities_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%densities,1)
    dim2 = size(wrk%densities,2)
  end subroutine
  
  subroutine photochemwrk_densities_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%densities
  end subroutine

  subroutine photochemwrk_rx_rates_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%rx_rates,1)
    dim2 = size(wrk%rx_rates,2)
  end subroutine
  
  subroutine photochemwrk_rx_rates_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%rx_rates
  end subroutine

  subroutine photochemwrk_mubar_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%mubar,1)
  end subroutine
  
  subroutine photochemwrk_mubar_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%mubar
  end subroutine
  
  subroutine photochemwrk_prates_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%prates,1)
    dim2 = size(wrk%prates,2)
  end subroutine
  
  subroutine photochemwrk_prates_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%prates
  end subroutine
  
  subroutine photochemwrk_amean_grd_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%amean_grd,1)
    dim2 = size(wrk%amean_grd,2)
  end subroutine
  
  subroutine photochemwrk_amean_grd_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%amean_grd
  end subroutine
  
  subroutine photochemwrk_optical_depth_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%optical_depth,1)
    dim2 = size(wrk%optical_depth,2)
  end subroutine
  
  subroutine photochemwrk_optical_depth_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%optical_depth
  end subroutine
  
  subroutine photochemwrk_surf_radiance_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%surf_radiance,1)
  end subroutine
  
  subroutine photochemwrk_surf_radiance_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrk), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%surf_radiance
  end subroutine

! PhotochemWrkEvo

  subroutine photochemwrkevo_pressure_hydro_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(PhotochemWrkEvo), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    dim1 = size(wrk%pressure_hydro,1)
  end subroutine
  
  subroutine photochemwrkevo_pressure_hydro_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(PhotochemWrkEvo), pointer :: wrk
    call c_f_pointer(ptr, wrk)
    arr = wrk%pressure_hydro
  end subroutine