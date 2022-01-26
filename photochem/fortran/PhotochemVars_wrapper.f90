module PhotochemVars_wrapper
  use photochem_const, only: s_str_len
  use photochem_types, only: PhotochemVars
  use iso_c_binding
  implicit none
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine allocate_photochemvars(ptr) bind(c)
    type(c_ptr), intent(out) :: ptr
    type(PhotochemVars), pointer :: dat
    allocate(dat)
    ptr = c_loc(dat)
  end subroutine
  
  subroutine deallocate_photochemvars(ptr) bind(c)
    type(c_ptr), intent(in) :: ptr
    type(PhotochemVars), pointer :: dat
    call c_f_pointer(ptr, dat)
    deallocate(dat)
  end subroutine
  
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
  
  subroutine photochemvars_at_photo_equilibrium_get(ptr, at_photo_equilibrium) bind(c)
    type(c_ptr), intent(in) :: ptr
    logical(4), intent(out) :: at_photo_equilibrium
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
  
end module