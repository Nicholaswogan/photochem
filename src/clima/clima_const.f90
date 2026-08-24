module clima_const
  use iso_fortran_env, only: dp => real64
  implicit none
  public
  
  integer, parameter :: s_str_len = 20
  
  ! physical constants 
  real(dp), parameter :: Rgas = 8.31446261815324e7_dp ! ideal gas constant (erg/(mol*K))
  real(dp), parameter :: k_boltz = 1.380649e-16_dp ! boltzmann's constant cgs units (egs/K)
  real(dp), parameter :: k_boltz_si = 1.380649e-23_dp ! boltzmann's constant si units (J/K)
  real(dp), parameter :: G_grav = 6.67430e-11_dp ! gravitational constant (N * m2 / kg)
  real(dp), parameter :: plank = 6.62607004e-34_dp ! planks constant (m2 kg / s)
  real(dp), parameter :: c_light = 299792458.0_dp ! Speed of light (m / s)
  real(dp), parameter :: N_avo = 6.02214076e23_dp ! avagadros number
  real(dp), parameter :: sigma_si = 5.670374419e-8_dp ! Stefan-Boltzmann constant (W / m^2 / K^4)
  real(dp), parameter :: pi = 3.14159265358979323846e0_dp
  real(dp), parameter :: von_karman_const = 0.41_dp ! from wikipedia
  
  ! useful
  real(dp), parameter :: log10tiny = log10(sqrt(tiny(1.0_dp)))
  
  ! other
  real(dp), parameter :: max_w0 = 0.99999_dp
  
end module