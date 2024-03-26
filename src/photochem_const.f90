

module photochem_const
  use, intrinsic :: iso_fortran_env, only : dp => real64
  implicit none
  ! public
  integer, parameter :: err_len = 1024
  integer, parameter :: str_len = 1024
  integer, parameter :: m_str_len = 100
  integer, parameter :: s_str_len = 20
  real(dp), parameter :: Rgas = 8.31446261815324e0_dp ! ideal gas constant (j/(mol*K))
  real(dp), parameter :: k_boltz = 1.380649e-16_dp ! boltzmann's constant cgs units (egs/K)
  real(dp), parameter :: G_grav = 6.67430e-11_dp ! gravitational constant (N * m2 / kg)
  real(dp), parameter :: plank = 6.62607004e-34_dp ! planks constant (m2 kg / s)
  real(dp), parameter :: c_light = 299792458.0_dp ! Speed of light (m / s)
  real(dp), parameter :: N_avo = 6.02214076e23_dp ! avagadros number
  real(dp), parameter :: T_crit_H2O = 647.0_dp ! critical point H2O (K)
  real(dp), parameter :: pi = 3.14159265358979323846e0_dp

  real(dp), parameter :: smallest_real = tiny(1.0_dp)
  real(dp), parameter :: smaller_real = tiny(1.0_dp)**0.5_dp
  real(dp), parameter :: small_real = tiny(1.0_dp)**0.25_dp
  real(dp), parameter :: ln_small_real = log(small_real)

  !> number of steps to save during integration, for steady-state checking
  integer, parameter :: nsteps_save = 500

end module
