

module photochem_const
  implicit none
  ! public
  integer, parameter :: err_len = 1024
  integer, parameter :: str_len = 1024
  integer, parameter :: m_str_len = 100
  integer, parameter :: s_str_len = 20
  integer, parameter :: real_kind = kind(1.0d0)
  real(real_kind), parameter :: Rgas = 8.31446261815324d0 ! ideal gas constant (j/(mol*K))
  real(real_kind), parameter :: k_boltz = 1.380649d-16 ! boltzmann's constant cgs units (egs/K)
  real(real_kind), parameter :: G_grav = 6.67430d-11 ! gravitational constant (N * m2 / kg)
  real(real_kind), parameter :: plank = 6.62607004d-34 ! planks constant (m2 kg / s)
  real(real_kind), parameter :: c_light = 299792458.d0 ! Speed of light (m / s)
  real(real_kind), parameter :: N_avo = 6.02214076d23 ! avagadros number
  real(real_kind), parameter :: pi = 3.14159265358979323846d0

  real(real_kind), parameter :: smallest_real = tiny(1.d0)
  real(real_kind), parameter :: smaller_real = tiny(1.d0)**0.5d0
  real(real_kind), parameter :: small_real = tiny(1.d0)**0.25d0
  real(real_kind), parameter :: ln_small_real = log(small_real)
end module
