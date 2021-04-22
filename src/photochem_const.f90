

module photochem_const
  implicit none
  public
  integer, private, parameter :: real_kind = kind(1.0d0)
  real(real_kind), parameter :: Rgas = 8.31446261815324d0 ! ideal gas constant (j/(mol*K))
  real(real_kind), parameter :: k_boltz = 1.3807d-16 ! boltzmann's constant cgs units (egs/K)
end module