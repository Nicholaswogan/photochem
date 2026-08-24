module clima_eqns_water
  use clima_const, only: dp
  implicit none
  private

  ! This module contains latent heat and saturation pressure of H2O
  ! fits for the following data:
  ! https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-d_599.html
  ! https://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html

  public :: latent_heat_H2O, latent_heat_H2O_vap, latent_heat_H2O_sub
  public :: sat_pressure_H2O, sat_pressure_H2O_vap, sat_pressure_H2O_sub
  public :: T_freeze, mu_H2O, Rgas

  ! Here are the exact Rgas and mu_H2O that I used to fit the data.
  real(dp), parameter :: Rgas = 8.31446261815324e7_dp ! Ideal gas constant (erg/(mol*K))
  real(dp), parameter :: mu_H2O = 18.01534_dp ! H2O molar mass (g/mol)

  ! Latent heat fit parameters (no units)
  real(dp), parameter :: A_v = -3413485157036.1396_dp
  real(dp), parameter :: B_v = 4.093669788667096e-06_dp
  real(dp), parameter :: C_v = 3441894705040.859_dp

  real(dp), parameter :: A_s = -208246976589.85126_dp
  real(dp), parameter :: B_s = -2.0162205697439128e-05_dp
  real(dp), parameter :: C_s = 235714178130.73007_dp

  ! We compute the saturation vapor pressure of H2O
  ! relative to its value at T = 373.15 K
  real(dp), parameter :: T0 = 373.15_dp ! K
  real(dp), parameter :: P0 = 1.0142e6_dp ! dynes/cm2

  real(dp), parameter :: T_freeze = 273.15_dp ! K, freezing temperature

contains

  !> Latent heat of vaporization of H2O
  function latent_heat_H2O_vap(T) result(L)
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    L = A_v*exp(B_v*T)+C_v
  end function

  !> Latent heat of sublimation of H2O
  function latent_heat_H2O_sub(T) result(L)
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    L = A_s*exp(B_s*T)+C_s
  end function

  !> Latent heat of vaporization or sublimation of H2O
  function latent_heat_H2O(T) result(L)
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    if (T > T_freeze) then
      L = latent_heat_H2O_vap(T)
    else ! (T <= T_freeze) then
      L = latent_heat_H2O_sub(T)
    endif
  end function

  !> This is $\int L/T^2 dT$
  function integral_fcn(A, B, C, T) result(out)
    use futils, only: expi
    real(dp), intent(in) :: A, B, C, T !! K
    real(dp) :: out
    out = (-A*B*T*expi(B*T) + A*exp(B*T) + C)/T
  end function

  !> Saturation pressure of H2O over water
  function sat_pressure_H2O_vap(T) result(p_H2O_sat)
    real(dp), intent(in) :: T !! K
    real(dp) :: p_H2O_sat !! dynes/cm2
    real(dp) :: tmp
    ! tmp = integral_fcn(A_v, B_v, C_v, T) - integral_fcn(A_v, B_v, C_v, T0)
    tmp = integral_fcn(A_v, B_v, C_v, T) - (-20369368.110596914_dp)
    p_H2O_sat = P0*exp((mu_H2O/Rgas)*(-tmp))
  end function

  !> Saturation pressure of H2O over ice
  function sat_pressure_H2O_sub(T) result(p_H2O_sat)
    real(dp), intent(in) :: T !! K
    real(dp) :: p_H2O_sat !! dynes/cm2
    real(dp) :: tmp
    ! tmp = (integral_fcn(A_v, B_v, C_v, T_freeze) - integral_fcn(A_v, B_v, C_v, T0)) + &
    !       (integral_fcn(A_s, B_s, C_s, T) - integral_fcn(A_s, B_s, C_s, T_freeze))
    tmp = (3141290.0653794562_dp - (-20369368.110596914_dp)) + &
          (integral_fcn(A_s, B_s, C_s, T) - 124184300.01342696_dp)
    p_H2O_sat = P0*exp((mu_H2O/Rgas)*(-tmp))
  end function

  !> Saturation pressure of H2O over ice or water
  function sat_pressure_H2O(T) result(p_H2O_sat)
    real(dp), intent(in) :: T !! K
    real(dp) :: p_H2O_sat !! dynes/cm2
    if (T > T_freeze) then
      p_H2O_sat = sat_pressure_H2O_vap(T)
    else ! (T <= T_freeze) then
      p_H2O_sat = sat_pressure_H2O_sub(T)
    endif
  end function

end module