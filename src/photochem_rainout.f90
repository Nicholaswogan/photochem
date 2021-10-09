
module photochem_rainout
  implicit none
  integer, private, parameter :: real_kind = kind(1.0d0)
contains
  
subroutine rainout(nq, nz, trop_ind, usol, T, den, rainout_rates)
  use photochem_const, only: k_boltz, N_avo, small_real
  use photochem_data, only: henry_data, LH2O, np, species_names
  use photochem_vars, only: edd, dz, z
  
  integer, intent(in) :: nq, nz, trop_ind
  real(real_kind), intent(in) :: usol(nq, nz)
  real(real_kind), intent(in) :: T(nz), den(nz)
  real(real_kind), intent(out) :: rainout_rates(nq, trop_ind)
  
  integer :: i, j
  
  real(real_kind) :: wH2O(trop_ind)
  real(real_kind) :: slope, intercept
  real(real_kind) :: denav_p, eddav_p, denav_m, eddav_m
  real(real_kind) :: scale_factor
  real(real_kind) :: k_bar, Q_i, H_coeff
  
  real(real_kind), parameter :: earth_rainfall_rate = 1.1d17 ! molecules/cm2/s
  real(real_kind), parameter :: fz = 0.05d0 ! fraction of the time it rains
  real(real_kind), parameter :: gamma = 4.d5 ! average time of storm cycles (s)
  real(real_kind), parameter :: LLL = 1.d0 ! g H2O/m3 of clouds
  real(real_kind), parameter :: C1 = 1.d-6 !dynes/bar
  real(real_kind), parameter :: C2 = 1.d-9 ![m3/cm3][L H2O/g H2O]
  real(real_kind), parameter :: MH2O = 18.d0 ! g H2O/mol H2O
  real(real_kind), parameter :: rho_H2O = 1000.d0 ! g H2O/L H2O
  
  !!!!!!! calculate raining rate !!!!!!!
  ! middle of atmosphere
  wH2O = 0.d0
  do i = 2,trop_ind-1
    denav_p = sqrt(den(i+1)*den(i))
    eddav_p = sqrt(edd(i+1)*edd(i))
    denav_m = sqrt(den(i-1)*den(i))
    eddav_m = sqrt(edd(i-1)*edd(i))
    wH2O(i) = (eddav_p*denav_p/dz(i)**2.d0) * usol(LH2O,i+1) &
            - (eddav_p*denav_p/dz(i)**2.d0 + eddav_m*denav_m/dz(i)**2.d0) * usol(lH2O,i) &
            + (eddav_m*denav_m/dz(i)**2.d0) * usol(lH2O,i-1)
    if (wH2O(i) < 0.d0) then
      wH2O(i) = 1.d-20
    endif
  enddo
  ! lets just linearly extrapolate wH2O to bottom and top grid cell
  !!! lower boundary !!!
  slope = (wH2O(3) - wH2O(2))/(dz(2))
  intercept = wH2O(2) - slope*z(2)
  wH2O(1) = slope*z(1) + intercept
  if (wH2O(1) < 0.d0) then
    wH2O(1) = 1.d-20
  endif
  !!! upper boundary !!!
  slope = (wH2O(trop_ind-1) - wH2O(trop_ind-2))/(dz(trop_ind-1))
  intercept = wH2O(trop_ind-1) - slope*z(trop_ind-1)
  wH2O(trop_ind) = slope*z(trop_ind) + intercept
  if (wH2O(trop_ind) < 0.d0) then
    wH2O(trop_ind) = 1.d-20
  endif
  ! Earth's globally averaged rainfall is 1.1e+17 molecules/cm2/s (Giorgi+1985)
  ! here we rescale the raining rate so that the total integrated value is
  ! the same as Earth's.
  scale_factor = earth_rainfall_rate/sum(wH2O*dz(1))
  wH2O = wH2O*scale_factor
  !!!!!!! end calculate raining rate !!!!!!!
  
  !!!!!!! dissolve gas in the rain !!!!!!!!!
  do j = 1,trop_ind
    do i = 1,nq
      H_coeff = henrys_law(max(T(j),273.15d0),henry_data(1,i),henry_data(2,i))*(1.d5)
      H_coeff = max(H_coeff, small_real)
      k_bar = (C1*k_boltz*T(j)*H_coeff/(1.d0+C1*C2*N_avo*LLL*k_boltz*T(j)*H_coeff)) &
               * (WH2O(j)*MH2O/rho_H2O) 
      Q_i = (1.d0-fz) + (fz/(gamma*k_bar))*(1.d0 - exp(-k_bar*gamma))
      rainout_rates(i,j) = (1.d0/(gamma*Q_i)) * (1.d0 - exp(-k_bar*gamma))
    enddo
  enddo
  !!!!!!! end dissolve gas in the rain !!!!!!!!!
  
end subroutine

function henrys_law(T, A, B) result(H)
  real(real_kind), intent(in) :: T
  real(real_kind), intent(in) :: A
  real(real_kind), intent(in) :: B
  real(real_kind) :: H ! mol/(kg * Pa)
  H = A*exp(B*(1.d0/298.15d0 - 1.d0/T))
end function
  
end module