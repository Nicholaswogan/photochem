module photochem_eqns
  use photochem_const, only: real_kind
  implicit none
  
contains
  
  subroutine gravity(radius, mass, nz, z, grav)
    use photochem_const, only: G_grav
    real(real_kind), intent(in) :: radius, mass ! radius in cm, mass in grams
    integer, intent(in) :: nz
    real(real_kind), intent(in) :: z(nz) ! cm
    real(real_kind), intent(out) :: grav(nz) ! cm/s2

    integer :: i
    
    do i = 1, nz              
      grav(i) = G_grav * (mass/1.d3) / ((radius + z(i))/1.d2)**2.d0
      grav(i) = grav(i)*1.d2 ! convert to cgs
    enddo 
    
  end subroutine
  
  subroutine vertical_grid(bottom, top, nz, z, dz)
    real(real_kind), intent(in) :: bottom, top
    integer, intent(in) :: nz
    real(real_kind), intent(out) :: z(nz), dz(nz)
  
    integer :: i
  
    dz = (top - bottom)/nz
    z(1) = dz(1)/2.d0
    do i = 2,nz
      z(i) = z(i-1) + dz(i)
    enddo
  end subroutine
  
  subroutine gibbs_energy_shomate(coeffs, T, gibbs)
    real(real_kind), intent(in) :: coeffs(7)
    real(real_kind), intent(in) :: T
    real(real_kind), intent(out) :: gibbs
    
    real(real_kind) :: enthalpy, entropy, TT
    
    TT = T/1000.d0
    enthalpy = coeffs(1)*TT + (coeffs(2)*TT**2)/2.d0 &
             + (coeffs(3)*TT**3)/3.d0  + (coeffs(4)*TT**4)/4.d0 &
             - coeffs(5)/TT + coeffs(6)
    entropy = coeffs(1)*log(TT) + coeffs(2)*TT &
            + (coeffs(3)*TT**2)/2.d0 + (coeffs(4)*TT**3)/3.d0 &
            - coeffs(5)/(2.d0 * TT**2) + coeffs(7)
    gibbs = enthalpy*1000.d0 - T*entropy
  end subroutine
  
  subroutine round(in,precision)
    implicit none
    real(real_kind), intent(inout) :: in
    integer, intent(in) :: precision
    integer :: order
    order = nint(log10(abs(in)))
    in = nint(in * 10.d0**(-precision-order),8)*10.d0**(precision+order)
  end subroutine
  
  subroutine press_and_den(nz, T, grav, Psurf, dz, &
                           mubar, pressure, density)
    use photochem_const, only: k_boltz, N_avo
  
    integer, intent(in) :: nz
    real(real_kind), intent(in) :: T(nz), grav(nz)
    real(real_kind), intent(in) :: Psurf, dz(nz), mubar(nz)
  
    real(real_kind), intent(out) :: pressure(nz)
    real(real_kind), intent(out) :: density(nz)
  
    real(real_kind) :: T_temp
    integer :: i
  
    ! first layer
    T_temp = T(1)
    pressure(1) = Psurf * exp(-((mubar(1) * grav(1))/(N_avo * k_boltz * T_temp)) * 0.5d0 * dz(1))
    density(1) = pressure(1)/(k_boltz * T(1))
    ! other layers
    do i = 2,nz
      T_temp = (T(i) + T(i-1))/2.d0
      pressure(i) = pressure(i-1) * exp(-((mubar(i) * grav(i))/(N_avo * k_boltz * T_temp))* dz(i))
      density(i) = pressure(i)/(k_boltz * T(i))
    enddo
  
  end subroutine
  
  subroutine molar_weight(nll, usol_layer, sum_usol_layer, masses, background_mu, mubar_layer)
    implicit none
    integer, intent(in) :: nll
    real(real_kind), intent(in) :: usol_layer(nll)
    real(real_kind), intent(in) :: sum_usol_layer
    real(real_kind), intent(in) :: masses(nll)
    real(real_kind), intent(in) :: background_mu
    real(real_kind), intent(out) :: mubar_layer
    integer :: j
    real(real_kind) :: f_background
  
    mubar_layer = 0.d0
    do j = 1, nll
      mubar_layer = mubar_layer + usol_layer(j) * masses(j)
    enddo
    f_background = 1.d0 - sum_usol_layer
    mubar_layer = mubar_layer + f_background * background_mu
  
  end subroutine
  
  function saturation_pressure(T, A, B, C) result(Psat)
      real(real_kind), intent(in) :: T ! Kelvin
      real(real_kind), intent(in) :: A, B, C ! parameters
      real(real_kind) :: Psat ! saturation pressure [bar]
      ! Simple exponential fit to saturation pressure data
      Psat = exp(A + B/T + C/T**2.d0)
  end function
  
  function saturation_density(T, A, B, C) result(nsat)
      use photochem_const, only: k_boltz
      real(real_kind), intent(in) :: T ! Kelvin
      real(real_kind), intent(in) :: A, B, C ! parameters
      real(real_kind) :: nsat ! saturation density [molecules/cm3]
      nsat = saturation_pressure(T, A, B, C)*1.d6/(k_boltz*T)
  end function
  
  function dynamic_viscosity_air(T) result(eta)
      real(real_kind), intent(in) :: T
      real(real_kind) :: eta ! dynamic viscosity [dynes s/cm^2]
      ! parameters speceific to Modern Earth air
      real(real_kind), parameter :: T0 = 273.d0 ! K
      real(real_kind), parameter :: eta0 = 1.716d-5 ! N s /m^2
      real(real_kind), parameter :: S = 111.d0 ! K
      real(real_kind), parameter :: unit_conversion = 10.d0 ! [dynes s/cm^2]/[N s/m^2]
      ! Dynamic viscosity of Air using the Sutherland relation.
      ! Reference: White (2006) "Viscous Fluid Flow"
      ! Equation 1-36 and Table 1-2 (air)
      eta = unit_conversion*eta0*(T/T0)**(3.d0/2.d0)*(T0 + S)/(T + S)
  end function

  function fall_velocity(gravity, partical_radius, particle_density, air_density, viscosity) result(wfall)
    real(real_kind), intent(in) :: gravity ! cm/s^2
    real(real_kind), intent(in) :: partical_radius ! cm
    real(real_kind), intent(in) :: particle_density ! g/cm^3
    real(real_kind), intent(in) :: air_density ! g/cm^3
    real(real_kind), intent(in) :: viscosity ! dynes s/cm^2
    real(real_kind) :: wfall ! fall velocity [cm/s]
    ! fall velocity from stokes law
    ! derived using Equation 9.29 in Seinfeld (2006) 
    ! title: "Atmospheric Chemistry and Physics"
    wfall = (2.d0/9.d0)*gravity*partical_radius**2.d0* &
            (particle_density - air_density)/(viscosity)
  end function
  
  function slip_correction_factor(partical_radius, density) result(correct_fac)
    real(real_kind), intent(in) :: partical_radius ! cm
    real(real_kind), intent(in) :: density ! molecules/cm3
    real(real_kind), parameter :: area_of_molecule = 6.0d-15 ! cm2
    real(real_kind) :: correct_fac
    real(real_kind) :: mean_free_path
    
    mean_free_path = 1.d0/(density*area_of_molecule)
    ! slip correction factor
    ! Equation 9.34 in Seinfeld (2006) 
    ! title: "Atmospheric Chemistry and Physics"
    correct_fac =  1.d0 + (mean_free_path/partical_radius)* &
                          (1.257d0 + 0.4d0*exp((-1.1d0*partical_radius)/(mean_free_path)))
  end function
  
  function binary_diffusion_param(mu_i, mubar, T) result(b)
    real(real_kind), intent(in) :: mu_i, mubar, T
    real(real_kind) :: b
    ! Banks and Kockarts 1973, Eq 15.29
    ! also Catling and Kasting 2017, Eq B.4 (although Catling has a typo,
    ! and is missing a power of 0.5)
    b = 1.52d18*((1.d0/mu_i+1.d0/mubar)**0.5d0)*(T**0.5d0)
  end function
  
  function arrhenius_rate(A, b, Ea, T) result(k)
    real(real_kind), intent(in) :: A, b, Ea, T
    real(real_kind) :: k
    k = A * T**b * exp(-Ea/T)
  end function

  function falloff_rate(kinf, Pr, F) result(k)
    real(real_kind), intent(in) :: kinf, Pr, F
    real(real_kind) :: k
    
    k = kinf * (Pr / (1.d0 + Pr)) * F
  end function

  function Troe_noT2(A, T1, T3, T, Pr) result(F)
    real(real_kind), intent(in) :: A, T1, T3, T, Pr
    real(real_kind) :: F
    
    real(real_kind) :: log10Fcent, f1, C, N
    
    log10Fcent = log10((1.d0-A)*exp(-T/T3) + A*exp(-T/T1))
    C = -0.4d0 - 0.67d0*log10Fcent
    N = 0.75d0 - 1.27d0*log10Fcent
    f1 = (log10(Pr) + C)/(N - 0.14d0*(log10(Pr + C)))
    F = 10.d0**((log10Fcent)/(1.d0 + f1**2.d0))
  end function

  function Troe_withT2(A, T1, T2, T3, T, Pr) result(F)
    real(real_kind), intent(in) :: A, T1, T2, T3, T, Pr
    real(real_kind) :: F
    
    real(real_kind) :: log10Fcent, f1, C, N
    
    log10Fcent = log10((1.d0-A)*exp(-T/T3) + A*exp(-T/T1) + exp(-T2/T))
    C = -0.4d0 - 0.67d0*log10Fcent
    N = 0.75d0 - 1.27d0*log10Fcent
    f1 = (log10(Pr) + C)/(N - 0.14d0*(log10(Pr + C)))
    F = 10.d0**((log10Fcent)/(1.d0 + f1**2.d0))
  end function
  
  function sat_pressure_H2O(T) result(p_H2O)
    real(real_kind), intent(in) :: T ! temperature in K
    real(real_kind) :: p_H2O
    real(real_kind), parameter :: lc = 2.5d6 ! specific enthalpy of H2O vaporization
    real(real_kind), parameter :: Rc = 461.d0 ! gas constant for water
    real(real_kind), parameter :: e0 = 611.d0 ! Pascals
    real(real_kind), parameter :: T0 = 273.15d0 ! K
    ! Catling and Kasting (Equation 1.49)
    p_H2O = 10.d0*e0*exp(lc/Rc*(1/T0 - 1/T))
    ! output is in dynes/cm2
  end function
  
  function damp_condensation_rate(A, rh0, rh) result(k)
    use photochem_const, only: pi
    real(real_kind), intent(in) :: A
    real(real_kind), intent(in) :: rh0
    real(real_kind), intent(in) :: rh ! the relative humidity
    real(real_kind) :: k
    ! k = 0 for rh = 0
    ! k = A for rh = infinity, approaches asymptotically
    ! k = 0.5*A for rh = rh0
    k = A*(2.d0/pi)*atan((rh - 1.d0)/(rh0 - 1.d0))
  end function
  
  
end module