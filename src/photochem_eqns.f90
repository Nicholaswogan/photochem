#:set TYPES = ['real(dp)', 'type(dual)']
#:set NAMES = ['real', 'dual']
#:set TYPES_NAMES = list(zip(TYPES, NAMES))

module photochem_eqns
  use photochem_const, only: dp
  implicit none

  interface damp_condensation_rate
    module procedure :: damp_condensation_rate_real, damp_condensation_rate_dual 
  end interface
  
contains
  
  pure subroutine gravity(radius, mass, nz, z, grav)
    use photochem_const, only: G_grav
    real(dp), intent(in) :: radius, mass ! radius in cm, mass in grams
    integer, intent(in) :: nz
    real(dp), intent(in) :: z(nz) ! cm
    real(dp), intent(out) :: grav(nz) ! cm/s2

    integer :: i
    
    do i = 1, nz              
      grav(i) = G_grav * (mass/1.e3_dp) / ((radius + z(i))/1.e2_dp)**2.0_dp
      grav(i) = grav(i)*1.e2_dp ! convert to cgs
    enddo 
    
  end subroutine
  
  pure subroutine vertical_grid(bottom, top, nz, z, dz)
    real(dp), intent(in) :: bottom, top
    integer, intent(in) :: nz
    real(dp), intent(out) :: z(nz), dz(nz)
  
    integer :: i
  
    dz = (top - bottom)/nz
    z(1) = dz(1)/2.0_dp
    do i = 2,nz
      z(i) = z(i-1) + dz(i)
    enddo
  end subroutine
  
  pure function gibbs_energy_shomate(coeffs, T) result(gibbs)
    real(dp), intent(in) :: coeffs(7)
    real(dp), intent(in) :: T
    real(dp) :: gibbs
    
    real(dp) :: enthalpy, entropy, TT
    
    TT = T/1000.0_dp
    enthalpy = coeffs(1)*TT + (coeffs(2)*TT**2)/2.0_dp &
             + (coeffs(3)*TT**3)/3.0_dp  + (coeffs(4)*TT**4)/4.0_dp &
             - coeffs(5)/TT + coeffs(6)
    entropy = coeffs(1)*log(TT) + coeffs(2)*TT &
            + (coeffs(3)*TT**2)/2.0_dp + (coeffs(4)*TT**3)/3.0_dp &
            - coeffs(5)/(2.0_dp * TT**2) + coeffs(7)
    gibbs = enthalpy*1000.0_dp - T*entropy
  end function
  
  pure function gibbs_energy_nasa9(coeffs, T) result(gibbs)
    use photochem_const, only: Rgas
    real(dp), intent(in) :: coeffs(9)
    real(dp), intent(in) :: T
    real(dp) :: gibbs
    
    real(dp) :: enthalpy, entropy
    
    enthalpy = (- coeffs(1)*T**(-2.0_dp) + coeffs(2)*log(T)/T &
                + coeffs(3) + coeffs(4)*T/2.0_dp + coeffs(5)*T**(2.0_dp)/3.0_dp &
                + coeffs(6)*T**(3.0_dp)/4.0_dp + coeffs(7)*T**(4.0_dp)/5.0_dp &
                + coeffs(8)/T)*T*Rgas
             
    entropy = (- coeffs(1)*T**(-2.0_dp)/2.0_dp - coeffs(2)*T**(-1.0_dp) &
               + coeffs(3)*log(T) + coeffs(4)*T + coeffs(5)*T**(2.0_dp)/2.0_dp &
               + coeffs(6)*T**(3.0_dp)/3.0_dp + coeffs(7)*T**(4.0_dp)/4.0_dp &
               + coeffs(9))*Rgas
               
    gibbs = enthalpy - T*entropy
  end function

  pure function gibbs_energy_nasa7(coeffs, T) result(gibbs)
    use photochem_const, only: Rgas
    real(dp), intent(in) :: coeffs(:)
    real(dp), intent(in) :: T
    real(dp) :: gibbs

    real(dp) :: enthalpy, entropy
    real(dp) :: a0, a1, a2, a3, a4, a5, a6

    a0 = coeffs(1)
    a1 = coeffs(2)
    a2 = coeffs(3)
    a3 = coeffs(4)
    a4 = coeffs(5)
    a5 = coeffs(6)
    a6 = coeffs(7)

    enthalpy = &
      a0 + &
      a1*T/2.0_dp + &
      a2*T**2.0_dp/3.0_dp + &
      a3*T**3.0_dp/4.0_dp + &
      a4*T**4.0_dp/5.0_dp + &
      a5/T
    enthalpy = enthalpy*Rgas*T
    
    entropy = &
      a0*log(T) + &
      a1*T + &
      a2*T**2.0_dp/2.0_dp + &
      a3*T**3.0_dp/3.0_dp + &
      a4*T**4.0_dp/4.0_dp + &
      a6
    entropy = entropy*Rgas

    gibbs = enthalpy - T*entropy

  end function
  
  pure subroutine gibbs_energy_eval(thermo, T, found, gibbs_energy)
    use photochem_enum, only: ShomatePolynomial, Nasa9Polynomial, Nasa7Polynomial
    use photochem_types, only: ThermodynamicData
    
    type(ThermodynamicData), intent(in) :: thermo
    real(dp), intent(in) :: T
    logical, intent(out) :: found
    real(dp), intent(out) :: gibbs_energy
    
    integer :: k
    
    found = .false.
    do k = 1,thermo%ntemps
      if (T >= thermo%temps(k) .and. &
          T <  thermo%temps(k+1)) then
          
        found = .true.
        if (thermo%dtype == ShomatePolynomial) then
          gibbs_energy = gibbs_energy_shomate(thermo%data(1:7,k), T)
        elseif (thermo%dtype == Nasa9Polynomial) then
          gibbs_energy = gibbs_energy_nasa9(thermo%data(1:9,k), T)
        elseif (thermo%dtype == Nasa7Polynomial) then
          gibbs_energy = gibbs_energy_nasa7(thermo%data(1:7,k), T)          
        endif
        
        exit
        
      endif
    enddo

  end subroutine

  pure function heat_capacity_shomate(coeffs, T) result(cp)
    real(dp), intent(in) :: coeffs(7)
    real(dp), intent(in) :: T !! K
    real(dp) :: cp !! J/(mol K)
    
    real(dp) :: TT
    
    TT = T/1000.0_dp
    cp = coeffs(1) + coeffs(2)*TT + coeffs(3)*TT**2 + &
         coeffs(4)*TT**3 + coeffs(5)/TT**2
  end function
  
  pure subroutine heat_capacity_eval(thermo, T, found, cp)
    use photochem_enum, only: ShomatePolynomial, Nasa9Polynomial
    use photochem_types, only: ThermodynamicData
  
    type(ThermodynamicData), intent(in) :: thermo
    real(dp), intent(in) :: T !! K
    logical, intent(out) :: found
    real(dp), intent(out) :: cp !! J/(mol*K)
  
    integer :: k

    found = .false.
    do k = 1,thermo%ntemps
      if (T >= thermo%temps(k) .and. &
          T <  thermo%temps(k+1)) then
  
        found = .true.
        if (thermo%dtype == ShomatePolynomial) then
          cp = heat_capacity_shomate(thermo%data(1:7,k), T)
        elseif (thermo%dtype == Nasa9Polynomial) then
          ! gibbs_energy = heat_capacity_nasa9(thermo%data(1:9,k), T)
          found = .false.         
        endif
  
        exit
  
      endif
    enddo

  end subroutine

  pure subroutine press_and_den(nz, T, grav, Psurf, dz, &
                           mubar, pressure, density)
    use photochem_const, only: k_boltz, N_avo
  
    integer, intent(in) :: nz
    real(dp), intent(in) :: T(nz), grav(nz)
    real(dp), intent(in) :: Psurf, dz(nz), mubar(nz)
  
    real(dp), intent(out) :: pressure(nz)
    real(dp), intent(out) :: density(nz)
  
    real(dp) :: T_temp
    integer :: i
  
    ! first layer
    T_temp = T(1)
    pressure(1) = Psurf * exp(-((mubar(1) * grav(1))/(N_avo * k_boltz * T_temp)) * 0.5e0_dp * dz(1))
    density(1) = pressure(1)/(k_boltz * T(1))
    ! other layers
    do i = 2,nz
      T_temp = (T(i) + T(i-1))/2.0_dp
      pressure(i) = pressure(i-1) * exp(-((mubar(i) * grav(i))/(N_avo * k_boltz * T_temp))* dz(i))
      density(i) = pressure(i)/(k_boltz * T(i))
    enddo
  
  end subroutine
  
  pure function dynamic_viscosity_air(T) result(eta)
      real(dp), intent(in) :: T
      real(dp) :: eta ! dynamic viscosity [dynes s/cm^2]
      ! parameters speceific to Modern Earth air
      real(dp), parameter :: T0 = 273.0_dp ! K
      real(dp), parameter :: eta0 = 1.716e-5_dp ! N s /m^2
      real(dp), parameter :: S = 111.0_dp ! K
      real(dp), parameter :: unit_conversion = 10.0_dp ! [dynes s/cm^2]/[N s/m^2]
      ! Dynamic viscosity of Air using the Sutherland relation.
      ! Reference: White (2006) "Viscous Fluid Flow"
      ! Equation 1-36 and Table 1-2 (air)
      eta = unit_conversion*eta0*(T/T0)**(3.0_dp/2.0_dp)*(T0 + S)/(T + S)
  end function

  pure function fall_velocity(grav, partical_radius, particle_density, air_density, viscosity) result(wfall)
    real(dp), intent(in) :: grav ! cm/s^2
    real(dp), intent(in) :: partical_radius ! cm
    real(dp), intent(in) :: particle_density ! g/cm^3
    real(dp), intent(in) :: air_density ! g/cm^3
    real(dp), intent(in) :: viscosity ! dynes s/cm^2
    real(dp) :: wfall ! fall velocity [cm/s]
    ! fall velocity from stokes law
    ! derived using Equation 9.29 in Seinfeld (2006) 
    ! title: "Atmospheric Chemistry and Physics"
    wfall = (2.0_dp/9.0_dp)*grav*partical_radius**2.0_dp* &
            (particle_density - air_density)/(viscosity)
  end function

  pure function mean_free_path(den, T, mubar) result(lambda)
    use photochem_const, only: pi, k_boltz, N_avo
    real(dp), intent(in) :: den ! molecules/cm3
    real(dp), intent(in) :: T ! Kelvin
    real(dp), intent(in) :: mubar ! g/mol
    real(dp) :: lambda ! cm
    
    real(dp) :: eta ! dynamic viscosity [dynes s/cm^2]

    eta = dynamic_viscosity_air(T)
    lambda = (2.0_dp*eta/den)*sqrt(N_avo*pi/(8.0_dp*k_boltz*T*mubar))
  end function
  
  pure function slip_correction_factor(partical_radius, density) result(correct_fac)
    real(dp), intent(in) :: partical_radius ! cm
    real(dp), intent(in) :: density ! molecules/cm3
    real(dp), parameter :: area_of_molecule = 6.0e-15_dp ! cm2
    real(dp) :: correct_fac
    real(dp) :: lambda
    
    ! The density dependence really is nuts. I'm going to
    ! dampen out dependence with a power.
    ! mean_free_path = 1.0_dp/(density*area_of_molecule)
    lambda = 1.0_dp/(density**(1.0e0_dp)*area_of_molecule)
    ! slip correction factor
    ! Equation 9.34 in Seinfeld (2006) 
    ! title: "Atmospheric Chemistry and Physics"
    correct_fac = 1.0_dp + (lambda/partical_radius)* &
                          (1.257e0_dp + 0.4e0_dp*exp((-1.1e0_dp*partical_radius)/(lambda)))
  end function
  
  function default_binary_diffusion_param(mu_i, mubar, T) result(b)
    use iso_c_binding, only: c_double
    real(c_double), value, intent(in) :: mu_i, mubar, T
    real(c_double) :: b
    ! Banks and Kockarts 1973, Eq 15.29
    ! also Catling and Kasting 2017, Eq B.4 (although Catling has a typo,
    ! and is missing a power of 0.5)
    b = 1.52e18_dp*((1.0_dp/mu_i+1.0_dp/mubar)**0.5e0_dp)*(T**0.5e0_dp)
  end function
  
  pure function arrhenius_rate(A, b, Ea, T) result(k)
    real(dp), intent(in) :: A, b, Ea, T
    real(dp) :: k
    k = A * T**b * exp(-Ea/T)
  end function

  pure function falloff_rate(kinf, Pr, F) result(k)
    real(dp), intent(in) :: kinf, Pr, F
    real(dp) :: k
    
    k = kinf * (Pr / (1.0_dp + Pr)) * F
  end function

  pure function Troe_noT2(A, T1, T3, T, Pr) result(F)
    real(dp), intent(in) :: A, T1, T3, T, Pr
    real(dp) :: F
    
    real(dp) :: log10Fcent, f1, C, N
    
    log10Fcent = log10((1.0_dp-A)*exp(-T/T3) + A*exp(-T/T1))
    C = -0.4e0_dp - 0.67e0_dp*log10Fcent
    N = 0.75e0_dp - 1.27e0_dp*log10Fcent
    f1 = (log10(Pr) + C)/(N - 0.14e0_dp*(log10(Pr + C)))
    F = 10.0_dp**((log10Fcent)/(1.0_dp + f1**2.0_dp))
  end function

  pure function Troe_withT2(A, T1, T2, T3, T, Pr) result(F)
    real(dp), intent(in) :: A, T1, T2, T3, T, Pr
    real(dp) :: F
    
    real(dp) :: log10Fcent, f1, C, N
    
    log10Fcent = log10((1.0_dp-A)*exp(-T/T3) + A*exp(-T/T1) + exp(-T2/T))
    C = -0.4e0_dp - 0.67e0_dp*log10Fcent
    N = 0.75e0_dp - 1.27e0_dp*log10Fcent
    f1 = (log10(Pr) + C)/(N - 0.14e0_dp*(log10(Pr + C)))
    F = 10.0_dp**((log10Fcent)/(1.0_dp + f1**2.0_dp))
  end function

  #:for TYPE1, NAME in TYPES_NAMES
  pure function damp_condensation_rate_${NAME}$(A, rhc, rh0, rh) result(k)
    #:if NAME == 'dual'
    use differentia
    #:endif
    use photochem_const, only: pi
    real(dp), intent(in) :: A
    real(dp), intent(in) :: rhc, rh0
    ${TYPE1}$, intent(in) :: rh ! the relative humidity
    ${TYPE1}$ :: k 
    ! k = 0 for rh = rhc
    ! k = A for rh = infinity, approaches asymptotically
    ! k = 0.5*A for rh = rh0
    k = A*(2.0_dp/pi)*atan((rh - rhc)/(rh0 - rhc))
  end function
  
  #:endfor
  pure function henrys_law(T, A, B) result(H)
    real(dp), intent(in) :: T
    real(dp), intent(in) :: A
    real(dp), intent(in) :: B
    real(dp) :: H ! mol/(kg * Pa)
    H = A*exp(B*(1.0_dp/298.15_dp - 1.0_dp/T))
  end function

  pure function zahnle_Hescape_coeff(S1) result(coeff)
    real(dp), intent(in) :: S1 ! EUV radiation relative to Modern Earth
    real(dp) :: coeff ! molecules/cm2/s

    real(dp), parameter :: A = 2.0e12_dp ! molecules/cm2/s
    real(dp), parameter :: B2 = 0.006_dp ! no units

    ! Equation 47 in Zahnle et al. (2020), PSJ.
    ! `coeff` is the hydrogen escape coefficient, such that
    ! [Hydrogen Escape] = coeff*[H2 mixing ratio] 

    coeff = (A*S1)/sqrt(1.0_dp + B2*S1**2.0_dp)

  end function

end module