submodule(photochem_evoatmosphere:photochem_evoatmosphere_rhs) photochem_evoatmosphere_rhs_climate
  implicit none

contains

  module subroutine temperature_objective(self, usol_den, T_surf, T_trop, temperature, err)
    use photochem_const, only: k_boltz
    use photochem_eqns, only: sat_pressure_H2O
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_den(:,:)
    real(dp), intent(in) :: T_surf, T_trop
    real(dp), intent(out) :: temperature(:)
    character(:), allocatable, intent(out) :: err

    real(dp) :: n_dry(self%var%nz), mu_dry(self%var%nz)
    real(dp) :: f_H2O_trop, P_H2O, n_H2O(self%var%nz), pressure(self%var%nz)

    real(dp), allocatable :: densities(:,:)

    integer :: i, j, k, trop_ind

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! compute dry atmospheric properties
    do j = 1,var%nz
      n_dry(j) = 0.0_dp
      do i = dat%ng_1,dat%nq
        if (i /= dat%LH2O) then
          n_dry(j) = n_dry(j) + usol_den(i,j)
        endif
      enddo
      mu_dry(j) = 0.0_dp
      do i = dat%ng_1,dat%nq
        if (i /= dat%LH2O) then
          mu_dry(j) = mu_dry(j) + (usol_den(i,j)/n_dry(j))*dat%species_mass(i)
        endif
      enddo
    enddo

    call draw_T_profile(self, usol_den, n_dry, mu_dry, T_surf, T_trop, temperature, trop_ind, err)
    if (allocated(err)) return

    ! For raditave transfer we use an estimate of the H2O density profile (NOT the real one)
    ! H2O is fixed to some relative humidity in the troposphere
    ! H2O is constant above the troposphere
    ! Guesstimate the H2O density
    do j = 1,trop_ind
      P_H2O = sat_pressure_H2O(temperature(j))
      n_H2O(j) = P_H2O/(k_boltz*temperature(j))
    enddo
    ! Guess mixing ratio of H2O at the tropopause
    f_H2O_trop = n_H2O(trop_ind)/(n_dry(trop_ind) + n_H2O(trop_ind))
    ! assume H2O mixing ratio is constant above the troposphere
    do j = trop_ind+1,var%nz
      n_H2O(j) = n_dry(j)*(f_H2O_trop/(1.0_dp - f_H2O_trop))
    enddo

    allocate(densities(var%nz,dat%nll))
    do i = 1,dat%nll
      do j = 1,var%nz
        k = i + dat%npq
        densities(j,i) = usol_den(k,j)
      enddo
    enddo
    k = dat%LH2O + dat%npq
    densities(:,k) = n_H2O

    do j = 1,var%nz
      ! pressure in bars
      pressure(j) = sum(densities(j,:))*k_boltz*temperature(j)/1.0e6_dp
    enddo

    ! do raditative transfer
    ! call self%rad%radiate(T_surf, temperature, pressure, densities, var%dz, err)
    ! if (allocated(err)) return

  end subroutine

  
  module subroutine draw_T_profile(self, usol_den, n_dry, mu_dry, T_surf, T_trop, temperature, trop_ind, err)
    use photochem_const, only: k_boltz
    use photochem_eqns, only: heat_capacity_eval
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_den(:,:)
    real(dp), intent(in) :: n_dry(:), mu_dry(:)
    real(dp), intent(in) :: T_surf, T_trop
    real(dp), intent(out) :: temperature(:)
    integer, intent(out) :: trop_ind
    character(:), allocatable, intent(out) :: err

    real(dp) :: cp_dry
    real(dp) :: T_edge, dT_dz
    integer :: j
    logical :: trop_hit

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! integrate from surface to center of first grid cell
    cp_dry = dry_heat_capacity_of_layer(dat, T_surf, usol_den(:,1), n_dry(1), err)
    if (allocated(err)) return
    dT_dz = dT_dz_moist(T_surf, n_dry(1), mu_dry(1), cp_dry, var%grav(1))    
    temperature(1) = T_surf + dT_dz*(0.5_dp*var%dz(1))

    if (temperature(1) < T_trop) then
      err = 'ground grid cell temperature is less than the troposphere temperature.'
      return
    endif

    trop_hit = .false. 

    do j = 1,var%nz-1

      ! Integrate to top edge of layer j
      cp_dry = dry_heat_capacity_of_layer(dat, temperature(j), usol_den(:,j), n_dry(j), err)
      if (allocated(err)) return
      dT_dz = dT_dz_moist(temperature(j), n_dry(j), mu_dry(j), cp_dry, var%grav(j))    
      T_edge = temperature(j) + dT_dz*(0.5_dp*var%dz(j))
      
      ! Integrate from top edge of j to middle of j+1
      cp_dry = dry_heat_capacity_of_layer(dat, T_edge, usol_den(:,j+1), n_dry(j+1), err)
      if (allocated(err)) return
      dT_dz = dT_dz_moist(T_edge, n_dry(j+1), mu_dry(j+1), cp_dry, var%grav(j+1))
      temperature(j+1) = T_edge + dT_dz*(0.5_dp*var%dz(j+1))

      if (temperature(j+1) <= T_trop) then
        trop_hit = .true.
        trop_ind = j
        temperature(j+1:) = temperature(j+1)
        exit
      endif
    enddo

    if (.not. trop_hit) then
      err = "Tropopause was never reached!"
      return
    endif

  end subroutine

  function dry_heat_capacity_of_layer(dat, T, usol_den_layer, n_dry, err) result(cp_dry)
    use photochem_eqns, only: heat_capacity_eval
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: T, usol_den_layer(:)
    real(dp), intent(in) :: n_dry
    character(:), allocatable, intent(out) :: err
    real(dp) :: cp_dry

    real(dp) :: cp
    logical :: found
    integer :: i

    cp_dry = 0.0_dp
    do i = dat%ng_1,dat%nq
      if (i /= dat%LH2O) then

        ! J/(mol*K)
        call heat_capacity_eval(dat%thermo_data(i-dat%npq), T, found, cp)
        if (.not. found) then
          err = 'Could not compute heat capacity for '//trim(dat%species_names(i))
          return
        endif

        ! J/(kg*K)
        cp_dry = cp_dry + (usol_den_layer(i)/n_dry)*cp*(1.0_dp/(dat%species_mass(i)*1.0e-3_dp))
      endif
    enddo
    ! convert to erg/(g*K)
    cp_dry = cp_dry*1.0e4_dp

  end function

  pure function dP_dry_dT(T, P_dry, mu_dry, cp_dry, P_H2O, mu_H2O, cp_H2O, L_H2O) result(dPdT)
    use photochem_const, only: Rgas_SI => Rgas
    real(dp), intent(in) :: T
    real(dp), intent(in) :: P_dry, mu_dry, cp_dry
    real(dp), intent(in) :: P_H2O, mu_H2O, cp_H2O, L_H2O
    real(dp) :: dPdT

    real(dp), parameter :: Rgas = Rgas_SI*1.0e7_dp

    dPdT = (P_dry*mu_dry*cp_dry)/(Rgas*T)*&
           (P_dry + (cp_H2O/cp_dry + (L_H2O*mu_H2O/(Rgas*T) - 1.0_dp)&
                  *(L_H2O/(cp_dry*T)))*((mu_H2O*P_H2O)/(mu_dry)))/ &
           (P_dry + (L_H2O*mu_dry)/(Rgas*T)*(mu_H2O*P_H2O)/(mu_dry))
  end function
  
  pure function dP_H2O_dT(T, P_H2O, mu_H2O, L_H2O) result(dPdT)
    use photochem_const, only: Rgas_SI => Rgas
    real(dp), intent(in) :: T
    real(dp), intent(in) :: P_H2O, mu_H2O, L_H2O
    real(dp) :: dPdT

    real(dp), parameter :: Rgas = Rgas_SI*1.0e7_dp

    dPdT = (L_H2O*mu_H2O*P_H2O)/(Rgas*T**2.0_dp)
  end function

  pure function dT_dz_moist(T, n_dry, mu_dry, cp_dry, grav) result(dTdz)
    use photochem_eqns, only: sat_pressure_H2O
    use photochem_const, only: k_boltz, N_avo
    real(dp), intent(in) :: T
    real(dp), intent(in) :: n_dry, mu_dry, cp_dry, grav
    real(dp) :: dTdz

    real(dp) :: dPdT_dry, dPdT_H2O
    real(dp) :: P_H2O, n_H2O
    real(dp) :: cp_H2O, L_H2O
    real(dp) :: P_dry
    real(dp) :: n_tot, mubar

    real(dp), parameter :: mu_H2O = 18.0_dp

    P_H2O = sat_pressure_H2O(T)
    n_H2O = P_H2O/(k_boltz*T)
    L_H2O = latent_heat_vap_H2O(T)
    cp_H2O = heat_capacity_H2O(T) ! J/(mol*K)
    cp_H2O = cp_H2O*(1.0_dp/(mu_H2O*1.0e-3_dp)) ! J/(kg*K)
    cp_H2O = cp_H2O*1.0e4_dp ! convert to erg/(g*K)

    ! dry pressure
    P_dry = n_dry*k_boltz*T

    ! combine dry and moist atmosphere
    n_tot = n_dry + n_H2O
    mubar = (n_dry/n_tot)*mu_dry + (n_H2O/n_tot)*mu_H2O

    dPdT_dry = dP_dry_dT(T, P_dry, mu_dry, cp_dry, P_H2O, mu_H2O, cp_H2O, L_H2O)
    dPdT_H2O = dP_H2O_dT(T, P_H2O, mu_H2O, L_H2O) 
    
    dTdz = (-grav*mubar*n_tot/N_avo) * (dPdT_dry + dPdT_H2O)**(-1.0_dp)

  end function

  pure function latent_heat_vap_H2O(T) result(L_H2O)
    real(dp), intent(in) :: T ! K
    real(dp) :: L_H2O ! erg/g
    L_H2O = 1.91846e6_dp*(T/(T - 33.91_dp))**2.0_dp ! J/kg
    L_H2O = L_H2O*1.0e4_dp ! convert to erg/g
  end function

  pure function heat_capacity_H2O(T) result(cp)
    use photochem_eqns, only: heat_capacity_shomate
    real(dp), intent(in) :: T !! K
    real(dp) :: cp !! J/(mol*K)
    
    real(dp), parameter :: coeffs_1(7) = &
          [30.092, 6.832514, 6.793435, -2.53448, 0.082139, -250.881, 223.3967]
    real(dp), parameter :: coeffs_2(7) = &
          [41.96426, 8.622053, -1.49978, 0.098119, -11.15764, -272.1797, 219.7809]
    
    if (T > 0.0_dp .and. T <= 1700.0_dp) then
      cp = heat_capacity_shomate(coeffs_1, T)
    elseif (T > 1700.0_dp .and. T <= 6000.0_dp) then
      cp = heat_capacity_shomate(coeffs_2, T)
    endif

  end function

end submodule