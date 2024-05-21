#:set TYPES = ['real(dp)', 'type(dual)']
#:set NAMES = ['real', 'dual']
#:set TYPES_NAMES = list(zip(TYPES, NAMES))

submodule(photochem_atmosphere) photochem_atmosphere_rhs
  implicit none
  
  interface dochem
    module procedure :: dochem_real, dochem_dual
  end interface
  
contains
      
  #:for TYPE1, NAME in TYPES_NAMES
  subroutine dochem_${NAME}$(self, neqs, nsp, np, nsl, nq, nz, trop_ind, nrT, usol, density, rx_rates, &
                    gas_sat_den, molecules_per_particle, &
                    H2O_sat_mix, H2O_rh, rainout_rates, scale_height, wfall, bg_density, &
                    densities, xp, xl, rhs)                 
    use photochem_enum, only: CondensingParticle
    use photochem_common, only: chempl, chempl_sl
    use photochem_eqns, only: damp_condensation_rate
    use photochem_const, only: N_avo, pi, small_real, T_crit_H2O
    #:if NAME == 'dual'
    use forwarddiff
    #:endif
    class(Atmosphere), target, intent(in) :: self
    integer, intent(in) :: neqs, nsp, np, nsl, nq, nz, trop_ind, nrT
    ${TYPE1}$, intent(in) :: usol(nq,nz)
    real(dp), intent(in) :: density(nz)
    real(dp), intent(in) :: rx_rates(nz,nrT)
    real(dp), intent(in) :: gas_sat_den(np,nz)
    real(dp), intent(in) :: molecules_per_particle(np,nz)
    real(dp), intent(in) :: H2O_sat_mix(nz), H2O_rh(trop_ind)
    real(dp), intent(in) :: rainout_rates(nq, trop_ind), scale_height(nz), wfall(np,nz)
    real(dp), intent(in) :: bg_density(nz)
    ${TYPE1}$, intent(inout) :: densities(nsp+1,nz), xp(nz), xl(nz)
    ${TYPE1}$, intent(inout) :: rhs(neqs)
    
    real(dp) :: cond_rate0
    ${TYPE1}$ :: rh, df_gas_dt, df_particle_dt, cond_rate
    integer :: i, ii, j, k, kk
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var
    
    #:if NAME == 'dual'
    rh = dual(size(usol(1,1)%der))
    df_gas_dt = dual(size(usol(1,1)%der))
    df_particle_dt = dual(size(usol(1,1)%der))
    cond_rate = dual(size(usol(1,1)%der))
    #:endif
    do j = 1,var%nz
      do k = 1,dat%np
        densities(k,j) = max(usol(k,j)*(density(j)/molecules_per_particle(k,j)),small_real)
      enddo
      do k = dat%ng_1,dat%nq
        densities(k,j) = usol(k,j)*density(j)
      enddo
      densities(nsp,j) = bg_density(j)
      densities(nsp+1,j) = 1.0_dp ! for hv
    enddo
    
    ! short lived
    do k = dat%nq+1,dat%nq+dat%nsl
      call chempl_sl(self%dat, self%var, densities, rx_rates, k, xp, xl) 
      densities(k,:) = xp/xl
    enddo
    
    rhs = 0.0_dp
    
    ! long lived              
    do i = dat%ng_1,dat%nq
      call chempl(self%dat, self%var, densities, rx_rates, i, xp, xl)
      do j = 1,var%nz
        k = i + (j - 1) * dat%nq
        rhs(k) = (xp(j) - xl(j))*(1.0_dp/density(j))
      enddo
    enddo
    
    if (dat%gas_rainout) then
      ! rainout rates
      do j = 1,trop_ind
        do i = 1,dat%nq
          k = i + (j - 1) * dat%nq
          rhs(k) = rhs(k) - rainout_rates(i,j)*usol(i,j)
        enddo
      enddo
    endif
    
    !!! Deal with H2O !!!

    ! to fix water in the troposphere, we produce or destroy water
    ! in the troposphere so that it reaches the target relative humidity
    if (dat%fix_water_in_trop) then
      do j = 1,trop_ind
        k = dat%LH2O + (j - 1) * dat%nq
        rhs(k) = rhs(k) + var%fast_arbitrary_rate*(H2O_sat_mix(j)*H2O_rh(j) - usol(dat%LH2O,j))
      enddo
    endif
    ! H2O condensation
    if (dat%water_cond) then
      if (dat%fix_water_in_trop) then
        ! need to start above tropopause if fixing
        ! water in troposphere
        i = trop_ind+1
      else
        i = 1
      endif
      do j = i,var%nz
        if (var%temperature(j) < T_crit_H2O) then
          ! water will condense if it is below the critical point.

          k = dat%LH2O + (j - 1) * dat%nq ! gas phase rhs index

          rh = max(usol(dat%LH2O,j)/H2O_sat_mix(j),small_real)

          if (rh > var%H2O_cond_params%RHc) then

            cond_rate0 = var%H2O_cond_params%k_cond*(var%edd(j)/scale_height(j)**2.0_dp)
            cond_rate = damp_condensation_rate(cond_rate0, &
                                               var%H2O_cond_params%RHc, &
                                               (1.0_dp + var%H2O_cond_params%smooth_factor)*var%H2O_cond_params%RHc, &
                                               rh)
                   
            ! Rate H2O gas is destroyed (1/s)
            df_gas_dt = - cond_rate*usol(dat%LH2O,j)
            rhs(k) = rhs(k) + df_gas_dt
            
          endif
        endif
      enddo
    endif
    
    if (dat%there_are_particles) then
      ! formation from reaction
      do i = 1,dat%np
        call chempl(self%dat, self%var, densities, rx_rates, i, xp, xl)
        do j = 1,var%nz
          k = i + (j - 1) * nq
          rhs(k) = rhs(k) + (xp(j) - xl(j))/density(j)
        enddo
      enddo
    
      ! particle condensation
      do j = 1,var%nz
        do i = 1,dat%np
          ! if this particle forms from condensation
          if (dat%particle_formation_method(i) == CondensingParticle) then
            ii = dat%particle_gas_phase_ind(i) ! index of gas phase
            kk = ii + (j - 1) * dat%nq ! gas phase rhs index
            k = i + (j - 1) * dat%nq ! particle rhs index

            ! compute the relative humidity
            rh = max(densities(ii,j)/gas_sat_den(i,j),small_real)
            
            ! if the gas phase is super-saturated
            if (rh > var%cond_params(i)%RHc) then
              ! Condensation occurs

              ! Condensation rate is based on how rapidly gases are mixing
              ! in the atmosphere. We also smooth the rate to prevent stiffness.
              cond_rate0 = var%cond_params(i)%k_cond*(var%edd(j)/scale_height(j)**2.0_dp)
              cond_rate = damp_condensation_rate(cond_rate0, &
                                                 var%cond_params(i)%RHc, &
                                                 (1.0_dp + var%cond_params(i)%smooth_factor)*var%cond_params(i)%RHc, &
                                                 rh)
            
              ! Rate gases are being destroyed (1/s)
              df_gas_dt = - cond_rate*usol(ii,j)
              rhs(kk) = rhs(kk) + df_gas_dt          
            
              ! Rate particles are being produced (1/s)
              df_particle_dt = - df_gas_dt
              rhs(k) = rhs(k) + df_particle_dt

              elseif (rh <= var%cond_params(i)%RHc .and. var%evaporation) then
              ! Evaporation occurs

              ! Evaporation rate is based on how rapidly particles are falling
              ! in the atmosphere. We also smooth the rate to prevent stiffness.
              cond_rate0 = var%cond_params(i)%k_evap*(wfall(i,j)/scale_height(j))
              cond_rate = damp_condensation_rate(cond_rate0, &
                                                 1.0_dp/var%cond_params(i)%RHc, &
                                                 (1.0_dp + var%cond_params(i)%smooth_factor)/var%cond_params(i)%RHc, &
                                                 1.0_dp/rh)

              ! Rate gases are being produced (1/s)                   
              df_gas_dt = cond_rate*usol(i,j)
              rhs(kk) = rhs(kk) + df_gas_dt

              ! Rate particles are being destroyed (1/s)
              df_particle_dt = - df_gas_dt
              rhs(k) = rhs(k) + df_particle_dt

            endif
          endif          
        enddo
      enddo
      
    endif
    
  end subroutine

  #:endfor
  subroutine diffusion_coefficients(dat, var, den, mubar, &
                                         DU, DD, DL, ADU, ADL, ADD, wfall, VH2_esc, VH_esc)
    use photochem_eqns, only: dynamic_viscosity_air, fall_velocity, slip_correction_factor, &
                              default_binary_diffusion_param
    use photochem_const, only: k_boltz, N_avo
    use photochem_enum, only: DiffusionLimHydrogenEscape
    use photochem_types, only: binary_diffusion_fcn

    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var

    real(dp), intent(in) :: den(:)
    real(dp), intent(in) :: mubar(:)
    
    real(dp), intent(out) :: DU(:,:), DL(:,:), DD(:,:) ! (nq,nz)
    real(dp), intent(out) :: ADU(:,:), ADL(:,:), ADD(:,:) ! (nq,nz)
    real(dp), intent(out) :: wfall(:,:) ! (npq,nz)
    real(dp), intent(out) :: VH2_esc, VH_esc

    real(dp) :: eddav_p, eddav_m, denav_p, denav_m, tav_p, tav_m
    real(dp) :: bx1x2_p, bx1x2_m, zeta_p, zeta_m
    real(dp) :: grav_p, grav_m, mubar_p, mubar_m
    real(dp) :: bx1x2
    
    ! for particles
    real(dp) :: air_density_pp, air_density
    real(dp) :: wfall_pp, wfall_i
    real(dp) :: viscosity_pp, viscosity
    real(dp) :: FF2, FF1

    ! molecular diffusion function
    procedure(binary_diffusion_fcn), pointer :: binary_diffusion_param

    integer :: i, j

    if (associated(var%custom_binary_diffusion_fcn)) then
      binary_diffusion_param => var%custom_binary_diffusion_fcn
    else
      binary_diffusion_param => default_binary_diffusion_param
    endif
    
    ! eddy diffusion. gases only
    ! middle
    do i = 2,var%nz-1
      eddav_p = sqrt(var%edd(i)*var%edd(i+1))
      eddav_m = sqrt(var%edd(i)*var%edd(i-1))
      denav_p = sqrt(den(i)*den(i+1))
      denav_m = sqrt(den(i)*den(i-1))
      do j = dat%ng_1,dat%nq
        DU(j,i) = (eddav_p*denav_p)/(den(i)*var%dz(i)**2.0_dp)
        DL(j,i) = (eddav_m*denav_m)/(den(i)*var%dz(i)**2.0_dp)
        DD(j,i) = - DU(j,i) - DL(j,i)
      enddo
    enddo
    ! top and bottom
    eddav_p = sqrt(var%edd(1)*var%edd(2))
    eddav_m = sqrt(var%edd(var%nz)*var%edd(var%nz-1))
    denav_p = sqrt(den(1)*den(2))
    denav_m = sqrt(den(var%nz)*den(var%nz-1))
    do j = dat%ng_1,dat%nq
      DU(j,1) = (eddav_p*denav_p)/(den(1)*var%dz(1)**2.0_dp)
      DD(j,1) = - DU(j,1)
      DL(j,var%nz) = (eddav_m*denav_m)/(den(var%nz)*var%dz(var%nz)**2.0_dp)
      DD(j,var%nz) = - DL(j,var%nz)
    enddo

    ! We do not include eddy diffusion for particles
    do j = 1,var%nz
      do i = 1,dat%np
        DU(i,j) = 0.0_dp
        DL(i,j) = 0.0_dp
        DD(i,j) = 0.0_dp
        ADU(i,j) = 0.0_dp
        ADL(i,j) = 0.0_dp
        ADD(i,j) = 0.0_dp
      enddo
    enddo
    
    ! Molecular diffusion. Only gas species
    ! middle
    do i = 2,var%nz-1
      tav_p = sqrt(var%temperature(i)*var%temperature(i+1))
      tav_m = sqrt(var%temperature(i)*var%temperature(i-1))
      grav_p = sqrt(var%grav(i)*var%grav(i+1))
      grav_m = sqrt(var%grav(i)*var%grav(i-1))
      mubar_p = sqrt(mubar(i)*mubar(i+1))
      mubar_m = sqrt(mubar(i)*mubar(i-1))
      do j = dat%ng_1, dat%nq
        bx1x2_p = binary_diffusion_param(dat%species_mass(j), mubar_p, tav_p)
        bx1x2_m = binary_diffusion_param(dat%species_mass(j), mubar_m, tav_m)
        
        DU(j,i) = DU(j,i) + bx1x2_p/(var%dz(i)**2.0_dp*den(i))
        DL(j,i) = DL(j,i) + bx1x2_m/(var%dz(i)**2.0_dp*den(i))
        DD(j,i) = - DU(j,i) - DL(j,i)
        
        zeta_p =  bx1x2_p*((dat%species_mass(j)*grav_p)/(k_boltz*tav_p*N_avo) &
                           - (mubar_p*grav_p)/(k_boltz*tav_p*N_avo) &
                           + 0.0_dp) ! zeroed out thermal diffusion   
        zeta_m =  bx1x2_m*((dat%species_mass(j)*grav_m)/(k_boltz*tav_m*N_avo) &
                          - (mubar_m*grav_m)/(k_boltz*tav_m*N_avo) &
                          + 0.0_dp) ! zeroed out thermal diffusion
                          
        ADU(j,i) = zeta_p/(2.0_dp*var%dz(i)*den(i)) 
        ADL(j,i) = - zeta_m/(2.0_dp*var%dz(i)*den(i))
        ADD(j,i) = ADU(j,i) + ADL(j,i)
      enddo
    enddo
    ! top and bottom
    tav_p = sqrt(var%temperature(1)*var%temperature(2))
    tav_m = sqrt(var%temperature(var%nz)*var%temperature(var%nz-1))
    grav_p = sqrt(var%grav(1)*var%grav(2))
    grav_m = sqrt(var%grav(var%nz)*var%grav(var%nz-1))
    mubar_p = sqrt(mubar(1)*mubar(2))
    mubar_m = sqrt(mubar(var%nz)*mubar(var%nz-1))
    ! lower boundary
    i = 1
    do j = dat%ng_1, dat%nq
      bx1x2_p = binary_diffusion_param(dat%species_mass(j), mubar_p, tav_p)
      DU(j,i) = DU(j,i) + bx1x2_p/(var%dz(i)**2.0_dp*den(i))
      DD(j,i) = - DU(j,i)

      zeta_p =  bx1x2_p*((dat%species_mass(j)*grav_p)/(k_boltz*tav_p*N_avo) &
                         - (mubar_p*grav_p)/(k_boltz*tav_p*N_avo) &
                         + 0.0_dp) ! zeroed out thermal diffusion   
      ADU(j,i) = zeta_p/(2.0_dp*var%dz(i)*den(i)) 
      ADD(j,i) = ADU(j,i)
    enddo
    ! upper boundary
    i = var%nz
    do j = dat%ng_1, dat%nq
      bx1x2_m = binary_diffusion_param(dat%species_mass(j), mubar_m, tav_m)
      
      DL(j,i) = DL(j,i) + bx1x2_m/(var%dz(i)**2.0_dp*den(i))
      DD(j,i) = - DL(j,i)
      
      zeta_m =  bx1x2_m*((dat%species_mass(j)*grav_m)/(k_boltz*tav_m*N_avo) &
                        - (mubar_m*grav_m)/(k_boltz*tav_m*N_avo) &
                        + 0.0_dp) ! zeroed out thermal diffusion
                        
      ADL(j,i) = - zeta_m/(2.0_dp*var%dz(i)*den(i))
      ADD(j,i) = ADL(j,i)
    enddo
    
    ! Falling particles. Only particles
    ! middle
    do i = 2,var%nz-1
      do j = 1,dat%npq
        
        air_density = (den(i)/N_avo)*mubar(i)
        viscosity = dynamic_viscosity_air(var%temperature(i))
        wfall_i = fall_velocity(var%grav(i), var%particle_radius(j,i), &
                              dat%particle_density(j), air_density, viscosity) &
                   *slip_correction_factor(var%particle_radius(j,i), den(i))

        air_density_pp = (den(i+1)/N_avo)*mubar(i+1)
        viscosity_pp = dynamic_viscosity_air(var%temperature(i+1))
        wfall_pp = fall_velocity(var%grav(i+1), var%particle_radius(j,i+1), &
                                 dat%particle_density(j), air_density_pp, viscosity_pp) &
                   *slip_correction_factor(var%particle_radius(j,i+1), den(i+1))

        FF2 = (wfall_pp*den(i+1))/(var%dz(i)*den(i))

        FF1 = (-wfall_i)/(var%dz(i))
      
        ADU(j,i) = FF2
        ADD(j,i) = FF1
        ADL(j,i) = 0.0_dp
      enddo
    enddo
    ! Lower boundary
    i = 1
    do j = 1,dat%npq
      air_density_pp = (den(i+1)/N_avo)*mubar(i+1)
      viscosity_pp = dynamic_viscosity_air(var%temperature(i+1))
      wfall_pp = fall_velocity(var%grav(i+1), var%particle_radius(j,i+1), &
                                dat%particle_density(j), air_density_pp, viscosity_pp) &
                  *slip_correction_factor(var%particle_radius(j,i+1), den(i+1))

      FF2 = (wfall_pp*den(i+1))/(var%dz(i)*den(i))

      ADU(j,i) = FF2
      ADD(j,i) = 0.0_dp
    enddo
    ! Upper boundary
    i = var%nz
    do j = 1,dat%npq
      air_density = (den(i)/N_avo)*mubar(i)
      viscosity = dynamic_viscosity_air(var%temperature(i))
      wfall_i = fall_velocity(var%grav(i), var%particle_radius(j,i), &
                              dat%particle_density(j), air_density, viscosity) &
                  *slip_correction_factor(var%particle_radius(j,i), den(i))
      FF1 = (-wfall_i)/(var%dz(i))

      ADD(j,i) = FF1
      ADL(j,i) = 0.0_dp
    enddo
    
    ! option to turn off everything but eddy diffusion
    do i = 1,dat%nq
      if (var%only_eddy(i)) then
        ADL(i,:) = 0.0_dp
        ADU(i,:) = 0.0_dp
        ADD(i,:) = 0.0_dp
      endif
    enddo
    
    ! H2 escape
    if (dat%H_escape_type == DiffusionLimHydrogenEscape) then
      if (dat%back_gas_name /= "H2") then
        bx1x2 = binary_diffusion_param(dat%species_mass(dat%LH2), mubar(var%nz), var%temperature(var%nz))
        VH2_esc = bx1x2/den(var%nz)*(-(dat%species_mass(dat%LH2)*var%grav(var%nz))/(k_boltz*var%temperature(var%nz)*N_avo) &
                                 + (mubar(var%nz)*var%grav(var%nz))/(k_boltz*var%temperature(var%nz)*N_avo))                     
      endif
      bx1x2 = binary_diffusion_param(dat%species_mass(dat%LH), mubar(var%nz), var%temperature(var%nz))
      VH_esc = bx1x2/den(var%nz)*(-(dat%species_mass(dat%LH)*var%grav(var%nz))/(k_boltz*var%temperature(var%nz)*N_avo) &
                              + (mubar(var%nz)*var%grav(var%nz))/(k_boltz*var%temperature(var%nz)*N_avo))
    endif
    
    ! wfall in center of grid cells. For boundary fluxes
    do i = 1,var%nz
      do j = 1,dat%npq
        air_density = (den(i)/N_avo)*mubar(i)
        viscosity = dynamic_viscosity_air(var%temperature(i))
        wfall(j,i) = fall_velocity(var%grav(i), var%particle_radius(j,i), &
                                   dat%particle_density(j), air_density, viscosity) &
                     *slip_correction_factor(var%particle_radius(j,i), den(i))
      enddo
    enddo
    
  end subroutine

  module subroutine prep_atm_background_gas(self, usol_in, usol, molecules_per_particle)
    use photochem_common, only: molec_per_particle
    use photochem_const, only: small_real
    use photochem_enum, only: MixingRatioBC
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    real(dp), intent(out) :: usol(:,:)
    real(dp), intent(out) :: molecules_per_particle(:,:)

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    integer :: i, j

    dat => self%dat
    var => self%var
    
    !!! alter input usol
    do j = 1,var%nz
      do i = 1,dat%nq
        if (usol_in(i,j) < 0.0_dp) then
          usol(i,j) = min(usol_in(i,j),-small_real)
        else
          usol(i,j) = max(usol_in(i,j), small_real)
        endif
      enddo
    enddo

    do i = 1,dat%nq
      if (var%lowerboundcond(i) == MixingRatioBC) then
        usol(i,1) = var%lower_fix_mr(i)
      endif
    enddo

    !!! molecules/particle
    if (dat%there_are_particles) then
      call molec_per_particle(dat, var, molecules_per_particle)
    endif

  end subroutine
  
  module subroutine prep_all_background_gas(self, usol_in, err)
    use photochem_enum, only: CondensingParticle
    use photochem_enum, only: ArrheniusSaturation, H2SO4Saturation
    use photochem_enum, only: DiffusionLimHydrogenEscape
    use photochem_common, only: reaction_rates, rainout, photorates
    use photochem_common, only: gas_saturation_density
    use photochem_eqns, only: saturation_density
    use clima_eqns_water, only: sat_pressure_H2O
    use photochem_eqns, only: molar_weight, press_and_den
    use photochem_const, only: pi, k_boltz, N_avo, small_real
  
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    character(:), allocatable, intent(out) :: err
  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    integer :: i, j, k
  
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    call prep_atm_background_gas(self, usol_in, wrk%usol, wrk%molecules_per_particle)

    !!! pressure, density and mean molcular weight
    do i = 1,var%nz
      wrk%sum_usol(i) = sum(wrk%usol(dat%ng_1:,i))
      if (wrk%sum_usol(i) > 1.0e0_dp) then
        err = 'Mixing ratios sum to >1.0 at some altitude (should be <=1).' // &
              ' The atmosphere is probably in a run-away state'
        return
      endif
    enddo

    do i = 1,var%nz
      call molar_weight(dat%nll, wrk%usol(dat%ng_1:,i), wrk%sum_usol(i), dat%species_mass(dat%ng_1:), dat%back_gas_mu, wrk%mubar(i))
    enddo
    
    call press_and_den(var%nz, var%temperature, var%grav, var%surface_pressure*1.e6_dp, var%dz, &
                       wrk%mubar, wrk%pressure, wrk%density)
    
    !!! H2O saturation
    if (dat%fix_water_in_trop) then
      do i = 1,var%trop_ind
        if (var%use_manabe) then
          ! manabe formula
          wrk%H2O_rh(i) = 0.77e0_dp*(wrk%pressure(i)/wrk%pressure(1)-0.02e0_dp)/0.98e0_dp
        else
          wrk%H2O_rh(i) = var%relative_humidity 
        endif
      enddo
    endif

    if (dat%water_cond .or. dat%fix_water_in_trop) then
      do i = 1,var%nz
        wrk%H2O_sat_mix(i) = sat_pressure_H2O(var%temperature(i))/wrk%pressure(i)
      enddo
    endif
  
    !!! diffusion and advection coefficients
    call diffusion_coefficients(dat, var, wrk%density, wrk%mubar, &
                                wrk%DU, wrk%DD, wrk%DL, wrk%ADU, wrk%ADL, wrk%ADD, &
                                wrk%wfall, wrk%VH2_esc, wrk%VH_esc)

    wrk%scale_height = (k_boltz*var%temperature(:)*N_avo)/(wrk%mubar(:)*var%grav(:))

    !!! H and H2 escape
    wrk%upper_veff_copy = var%upper_veff
    wrk%lower_vdep_copy = var%lower_vdep
    if (dat%H_escape_type == DiffusionLimHydrogenEscape) then
      if (dat%back_gas_name /= "H2") then
        wrk%upper_veff_copy(dat%LH2) = wrk%VH2_esc                     
      endif
      wrk%upper_veff_copy(dat%LH) = wrk%VH_esc 
    endif

    !!! Particle lower boundary, and saturation properties
    if (dat%there_are_particles) then
      do i = 1,dat%np
        ! Here we impose a lower boundary condition for particles. They fall out
        ! of the model according to the fall velocity.
        wrk%lower_vdep_copy(i) = wrk%lower_vdep_copy(i) + wrk%wfall(i,1)
      enddo

      call gas_saturation_density(dat, var, wrk%usol(dat%LH2O,:), wrk%pressure, &
                                  wrk%gas_sat_den)
    endif
  
    !!! densities
    do j = 1,var%nz
      do i = 1,dat%npq
        wrk%densities(i,j) = max(wrk%usol(i,j)*(wrk%density(j)/wrk%molecules_per_particle(i,j)), small_real)
      enddo
      do i = dat%ng_1,dat%nq
        wrk%densities(i,j) = wrk%usol(i,j)*wrk%density(j)
      enddo
      wrk%densities(dat%nsp,j) = (1.0_dp-wrk%sum_usol(j))*wrk%density(j) ! background gas
      wrk%densities(dat%nsp+1,j) = 1.0_dp ! for hv
    enddo
    
    !!! reaction rates
    call reaction_rates(self%dat, self%var, wrk%pressure, wrk%density, wrk%densities, wrk%rx_rates)

    ! Update the photon_flux if the function is associated.
    ! we use time wrk%tn, which MUST be updated.
    if (associated(var%photon_flux_fcn)) then
      call var%photon_flux_fcn(wrk%tn, dat%nw, var%photon_flux)
    endif
    call photorates(dat, var, wrk%densities, &
                    wrk%prates, wrk%surf_radiance, wrk%amean_grd, wrk%optical_depth, err)
    if (allocated(err)) return
  
    do i = 1,dat%kj
      k = dat%photonums(i)
      wrk%rx_rates(:,k) = wrk%prates(:,i) 
    enddo 
  
    !!! rainout rates
    if (dat%gas_rainout) then
      call rainout(self%dat, self%var, &
                   wrk%usol(dat%LH2O,:), wrk%density, wrk%rainout_rates)
    endif
  
  end subroutine
  
  module subroutine rhs_background_gas(self, neqs, tn, usol_flat, rhs, err)
    use photochem_enum, only: MosesBC, VelocityBC, MixingRatioBC, FluxBC, VelocityDistributedFluxBC
    use photochem_enum, only: ZahnleHydrogenEscape
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real  
    
    class(Atmosphere), target, intent(inout) :: self
    integer, intent(in) :: neqs
    real(dp), intent(in) :: tn
    real(dp), target, intent(in) :: usol_flat(neqs)
    real(dp), intent(out) :: rhs(neqs)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: disth, ztop, ztop1    
    integer :: i, k, j, jdisth
    
    real(dp), pointer :: usol_in(:,:)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    ! reshape
    usol_in(1:dat%nq,1:var%nz) => usol_flat(1:neqs)
    
    if (any(usol_flat /= usol_flat)) then
      err = 'Input mixing ratios to the rhs contains NaNs. This is typically '//&
            'related to some mixing ratios getting too negative.'
      return 
    endif

    ! time
    wrk%tn = tn
    
    ! fills self%wrk with data
    call prep_all_background_gas(self, usol_in, err)
    if (allocated(err)) return
    
    call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, &
                var%trop_ind, dat%nrT, wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, wrk%densities(dat%nsp,:), &
                wrk%densities, wrk%xp, wrk%xl, rhs) 

    ! Extra functions specifying production or destruction
    do i = 1,dat%nq
      if (associated(var%rate_fcns(i)%fcn)) then
        call var%rate_fcns(i)%fcn(tn, var%nz, wrk%xp) ! using wrk%xp space.
        do j = 1,var%nz
          k = i + (j-1)*dat%nq
          ! (moelcules/cm^3/s) (cm^3/molecules) = (1/s)
          rhs(k) = rhs(k) + wrk%xp(j)/wrk%density(j)
        enddo
      endif
    enddo

    ! diffusion (interior grid points)
    do j = 2,var%nz-1
      do i = 1,dat%nq
        k = i + (j-1)*dat%nq
        rhs(k) = rhs(k) + wrk%DU(i,j)*wrk%usol(i,j+1) + wrk%ADU(i,j)*wrk%usol(i,j+1) &
                        + wrk%DD(i,j)*wrk%usol(i,j) + wrk%ADD(i,j)*wrk%usol(i,j) &
                        + wrk%DL(i,j)*wrk%usol(i,j-1) + wrk%ADL(i,j)*wrk%usol(i,j-1)
      enddo
    enddo
    
    ! Lower boundary
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityBC .or. &
          var%lowerboundcond(i) == VelocityDistributedFluxBC) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        + wrk%DD(i,1)*wrk%usol(i,1) + wrk%ADD(i,1)*wrk%usol(i,1) &
                        - wrk%lower_vdep_copy(i)*wrk%usol(i,1)/var%dz(1)
      elseif (var%lowerboundcond(i) == MixingRatioBC) then
        rhs(i) = 0.0_dp
      elseif (var%lowerboundcond(i) == FluxBC) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        + wrk%DD(i,1)*wrk%usol(i,1) + wrk%ADD(i,1)*wrk%usol(i,1) &
                        + var%lower_flux(i)/(wrk%density(1)*var%dz(1))
      ! Moses (2001) boundary condition for gas giants
      ! A deposition velocity controled by how quickly gases
      ! turbulantly mix vertically
      elseif (var%lowerboundcond(i) == MosesBC) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        + wrk%DD(i,1)*wrk%usol(i,1) + wrk%ADD(i,1)*wrk%usol(i,1) &
                        - (var%edd(1)/wrk%scale_height(1))*wrk%usol(i,1)/var%dz(1)
      endif
    enddo

    ! Upper boundary
    do i = 1,dat%nq
      k = i + (var%nz-1)*dat%nq
      if (var%upperboundcond(i) == VelocityBC) then
        rhs(k) = rhs(k) + wrk%DD(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADD(i,var%nz)*wrk%usol(i,var%nz) &
                        + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1) &
                        - wrk%upper_veff_copy(i)*wrk%usol(i,var%nz)/var%dz(var%nz)    
      elseif (var%upperboundcond(i) == FluxBC) then
        rhs(k) = rhs(k) + wrk%DD(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADD(i,var%nz)*wrk%usol(i,var%nz) &
                        + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1) &
                        - var%upper_flux(i)/(wrk%density(var%nz)*var%dz(var%nz))
      endif
    enddo
    
    ! Distributed (volcanic) sources
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityDistributedFluxBC) then
        disth = var%lower_dist_height(i)*1.e5_dp        
        jdisth = minloc(var%Z,1, var%Z >= disth) - 1
        jdisth = max(jdisth,2)
        ztop = var%z(jdisth)-var%z(1)
        ztop1 = var%z(jdisth) + 0.5e0_dp*var%dz(jdisth)
        do j = 2,jdisth
          k = i + (j-1)*dat%nq
          rhs(k) = rhs(k) + 2.0_dp*var%lower_flux(i)*(ztop1-var%z(j))/(wrk%density(j)*ztop**2.0_dp)
        enddo
      endif
    enddo 

    ! zahnle hydrogen escape
    if (dat%H_escape_type == ZahnleHydrogenEscape) then

      ! for Zahnle hydrogen escape, we pull H2 out of 
      ! the bottom grid cell of the model.

      rhs(dat%LH2) = rhs(dat%LH2) &
      - dat%H_escape_coeff*wrk%usol(dat%LH2,1)/(wrk%density(1)*var%dz(1))
    endif
    
  end subroutine

  subroutine autodiff_chemistry_jacobian(self, usol, rhs, djac, err)
    use forwarddiff, only: jacobian, dual, BlockDiagonalJacobian, initialize_dual_array
    class(Atmosphere), target, intent(inout) :: self
    real(dp), target, contiguous, intent(in) :: usol(:,:)
    real(dp), intent(out) :: rhs(:), djac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), pointer :: usol_flat(:)
    type(dual), allocatable :: densities(:,:), xp(:), xl(:)
    integer :: blocksize

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    dat => self%dat
    var => self%var
    wrk => self%wrk

    usol_flat(1:var%neqs) => usol(:,:) ! reshape input
    blocksize = dat%nq ! blocksize is number of integrated species
    ! allocate some work memory
    allocate(densities(size(wrk%densities,1),size(wrk%densities,2)))
    call initialize_dual_array(densities, blocksize)
    allocate(xp(size(wrk%xp)))
    call initialize_dual_array(xp, blocksize)
    allocate(xl(size(wrk%xl)))
    call initialize_dual_array(xl, blocksize)

    ! Compute the jacobian
    call jacobian(fcn, usol_flat, rhs, djac, jt=BlockDiagonalJacobian, blocksize=blocksize, err=err)
    if (allocated(err)) return

  contains
    subroutine fcn(x_, f_)
      type(dual), target, intent(in) :: x_(:)
      type(dual), target, intent(out) :: f_(:)
      type(dual), pointer :: usol_(:,:)
      usol_(1:dat%nq,1:var%nz) => x_(:)
      call initialize_dual_array(f_, blocksize)
      call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, var%trop_ind, dat%nrT, &
                  usol_, wrk%density, wrk%rx_rates, &
                  wrk%gas_sat_den, wrk%molecules_per_particle, &
                  wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, wrk%densities(dat%nsp,:), &
                  densities, xp, xl, f_) 
    end subroutine
  end subroutine
  
  module subroutine jac_background_gas(self, lda_neqs, neqs, tn, usol_flat, jac, err)
    use photochem_enum, only: MosesBC, VelocityBC, MixingRatioBC, FluxBC, VelocityDistributedFluxBC
    use photochem_enum, only: ZahnleHydrogenEscape
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real
    
    class(Atmosphere), target, intent(inout) :: self
    integer, intent(in) :: lda_neqs, neqs
    real(dp), intent(in) :: tn
    real(dp), target, intent(in) :: usol_flat(neqs)
    real(dp), intent(out), target :: jac(lda_neqs)
    character(:), allocatable, intent(out) :: err
    
    real(dp), pointer :: usol_in(:,:)
    real(dp), pointer :: djac(:,:)
    real(dp) :: rhs(self%var%neqs)
  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    integer :: i, k, j, m, mm
  
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    ! reshape jac and input mixing ratios
    usol_in(1:dat%nq,1:var%nz) => usol_flat(1:neqs)
    djac(1:dat%lda,1:var%neqs) => jac
    
    if (any(usol_flat /= usol_flat)) then
      err = 'Input mixing ratios to the rhs contains NaNs. This is typically '//&
            'related to some mixing ratios getting too negative.'
      return 
    endif

    ! time
    wrk%tn = tn
  
    call prep_all_background_gas(self, usol_in, err)
    if (allocated(err)) return

    jac = 0.0_dp

    if (.not. var%autodiff) then; block
    real(dp) :: usol_perturb(dat%nq,var%nz)
    real(dp) :: R(var%nz)
    real(dp) :: rhs_perturb(var%neqs)
    real(dp) :: densities(dat%nsp+1,var%nz), xl(var%nz), xp(var%nz)

    ! Finite differenced Jacobian
  
    ! compute chemistry contribution to jacobian using forward differences
    call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, var%trop_ind, dat%nrT, &
                wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, wrk%densities(dat%nsp,:), &
                densities, xp, xl, rhs) 
    !$omp parallel private(i, j, k, m, mm, usol_perturb, R, densities, xl, xp, rhs_perturb)
    usol_perturb = wrk%usol
    !$omp do
    do i = 1,dat%nq
      do j = 1,var%nz
        R(j) = var%epsj*abs(wrk%usol(i,j))
        usol_perturb(i,j) = wrk%usol(i,j) + R(j)
      enddo
  
      call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, var%trop_ind, dat%nrT, &
                  usol_perturb, wrk%density, wrk%rx_rates, &
                  wrk%gas_sat_den, wrk%molecules_per_particle, &
                  wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, wrk%densities(dat%nsp,:), &
                  densities, xp, xl, rhs_perturb) 
  
      do m = 1,dat%nq
        mm = m - i + dat%kd
        do j = 1,var%nz
          k = i + (j-1)*dat%nq
          djac(mm,k) = (rhs_perturb(m + (j-1)*dat%nq) - rhs(m + (j-1)*dat%nq))/R(j)
        enddo
      enddo
  
      do j= 1,var%nz
        usol_perturb(i,j) = wrk%usol(i,j)
      enddo
    enddo
    !$omp enddo
    !$omp end parallel

    endblock; else

    ! Forward mode autodiff Jacobian
    call autodiff_chemistry_jacobian(self, wrk%usol, rhs, wrk%djac_chem, err)
    if (allocated(err)) return

    do mm = 1,dat%nq
      do m = 1,var%nz
        j = mm + dat%nq*(m-1)
        do i = 1,dat%nq
          djac(i + 2*dat%nq - (mm-1),j) = wrk%djac_chem(i,j)
        enddo
      enddo
    enddo

    endif
  
    ! diffusion (interior grid points)
    do j = 2,var%nz-1
      do i = 1,dat%nq
        k = i + (j-1)*dat%nq      
        djac(dat%ku,k+dat%nq) = wrk%DU(i,j) + wrk%ADU(i,j)
        djac(dat%kd,k) = djac(dat%kd,k) + wrk%DD(i,j) + wrk%ADD(i,j)     
        djac(dat%kl,k-dat%nq) = wrk%DL(i,j) + wrk%ADL(i,j)
      enddo
    enddo
  
    ! Lower boundary
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityBC .or. var%lowerboundcond(i) == VelocityDistributedFluxBC) then

        djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(dat%kd,i) = djac(dat%kd,i) + wrk%DD(i,1) + wrk%ADD(i,1) - wrk%lower_vdep_copy(i)/var%dz(1)
      elseif (var%lowerboundcond(i) == MixingRatioBC) then

        do m=1,dat%nq
          mm = dat%kd + i - m
          djac(mm,m) = 0.0_dp
        enddo
        djac(dat%ku,i+dat%nq) = 0.0_dp
        ! For some reason this term makes the integration
        ! much happier. I will keep it. Jacobians don't need to be perfect.
        djac(dat%kd,i) = - wrk%DU(i,1)
  
      elseif (var%lowerboundcond(i) == FluxBC) then
        djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(dat%kd,i) = djac(dat%kd,i) + wrk%DD(i,1) + wrk%ADD(i,1)
      elseif (var%lowerboundcond(i) == MosesBC) then
        djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(dat%kd,i) = djac(dat%kd,i) + wrk%DD(i,1) + wrk%ADD(i,1) - &
                         (var%edd(1)/wrk%scale_height(1))/var%dz(1)
      ! elseif (var%lowerboundcond(i) == MixDependentFluxBC) then
      !   djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
      !   djac(dat%kd,i) = djac(dat%kd,i) + wrk%DD(i,1) + wrk%ADD(i,1) &
      !                    + var%lower_mix_dep_flux(i)/(wrk%density(1)*var%dz(1))
      endif
    enddo
  
    ! Upper boundary
    do i = 1,dat%nq
      k = i + (var%nz-1)*dat%nq
      if (var%upperboundcond(i) == VelocityBC) then
        ! rhs(k) = rhs(k) - DL(i,nz)*usol(i,nz) &
        !                 + DL(i,nz)*usol(i,nz-1) + ADL(i,nz)*usol(i,nz-1) &
        !                 - upper_veff(i)*usol(i,nz)/dz(nz)    
  
        djac(dat%kd,k) = djac(dat%kd,k) + wrk%DD(i,var%nz) + wrk%ADD(i,var%nz) &
                        - wrk%upper_veff_copy(i)/var%dz(var%nz) 
        djac(dat%kl,k-dat%nq) = wrk%DL(i,var%nz) + wrk%ADL(i,var%nz)
      elseif (var%upperboundcond(i) == FluxBC) then
        djac(dat%kd,k) = djac(dat%kd,k) + wrk%DD(i,var%nz) + wrk%ADD(i,var%nz)
        djac(dat%kl,k-dat%nq) = wrk%DL(i,var%nz) + wrk%ADL(i,var%nz)
      endif
    enddo

    ! zahnle hydrogen escape
    if (dat%H_escape_type == ZahnleHydrogenEscape) then
      djac(dat%kd,dat%LH2) = djac(dat%kd,dat%LH2) + &
      - dat%H_escape_coeff/(wrk%density(1)*var%dz(1))
    endif
  
  end subroutine
  
  module subroutine right_hand_side_chem(self, usol, rhs, err)
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(out) :: rhs(:)
    character(:), allocatable, intent(out) :: err
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz .or. size(rhs) /= var%neqs) then
      err = "Input usol or rhs to dochem_implicit has the wrong dimensions"
      return
    endif

    call self%prep_atmosphere(usol, err)
    if (allocated(err)) return
    
    call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, &
                var%trop_ind, dat%nrT, wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, wrk%densities(dat%nsp,:), &
                wrk%densities, wrk%xp, wrk%xl, rhs)
                              
  end subroutine
  
  module subroutine production_and_loss(self, species, usol, pl, err)     
    use futils, only: argsort            
    use photochem_common, only: chempl_sl, chempl_t
    use photochem_types, only: ProductionLoss
    use photochem_const, only: small_real
  
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: species
    real(dp), intent(in) :: usol(:,:)
    type(ProductionLoss), intent(out) :: pl
    character(:), allocatable, intent(out) :: err
  
    real(dp) :: xl(self%var%nz), xp(self%var%nz)
    integer, allocatable :: prod_inds(:), loss_inds(:)
    integer :: sp_ind
    integer :: i, j, k, np, nl, nlT
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
  
    dat => self%dat
    var => self%var
    wrk => self%wrk
  
    
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = "Input usol to production_and_loss has the wrong dimensions"
      return
    endif
  
    sp_ind = findloc(dat%species_names(1:dat%nsp),species,1)
    if (sp_ind == 0) then
      err = "Species "//trim(species)//" is not in the list of species."
      return
    endif
    
    call self%prep_atmosphere(usol, err)
    if (allocated(err)) return
  
    np = dat%pl(sp_ind)%nump
    nl = dat%pl(sp_ind)%numl
    nlT = nl + 1 ! + 1 for rainout
    
    allocate(pl%production(var%nz,np))
    allocate(pl%loss(var%nz,nlT))
    allocate(pl%integrated_production(np), pl%integrated_loss(nlT))
    allocate(pl%loss_rx(nlT),pl%production_rx(np))
    allocate(prod_inds(np), loss_inds(nlT))
  
    do j = 1,var%nz
      do k = 1,dat%np
        wrk%densities(k,j) = max(wrk%usol(k,j)* &
                            (wrk%density(j)/wrk%molecules_per_particle(k,j)),small_real)
      enddo
      do k = dat%ng_1,dat%nq
        wrk%densities(k,j) = wrk%usol(k,j)*wrk%density(j)
      enddo
      wrk%densities(dat%nsp,j) = (1.0_dp-sum(wrk%usol(dat%ng_1:,j)))*wrk%density(j) ! background gas
      wrk%densities(dat%nsp+1,j) = 1.0_dp ! for hv
    enddo
  
    if (sp_ind <= dat%nq) then ! long lived or particle
      do k = dat%nq+1,dat%nq+dat%nsl
        call chempl_sl(self%dat, self%var, wrk%densities, wrk%rx_rates, &
                       k, xp, xl) 
        wrk%densities(k,:) = xp/xl
      enddo
    endif
    
    call chempl_t(self%dat, self%var, &
                  wrk%densities, wrk%rx_rates, sp_ind, pl%production, pl%loss)
  
    do i = 1,np
      pl%integrated_production(i) = sum(pl%production(:,i)*var%dz)
      k = dat%pl(sp_ind)%iprod(i) ! reaction number
      pl%production_rx(i) = dat%reaction_equations(k)
    enddo
    do i = 1,nl
      pl%integrated_loss(i) = sum(pl%loss(:,i)*var%dz)
      k = dat%pl(sp_ind)%iloss(i) ! reaction number
      pl%loss_rx(i) = dat%reaction_equations(k)
    enddo
    
    ! rainout
    pl%loss_rx(nl+1) = "rainout"
    pl%loss(:,nl+1) = 0.0_dp
    pl%integrated_loss(nl+1) = 0.0_dp
    if (dat%gas_rainout .and. sp_ind <= dat%nq) then
      pl%loss(1:var%trop_ind,nl+1) = &
          wrk%rainout_rates(sp_ind,1:var%trop_ind)*wrk%usol(sp_ind,1:var%trop_ind)*wrk%density(1:var%trop_ind)
      pl%integrated_loss(nl+1) = sum(pl%loss(:,nl+1)*var%dz)
    endif
    
    ! ignoring condensation, and fluxes
    
    ! sort 
    prod_inds = argsort(pl%integrated_production)
    loss_inds = argsort(pl%integrated_loss)
    prod_inds = prod_inds(np:1:-1)
    loss_inds = loss_inds(nlT:1:-1)
  
    pl%integrated_production = pl%integrated_production(prod_inds)
    pl%integrated_loss = pl%integrated_loss(loss_inds)
    
    pl%production = pl%production(:,prod_inds)
    pl%loss = pl%loss(:,loss_inds)
    
    pl%production_rx = pl%production_rx(prod_inds)
    pl%loss_rx = pl%loss_rx(loss_inds)
    
  end subroutine
  
end submodule