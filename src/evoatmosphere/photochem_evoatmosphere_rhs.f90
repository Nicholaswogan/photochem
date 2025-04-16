#:set TYPES = ['real(dp)', 'type(dual)']
#:set NAMES = ['real', 'dual']
#:set TYPES_NAMES = list(zip(TYPES, NAMES))

submodule(photochem_evoatmosphere) photochem_evoatmosphere_rhs
  implicit none
  
  interface dochem
    module procedure :: dochem_real, dochem_dual
  end interface

  interface
    module subroutine equilibrium_climate(self, usol_den, molecules_per_particle, T_trop, T_surf_guess, &
                                          T_surf, T, z_trop, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), target, intent(in) :: usol_den(:,:)
      real(dp), intent(in) :: molecules_per_particle(:,:)
      real(dp), target, intent(in) :: T_trop, T_surf_guess
      real(dp), target, intent(out) :: T_surf, T(:)
      real(dp), target, intent(out) :: z_trop
      character(:), allocatable, intent(out) :: err
    end subroutine
  end interface

contains

  #:for TYPE1, NAME in TYPES_NAMES
  subroutine dochem_${NAME}$(self, usol, rx_rates, &
                    gas_sat_den, molecules_per_particle, &
                    H2O_sat_mix, H2O_rh, rainout_rates, scale_height, wfall, &
                    density, mix, densities, xp, xl, rhs)                 
    use photochem_enum, only: CondensingParticle
    use photochem_common, only: chempl, chempl_sl
    use photochem_eqns, only: damp_condensation_rate
    use photochem_const, only: N_avo, pi, small_real, T_crit_H2O
    #:if NAME == 'dual'
    use differentia
    #:endif
    class(EvoAtmosphere), target, intent(in) :: self
    ${TYPE1}$, intent(in) :: usol(:,:)
    real(dp), intent(in) :: rx_rates(:,:)
    real(dp), intent(in) :: gas_sat_den(:,:)
    real(dp), intent(in) :: molecules_per_particle(:,:)
    real(dp), intent(in) :: H2O_sat_mix(:), H2O_rh(:)
    real(dp), intent(in) :: rainout_rates(:,:), scale_height(:), wfall(:,:)
    real(dp), intent(in) :: density(:)
    ${TYPE1}$, intent(inout) :: mix(:,:)
    ${TYPE1}$, intent(inout) :: densities(:,:), xp(:), xl(:)
    ${TYPE1}$, intent(inout) :: rhs(:) ! neqs

    real(dp) :: cond_rate0
    ${TYPE1}$ :: rh, dn_gas_dt, dn_particle_dt, cond_rate
    integer :: i, ii, j, k, kk

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var

    #:if NAME == 'dual'
    rh = dual(size(usol(1,1)%der))
    dn_gas_dt = dual(size(usol(1,1)%der))
    dn_particle_dt = dual(size(usol(1,1)%der))
    cond_rate = dual(size(usol(1,1)%der))
    #:endif
    do j = 1,var%nz
      mix(:,j) = usol(:,j)/density(j)
    enddo

    do j = 1,var%nz
      do i = 1,dat%npq
        densities(i,j) = max(usol(i,j)*(1.0_dp/molecules_per_particle(i,j)), small_real)
      enddo
      do i = dat%ng_1,dat%nq
        densities(i,j) = usol(i,j)
      enddo
      densities(dat%nsp+1,j) = 1.0_dp ! for hv
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
        rhs(k) = xp(j) - xl(j)
      enddo
    enddo

    if (dat%gas_rainout) then
      ! rainout rates
      do j = 1,var%trop_ind
        do i = 1,dat%nq
          k = i + (j - 1) * dat%nq
          rhs(k) = rhs(k) - rainout_rates(i,j)*usol(i,j)
        enddo
      enddo
    endif

    if (dat%fix_water_in_trop) then
      do j = 1,var%trop_ind
        k = dat%LH2O + (j - 1) * dat%nq
        rhs(k) = rhs(k) + var%fast_arbitrary_rate*(density(j)*H2O_sat_mix(j)*H2O_rh(j) - usol(dat%LH2O,j))
      enddo
    endif
    if (dat%water_cond) then
      if (dat%fix_water_in_trop) then
        i = var%trop_ind+1
      else
        i = 1
      endif
      do j = i,var%nz
        if (var%temperature(j) < T_crit_H2O) then
          ! water will condense if it is below the critical point.

          k = dat%LH2O + (j - 1) * dat%nq ! gas phase rhs index

          ! compute the relative humidity
          rh = max(mix(dat%LH2O,j)/H2O_sat_mix(j),small_real)

          if (rh > var%H2O_cond_params%RHc) then
            
            cond_rate0 = var%H2O_cond_params%k_cond*(var%edd(j)/scale_height(j)**2.0_dp)
            cond_rate = damp_condensation_rate(cond_rate0, &
                                               var%H2O_cond_params%RHc, &
                                               (1.0_dp + var%H2O_cond_params%smooth_factor)*var%H2O_cond_params%RHc, &
                                               rh)
            
            ! Rate H2O gas is destroyed (molecules/cm^3/s)
            dn_gas_dt = - cond_rate*usol(dat%LH2O,j)
            rhs(k) = rhs(k) + dn_gas_dt
            
          endif
        endif
      enddo
    endif

    if (dat%there_are_particles) then
      ! formation from reaction
      do i = 1,dat%np
        call chempl(self%dat, self%var, densities, rx_rates, i, xp, xl)
        do j = 1,var%nz
          k = i + (j - 1) * dat%nq
          rhs(k) = rhs(k) + (xp(j) - xl(j))
        enddo
      enddo

      ! 4 numbers for each condensing species.
      ! condensation rate coefficient ~ 100.0
      ! evaporation rate coefficient ~ 10.0
      ! RH of condensation ~ 1.0
      ! RH smoothing factor ~ 0.1
    
      ! particle condensation
      do j = 1,var%nz
        do i = 1,dat%np
          ! if this particle forms from condensation
          if (dat%particle_formation_method(i) == CondensingParticle) then
            ii = dat%particle_gas_phase_ind(i) ! index of gas phase
            kk = ii + (j - 1) * dat%nq ! gas phase rhs index
            k = i + (j - 1) * dat%nq ! particle rhs index

            ! compute the relative humidity
            rh = max(usol(ii,j)/gas_sat_den(i,j),small_real)

            if (rh > var%cond_params(i)%RHc) then
              ! Condensation occurs

              ! Condensation rate is based on how rapidly gases are mixing
              ! in the atmosphere. We also smooth the rate to prevent stiffness.
              cond_rate0 = var%cond_params(i)%k_cond*(var%edd(j)/scale_height(j)**2.0_dp)
              cond_rate = damp_condensation_rate(cond_rate0, &
                                                 var%cond_params(i)%RHc, &
                                                 (1.0_dp + var%cond_params(i)%smooth_factor)*var%cond_params(i)%RHc, &
                                                 rh)
            
              ! Rate gases are being destroyed (molecules/cm^3/s)
              dn_gas_dt = - cond_rate*usol(ii,j)
              rhs(kk) = rhs(kk) + dn_gas_dt          
            
              ! Rate particles are being produced (molecules/cm^3/s)
              dn_particle_dt = - dn_gas_dt
              rhs(k) = rhs(k) + dn_particle_dt

            elseif (rh <= var%cond_params(i)%RHc .and. var%evaporation) then
              ! Evaporation occurs

              ! Evaporation rate is based on how rapidly particles are falling
              ! in the atmosphere. We also smooth the rate to prevent stiffness.
              cond_rate0 = var%cond_params(i)%k_evap*(wfall(i,j)/scale_height(j))
              cond_rate = damp_condensation_rate(cond_rate0, &
                                                 1.0_dp/var%cond_params(i)%RHc, &
                                                 (1.0_dp + var%cond_params(i)%smooth_factor)/var%cond_params(i)%RHc, &
                                                 1.0_dp/rh)

              ! Rate gases are being produced (molecules/cm^3/s)                   
              dn_gas_dt = cond_rate*usol(i,j)
              rhs(kk) = rhs(kk) + dn_gas_dt

              ! Rate particles are being destroyed (molecules/cm^3/s)
              dn_particle_dt = - dn_gas_dt
              rhs(k) = rhs(k) + dn_particle_dt

            endif
          endif          
        enddo
      enddo
      
    endif

  end subroutine

  #:endfor
  subroutine diffusion_coefficients_evo(dat, var, den, mubar, &
                                    DU, DD, DL, ADU, ADL, ADD, wfall, VH2_esc, VH_esc)
    use photochem_eqns, only: dynamic_viscosity_air, fall_velocity, slip_correction_factor, &
                              default_binary_diffusion_param
    use photochem_const, only: k_boltz, N_avo
    use photochem_enum, only: DiffusionLimHydrogenEscape
    use photochem_types, only: binary_diffusion_fcn

    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var

    real(dp), intent(in) :: den(:), mubar(:)

    real(dp), intent(out) :: DU(:,:), DL(:,:), DD(:,:) ! (nq,nz)
    real(dp), intent(out) :: ADU(:,:), ADL(:,:), ADD(:,:) ! (nq,nz)
    real(dp), intent(out) :: wfall(:,:) ! (npq,nz)
    real(dp), intent(out) :: VH2_esc, VH_esc

    real(dp) :: eddav(var%nz-1), denav(var%nz-1)
    real(dp) :: tav(var%nz-1), dTdz(var%nz-1)
    real(dp) :: scale_height_av(var%nz-1),scale_height_i_av
    real(dp) :: b1x2av(dat%nll,var%nz-1)
    real(dp) :: gamma_i_gas_av(dat%nll,var%nz-1), gamma_i_part_av(dat%np,var%nz-1)
    real(dp) :: gamma_i_gas_D(dat%nll,var%nz)
    real(dp) :: grav_av, mubar_av
    real(dp) :: bx1x2

    ! for particles
    real(dp) :: air_density_pp, air_density
    real(dp) :: wfall_pp, wfall_i
    real(dp) :: viscosity_pp, viscosity
    real(dp) :: FF2, FF1

    ! molecular diffusion function
    procedure(binary_diffusion_fcn), pointer :: binary_diffusion_param
    
    integer :: j, i, k

    if (associated(var%custom_binary_diffusion_fcn)) then
      binary_diffusion_param => var%custom_binary_diffusion_fcn
    else
      binary_diffusion_param => default_binary_diffusion_param
    endif

    if (var%upwind_molec_diff) then; block
    real(dp) :: dTdz_tmp, scale_height_i_tmp, b1x2_tmp, grav_tmp, mubar_tmp
    ! If using upwind scheme for molecular diffusion, then compute the needed
    ! advective velocities.
    do j = 1,var%nz-1
      dTdz_tmp = (var%temperature(j+1) - var%temperature(j))/var%dz(j)
      grav_tmp = var%grav(j)
      mubar_tmp = mubar(j)
      do i = dat%ng_1,dat%nq
        k = i - dat%np
        scale_height_i_tmp = (N_avo*k_boltz*var%temperature(j))/(grav_tmp*dat%species_mass(i))
        b1x2_tmp = binary_diffusion_param(dat%species_mass(i), mubar_tmp, var%temperature(j))
        gamma_i_gas_D(k,j) = (b1x2_tmp/den(j))*(1.0_dp/scale_height_i_tmp + (1.0_dp/var%temperature(j))*dTdz_tmp)
      enddo
    enddo
    j = var%nz
    dTdz_tmp = (var%temperature(var%nz) - var%temperature(var%nz-1))/var%dz(var%nz)
    grav_tmp = var%grav(j)
    mubar_tmp = mubar(j)
    do i = dat%ng_1,dat%nq
      k = i - dat%np
      scale_height_i_tmp = (N_avo*k_boltz*var%temperature(j))/(grav_tmp*dat%species_mass(i))
      b1x2_tmp = binary_diffusion_param(dat%species_mass(i), mubar_tmp, var%temperature(j))
      gamma_i_gas_D(k,j) = (b1x2_tmp/den(j))*(1.0_dp/scale_height_i_tmp + (1.0_dp/var%temperature(j))*dTdz_tmp)     
    enddo
    endblock; endif

    ! compute relevant parameters at the edges of the grid cells
    do j = 1,var%nz-1
      eddav(j) = sqrt(var%edd(j)*var%edd(j+1))
      denav(j) = sqrt(den(j)*den(j+1))
      tav(j) = sqrt(var%temperature(j)*var%temperature(j+1))
      dTdz(j) = (var%temperature(j+1) - var%temperature(j))/var%dz(j)
      grav_av = sqrt(var%grav(j)*var%grav(j+1))
      mubar_av = sqrt(mubar(j)*mubar(j+1))
      scale_height_av(j) = (N_avo*k_boltz*tav(j))/(grav_av*mubar_av)
      do i = 1,dat%np
        gamma_i_part_av(i,j) = eddav(j)*(1.0_dp/scale_height_av(j) + (1.0_dp/tav(j))*dTdz(j))
      enddo
      do i = dat%ng_1,dat%nq
        k = i - dat%np
        scale_height_i_av = (N_avo*k_boltz*tav(j))/(grav_av*dat%species_mass(i))
        b1x2av(k,j) = binary_diffusion_param(dat%species_mass(i), mubar_av, tav(j))

        gamma_i_gas_av(k,j) = eddav(j)*(1.0_dp/scale_height_av(j) + (1.0_dp/tav(j))*dTdz(j))

        ! If we don't use upwind scheme for molecular diffusion, then we add the centered
        ! molecular diffusion terms here.
        if (.not.var%upwind_molec_diff) then
          gamma_i_gas_av(k,j) = gamma_i_gas_av(k,j) + &
            (b1x2av(k,j)/denav(j))*(1.0_dp/scale_height_i_av + (1.0_dp/tav(j))*dTdz(j))
        endif

      enddo
    enddo

    ! gases
    ! middle
    do j = 2,var%nz-1
      do i = dat%ng_1,dat%nq
        k = i - dat%np
        ! diffusion
        DU(i,j) = (eddav(j) + (b1x2av(k,j)/denav(j)))/(var%dz(j)**2.0_dp)
        DL(i,j) = (eddav(j-1) + (b1x2av(k,j-1)/denav(j-1)))/(var%dz(j)**2.0_dp)
        DD(i,j) = - DU(i,j) - DL(i,j)

        ! advection
        ADU(i,j) = gamma_i_gas_av(k,j)/(2.0_dp*var%dz(j))
        ADL(i,j) = - gamma_i_gas_av(k,j-1)/(2.0_dp*var%dz(j))
        ADD(i,j) = ADU(i,j) + ADL(i,j)

        if (var%upwind_molec_diff) then
          ADU(i,j) = ADU(i,j) + gamma_i_gas_D(k,j+1)/var%dz(j)
          ADL(i,j) = ADL(i,j) + 0.0_dp
          ADD(i,j) = ADD(i,j) + (- gamma_i_gas_D(k,j)/var%dz(j))
        endif
      enddo
    enddo
    ! lower boundary
    j = 1
    do i = dat%ng_1,dat%nq
      k = i - dat%np
      DU(i,j) = (eddav(j) + (b1x2av(k,j)/denav(j)))/(var%dz(j)**2.0_dp)
      DD(i,j) = - DU(i,j)

      ADU(i,j) = gamma_i_gas_av(k,j)/(2.0_dp*var%dz(j))
      ADD(i,j) = ADU(i,j)

      if (var%upwind_molec_diff) then
        ADU(i,j) = ADU(i,j) + gamma_i_gas_D(k,j+1)/var%dz(j)
        ADD(i,j) = ADD(i,j) + 0.0_dp
      endif
    enddo
    ! upper boundary
    j = var%nz
    do i = dat%ng_1,dat%nq
      k = i - dat%np
      DL(i,j) = (eddav(j-1) + (b1x2av(k,j-1)/denav(j-1)))/(var%dz(j)**2.0_dp)
      DD(i,j) = - DL(i,j)

      ADL(i,j) = - gamma_i_gas_av(k,j-1)/(2.0_dp*var%dz(j))
      ADD(i,j) = ADL(i,j)

      if (var%upwind_molec_diff) then
        ADL(i,j) = ADL(i,j) + 0.0_dp
        ADD(i,j) = ADD(i,j) + (- gamma_i_gas_D(k,j)/var%dz(j))
      endif
    enddo
    
    ! ! particles (eddy diffusion)
    ! ! middle
    ! do j = 2,var%nz-1
    !   do i = 1,dat%np
    !     ! diffusion
    !     DU(i,j) = eddav(j)/(var%dz(j)**2.0_dp)
    !     DL(i,j) = eddav(j-1)/(var%dz(j)**2.0_dp)
    !     DD(i,j) = - DU(i,j) - DL(i,j)

    !     ! advection
    !     ADU(i,j) = gamma_i_part_av(i,j)/(2.0_dp*var%dz(j))
    !     ADL(i,j) = - gamma_i_part_av(i,j-1)/(2.0_dp*var%dz(j))
    !     ADD(i,j) = ADU(i,j) + ADL(i,j)
    !   enddo
    ! enddo
    ! ! lower boundary
    ! j = 1
    ! do i = 1,dat%np
    !   DU(i,j) = eddav(j)/(var%dz(j)**2.0_dp)
    !   DD(i,j) = - DU(i,j)

    !   ADU(i,j) = gamma_i_part_av(i,j)/(2.0_dp*var%dz(j))
    !   ADD(i,j) = ADU(i,j)
    ! enddo
    ! ! upper boundary
    ! j = var%nz
    ! do i = 1,dat%np
    !   DL(i,j) = eddav(j-1)/(var%dz(j)**2.0_dp)
    !   DD(i,j) = - DL(i,j)

    !   ADL(i,j) = - gamma_i_part_av(i,j-1)/(2.0_dp*var%dz(j))
    !   ADD(i,j) = ADL(i,j)
    ! enddo

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

    ! particles (falling)
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

        FF2 = wfall_pp/var%dz(i)

        FF1 = -wfall_i/var%dz(i)
      
        ADU(j,i) = ADU(j,i) + FF2
        ADD(j,i) = ADD(j,i) + FF1
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

      FF2 = wfall_pp/var%dz(i)

      ADU(j,i) = ADU(j,i) + FF2
    enddo
    ! Upper boundary
    i = var%nz
    do j = 1,dat%npq
      air_density = (den(i)/N_avo)*mubar(i)
      viscosity = dynamic_viscosity_air(var%temperature(i))
      wfall_i = fall_velocity(var%grav(i), var%particle_radius(j,i), &
                              dat%particle_density(j), air_density, viscosity) &
                  *slip_correction_factor(var%particle_radius(j,i), den(i))
      FF1 = -wfall_i/var%dz(i)

      ADD(j,i) = ADD(j,i) + FF1
    enddo

    ! not going to work. Advection is required to 
    ! maintain hydrostatic equilibrium
    ! do i = 1,dat%nq
    !   if (var%only_eddy(i)) then
    !     ADL(i,:) = 0.0_dp
    !     ADU(i,:) = 0.0_dp
    !     ADD(i,:) = 0.0_dp
    !   endif
    ! enddo
    
    ! H2 escape
    if (dat%H_escape_type == DiffusionLimHydrogenEscape) then
      bx1x2 = binary_diffusion_param(dat%species_mass(dat%LH2), mubar(var%nz), var%temperature(var%nz))
      VH2_esc = bx1x2/den(var%nz)*(-(dat%species_mass(dat%LH2)*var%grav(var%nz))/(k_boltz*var%temperature(var%nz)*N_avo) &
                                + (mubar(var%nz)*var%grav(var%nz))/(k_boltz*var%temperature(var%nz)*N_avo))                     

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

  module subroutine prep_atm_evo_gas(self, usol_in, usol, &
                                     molecules_per_particle, pressure, density, mix, mubar, &
                                     pressure_hydro, density_hydro, err)
    use photochem_eqns, only: press_and_den
    use photochem_common, only: molec_per_particle
    use photochem_const, only: small_real, N_avo, k_boltz
    use photochem_enum, only: DensityBC, PressureBC
    use photochem_input, only: compute_gibbs_energy, interp2xsdata
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    real(dp), intent(out) :: usol(:,:)
    real(dp), intent(out) :: molecules_per_particle(:,:)
    real(dp), intent(out) :: pressure(:), density(:), mix(:,:), mubar(:)
    real(dp), intent(out) :: pressure_hydro(:), density_hydro(:)
    character(:), allocatable, intent(out) :: err

    real(dp) :: T_surf_guess, Psat
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
      if (var%lowerboundcond(i) == DensityBC) then
        usol(i,1) = var%lower_fix_den(i)
      elseif (var%lowerboundcond(i) == PressureBC) then
        Psat = huge(1.0_dp)
        if (dat%gas_particle_ind(i) /= 0) then
          j = dat%gas_particle_ind(i)
          Psat = dat%particle_sat(j)%sat_pressure(var%temperature(1))*var%cond_params(j)%RHc
        endif
        usol(i,1) = min(var%lower_fix_press(i), Psat)/(k_boltz*var%temperature(1))
      endif
    enddo

    !!! molecules/particle
    if (dat%there_are_particles) then
      call molec_per_particle(dat, var, molecules_per_particle)
    endif

    !!! climate model
    if (self%evolve_climate) then
      ! This block changes the following variables, which normally do not change
      ! during integration.
      ! self%T_surf
      ! var%temperature
      ! var%trop_alt
      ! var%xs_x_qy
      ! var%gibbs_energy

      T_surf_guess = self%T_surf
      ! update the temperature profile
      call equilibrium_climate(self, usol, molecules_per_particle, self%T_trop, T_surf_guess, &
                               self%T_surf, var%temperature, var%trop_alt, err)
      if (allocated(err)) return

      ! Update all things that depend on temperature

      ! Need to interpolate xsections, but more work is needed
      call interp2xsdata(dat, var, err)
      if (allocated(err)) return

      if (dat%reverse) then
        call compute_gibbs_energy(dat, var, err)
        if (allocated(err)) return
      endif

    endif

    !!! pressure, density and mean molcular weight
    do j = 1,var%nz
      density(j) = sum(usol(dat%ng_1:,j))
      mix(:,j) = usol(:,j)/density(j) ! mixing ratios
      mubar(j) = sum(dat%species_mass(dat%ng_1:dat%nq)*mix(dat%ng_1:,j))
    enddo
    
    ! surface pressure by adding up all the mass in the atmosphere (bars)
    var%surface_pressure = sum(density(:)*mubar(:)*var%grav(:)*var%dz(:))/N_avo/1.0e6_dp
    call press_and_den(var%nz, var%temperature, var%grav, var%surface_pressure*1.0e6_dp, var%dz, &
                       mubar, pressure_hydro, density_hydro)
    pressure(:) = density(:)*k_boltz*var%temperature(:)

  end subroutine

  module subroutine set_trop_ind(self, usol_in, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    character(:), allocatable, intent(out) :: err

    type(PhotochemWrkEvo), pointer :: wrk

    wrk => self%wrk

    call prep_atm_evo_gas(self, usol_in, wrk%usol, &
                          wrk%molecules_per_particle, wrk%pressure, wrk%density, wrk%mix, wrk%mubar, &
                          wrk%pressure_hydro, wrk%density_hydro, err)
    if (allocated(err)) return

    if (self%dat%fix_water_in_trop .or. self%dat%gas_rainout) then
      self%var%trop_ind = max(minloc(abs(self%var%z - self%var%trop_alt), 1) - 1, 1)
      
      if (self%var%trop_ind < 3) then
        err = 'Tropopause is too low.'
        return
      elseif (self%var%trop_ind > self%var%nz-2) then
        err = 'Tropopause is too high.'
        return
      endif
    else
      self%var%trop_ind = 1
    endif

  end subroutine

  module subroutine prep_all_evo_gas(self, usol_in, err)

    use photochem_common, only: reaction_rates, rainout, photorates
    use photochem_common, only: gas_saturation_density
    use clima_eqns_water, only: sat_pressure_H2O
    use photochem_const, only: pi, N_avo, small_real, k_boltz
    use photochem_enum, only: DiffusionLimHydrogenEscape

    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    integer :: i, j, k

    dat => self%dat
    var => self%var
    wrk => self%wrk

    call prep_atm_evo_gas(self, usol_in, wrk%usol, &
                          wrk%molecules_per_particle, wrk%pressure, wrk%density, wrk%mix, wrk%mubar, &
                          wrk%pressure_hydro, wrk%density_hydro, err)
    if (allocated(err)) return

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
    call diffusion_coefficients_evo(dat, var, wrk%density, wrk%mubar, &
    wrk%DU, wrk%DD, wrk%DL, wrk%ADU, wrk%ADL, wrk%ADD, wrk%wfall, wrk%VH2_esc, wrk%VH_esc)
    
    wrk%scale_height = (k_boltz*var%temperature(:)*N_avo)/(wrk%mubar(:)*var%grav(:))

    !!! H and H2 escape
    wrk%upper_veff_copy = var%upper_veff
    wrk%lower_vdep_copy = var%lower_vdep
    if (dat%H_escape_type == DiffusionLimHydrogenEscape) then
      wrk%upper_veff_copy(dat%LH2) = wrk%VH2_esc                     
      wrk%upper_veff_copy(dat%LH) = wrk%VH_esc 
    endif

    !!! Particle lower boundary, and saturation properties
    if (dat%there_are_particles) then
      do i = 1,dat%np
        ! Here we impose a lower boundary condition for particles. They fall out
        ! of the model according to the fall velocity.
        wrk%lower_vdep_copy(i) = wrk%lower_vdep_copy(i) + wrk%wfall(i,1)
      enddo

      call gas_saturation_density(dat, var, wrk%gas_sat_den)
    endif

    !!! densities
    do j = 1,var%nz
      do i = 1,dat%npq
        wrk%densities(i,j) = max(wrk%usol(i,j)*(1.0_dp/wrk%molecules_per_particle(i,j)), small_real)
      enddo
      do i = dat%ng_1,dat%nq
        wrk%densities(i,j) = wrk%usol(i,j)
      enddo
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
                   wrk%mix(dat%LH2O,:), wrk%density, wrk%rainout_rates)
    endif

  end subroutine

  module subroutine rhs_evo_gas(self, neqs, tn, usol_flat, rhs, err)
    use photochem_enum, only: MosesBC, VelocityBC, DensityBC, PressureBC, FluxBC, VelocityDistributedFluxBC
    use photochem_enum, only: ZahnleHydrogenEscape
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real  
    
    class(EvoAtmosphere), target, intent(inout) :: self
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
    type(PhotochemWrkEvo), pointer :: wrk
    
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
    call prep_all_evo_gas(self, usol_in, err)
    if (allocated(err)) return

    call dochem(self, wrk%usol, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                wrk%density, wrk%mix, wrk%densities, wrk%xp, wrk%xl, rhs)  

    ! Extra functions specifying production or destruction
    do i = 1,dat%nq
      if (associated(var%rate_fcns(i)%fcn)) then
        call var%rate_fcns(i)%fcn(tn, var%nz, wrk%xp) ! using wrk%xp space.
        do j = 1,var%nz
          k = i + (j-1)*dat%nq
          rhs(k) = rhs(k) + wrk%xp(j) ! (molecules/cm^3/s)
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
      elseif (var%lowerboundcond(i) == DensityBC .or. &
              var%lowerboundcond(i) == PressureBC) then
        rhs(i) = 0.0_dp
      elseif (var%lowerboundcond(i) == FluxBC) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        + wrk%DD(i,1)*wrk%usol(i,1) + wrk%ADD(i,1)*wrk%usol(i,1) &
                        + var%lower_flux(i)/var%dz(1)
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
                        - var%upper_flux(i)/var%dz(var%nz)
      endif
    enddo

    ! Distributed (volcanic) sources
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityDistributedFluxBC) then
        disth = var%lower_dist_height(i)*1.e5_dp        
        if (disth < var%z(1) - 0.5_dp*var%dz(1)) then
        ! If the height is below the model domain, then we will put all flux into
        ! lowest layer.
        rhs(i) = rhs(i) + var%lower_flux(i)/var%dz(1)
        else
        ! If the height is within the model domain, then we will distribute the flux
        ! throught the model.
        jdisth = minloc(var%Z,1, var%Z >= disth) - 1
        jdisth = max(jdisth,2)
        ztop = var%z(jdisth)-var%z(1)
        ztop1 = var%z(jdisth) + 0.5e0_dp*var%dz(jdisth)
        do j = 2,jdisth
          k = i + (j-1)*dat%nq
          rhs(k) = rhs(k) + 2.0_dp*var%lower_flux(i)*(ztop1-var%z(j))/(ztop**2.0_dp)
        enddo
        endif
      endif
    enddo 

    ! zahnle hydrogen escape
    if (dat%H_escape_type == ZahnleHydrogenEscape) then

      ! for Zahnle hydrogen escape, we pull H2 out of 
      ! the bottom grid cell of the model.

      rhs(dat%LH2) = rhs(dat%LH2) &
      - dat%H_escape_coeff*wrk%mix(dat%LH2,1)/var%dz(1)
      
    endif

  end subroutine

  subroutine autodiff_chemistry_jacobian(self, usol, rhs, djac, err)
    use differentia, only: jacobian, dual, BlockDiagonalJacobian, initialize_dual_array
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), target, contiguous, intent(in) :: usol(:,:)
    real(dp), intent(out) :: rhs(:), djac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), pointer :: usol_flat(:)
    type(dual), allocatable :: mix(:,:), densities(:,:), xp(:), xl(:)
    integer :: blocksize

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk

    dat => self%dat
    var => self%var
    wrk => self%wrk

    usol_flat(1:var%neqs) => usol(:,:) ! reshape input
    blocksize = dat%nq
    ! allocate some work memory
    allocate(mix(dat%nq,var%nz))
    call initialize_dual_array(mix, blocksize)
    allocate(densities(size(wrk%densities,1),size(wrk%densities,2)))
    call initialize_dual_array(densities, blocksize)
    allocate(xp(size(wrk%xp)))
    call initialize_dual_array(xp, blocksize)
    allocate(xl(size(wrk%xl)))
    call initialize_dual_array(xl, blocksize)

    call jacobian(fcn, usol_flat, rhs, djac, jt=BlockDiagonalJacobian, blocksize=blocksize, err=err)
    if (allocated(err)) return

  contains
    subroutine fcn(x_, f_)
      type(dual), target, intent(in) :: x_(:)
      type(dual), target, intent(out) :: f_(:)
      type(dual), pointer :: usol_(:,:)
      usol_(1:dat%nq,1:var%nz) => x_(:)
      call initialize_dual_array(f_, blocksize)
      call dochem(self, usol_, wrk%rx_rates, &
                  wrk%gas_sat_den, wrk%molecules_per_particle, &
                  wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                  wrk%density, mix, densities, xp, xl, f_) 
    end subroutine
  end subroutine

  module subroutine jac_evo_gas(self, lda_neqs, neqs, usol_flat, jac, err)
    use photochem_enum, only: MosesBC, VelocityBC, DensityBC, PressureBC, FluxBC, VelocityDistributedFluxBC
    use photochem_enum, only: ZahnleHydrogenEscape
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real
    
    class(EvoAtmosphere), target, intent(inout) :: self
    integer, intent(in) :: lda_neqs, neqs
    real(dp), target, intent(in) :: usol_flat(neqs)
    real(dp), intent(out), target :: jac(lda_neqs)
    character(:), allocatable, intent(out) :: err
    
    real(dp), pointer :: usol_in(:,:)
    real(dp), pointer :: djac(:,:)
    real(dp) :: rhs(self%var%neqs)

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk

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
  
    call prep_all_evo_gas(self, usol_in, err)
    if (allocated(err)) return
  
    jac = 0.0_dp
  
    if (.not. var%autodiff) then; block
    real(dp) :: usol_perturb(dat%nq,var%nz)
    real(dp) :: R(var%nz)
    real(dp) :: rhs_perturb(var%neqs)
    real(dp) :: mix(dat%nq,var%nz)
    real(dp) :: densities(dat%nsp+1,self%var%nz), xl(var%nz), xp(var%nz)
  
    ! Finite differenced Jacobian

    ! compute chemistry contribution to jacobian using forward differences
    call dochem(self, wrk%usol, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                wrk%density, mix, densities, xp, xl, rhs) 

    !$omp parallel private(i, j, k, m, mm, usol_perturb, R, mix, densities, xl, xp, rhs_perturb)
    usol_perturb = wrk%usol
    !$omp do
    do i = 1,dat%nq
      do j = 1,var%nz
        R(j) = var%epsj*abs(wrk%usol(i,j))
        usol_perturb(i,j) = wrk%usol(i,j) + R(j)
      enddo
      
      call dochem(self, usol_perturb, wrk%rx_rates, &
                  wrk%gas_sat_den, wrk%molecules_per_particle, &
                  wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                  wrk%density, mix, densities, xp, xl, rhs_perturb) 
  
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
      elseif (var%lowerboundcond(i) == DensityBC .or. &
              var%lowerboundcond(i) == PressureBC) then

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
      endif
    enddo
  
    ! Upper boundary
    do i = 1,dat%nq
      k = i + (var%nz-1)*dat%nq
      if (var%upperboundcond(i) == VelocityBC) then
  
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

      djac(dat%kd,dat%LH2) = djac(dat%kd,dat%LH2) &
      - (dat%H_escape_coeff/var%dz(1)) &
        *((wrk%density(1) - wrk%usol(dat%LH2,1))/wrk%density(1)**2.0_dp)

      do m = dat%ng_1,dat%nq 
        if (m /= dat%LH2) then
          mm = dat%kd + dat%LH2 - m
          djac(mm,m) = djac(mm,m) &
          - (dat%H_escape_coeff/var%dz(1)) &
            *(-wrk%usol(dat%LH2,1)/wrk%density(1)**2.0_dp)
        endif
      enddo
    endif
  
  end subroutine

  module subroutine right_hand_side_chem(self, usol, rhs, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(out) :: rhs(:)
    character(:), allocatable, intent(out) :: err
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz .or. size(rhs) /= var%neqs) then
      err = "Input usol or rhs has the wrong dimensions"
      return
    endif

    call self%prep_atmosphere(usol, err)
    if (allocated(err)) return
    
    call dochem(self, wrk%usol, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                wrk%density, wrk%mix, wrk%densities, wrk%xp, wrk%xl, rhs) 
                              
  end subroutine

  module subroutine production_and_loss(self, species, usol, pl, err)     
    use futils, only: argsort            
    use photochem_common, only: chempl_sl, chempl_t
    use photochem_types, only: ProductionLoss
    use photochem_const, only: small_real
  
    class(EvoAtmosphere), target, intent(inout) :: self
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
    type(PhotochemWrkEvo), pointer :: wrk
  
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
      do i = 1,dat%npq
        wrk%densities(i,j) = max(wrk%usol(i,j)*(1.0_dp/wrk%molecules_per_particle(i,j)), small_real)
      enddo
      do i = dat%ng_1,dat%nq
        wrk%densities(i,j) = wrk%usol(i,j)
      enddo
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
          wrk%rainout_rates(sp_ind,1:var%trop_ind)*wrk%usol(sp_ind,1:var%trop_ind)
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


