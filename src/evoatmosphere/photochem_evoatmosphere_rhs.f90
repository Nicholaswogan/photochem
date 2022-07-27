
submodule(photochem_evoatmosphere) photochem_evoatmosphere_rhs
  implicit none

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
  subroutine dochem(self, usol, rx_rates, &
                    gas_sat_den, molecules_per_particle, &
                    H2O_sat_mix, H2O_rh, rainout_rates, &
                    density, mix, densities, xp, xl, rhs)                 
    use photochem_enum, only: CondensingParticle
    use photochem_common, only: chempl, chempl_sl
    use photochem_eqns, only: damp_condensation_rate
    use photochem_const, only: N_avo, pi, small_real, T_crit_H2O, fast_arbitrary_rate

    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: rx_rates(:,:)
    real(dp), intent(in) :: gas_sat_den(:,:)
    real(dp), intent(in) :: molecules_per_particle(:,:)
    real(dp), intent(in) :: H2O_sat_mix(:), H2O_rh(:)
    real(dp), intent(in) :: rainout_rates(:,:)
    real(dp), intent(inout) :: density(:), mix(:,:), densities(:,:), xp(:), xl(:)
    real(dp), intent(out) :: rhs(:) ! neqs

    real(dp) :: dn_gas_dt, dn_particle_dt, H2O_cold_trap, cond_rate
    integer :: i, ii, j, k, kk

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var

    do j = 1,var%nz
      density(j) = sum(usol(dat%ng_1:,j)) 
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
        rhs(k) = rhs(k) + fast_arbitrary_rate*(density(j)*H2O_sat_mix(j)*H2O_rh(j) - usol(dat%LH2O,j))
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

          k = dat%LH2O + (j - 1) * dat%nq
          H2O_cold_trap = var%H2O_condensation_rate(2)*H2O_sat_mix(j)
          if (mix(dat%LH2O,j) >= H2O_cold_trap) then
            
            cond_rate = damp_condensation_rate(var%H2O_condensation_rate(1), &
                                              var%H2O_condensation_rate(2), &
                                              var%H2O_condensation_rate(3), &
                                              mix(dat%LH2O,j)/H2O_sat_mix(j))
            rhs(k) = rhs(k) - cond_rate*(mix(dat%LH2O,j) - H2O_cold_trap)*density(j)
            
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
    
      ! particle condensation
      do j = 1,var%nz
        do i = 1,dat%np
          ! if this particle forms from condensation
          if (dat%particle_formation_method(i) == CondensingParticle) then
            ii = dat%particle_gas_phase_ind(i) ! index of gas phase
            kk = ii + (j - 1) * dat%nq ! gas phase rhs index
            k = i + (j - 1) * dat%nq ! particle rhs index
            
            ! if the gas phase is super-saturated
            if (densities(ii,j) >= gas_sat_den(i,j)*var%condensation_rate(2,i)) then
              ! compute condensation rate
              cond_rate = damp_condensation_rate(var%condensation_rate(1,i), &
                                                 var%condensation_rate(2,i), &
                                                 var%condensation_rate(3,i), &
                                                 densities(ii,j)/gas_sat_den(i,j))
            
              ! rate the gas molecules are going to particles
              ! in molecules/cm3/s
              dn_gas_dt = - cond_rate* &
                            (densities(ii,j) - gas_sat_den(i,j)*var%condensation_rate(2,i))
              ! add to rhs vector, convert to change in mixing ratio/s
              rhs(kk) = rhs(kk) + dn_gas_dt          
            
              ! rate of particle production from gas condensing
              ! in particles/cm3/s
              dn_particle_dt = - dn_gas_dt
              ! add to rhs vector, convert to moles/cm3/s
              rhs(k) = rhs(k) + dn_particle_dt
            else
              ! particles don't change!
              rhs(k) = rhs(k) + 0.0_dp
            endif
          endif          
        enddo
      enddo
      
    endif

  end subroutine

  subroutine diffusion_coefficients_evo(dat, var, den, mubar, &
                                    DU, DD, DL, ADU, ADL, ADD, wfall, VH2_esc, VH_esc)
    use photochem_eqns, only: dynamic_viscosity_air, fall_velocity, slip_correction_factor, &
    binary_diffusion_param
    use photochem_const, only: k_boltz, N_avo

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
    real(dp) :: grav_av, mubar_av
    real(dp) :: bx1x2

    ! for particles
    real(dp) :: air_density_pp, air_density
    real(dp) :: wfall_pp, wfall_i
    real(dp) :: viscosity_pp, viscosity
    real(dp) :: FF2, FF1

    
    integer :: j, i, k

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

        gamma_i_gas_av(k,j) = eddav(j)*(1.0_dp/scale_height_av(j) + (1.0_dp/tav(j))*dTdz(j)) + &
                           (b1x2av(k,j)/denav(j))*(1.0_dp/scale_height_i_av + (1.0_dp/tav(j))*dTdz(j))
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
    enddo
    ! upper boundary
    j = var%nz
    do i = dat%ng_1,dat%nq
      k = i - dat%np
      DL(i,j) = (eddav(j-1) + (b1x2av(k,j-1)/denav(j-1)))/(var%dz(j)**2.0_dp)
      DD(i,j) = - DL(i,j)

      ADL(i,j) = - gamma_i_gas_av(k,j-1)/(2.0_dp*var%dz(j))
      ADD(i,j) = ADL(i,j)
    enddo
    
    ! particles (eddy diffusion)
    ! middle
    do j = 2,var%nz-1
      do i = 1,dat%np
        ! diffusion
        DU(i,j) = eddav(j)/(var%dz(j)**2.0_dp)
        DL(i,j) = eddav(j-1)/(var%dz(j)**2.0_dp)
        DD(i,j) = - DU(i,j) - DL(i,j)

        ! advection
        ADU(i,j) = gamma_i_part_av(i,j)/(2.0_dp*var%dz(j))
        ADL(i,j) = - gamma_i_part_av(i,j-1)/(2.0_dp*var%dz(j))
        ADD(i,j) = ADU(i,j) + ADL(i,j)
      enddo
    enddo
    ! lower boundary
    j = 1
    do i = 1,dat%np
      DU(i,j) = eddav(j)/(var%dz(j)**2.0_dp)
      DD(i,j) = - DU(i,j)

      ADU(i,j) = gamma_i_part_av(i,j)/(2.0_dp*var%dz(j))
      ADD(i,j) = ADU(i,j)
    enddo
    ! upper boundary
    j = var%nz
    do i = 1,dat%np
      DL(i,j) = eddav(j-1)/(var%dz(j)**2.0_dp)
      DD(i,j) = - DL(i,j)

      ADL(i,j) = - gamma_i_part_av(i,j-1)/(2.0_dp*var%dz(j))
      ADD(i,j) = ADL(i,j)
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
    if (dat%diff_H_escape) then
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
    use photochem_enum, only: DensityBC
    use photochem_input, only: compute_gibbs_energy, interp2xsdata
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    real(dp), intent(out) :: usol(:,:)
    real(dp), intent(out) :: molecules_per_particle(:,:)
    real(dp), intent(out) :: pressure(:), density(:), mix(:,:), mubar(:)
    real(dp), intent(out) :: pressure_hydro(:), density_hydro(:)
    character(:), allocatable, intent(out) :: err

    real(dp) :: T_surf_guess
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
      ! var%xs_x_qy (NOT YET BUT WILL SOON)
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

    self%var%trop_ind = minloc(abs(self%var%z - self%var%trop_alt), 1) - 1

    if (self%var%trop_ind < 3) then
      err = 'Tropopause is too low.'
    elseif (self%var%trop_ind > self%var%nz-2) then
      err = 'Tropopause is too high.'
    endif

  end subroutine

  module subroutine prep_all_evo_gas(self, usol_in, err)

    use photochem_common, only: reaction_rates, rainout, photorates
    use photochem_common, only: gas_saturation_density
    use photochem_eqns, only: sat_pressure_H2O
    use photochem_const, only: pi, N_avo, small_real, k_boltz

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
    
    if (dat%water_cond) then
      do i = 1,var%nz
        wrk%H2O_sat_mix(i) = sat_pressure_H2O(var%temperature(i))/wrk%pressure(i)
      enddo
    endif

    !!! diffusion and advection coefficients
    call diffusion_coefficients_evo(dat, var, wrk%density, wrk%mubar, &
    wrk%DU, wrk%DD, wrk%DL, wrk%ADU, wrk%ADL, wrk%ADD, wrk%wfall, wrk%VH2_esc, wrk%VH_esc)
    
    wrk%surface_scale_height = (k_boltz*var%temperature(1)*N_avo)/(wrk%mubar(1)*var%grav(1))
    
    !!! H and H2 escape
    wrk%upper_veff_copy = var%upper_veff
    wrk%lower_vdep_copy = var%lower_vdep
    if (dat%diff_H_escape) then
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

      call gas_saturation_density(dat, var, wrk%mix(dat%LH2O,:), wrk%pressure, &
                                  wrk%gas_sat_den)
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
    call reaction_rates(self%dat, self%var, wrk%density, wrk%densities, wrk%rx_rates)
    
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

  module subroutine rhs_evo_gas(self, neqs, usol_flat, rhs, err)
    use photochem_enum, only: MosesBC, VelocityBC, DensityBC, FluxBC, VelocityDistributedFluxBC
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real  
    
    class(EvoAtmosphere), target, intent(inout) :: self
    integer, intent(in) :: neqs
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
    
    ! fills self%wrk with data
    call prep_all_evo_gas(self, usol_in, err)
    if (allocated(err)) return

    call dochem(self, wrk%usol, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, &
                wrk%density, wrk%mix, wrk%densities, wrk%xp, wrk%xl, rhs)  

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
      elseif (var%lowerboundcond(i) == DensityBC) then
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
                        - (var%edd(1)/wrk%surface_scale_height)*wrk%usol(i,1)/var%dz(1)
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
        jdisth = minloc(var%Z,1, var%Z >= disth) - 1
        jdisth = max(jdisth,2)
        ztop = var%z(jdisth)-var%z(1)
        ztop1 = var%z(jdisth) + 0.5e0_dp*var%dz(jdisth)
        do j = 2,jdisth
          k = i + (j-1)*dat%nq
          rhs(k) = rhs(k) + 2.0_dp*var%lower_flux(i)*(ztop1-var%z(j))/(ztop**2.0_dp)
        enddo
      endif
    enddo 

  end subroutine

  module subroutine jac_evo_gas(self, lda_neqs, neqs, usol_flat, jac, err)
    use photochem_enum, only: MosesBC, VelocityBC, DensityBC, FluxBC, VelocityDistributedFluxBC
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real
    
    class(EvoAtmosphere), target, intent(inout) :: self
    integer, intent(in) :: lda_neqs, neqs
    real(dp), target, intent(in) :: usol_flat(neqs)
    real(dp), intent(out), target :: jac(lda_neqs)
    character(:), allocatable, intent(out) :: err
    
    real(dp), pointer :: usol_in(:,:)
    real(dp), pointer :: djac(:,:)
    real(dp) :: usol_perturb(self%dat%nq,self%var%nz)
    real(dp) :: R(self%var%nz)
    real(dp) :: rhs(self%var%neqs)
    real(dp) :: rhs_perturb(self%var%neqs)
    real(dp) :: density(self%var%nz), mix(self%dat%nq,self%var%nz)
    real(dp) :: densities(self%dat%nsp+1,self%var%nz), xl(self%var%nz), xp(self%var%nz)
  
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
  
    ! compute chemistry contribution to jacobian using forward differences
    jac = 0.0_dp
    call dochem(self, wrk%usol, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, &
                density, mix, densities, xp, xl, rhs) 

    !$omp parallel private(i, j, k, m, mm, usol_perturb, R, density, mix, densities, xl, xp, rhs_perturb)
    usol_perturb = wrk%usol
    !$omp do
    do i = 1,dat%nq
      do j = 1,var%nz
        R(j) = var%epsj*abs(wrk%usol(i,j))
        usol_perturb(i,j) = wrk%usol(i,j) + R(j)
      enddo
      
      call dochem(self, usol_perturb, wrk%rx_rates, &
                  wrk%gas_sat_den, wrk%molecules_per_particle, &
                  wrk%H2O_sat_mix, wrk%H2O_rh, wrk%rainout_rates, &
                  density, mix, densities, xp, xl, rhs_perturb) 
  
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
      elseif (var%lowerboundcond(i) == DensityBC) then

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
                         (var%edd(1)/wrk%surface_scale_height)/var%dz(1)
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
  
  end subroutine


end submodule


