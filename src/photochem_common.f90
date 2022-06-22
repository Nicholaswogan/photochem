module photochem_common
  use photochem_const, only: dp
  use photochem_types, only : PhotochemData, PhotochemVars
  implicit none
  
contains
  
  pure subroutine reaction_rates(dat, var, density, densities, rx_rates)
    use photochem_enum, only: NoFalloff, TroeWithoutT2Falloff, TroeWithT2Falloff, JPLFalloff
    use photochem_types, only: ElementaryRate, ThreeBodyRate, FalloffRate
    use photochem_eqns, only: arrhenius_rate, Troe_noT2, Troe_withT2, falloff_rate
    use photochem_const, only: Rgas, k_boltz, smallest_real ! constants
    
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: density(:) ! (nz)
    real(dp), intent(in) :: densities(:,:) ! (nsp+1,nz)
    real(dp), intent(out) :: rx_rates(:,:) ! (nz,nrT) 
    
    integer :: i, j, k, n, l, m
    real(dp) :: eff_den(var%nz), F(var%nz)
    real(dp) :: k0, kinf(var%nz), Pr(var%nz)
    real(dp) :: gibbR_forward, gibbP_forward
    real(dp) :: Dg_forward
    
    
    do i = 1,dat%nrF
      select type(rp => dat%rx(i)%rp)
      class is (ElementaryRate)
        do j = 1,var%nz 
          rx_rates(j,i) = arrhenius_rate(rp%A, rp%b, &
                                         rp%Ea, var%temperature(j))
          
        enddo       
      class is (ThreeBodyRate)
        n = rp%eff%n_eff
        do j = 1,var%nz
          eff_den(j) = density(j)*rp%eff%def_eff
          ! subtract the default efficiency, then add the perscribed one
          do k = 1,n ! if n is 0 then it will be skipped
            l = rp%eff%eff_sp_inds(k)
            eff_den(j) = eff_den(j) - rp%eff%def_eff*densities(l,j) &
                                    + rp%eff%efficiencies(k)*densities(l,j)
          enddo
          rx_rates(j,i) = arrhenius_rate(rp%A, rp%b, &
                                         rp%Ea, var%temperature(j)) &
                                         * eff_den(j) ! we multiply density here
          
        enddo
        
      class is (FalloffRate)
        n = rp%eff%n_eff
        do j = 1,var%nz
          eff_den(j) = density(j)*rp%eff%def_eff
          ! subtract the default efficiency, then add the perscribed one
          do k = 1,n ! if n is 0 then it will be skipped
            l = rp%eff%eff_sp_inds(k)
            eff_den(j) = eff_den(j) - rp%eff%def_eff*densities(l,j) &
                                    + rp%eff%efficiencies(k)*densities(l,j)
          enddo
          
          k0  = arrhenius_rate(rp%A0, rp%b0, &
                               rp%Ea0, var%temperature(j))         
          kinf(j) = arrhenius_rate(rp%Ainf, rp%binf, &
                                   rp%Eainf, var%temperature(j))
          kinf(j) = max(kinf(j),smallest_real)
          Pr(j) = k0*eff_den(j)/kinf(j)                 
        enddo
        
        ! compute falloff function
        if (rp%falloff_type == NoFalloff) then ! no falloff function
          F = 1.0_dp
        elseif (rp%falloff_type == TroeWithoutT2Falloff) then ! Troe falloff without T2
          do j = 1,var%nz
            F(j) = Troe_noT2(rp%A_t, rp%T1, rp%T3, var%temperature(j), Pr(j))
          enddo
        elseif (rp%falloff_type == TroeWithT2Falloff) then ! Troe falloff with T2
          do j = 1,var%nz
            F(j) = Troe_withT2(rp%A_t, rp%T1, rp%T2, rp%T3, var%temperature(j), Pr(j))
          enddo
        elseif (rp%falloff_type == JPLFalloff) then ! JPL falloff function
          do j = 1,var%nz
            F(j) = 0.6_dp**(1.0_dp/(1.0_dp + (log10(Pr(j)))**2.0_dp ))
          enddo
        endif
        
        ! compute rate
        do j = 1,var%nz
          rx_rates(j,i) = falloff_rate(kinf(j), Pr(j), F(j))
        enddo

      end select
      
    enddo ! end loop over forward reactions
    
    if (dat%reverse) then ! if there are reverse reactions
      ! compute reverse rate
      do i = dat%nrF+1,dat%nrT
        
        n = dat%rx(i)%reverse_info ! Reaction number of the forward
        l = dat%rx(n)%nreact ! number of reactants for the forward reaction
        m = dat%rx(n)%nprod ! number of products for the forward reaction
        do j = 1,var%nz
          gibbR_forward = 0.0_dp
          do k = 1,l
            gibbR_forward = gibbR_forward + &
                            var%gibbs_energy(j,dat%rx(n)%react_sp_inds(k)-dat%npq)
          enddo
          gibbP_forward = 0.0_dp
          do k = 1,m
            gibbP_forward = gibbP_forward +  &
                            var%gibbs_energy(j,dat%rx(n)%prod_sp_inds(k)-dat%npq)
          enddo
          Dg_forward = gibbP_forward - gibbR_forward ! DG of the forward reaction (J/mol)
          ! compute the reverse rate
          rx_rates(j,i) = rx_rates(j,n) * &
                          (1.0_dp/exp(-Dg_forward/(Rgas * var%temperature(j)))) * &
                          (k_boltz*var%temperature(j)/1.e6_dp)**(m-l)
        enddo
      enddo
      
    endif

  end subroutine
  
  subroutine photorates(dat, var, densities, &
                        prates, surf_radiance, amean_grd, optical_depth, err)
    use photochem_radtran, only: two_stream
    use photochem_const, only: pi
    ! input
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: densities(:,:) ! (nsp+1,nz)
    
    ! output
    real(dp), intent(out) :: prates(:,:) ! (nz,kj)
    real(dp), intent(out) :: surf_radiance(:) ! (nw)
    real(dp), intent(out) :: amean_grd(:,:) ! (nz,nw)
    real(dp), intent(out) :: optical_depth(:,:) ! (nz,nw)
    character(:), allocatable, intent(out) :: err
    
    ! local
    real(dp) :: partial_prates(var%nz,dat%kj)
    real(dp) :: tausg(var%nz), taua(var%nz), tau(var%nz), w0(var%nz), gt(var%nz)
    real(dp) :: taup(var%nz), tausp(var%nz)
    real(dp) :: amean(var%nz+1), surf_rad, flx
    real(dp) :: taup_1, gt_1, tausp_1(dat%np,var%nz)
    real(dp) :: u0
    integer :: l, i, j, jj, k, n, ie, ierr, ierrs
    
    u0 = cos(var%solar_zenith*pi/180.0_dp)

    ierrs = 0
    prates = 0.0_dp
    !$omp parallel private(l, i, j, jj, k, n, ie, ierr, partial_prates, &
    !$omp& taup, taup_1, tausp, tausp_1, tausg, taua, tau, w0, gt, gt_1, &
    !$omp& amean, surf_rad, &
    !$omp& flx)
    ierr = 0
    partial_prates = 0.0_dp
    !$omp do
    do l = 1,dat%nw
      
      ! rayleigh scattering
      tausg = 0.0_dp
      do i = 1,dat%nray
        j = dat%raynums(i)
        do k = 1,var%nz
          n = var%nz+1-k
          tausg(n) = tausg(n) + dat%sigray(i,l)*densities(j,k)*var%dz(k)
        enddo
      enddo
      
      ! photolysis
      taua = 0.0_dp
      do i = 1,dat%kj
        jj = dat%photonums(i)
        j = dat%rx(jj)%react_sp_inds(1)
        do k = 1,var%nz
          n = var%nz+1-k
          taua(n) = taua(n) + var%xs_x_qy(k,i,l)*densities(j,k)*var%dz(k)
        enddo
      enddo
      
      ! particles
      taup = 0.0_dp
      tausp = 0.0_dp
      tausp_1 = 0.0_dp
      do k = 1,var%nz
        n = var%nz+1-k
        do i = 1,dat%np
          if (var%particle_xs(i)%ThereIsData) then
            taup_1 = var%particle_xs(i)%qext(k,l)*pi*var%particle_radius(i,j)**2 &
                     *densities(i,k)*var%dz(k)
            taup(n) = taup(n) + taup_1
            tausp_1(i,n) = var%particle_xs(i)%w0(k,l)*taup_1
            tausp(n) = tausp(n) + tausp_1(i,n)
          endif
        enddo
      enddo
      gt = 0.0_dp
      do k = 1,var%nz
        n = var%nz+1-k
        gt_1 = 0.0_dp
        do i = 1,dat%np    
          if (var%particle_xs(i)%ThereIsData) then
            gt_1 = gt_1 + var%particle_xs(i)%gt(k,l)*tausp_1(i,n) &
                    /(tausp(n)+tausg(n))
          endif
        enddo
        gt(n) = min(gt_1,0.999999e0_dp)
      enddo
      
      ! sum of all contributions
      tau = tausg + taua + taup + tausp
      optical_depth(:,l) = tau
      do i = 1,var%nz
        w0(i) = min(0.99999e0_dp,(tausg(i) + tausp(i))/tau(i))
      enddo
      
      call two_stream(var%nz, tau, w0, gt, u0, var%surface_albedo, amean, surf_rad, ie)
      surf_radiance(l) = surf_rad
      ierr = ierr + ie
      do i = 1, var%nz+1
        if (amean(i) < -1.d-5) then
          ierr = ierr + 1
        endif
        amean(i) = abs(amean(i))
        ! amean(i) = max(amean(i),0.0_dp)
      enddo
      
      ! convert amean to photolysis grid
      do i = 1,var%nz
        n = var%nz+1-i
        amean_grd(i,l) = sqrt(amean(n)*amean(n+1))        
      enddo
      
      flx = var%photon_flux(l)*var%diurnal_fac*var%photon_scale_factor ! photons/cm2/s
      
      do i=1,dat%kj
        do j=1,var%nz
          partial_prates(j,i) = partial_prates(j,i) + flx*amean_grd(j,l)*var%xs_x_qy(j,i,l)
        enddo
      enddo
      
    enddo
    !$omp enddo
    !$omp critical
    prates = prates + partial_prates
    ierrs = ierrs + ierr
    !$omp end critical
    !$omp end parallel

    if (ierrs /= 0) then
      err = 'Tridagiagonal linear solve in two stream radiative transfer failed.'
      return
    endif
    
  end subroutine
  
  pure subroutine diffusion_coefficients(dat, var, den, mubar, &
                                         DU, DD, DL, ADU, ADL, wfall, VH2_esc, VH_esc)
    use photochem_eqns, only: dynamic_viscosity_air, fall_velocity, slip_correction_factor, &
                              binary_diffusion_param
    use photochem_const, only: k_boltz, N_avo
    
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var

    real(dp), intent(in) :: den(:)
    real(dp), intent(in) :: mubar(:)
    
    real(dp), intent(out) :: DU(:,:), DL(:,:), DD(:,:) ! (nq,nz)
    real(dp), intent(out) :: ADU(:,:), ADL(:,:) ! (nq,nz)
    real(dp), intent(out) :: wfall(:,:) ! (npq,nz)
    real(dp), intent(out) :: VH2_esc, VH_esc
    
    real(dp) :: eddav_p, eddav_m, denav_p, denav_m, tav_p, tav_m
    real(dp) :: bx1x2_p, bx1x2_m, zeta_p, zeta_m
    real(dp) :: grav_p, grav_m, mubar_p, mubar_m
    real(dp) :: bx1x2
    
    ! for particles
    real(dp) :: air_density_p, air_density_m
    real(dp) :: wfall_p, wfall_m
    real(dp) :: viscosity_p, viscosity_m
    real(dp) :: radius_p, radius_m
    
    integer :: i, j
    
    ! eddy diffusion. particles and gases
    ! middle
    do i = 2,var%nz-1
      eddav_p = sqrt(var%edd(i)*var%edd(i+1))
      eddav_m = sqrt(var%edd(i)*var%edd(i-1))
      denav_p = sqrt(den(i)*den(i+1))
      denav_m = sqrt(den(i)*den(i-1))
      do j = 1,dat%nq
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
    do j = 1,dat%nq
      DU(j,1) = (eddav_p*denav_p)/(den(1)*var%dz(1)**2.0_dp)
      DL(j,var%nz) = (eddav_m*denav_m)/(den(var%nz)*var%dz(var%nz)**2.0_dp)
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
      zeta_p =  bx1x2_p*((dat%species_mass(j)*grav_p)/(k_boltz*tav_p*N_avo) &
                         - (mubar_p*grav_p)/(k_boltz*tav_p*N_avo) &
                         + 0.0_dp) ! zeroed out thermal diffusion   
      ADU(j,i) = zeta_p/(2.0_dp*var%dz(i)*den(i)) 
    enddo
    ! upper boundary
    i = var%nz
    do j = dat%ng_1, dat%nq
      bx1x2_m = binary_diffusion_param(dat%species_mass(j), mubar_m, tav_m)
      
      DL(j,i) = DL(j,i) + bx1x2_m/(var%dz(i)**2.0_dp*den(i))
      
      zeta_m =  bx1x2_m*((dat%species_mass(j)*grav_m)/(k_boltz*tav_m*N_avo) &
                        - (mubar_m*grav_m)/(k_boltz*tav_m*N_avo) &
                        + 0.0_dp) ! zeroed out thermal diffusion
                        
      ADL(j,i) = - zeta_m/(2.0_dp*var%dz(i)*den(i))
    enddo
    
    ! Falling particles. Only particles
    ! middle
    do i = 2,var%nz-1
      denav_p = sqrt(den(i)*den(i+1))
      denav_m = sqrt(den(i)*den(i-1))
      tav_p = sqrt(var%temperature(i)*var%temperature(i+1))
      tav_m = sqrt(var%temperature(i)*var%temperature(i-1))
      grav_p = sqrt(var%grav(i)*var%grav(i+1))
      grav_m = sqrt(var%grav(i)*var%grav(i-1))
      mubar_p = sqrt(mubar(i)*mubar(i+1))
      mubar_m = sqrt(mubar(i)*mubar(i-1))
      do j = 1,dat%npq
        
        radius_p = sqrt(var%particle_radius(j,i)*var%particle_radius(j,i+1))
        radius_m = sqrt(var%particle_radius(j,i)*var%particle_radius(j,i-1))
        
        air_density_p = (denav_p/N_avo)*mubar_p
        viscosity_p = dynamic_viscosity_air(tav_p)
        
        wfall_p = fall_velocity(grav_p, radius_p, &
                                dat%particle_density(j), air_density_p, viscosity_p) &
                   *slip_correction_factor(radius_p, denav_p)

        air_density_m = (denav_m/N_avo)*mubar_m
        viscosity_m = dynamic_viscosity_air(tav_m)
        wfall_m = fall_velocity(grav_m, radius_m, &
                                 dat%particle_density(j), air_density_m, viscosity_m) &
                   *slip_correction_factor(radius_m, denav_m)
      
        ADU(j,i) = wfall_p*denav_p/(2.0_dp*var%dz(i)*den(i))
        ADL(j,i) = - wfall_m*denav_m/(2.0_dp*var%dz(i)*den(i))
      enddo
    enddo
    ! top and bottom
    denav_p = sqrt(den(1)*den(2))
    denav_m = sqrt(den(var%nz)*den(var%nz-1))
    tav_p = sqrt(var%temperature(1)*var%temperature(2))
    tav_m = sqrt(var%temperature(var%nz)*var%temperature(var%nz-1))
    grav_p = sqrt(var%grav(1)*var%grav(2))
    grav_m = sqrt(var%grav(var%nz)*var%grav(var%nz-1))
    mubar_p = sqrt(mubar(1)*mubar(2))
    mubar_m = sqrt(mubar(var%nz)*mubar(var%nz-1))
    ! lower boundary
    i = 1
    do j = 1,dat%npq
      radius_p = sqrt(var%particle_radius(j,i)*var%particle_radius(j,i+1))
      
      air_density_p = (denav_p/N_avo)*mubar_p
      viscosity_p = dynamic_viscosity_air(tav_p)
      
      wfall_p = fall_velocity(grav_p, radius_p, &
                              dat%particle_density(j), air_density_p, viscosity_p) &
                 *slip_correction_factor(radius_p, denav_p)
    
      ADU(j,i) = wfall_p*denav_p/(2.0_dp*var%dz(i)*den(i))
    enddo
    ! Upper boundary
    i = var%nz
    do j = 1,dat%npq
      
      radius_m = sqrt(var%particle_radius(j,i)*var%particle_radius(j,i-1))

      air_density_m = (denav_m/N_avo)*mubar_m
      viscosity_m = dynamic_viscosity_air(tav_m)
      wfall_m = fall_velocity(grav_m, radius_m, &
                               dat%particle_density(j), air_density_m, viscosity_m) &
                 *slip_correction_factor(radius_m, denav_m)
    
      ADL(j,i) = - wfall_m*denav_m/(2.0_dp*var%dz(i)*den(i))
    enddo
    
    ! option to turn off everything but eddy diffusion
    do i = 1,dat%nq
      if (var%only_eddy(i)) then
        ADL(i,:) = 0.0_dp
        ADU(i,:) = 0.0_dp
      endif
    enddo
    
    ! H2 escape
    if (dat%diff_H_escape) then
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
        air_density_m = (den(i)/N_avo)*mubar(i)
        viscosity_m = dynamic_viscosity_air(var%temperature(i))
        wfall(j,i) = fall_velocity(var%grav(i), var%particle_radius(j,i), &
                                   dat%particle_density(j), air_density_m, viscosity_m) &
                     *slip_correction_factor(var%particle_radius(j,i), den(i))
      enddo
    enddo
    
  end subroutine
  
  pure subroutine rainout(dat, var, usol, den, rainout_rates)
    use photochem_const, only: k_boltz, N_avo, small_real
    use photochem_eqns, only: henrys_law
    
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: usol(:,:) ! (nq,nz)
    real(dp), intent(in) :: den(:) ! (nz)
    real(dp), intent(out) :: rainout_rates(:,:) ! (nq,trop_ind)
    
    integer :: i, j
    
    real(dp) :: wH2O(var%trop_ind)
    real(dp) :: total_rainfall ! molecules/cm2/s
    real(dp) :: slope, intercept
    real(dp) :: denav_p, eddav_p, denav_m, eddav_m
    real(dp) :: scale_factor
    real(dp) :: k_bar, Q_i, H_coeff
    
    real(dp), parameter :: earth_rainfall_rate = 1.1e17_dp ! molecules/cm2/s
    real(dp), parameter :: fz = 0.05e0_dp ! fraction of the time it rains
    real(dp), parameter :: gamma = 4.e5_dp ! average time of storm cycles (s)
    real(dp), parameter :: LLL = 1.0_dp ! g H2O/m3 of clouds
    real(dp), parameter :: C1 = 1.d-6 !dynes/bar
    real(dp), parameter :: C2 = 1.d-9 ![m3/cm3][L H2O/g H2O]
    real(dp), parameter :: MH2O = 18.0_dp ! g H2O/mol H2O
    real(dp), parameter :: rho_H2O = 1000.0_dp ! g H2O/L H2O
    
    !!!!!!! calculate raining rate !!!!!!!
    ! middle of atmosphere
    wH2O = 0.0_dp
    do i = 2,var%trop_ind-1
      denav_p = sqrt(den(i+1)*den(i))
      eddav_p = sqrt(var%edd(i+1)*var%edd(i))
      denav_m = sqrt(den(i-1)*den(i))
      eddav_m = sqrt(var%edd(i-1)*var%edd(i))
      wH2O(i) = (eddav_p*denav_p/var%dz(i)**2.0_dp) * usol(dat%LH2O,i+1) &
              - (eddav_p*denav_p/var%dz(i)**2.0_dp + eddav_m*denav_m/var%dz(i)**2.0_dp) * usol(dat%lH2O,i) &
              + (eddav_m*denav_m/var%dz(i)**2.0_dp) * usol(dat%lH2O,i-1)
      if (wH2O(i) < 0.0_dp) then
        wH2O(i) = 1.d-20
      endif
    enddo
    ! lets just linearly extrapolate wH2O to bottom and top grid cell
    !!! lower boundary !!!
    slope = (wH2O(3) - wH2O(2))/(var%dz(2))
    intercept = wH2O(2) - slope*var%z(2)
    wH2O(1) = slope*var%z(1) + intercept
    if (wH2O(1) < 0.0_dp) then
      wH2O(1) = 1.d-20
    endif
    !!! upper boundary !!!
    slope = (wH2O(var%trop_ind-1) - wH2O(var%trop_ind-2))/(var%dz(var%trop_ind-1))
    intercept = wH2O(var%trop_ind-1) - slope*var%z(var%trop_ind-1)
    wH2O(var%trop_ind) = slope*var%z(var%trop_ind) + intercept
    if (wH2O(var%trop_ind) < 0.0_dp) then
      wH2O(var%trop_ind) = 1.d-20
    endif
    ! Here we re-scale the rainfall rate so that is the the same as what
    ! is perscribed in the settings file. This means that distribution of
    ! raining is controlled by H2O vs z, but magnitude is fixed.
    total_rainfall = var%rainfall_rate*earth_rainfall_rate
    scale_factor = total_rainfall/sum(wH2O*var%dz(1))
    wH2O = wH2O*scale_factor
    !!!!!!! end calculate raining rate !!!!!!!
    
    !!!!!!! dissolve gas in the rain !!!!!!!!!
    do j = 1,var%trop_ind
      do i = 1,dat%nq
        H_coeff = henrys_law(max(var%temperature(j),273.15e0_dp),dat%henry_data(1,i),dat%henry_data(2,i))*(1.e5_dp)
        H_coeff = max(H_coeff, small_real)
        k_bar = (C1*k_boltz*var%temperature(j)*H_coeff/ &
                (1.0_dp+C1*C2*N_avo*LLL*k_boltz*var%temperature(j)*H_coeff)) &
                * (WH2O(j)*MH2O/rho_H2O) 
        Q_i = (1.0_dp-fz) + (fz/(gamma*k_bar))*(1.0_dp - exp(-k_bar*gamma))
        rainout_rates(i,j) = (1.0_dp/(gamma*Q_i)) * (1.0_dp - exp(-k_bar*gamma))
      enddo
    enddo
    !!!!!!! end dissolve gas in the rain !!!!!!!!!
    
  end subroutine
  
  pure subroutine chempl(dat, var, densities, rx_rates, k, xp, xl)
    
    ! input
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: densities(:,:) ! (nsp+1, nz) molecules/cm3 of each species
    real(dp), intent(in) :: rx_rates(:,:) ! (nz,nrT) reaction rates (various units)
    integer, intent(in) :: k ! species number
    
    ! output
    real(dp), intent(out) :: xp(:), xl(:) ! (nz) molecules/cm3/s. if loss_start_ind = 2, then xl is in units of 1/s
    
    ! local
    real(dp) :: DD
    integer :: i, ii, iii, m, l, j
    
    xp = 0.0_dp
    xl = 0.0_dp
    
    ! k is a species
    ! nump is number of reactions that produce species k
    do i = 1,dat%pl(k)%nump
      m = dat%pl(k)%iprod(i) ! m is reaction number
      l = dat%rx(m)%nreact ! l is the number of reactants
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD * densities(iii,j)
        enddo
        xp(j) = xp(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
    ! k is a species
    ! numl is number of reactions that destroy species k
    do i=1,dat%pl(k)%numl
      m = dat%pl(k)%iloss(i) ! This will JUST be reaction number
      l = dat%rx(m)%nreact ! number of reactants
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD * densities(iii,j)
        enddo
        xl(j) = xl(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
  pure subroutine chempl_sl(dat, var, densities, rx_rates, k, xp, xl)
    
    ! input
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: densities(:,:) ! (nsp+1, nz) molecules/cm3 of each species
    real(dp), intent(in) :: rx_rates(:,:) ! (nz,nrT) reaction rates (various units)
    integer, intent(in) :: k ! species number
    
    ! output
    real(dp), intent(out) :: xp(:), xl(:) ! (nz) molecules/cm3/s. if loss_start_ind = 2, then xl is in units of 1/s
    
    ! local
    real(dp) :: DD
    integer :: i, ii, iii, m, l, j
    
    xp = 0.0_dp
    xl = 0.0_dp
    
    ! k is a species
    ! nump is number of reactions that produce species k
    do i = 1,dat%pl(k)%nump
      m = dat%pl(k)%iprod(i) ! m is reaction number
      l = dat%rx(m)%nreact ! l is the number of reactants
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD * densities(iii,j)
        enddo
        xp(j) = xp(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
    ! k is a species
    ! numl is number of reactions that destroy species k
    do i=1,dat%pl(k)%numl
      m = dat%pl(k)%iloss(i) ! This will JUST be reaction number
      l = dat%rx(m)%nreact ! number of reactants
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          ! We skip the short-lived species.
          if (iii /= k) then
            DD = DD * densities(iii,j)
          endif
        enddo
        xl(j) = xl(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
  pure subroutine chempl_t(dat, var, densities, rx_rates, k, xpT, xlT)
    
    ! input
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: densities(:,:) ! (nsp+1, nz) molecules/cm3 of each species
    real(dp), intent(in) :: rx_rates(:,:) ! (nz,nrT) reaction rates (various units)
    integer, intent(in) :: k ! species number
    
    ! output
    real(dp), intent(out) :: xpT(:,:) ! (nz,nprod) molecules/cm3/s.
    real(dp), intent(out) :: xlT(:,:) ! (nz,nloss)
      
    ! local
    real(dp) :: DD
    integer :: i, ii, iii, m, l, j
    
    xpT = 0.0_dp
    xlT = 0.0_dp
    
    ! k is a species
    ! nump is number of reactions that produce species k
    do i = 1,dat%pl(k)%nump
      m = dat%pl(k)%iprod(i) ! m is reaction number
      l = dat%rx(m)%nreact ! l is the number of reactants
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD * densities(iii,j)
        enddo
        xpT(j,i) = rx_rates(j,m) * DD
      enddo
    enddo
    
    ! k is a species
    ! numl is number of reactions that destroy species k
    do i=1,dat%pl(k)%numl
      m = dat%pl(k)%iloss(i) ! This will JUST be reaction number
      l = dat%rx(m)%nreact ! number of reactants
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD * densities(iii,j)
        enddo
        xlT(j,i) = rx_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
end module