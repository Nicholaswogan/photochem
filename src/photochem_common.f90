#:set TYPES = ['real(dp)', 'type(dual)']
#:set NAMES = ['real', 'dual']
#:set TYPES_NAMES = list(zip(TYPES, NAMES))

module photochem_common
  use photochem_const, only: dp
  use photochem_types, only : PhotochemData, PhotochemVars
  implicit none

  interface chempl
    module procedure :: chempl_real, chempl_dual
  end interface

  interface chempl_sl
    module procedure :: chempl_sl_real, chempl_sl_dual
  end interface
  
contains
  
  subroutine reaction_rates(dat, var, pressure, density, densities, rx_rates)
    use futils, only: searchsorted
    use photochem_enum, only: NoFalloff, TroeWithoutT2Falloff, TroeWithT2Falloff, JPLFalloff
    use photochem_types, only: ElementaryRate, ThreeBodyRate, FalloffRate, PressDependentRate
    use photochem_eqns, only: arrhenius_rate, Troe_noT2, Troe_withT2, falloff_rate
    use photochem_const, only: Rgas, k_boltz, smallest_real ! constants
    
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: pressure(:) ! (nz)
    real(dp), intent(in) :: density(:) ! (nz)
    real(dp), intent(in) :: densities(:,:) ! (nsp+1,nz)
    real(dp), intent(out) :: rx_rates(:,:) ! (nz,nrT) 
    
    integer :: i, j, k, n, l, m
    real(dp) :: eff_den(var%nz), F(var%nz)
    real(dp) :: k0, kinf(var%nz), Pr(var%nz)
    real(dp) :: logP, k1, k2, logk
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

      class is (PressDependentRate)
        n = size(rp%logP) ! number of pressures in grid

        do j = 1,var%nz
          logP = log(pressure(j))

          if (logP <= rp%logP(1)) then
            ! If P below grid, then constant extrapolation below.
            k1 = 0.0_dp
            do k = 1,size(rp%rate(1)%A)
              k1 = k1 + &
                   arrhenius_rate(rp%rate(1)%A(k), rp%rate(1)%b(k), &
                                  rp%rate(1)%Ea(k), var%temperature(j))
            enddo
            rx_rates(j,i) = k1
          elseif (logP >= rp%logP(n)) then
            ! If P above grid, then constant extrapolation above.
            k1 = 0.0_dp
            do k = 1,size(rp%rate(n)%A)
              k1 = k1 + &
                   arrhenius_rate(rp%rate(n)%A(k), rp%rate(n)%b(k), &
                                  rp%rate(n)%Ea(k), var%temperature(j))
            enddo
            rx_rates(j,i) = k1
          else
            ! P lies within grid
            l = searchsorted(rp%logP, logP) - 1

            ! Compute rate at lower pressure
            k1 = 0.0_dp
            do k = 1,size(rp%rate(l)%A)
              k1 = k1 + &
                   arrhenius_rate(rp%rate(l)%A(k), rp%rate(l)%b(k), &
                                  rp%rate(l)%Ea(k), var%temperature(j))
            enddo
            k1 = log(k1)

            ! Compute rate at higher pressure
            k2 = 0.0_dp
            do k = 1,size(rp%rate(l+1)%A)
              k2 = k2 + &
                   arrhenius_rate(rp%rate(l+1)%A(k), rp%rate(l+1)%b(k), &
                                  rp%rate(l+1)%Ea(k), var%temperature(j))
            enddo
            k2 = log(k2)

            ! log-interpolate to get rate at intermediate pressure.
            logk = k1 + (k2 - k1)*((logP - rp%logP(l))/(rp%logP(l+1) - rp%logP(l)))
            rx_rates(j,i) = exp(logk)
          endif
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
            taup_1 = var%particle_xs(i)%qext(k,l)*pi*var%particle_radius(i,k)**2 &
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
  
  pure subroutine rainout(dat, var, fH2O, den, rainout_rates)
    use photochem_const, only: k_boltz, N_avo, small_real
    use photochem_eqns, only: henrys_law
    
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: fH2O(:) ! (nz)
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
      wH2O(i) = (eddav_p*denav_p/var%dz(i)**2.0_dp) * fH2O(i+1) &
              - (eddav_p*denav_p/var%dz(i)**2.0_dp + eddav_m*denav_m/var%dz(i)**2.0_dp) * fH2O(i) &
              + (eddav_m*denav_m/var%dz(i)**2.0_dp) * fH2O(i-1)
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
  
  #:for TYPE1, NAME in TYPES_NAMES
  subroutine chempl_${NAME}$(dat, var, densities, rx_rates, k, xp, xl)
    #:if NAME == 'dual'
    use forwarddiff
    #:endif
    
    ! input
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: densities(:,:) ! (nsp+1, nz) molecules/cm3 of each species
    real(dp), intent(in) :: rx_rates(:,:) ! (nz,nrT) reaction rates (various units)
    integer, intent(in) :: k ! species number
    
    ! output
    ${TYPE1}$, intent(inout) :: xp(:), xl(:) ! (nz) molecules/cm3/s. if loss_start_ind = 2, then xl is in units of 1/s
    
    ! local
    ${TYPE1}$ :: DD
    integer :: i, ii, iii, m, l, j

    #:if NAME == 'dual'
    DD = dual(size(densities(1,1)%der))
    #:endif
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
  
  #:endfor
  #:for TYPE1, NAME in TYPES_NAMES
  subroutine chempl_sl_${NAME}$(dat, var, densities, rx_rates, k, xp, xl)
    #:if NAME == 'dual'
    use forwarddiff
    #:endif
    ! input
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: densities(:,:) ! (nsp+1, nz) molecules/cm3 of each species
    real(dp), intent(in) :: rx_rates(:,:) ! (nz,nrT) reaction rates (various units)
    integer, intent(in) :: k ! species number
    
    ! output
    ${TYPE1}$, intent(inout) :: xp(:), xl(:) ! (nz) molecules/cm3/s. if loss_start_ind = 2, then xl is in units of 1/s
    
    ! local
    ${TYPE1}$ :: DD
    integer :: i, ii, iii, m, l, j

    #:if NAME == 'dual'
    DD = dual(size(densities(1,1)%der))
    #:endif
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
  
  #:endfor
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

  subroutine gas_saturation_density(dat, var, fH2O, pressure, &
                                    gas_sat_den)
    use photochem_const, only: k_boltz
    use photochem_enum, only: CondensingParticle
    use photochem_enum, only: ArrheniusSaturation, H2SO4Saturation
    use photochem_eqns, only: saturation_density

    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: fH2O(:)
    real(dp), intent(in) :: pressure(:)
    real(dp), intent(out) :: gas_sat_den(:,:)

    real(dp) :: P_H2SO4
    integer :: i, j

    ! compute The saturation density
    do j = 1,var%nz
      do i = 1,dat%np
        if (dat%particle_formation_method(i) == CondensingParticle) then
          if (dat%particle_sat_type(i) == ArrheniusSaturation) then ! arrhenius
            gas_sat_den(i,j) = saturation_density(var%temperature(j), &
                                                  dat%particle_sat_params(1,i), &
                                                  dat%particle_sat_params(2,i), &
                                                  dat%particle_sat_params(3,i))
          elseif (dat%particle_sat_type(i) == H2SO4Saturation) then ! interpolate to H2SO4 data
            call dat%H2SO4_sat%evaluate(var%temperature(j), &
                                        log10(fH2O(j)*pressure(j)/1.e6_dp),P_H2SO4)
            P_H2SO4 = 10.0_dp**(P_H2SO4)
            gas_sat_den(i,j) = (P_H2SO4*1.e6_dp)/(k_boltz*var%temperature(j))
          endif
        endif
      enddo
    enddo

  end subroutine

  subroutine molec_per_particle(dat, var, molecules_per_particle)
    use photochem_const, only: pi, N_avo
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(out) :: molecules_per_particle(:,:)

    integer :: i, j

    do j = 1,var%nz
      do i = 1,dat%np
        molecules_per_particle(i,j) = (4.0_dp/3.0_dp)*pi*var%particle_radius(i,j)**3.0_dp* &
                                       dat%particle_density(i)*(1/dat%species_mass(i))*N_avo
      enddo
    enddo

  end subroutine

  subroutine out2atmosphere_txt_base(dat, var, &
                                     pressure, density, densities, molecules_per_particle, &
                                     filename, overwrite, clip, err)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    real(dp), intent(in) :: pressure(:), density(:), densities(:,:), molecules_per_particle(:,:)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: overwrite, clip
    character(:), allocatable, intent(out) :: err
    
    character(len=100) :: tmp
    integer :: io, i, j
    real(dp) :: val, clip_value

    if (clip) then
      clip_value = 1.0e-40_dp
    else
      clip_value = - huge(1.0_dp)
    endif
    
    if (overwrite) then
      open(1, file=filename, form='formatted', status='replace', iostat=io)
      if (io /= 0) then
        err = "Unable to overwrite file "//trim(filename)
        return
      endif
    else
      open(1, file=filename, form='formatted', status='new', iostat=io)
      if (io /= 0) then
        err = "Unable to create file "//trim(filename)//" because it already exists"
        return
      endif
    endif
    
    tmp = 'alt'
    write(unit=1,fmt="(3x,a27)",advance='no') tmp
    tmp = 'press'
    write(unit=1,fmt="(a27)",advance='no') tmp
    tmp = 'den'
    write(unit=1,fmt="(a27)",advance='no') tmp
    tmp = 'temp'
    write(unit=1,fmt="(a27)",advance='no') tmp
    tmp = 'eddy'
    write(unit=1,fmt="(a27)",advance='no') tmp
    do j = 1,dat%nsp
      tmp = dat%species_names(j)
      write(unit=1,fmt="(a27)",advance='no') tmp
    enddo
    if (dat%there_are_particles) then
      do j = 1,dat%npq
        tmp = trim(dat%species_names(j))//"_r"
        write(unit=1,fmt="(a27)",advance='no') tmp
      enddo
    endif
    
    do i = 1,var%nz
      write(1,*)
      write(unit=1,fmt="(es27.17e3)",advance='no') var%z(i)/1.e5_dp
      write(unit=1,fmt="(es27.17e3)",advance='no') pressure(i)/1.e6_dp
      write(unit=1,fmt="(es27.17e3)",advance='no') density(i)
      write(unit=1,fmt="(es27.17e3)",advance='no') var%temperature(i)
      write(unit=1,fmt="(es27.17e3)",advance='no') var%edd(i)
      ! particles
      if (dat%there_are_particles) then
        do j = 1,dat%npq
          val = densities(j,i)*(molecules_per_particle(j,i)/density(i)) ! mixing ratio of particle
          write(unit=1,fmt="(es27.17e3)",advance='no') max(val, clip_value)
        enddo
      endif
      ! gases
      do j = dat%ng_1,dat%nsp
        val = densities(j,i)/density(i)
        write(unit=1,fmt="(es27.17e3)",advance='no') max(val, clip_value)
      enddo
      ! particle radii
      if (dat%there_are_particles) then
        do j = 1,dat%npq
          write(unit=1,fmt="(es27.17e3)",advance='no') var%particle_radius(j,i)
        enddo
      endif
    enddo
    
    close(1)

  end subroutine
  
end module