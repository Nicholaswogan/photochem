

module photochem
  implicit none
  
  ! public :: photochem_equil_backrnd_gas, photorates
  ! public :: rhs_backrnd_gas, jac_backrnd_gas
  
  ! private :: reaction_rates, compute_gibbs_energy
  ! private :: gibbs_energy_shomate, chempl, molar_weight
  ! private :: press_and_den, arrhenius_rate, falloff_rate
  ! private :: Troe_noT2, Troe_withT2
  
  integer, private, parameter :: real_kind = kind(1.0d0)
  integer, private, parameter :: err_len = 1000
  
contains
  
  subroutine reaction_rates(nsp, nz, nrT, temperature, density, densities, rx_rates, err)
    use photochem_data, only: rateparams, rxtypes, nreactants, &
                              nproducts, reactants_sp_inds, products_sp_inds, &
                              reverse_info, nrF, reverse, num_efficient, def_eff, &
                              eff_sp_inds, efficiencies, falloff_type
    use photochem_const, only: Rgas, k_boltz, smallest_real ! constants
                              
    integer, intent(in) :: nsp, nz, nrT
    real(real_kind), intent(in) :: temperature(nz)
    real(real_kind), intent(in) :: density(nz)
    real(real_kind), intent(in) :: densities(nsp+1,nz)
    real(real_kind), intent(out) :: rx_rates(nz, nrT)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k, n, l, m
    real(real_kind) :: gibbs_energy(nz,nsp)
    real(real_kind) :: eff_den(nz), F(nz), k0, kinf(nz), Pr(nz)
    real(real_kind) :: gibbR_forward, gibbP_forward
    real(real_kind) :: Dg_forward
    err = ''
    
    do i = 1,nrF
      if (rxtypes(i) == 1) then ! elementary        
        do j = 1,nz 
          rx_rates(j,i) = arrhenius_rate(rateparams(1,i), rateparams(2,i), &
                                         rateparams(3,i), temperature(j))
          
        enddo        
      elseif (rxtypes(i) == 2) then ! three-body
        n = num_efficient(i)
        do j = 1,nz
          eff_den(j) = density(j) * def_eff(i)
          ! subtract the default efficiency, then add the perscribed one
          do k = 1,n ! if n is 0 then it will be skipped
            l = eff_sp_inds(k,i)
            eff_den(j) = eff_den(j) - def_eff(i)*densities(l,j) &
                                    + efficiencies(k,i)*densities(l,j)
          enddo
          rx_rates(j,i) = arrhenius_rate(rateparams(1,i), rateparams(2,i), &
                                         rateparams(3,i), temperature(j)) &
                                         * eff_den(j) ! we multiply density here
        enddo
        
      
      elseif (rxtypes(i) == 3) then ! falloff
        ! compute eff_den, kinf, and Pr at all altitudes
        n = num_efficient(i)
        do j = 1,nz
          eff_den(j) = density(j) * def_eff(i)
          ! subtract the default efficiency, then add the perscribed one
          do k = 1,n
            l = eff_sp_inds(k,i)
            eff_den(j) = eff_den(j) - def_eff(i)*densities(l,j) &
                                    + efficiencies(k,i)*densities(l,j)
          enddo
          k0  = arrhenius_rate(rateparams(1,i), rateparams(2,i), &
                               rateparams(3,i), temperature(j))
          kinf(j) = arrhenius_rate(rateparams(4,i), rateparams(5,i), &
                                   rateparams(6,i), temperature(j))
          kinf(j) = max(kinf(j),smallest_real)
          Pr(j) = k0*eff_den(j)/kinf(j)                        
        enddo
        
        ! compute falloff function
        if (falloff_type(i) == 0) then ! no falloff function
          F = 1.d0
        elseif (falloff_type(i) == 1) then ! Troe falloff without T2
          do j = 1,nz
            F(j) = Troe_noT2(rateparams(7,i), rateparams(8,i), &
                             rateparams(10,i), temperature(j), Pr(j))
          enddo
        elseif (falloff_type(i) == 2) then ! Troe falloff with T2
          do j = 1,nz
            F(j) = Troe_withT2(rateparams(7,i), rateparams(8,i), rateparams(9,i), &
                               rateparams(10,i), temperature(j), Pr(j))
          enddo
        endif
        
        ! compute rate
        do j = 1,nz
          rx_rates(j,i) = falloff_rate(kinf(j), Pr(j), F(j))
        enddo
    
      endif
    enddo ! end loop over forward reactions

    if (reverse) then ! if there are reverse reactions
      ! compute gibbs energy at all altitudes
      call compute_gibbs_energy(nz, nsp, temperature, gibbs_energy, err)
      if (len_trim(err) /= 0) return
      ! compute reverse rate
      do i = nrF+1,nrT
        
        n = reverse_info(i) ! Reaction number of the forward
        l = nreactants(n) ! number of reactants for the forward reaction
        m = nproducts(n) ! number of products for the forward reaction
        do j = 1,nz
          gibbR_forward = 0.d0
          do k = 1,l
            gibbR_forward = gibbR_forward + gibbs_energy(j,reactants_sp_inds(k,n))
          enddo
          gibbP_forward = 0.d0
          do k = 1,m
            gibbP_forward = gibbP_forward + gibbs_energy(j,products_sp_inds(k,n))
          enddo
          Dg_forward = gibbP_forward - gibbR_forward ! DG of the forward reaction (J/mol)
          ! compute the reverse rate
          rx_rates(j,i) = rx_rates(j,n) * &
                          (1.d0/exp(-Dg_forward/(Rgas * temperature(j)))) * &
                          (k_boltz*temperature(j)/1.d6)**(m-l)
        enddo
      enddo
    endif

  end subroutine
  
  subroutine compute_gibbs_energy(nz, nsp, temperature, gibbs_energy, err)
    use photochem_data, only: thermo_data, thermo_temps, species_names
    
    integer, intent(in) :: nz, nsp
    real(real_kind), intent(in) :: temperature(nz)
    real(real_kind), intent(out) :: gibbs_energy(nz,nsp)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k

    err = ''
    
    do i = 1,nsp
      do j = 1,nz
        if ((temperature(j) >= thermo_temps(1,i)) .and. &
           (temperature(j) < thermo_temps(2,i))) then
           k = 1
        elseif ((temperature(j) >= thermo_temps(2,i)) .and. &
              (temperature(j) <= thermo_temps(3,i))) then
           k = 2
        else
          err = 'The temperature is not within the ranges '// &
                'given for the thermodynamic data for '//trim(species_names(i))
          return
        endif
        ! j/mol
        call gibbs_energy_shomate(thermo_data(1:7,k,i), temperature(j), gibbs_energy(j,i))
      enddo
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
  
  subroutine chempl(nz, nsp, nrT, densities, rx_rates, k, loss_start_ind, xp, xl)
    use photochem_data, only: nump, numl, iprod, iloss, &
                              reactants_sp_inds, nreactants
    
    ! input
    integer, intent(in) :: nz, nsp, nrT
    real(real_kind), intent(in) :: densities(nsp+1, nz) ! molecules/cm3 of each species
    real(real_kind), intent(in) :: rx_rates(nz,nrT) ! reaction rates (various units)
    integer, intent(in) :: k ! species number
    integer, intent(in) :: loss_start_ind ! = 1 if LL, and = 2 for short-lived
    
    ! output
    real(real_kind), intent(out) :: xp(nz), xl(nz) ! molecules/cm3/s. if loss_start_ind = 2, then xl is in units of 1/s
    
    ! local
    real(real_kind) :: DD
    integer :: np, nl
    integer :: i, ii, iii, m, l, j
    xp = 0.d0
    xl = 0.d0
    
    np = nump(k) ! k is a species
    ! np is number of reactions that produce species k
    do i=1,np
      m = iprod(i,k) ! m is reaction number
      l = nreactants(m) ! l is the number of reactants
      do j = 1,nz
        DD = 1.d0
        do ii = 1,l
          iii = reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xp(j) = xp(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
    nl = numl(k) ! k is a species
    ! nl is number of reactions that destroy species k
    do i=1,nl
      m = iloss(i,k) ! This will JUST be reaction number
      l = nreactants(m) ! number of reactants
      do j = 1,nz
        DD = 1.d0
        do ii = loss_start_ind,l
          iii = reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xl(j) = xl(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
  subroutine dochem(neqs, nsp, nsl, nq, nz, nrT, usol, density, rx_rates, &
                    densities, xp, xl, rhs)                    
    integer, intent(in) :: neqs, nsp, nsl, nq, nz, nrT
    real(real_kind), intent(in) :: usol(nq,nz), density(nz)
    real(real_kind), intent(in) :: rx_rates(nz,nrT)
    real(real_kind), intent(inout) :: densities(nsp+1,nz), xp(nz), xl(nz)
    real(real_kind), intent(out) :: rhs(neqs)
    
    integer :: i, j, k
    
    do j = 1,nz
      do k = 1,nq
        densities(k,j) = usol(k,j)*density(j)
      enddo
      densities(nsp,j) = (1.d0-sum(usol(:,j)))*density(j) ! background gas
      densities(nsp+1,j) = 1.d0 ! for hv
    enddo
    
    ! short lived
    do k = nq+1,nq+nsl
      call chempl(nz, nsp, nrT, densities, rx_rates, k, 2, xp, xl) 
      densities(k,:) = xp/xl
    enddo
    
    ! long lived              
    do i = 1,nq
      call chempl(nz, nsp, nrT, densities, rx_rates, i, 1, xp, xl)
      do j = 1,nz
        k = i + (j - 1) * nq
        rhs(k) = xp(j)/density(j) - xl(j)/density(j)
      enddo
    enddo
    
  end subroutine
  

  subroutine photorates(nz, nsp, kj, nw, dz, densities, xs_x_qy, &
                        flux, diurnal_fac, u0, Rsfc, &
                        prates, surf_radiance,err)
    use photochem_radtran, only: two_stream
    use photochem_data, only: photonums, reactants_sp_inds, nray, sigray, raynums

    ! input
    integer, intent(in) :: nz, nsp, kj, nw
    real(real_kind), intent(in) :: dz(nz)
    real(real_kind), intent(in) :: densities(nsp+1, nz)
    real(real_kind), intent(in) :: xs_x_qy(nz,kj,nw)
    real(real_kind), intent(in) :: flux(nw)
    real(real_kind), intent(in) :: diurnal_fac, u0, Rsfc
    
    ! output
    real(real_kind), intent(out) :: prates(nz,kj)
    real(real_kind), intent(out) :: surf_radiance(nw)
    character(len=err_len), intent(out) :: err
    
    ! local
    real(real_kind) :: partial_prates(nz,kj)
    real(real_kind) :: tausg(nz), taua(nz), tau(nz), w0(nz)
    real(real_kind) :: amean(nz+1), surf_rad, flx
    real(real_kind) :: amean_grd(nz)
    integer :: l, i, j, jj, k, n, ie, ierr, ierrs
    
    ierrs = 0
    prates = 0.d0
    err = ''
    !$omp parallel private(l, i, j, jj, k, n, ie, ierr, partial_prates, &
    !$omp& tausg, taua, tau, w0, amean, surf_rad, &
    !$omp& amean_grd, flx)
    ierr = 0
    partial_prates = 0.d0
    !$omp do
    do l = 1,nw
      
      tausg = 0.d0
      do i = 1,nray
        j = raynums(i)
        do k = 1,nz
          n = nz+1-k
          tausg(n) = tausg(n) + sigray(i,l)*densities(j,k)*dz(k)
        enddo
      enddo

      taua = 0.d0
      do i = 1,kj
        jj = photonums(i)
        j = reactants_sp_inds(1,jj)
        do k = 1,nz
          n = nz+1-k
          taua(n) = taua(n) + xs_x_qy(k,i,l)*densities(j,k)*dz(k)
        enddo
      enddo

      tau = tausg + taua
      do i = 1,nz
        w0(i) = min(0.99999d0,tausg(i)/tau(i))
      enddo
      
      call two_stream(nz, tau, w0, u0, Rsfc, amean, surf_rad, ie)
      surf_radiance(l) = surf_rad
      ierr = ierr + ie
      do i = 1, nz+1
        amean(i) = abs(amean(i))
        ! amean(i) = max(amean(i),0.d0)
      enddo
      
      ! convert amean to photolysis grid
      do i = 1,nz
        n = nz+1-i
        amean_grd(i) = sqrt(amean(n)*amean(n+1))        
      enddo
      
      flx = flux(l)*diurnal_fac ! photons/cm2/s

      do i=1,kj
        do j=1,nz
          partial_prates(j,i) = partial_prates(j,i) + flx*amean_grd(j)*xs_x_qy(j,i,l)
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
    
    do i=1,kj
      do j=1,nz
        call round(prates(j,i),-8)
      enddo
    enddo
    ! 
    ! do i=1,kj
    !   print*,prates(nz,i)
    ! enddo
    ! stop

  end subroutine
  
  subroutine round(in,precision)
    implicit none
    real(real_kind), intent(inout) :: in
    integer, intent(in) :: precision
    integer :: order
    order = nint(log10(abs(in)))
    in = nint(in * 10.d0**(-precision-order),8)*10.d0**(precision+order)
  end subroutine
  
  subroutine prep_atm_background_gas(nq, nz, trop_ind, sum_usol, usol, &
                                     density, mubar, pressure, fH2O, err)
    
    use photochem_data, only: species_mass, back_gas_mu, water_sat_trop             
    use photochem_vars, only: temperature, grav, dz, surface_pressure
    
    integer, intent(in) :: nq, nz, trop_ind
    real(real_kind), intent(in), target :: usol(nq,nz)
    real(real_kind), intent(inout) :: sum_usol(nz)
    real(real_kind), intent(out) :: density(nz)
    real(real_kind), intent(out) :: mubar(nz), pressure(nz), fH2O(trop_ind)
    character(len=err_len), intent(out) :: err
    
    integer :: i
    
    err = ''
    
    do i = 1,nz
      sum_usol(i) = sum(usol(:,i))
      if (sum_usol(i) > 1.0d0) then
        err = 'Mixing ratios sum to >1.0 at some altitude (should be <=1).' // &
              ' The atmosphere is probably in a run-away state'
        return
      endif
    enddo
    
    do i = 1,nz
      call molar_weight(nq, usol(:,i), sum_usol(i), species_mass, back_gas_mu, mubar(i))
    enddo
    
    call press_and_den(nz, temperature, grav, surface_pressure*1.d6, dz, &
                       mubar, pressure, density)
                       
    if (water_sat_trop) then
      do i = 1,trop_ind
        fH2O(i) = sat_pressure_H2O(temperature(i))/pressure(i)
      enddo
    endif
    
  end subroutine
  
  subroutine prep_all_background_gas(nsp, nq, nz, nrT, kj, nw, trop_ind, usol, densities, &
                                     density, rx_rates, mubar, pressure, fH2O, fH2O_save, &
                                     prates, surf_radiance, &
                                     DU, DD, DL, ADU, ADL, err)
    use photochem_const, only: pi
    use photochem_data, only: photonums, water_sat_trop, LH2O
    use photochem_vars, only: temperature, grav, dz, edd, &
                              xs_x_qy, photon_flux, diurnal_fac, solar_zenith, &
                              surface_albedo
    
    integer, intent(in) :: nsp, nq, nz, nrT, kj, nw, trop_ind
    real(real_kind), intent(inout) :: usol(nq,nz)
    
    real(real_kind), intent(inout) :: densities(nsp+1,nz)
    
    real(real_kind), intent(out) :: density(nz), rx_rates(nz,nrT)
    real(real_kind), intent(out) :: mubar(nz), pressure(nz)
    real(real_kind), intent(out) :: fH2O(trop_ind), fH2O_save(trop_ind)
    real(real_kind), intent(out) :: prates(nz,kj), surf_radiance(nw)
    real(real_kind), intent(out) :: DU(nq,nz), DD(nq,nz), DL(nq,nz)
    real(real_kind), intent(out) :: ADU(nq,nz), ADL(nq,nz)
    
    character(len=err_len), intent(out) :: err
    
    real(real_kind) :: sum_usol(nz)
    real(real_kind) :: u0
    integer :: i, j, k
    

    err = ''
    
    call prep_atm_background_gas(nq, nz, trop_ind, sum_usol, usol, &
                                 density, mubar, pressure, fH2O, err)  
    if (len_trim(err) /= 0) return  
    
    if (water_sat_trop) then
      do i = 1,trop_ind
        fH2O_save(i) = usol(LH2O,i)
        usol(LH2O,i) = fH2O(i)
      enddo
    endif 
        
    ! diffusion coefficients
    call diffusion_coefficients(nq, nz, dz, edd, temperature, density, grav, mubar, &
                                DU, DD, DL, ADU, ADL)
    
    do j = 1,nz
      do i = 1,nq
        densities(i,j) = usol(i,j)*density(j)
      enddo
      densities(nsp,j) = (1.d0-sum_usol(j))*density(j) ! background gas
      densities(nsp+1,j) = 1.d0 ! for hv
    enddo

    call reaction_rates(nsp, nz, nrT, temperature, density, &
                        densities, rx_rates, err)
    if (len_trim(err) /= 0) return
    
    u0 = cos(solar_zenith*pi/180.d0)
    call photorates(nz, nsp, kj, nw, dz, densities, xs_x_qy, &
                    photon_flux, diurnal_fac, u0, surface_albedo, &
                    prates, surf_radiance, err)
    if (len_trim(err) /= 0) return
    
    do i = 1,kj
      k = photonums(i)
      rx_rates(:,k) = prates(:,i) 
    enddo 
    
  end subroutine
  
  
  subroutine rhs_background_gas(neqs, usol_flat, rhs, err)
    use photochem_const, only: pi
    
    use photochem_data, only: nq, nsp, nsl, nrT, kj, nw, LH2O, water_sat_trop
    use photochem_vars, only: nz, z, dz, trop_ind, &
                              lowerboundcond, upperboundcond, lower_vdep, lower_flux, &
                              lower_dist_height, upper_veff, upper_flux, lower_fix_mr
  
    integer, intent(in) :: neqs
    real(real_kind), intent(in), target :: usol_flat(neqs)
    real(real_kind), intent(out) :: rhs(neqs)
    character(len=err_len), intent(out) :: err
    
    real(real_kind), pointer :: usol(:,:)
    
    real(real_kind) :: mubar(nz), pressure(nz)
    real(real_kind) :: density(nz), fH2O(nz)
    real(real_kind) :: densities(nsp+1,nz)
    real(real_kind) :: rx_rates(nz,nrT)
    real(real_kind) :: prates(nz,kj), surf_radiance(nw)
    real(real_kind) :: xp(nz), xl(nz)
    real(real_kind) :: DU(nq,nz), DD(nq,nz), DL(nq,nz), ADU(nq,nz), ADL(nq,nz)
    
    real(real_kind) :: fH2O_save(trop_ind), lower_fix_mr_save(nz)
    
    real(real_kind) :: disth, ztop, ztop1    
    integer :: i, k, j, jdisth
    
    err = ''
    ! reshape usol_flat with a pointer (no copying; same memory)
    usol(1:nq,1:nz) => usol_flat
    
    do i = 1,nq
      if (lowerboundcond(i) == 1) then
        lower_fix_mr_save(i) = usol(i,1)
        usol(i,1) = lower_fix_mr(i)
      endif
    enddo
    
    call prep_all_background_gas(nsp, nq, nz, nrT, kj, nw, trop_ind, usol, densities, &
                             density, rx_rates, mubar, pressure, fH2O, fH2O_save, &
                             prates, surf_radiance, &
                             DU, DD, DL, ADU, ADL, err)
    if (len_trim(err) /= 0) return
    
    call dochem(neqs, nsp, nsl, nq, nz, nrT, usol, density, rx_rates, &
                densities, xp, xl, rhs) 
    
    ! diffusion (interior grid points)
    do j = 2,nz-1
      do i = 1,nq
        k = i + (j-1)*nq
        rhs(k) = rhs(k) + DU(i,j)*usol(i,j+1) + ADU(i,j)*usol(i,j+1) &
                        + DD(i,j)*usol(i,j) &
                        + DL(i,j)*usol(i,j-1) + ADL(i,j)*usol(i,j-1)
      enddo
    enddo
    
    ! Lower boundary
    do i = 1,nq
      if (lowerboundcond(i) == 0 .or. lowerboundcond(i) == 3) then
        rhs(i) = rhs(i) + DU(i,1)*usol(i,2) + ADU(i,1)*usol(i,2) &
                        - DU(i,1)*usol(i,1) &
                        - lower_vdep(i)*usol(i,1)/dz(1)
      elseif (lowerboundcond(i) == 1) then
        ! rhs(i) = 0.d0
        rhs(i) = 0.d0
        ! rhs(i) = - 1.d-5*atan((usol(i,1) - lower_fix_mr(i)))
      else ! (lowerboundcond(i) == 2) then
        rhs(i) = rhs(i) + DU(i,1)*usol(i,2) + ADU(i,1)*usol(i,2) &
                        - DU(i,1)*usol(i,1) &
                        + lower_flux(i)/(density(1)*dz(1))
      endif
    enddo

    ! Upper boundary
    do i = 1,nq
      k = i + (nz-1)*nq
      if (upperboundcond(i) == 0) then
        rhs(k) = rhs(k) - DL(i,nz)*usol(i,nz) &
                        + DL(i,nz)*usol(i,nz-1) + ADL(i,nz)*usol(i,nz-1) &
                        - upper_veff(i)*usol(i,nz)/dz(nz)    
      elseif (upperboundcond(i) == 2) then
        rhs(k) = rhs(k) - DL(i,nz)*usol(i,nz) &
                        + DL(i,nz)*usol(i,nz-1) + ADL(i,nz)*usol(i,nz-1) &
                        - upper_flux(i)/(density(nz)*dz(nz))
      endif
    enddo
    
    ! Distributed (volcanic) sources
    do i = 1,nq
      if (lowerboundcond(i) == 3) then
        disth = lower_dist_height(i)*1.d5        
        jdisth = minloc(Z,1, Z >= disth) - 1
        jdisth = max(jdisth,2)
        ztop = z(jdisth)-z(1)
        ztop1 = z(jdisth) + 0.5d0*dz(jdisth)
        do j = 2,jdisth
          k = i + (j-1)*nq
          rhs(k) = rhs(k) + 2.d0*lower_flux(i)*(ztop1-z(j))/(density(j)*ztop**2.d0)
        enddo
      endif
    enddo 
    
    if (water_sat_trop) then
      do j = 1,trop_ind
        k = LH2O + (j-1)*nq
        rhs(k) = 0.d0
        usol(LH2O,j) = fH2O_save(j)
      enddo
    endif
    
    do i = 1,nq
      if (lowerboundcond(i) == 1) then
        usol(i,1) = lower_fix_mr_save(i)
      endif
    enddo
    
  end subroutine
  
  subroutine jac_background_gas(lda_neqs, neqs, usol_flat, jac, err)
    use photochem_const, only: pi
    
    use photochem_data, only: lda, kd, ku, kl, nq, nsp, nsl, nrT, kj, nw,  &
                              water_sat_trop, LH2O
    use photochem_vars, only: nz, dz, epsj, trop_ind, &
                              lowerboundcond, upperboundcond, lower_vdep, &
                              upper_veff, lower_fix_mr
  
    integer, intent(in) :: lda_neqs, neqs
    real(real_kind), intent(in), target :: usol_flat(neqs)
    real(real_kind), intent(out), target :: jac(lda_neqs)
    character(len=err_len), intent(out) :: err
    
    real(real_kind), pointer :: usol(:,:)
    real(real_kind), pointer :: djac(:,:)
    real(real_kind) :: usol_perturb(nq,nz)
    real(real_kind) :: R(nz)
    real(real_kind) :: rhs(neqs)
    real(real_kind) :: rhs_perturb(neqs)
    
    real(real_kind) :: mubar(nz), pressure(nz), density(nz), fH2O(trop_ind)
    real(real_kind) :: densities(nsp+1,nz)
    real(real_kind) :: rx_rates(nz,nrT)
    real(real_kind) :: prates(nz,kj), surf_radiance(nw)
    real(real_kind) :: xp(nz), xl(nz)
    real(real_kind) :: DU(nq,nz), DD(nq,nz), DL(nq,nz), ADU(nq,nz), ADL(nq,nz)
    real(real_kind) :: fH2O_save(trop_ind)
    
    integer :: i, k, j, m, mm
    
    err = ''
    ! reshape usol_flat with a pointer (no copying; same memory)
    usol(1:nq,1:nz) => usol_flat
    djac(1:lda,1:neqs) => jac
    
    do i = 1,nq
      if (lowerboundcond(i) == 1) then
        usol(i,1) = lower_fix_mr(i)
      endif
    enddo
    
    call prep_all_background_gas(nsp, nq, nz, nrT, kj, nw, trop_ind, usol, densities, &
                             density, rx_rates, mubar, pressure, fH2O, fH2O_save, &
                             prates, surf_radiance, &
                             DU, DD, DL, ADU, ADL, err)
    if (len_trim(err) /= 0) return
    
    ! compute chemistry contribution to jacobian using forward differences
    jac = 0.d0
    call dochem(neqs, nsp, nsl, nq, nz, nrT, usol, density, rx_rates, &
                densities, xp, xl, rhs) 
    !$omp parallel private(i,j,k,m,mm,usol_perturb, R, densities, xp, xl, rhs_perturb)
    usol_perturb = usol
    !$omp do
    do i = 1,nq
      do j = 1,nz
        R(j) = epsj*abs(usol(i,j))
        usol_perturb(i,j) = usol(i,j) + R(j)
      enddo

      call dochem(neqs, nsp, nsl, nq, nz, nrT, usol_perturb, density, rx_rates, &
                  densities, xp, xl, rhs_perturb) 

      do m = 1,nq
        mm = m - i + kd
        do j = 1,nz
          k = i + (j-1)*nq
          djac(mm,k) = (rhs_perturb(m + (j-1)*nq) - rhs(m + (j-1)*nq))/R(j)
        enddo
      enddo
      
      do j= 1,nz
        usol_perturb(i,j) = usol(i,j)
      enddo
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! diffusion (interior grid points)
    do j = 2,nz-1
      do i = 1,nq
        k = i + (j-1)*nq      
        djac(ku,k+nq) = djac(ku,k+nq) + DU(i,j) + ADU(i,j)
        djac(kd,k)    = djac(kd,k)    + DD(i,j)        
        djac(kl,k-nq) = djac(kl,k-nq) + DL(i,j) + ADL(i,j)
      enddo
    enddo
    
    ! Lower boundary
    do i = 1,nq
      if (lowerboundcond(i) == 0 .or. lowerboundcond(i) == 3) then
        ! rhs(i) = rhs(i) + DU(i,1)*usol(i,2) + ADU(i,1)*usol(i,2) &
        !                 - DU(i,1)*usol(i,1) &
        !                 - lower_vdep(i)*usol(i,1)/dz(1)
                        
        djac(ku,i+nq) = djac(ku,i+nq) + DU(i,1) + ADU(i,1)
        djac(kd,i)    = djac(kd,i)    - DU(i,1) - lower_vdep(i)/dz(1)
      elseif (lowerboundcond(i) == 1) then
        ! rhs(i) = - atan((usol(i,1) - lower_fix_mr(i)))
        do m=1,nq
          mm = kd + i - m
          djac(mm,m) = 0.d0
        enddo
        djac(ku,i+nq) = 0.d0
        djac(kd,i) = - DU(i,1)

        ! djac(kd,i) = 1.d0/((lower_fix_mr(i) - usol(i,1))**2.d0 + 1.d0)

      else ! (lowerboundcond(i) == 2) then
        ! rhs(i) = rhs(i) + DU(i,1)*usol(i,2) + ADU(i,1)*usol(i,2) &
        !                 - DU(i,1)*usol(i,1) &
        !                 + lower_flux(i)/(density(1)*dz(1))
                                              
        djac(ku,i+nq) = djac(ku,i+nq) + DU(i,1) + ADU(i,1)
        djac(kd,i)    = djac(kd,i)    - DU(i,1)         
      endif
    enddo

    ! Upper boundary
    do i = 1,nq
      k = i + (nz-1)*nq
      if (upperboundcond(i) == 0) then
        ! rhs(k) = rhs(k) - DL(i,nz)*usol(i,nz) &
        !                 + DL(i,nz)*usol(i,nz-1) + ADL(i,nz)*usol(i,nz-1) &
        !                 - upper_veff(i)*usol(i,nz)/dz(nz)    
        
        djac(kd,k) = djac(kd,k) - DL(i,nz) - upper_veff(i)/dz(nz) 
        djac(kl,k-nq) = djac(kl,k-nq) + DL(i,nz) + ADL(i,nz)
      elseif (upperboundcond(i) == 2) then
        djac(kd,k) = djac(kd,k) - DL(i,nz)
        djac(kl,k-nq) = djac(kl,k-nq) + DL(i,nz) + ADL(i,nz)
      endif
    enddo
    
    if (water_sat_trop) then
      do j = 1,trop_ind
        k = LH2O + (j-1)*nq
        do m = 1,nq
          mm = m - LH2O + kd
          djac(mm,k) = 0.d0
        enddo
        djac(kd,k) = 0.d0
        djac(ku,k+nq) = 0.d0
        usol(LH2O,j) = fH2O_save(j)
      enddo
      do j = 2,trop_ind
        k = LH2O + (j-1)*nq
        djac(kl,k-nq) = 0.d0
      enddo
    endif
    
  end subroutine
  
  integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
                                result(ierr) bind(c, name='RhsFn')
    use, intrinsic :: iso_c_binding
    use fcvode_mod
    use fsundials_nvector_mod
    use photochem_data, only: nq, species_names
    use photochem_vars, only: neqs, verbose, z
    use photochem_wrk, only: nsteps_previous, cvode_mem
    
    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined dat
    
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)
    character(len=err_len) :: err
    integer(c_long) :: nsteps(1)
    integer(c_int) :: loc_ierr
    real(c_double) :: hcur(1)
    real(real_kind) :: tmp, mx
    integer :: k(1),i,j,ii
    
    ierr = 0
    
    ! get data arrays from SUNDIALS vectors
    yvec(1:neqs) => FN_VGetArrayPointer(sunvec_y)
    fvec(1:neqs) => FN_VGetArrayPointer(sunvec_f)
    
    ! fill RHS vector
    call rhs_background_gas(neqs, yvec, fvec, err)
    loc_ierr = FCVodeGetNumSteps(cvode_mem, nsteps)
    
    if (nsteps(1) /= nsteps_previous .and. verbose) then
      loc_ierr = FCVodeGetCurrentStep(cvode_mem, hcur)
      
      tmp = 0.d0
      mx = tmp
      k(1) = 1
      do ii = 1,neqs
        tmp = abs(fvec(ii)/yvec(ii))
        if (tmp > mx .and. yvec(ii) > 1.d-30) then
          mx = tmp
          k(1) = ii
        endif
      enddo
      
      ! k = maxloc(fvec)
      j = k(1)/nq
      i = k(1)-j*nq

      ! print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5)", nsteps, tn, hcur(1)
      
      print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,a10,3x,es12.5,3x,es12.5,3x,es12.5)", &
              nsteps, tn, hcur(1),trim(species_names(i)),fvec(k(1)),yvec(k(1)),z(j+1)/1.d5
      ! print*,nsteps, tn, hcur(1), maxval(fvec),yvec(k),species_names(i),z(j+1)/1.d5
      ! k = 3 + (trop_ind)*nq
      
      ! print*,yvec(3),yvec(3 + (trop_ind-1)*nq),yvec(k)
      
      nsteps_previous = nsteps(1)
    endif
    
    if (len_trim(err) /= 0) then
      print*,trim(err)//". CVODE will attempt to correct the error."
      ierr = 1
    endif
    return
  end function
  
  integer(c_int) function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                                tmp1, tmp2, tmp3) &
                                result(ierr) bind(C,name='JacFn')
    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fcvode_mod
    use fsundials_nvector_mod
    use fnvector_serial_mod
    use fsunmatrix_band_mod
    use fsundials_matrix_mod
    
    use photochem_data, only: lda
    use photochem_vars, only: neqs
    ! use photochem_wrk, only: cvode_mem

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f
    type(SUNMatrix)        :: sunmat_J  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined data
    type(N_Vector)        :: tmp1, tmp2, tmp3

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: Jmat(:)
    character(len=err_len) :: err
    
    ! type(N_Vector)        :: sunvec_errwts
    ! real(c_double), pointer :: errwts(:)
    ! real(c_double) :: fnorm, hcur(1), sig0
    ! integer(c_int) :: ierr1
    ! integer(c_long) :: neqs_cint
    ! 
    ! neqs_cint = neqs
    ! sunvec_errwts = FN_VNew_Serial(neqs_cint)
    ! ierr1 = FCVodeGetErrWeights(cvode_mem,sunvec_errwts)
    ! errwts(1:neqs) => FN_VGetArrayPointer(sunvec_errwts)
    ! fnorm = FN_VWrmsNorm(sunvec_f,sunvec_errwts)
    ! ierr1 = FCVodeGetCurrentStep(cvode_mem, hcur)
    ! sig0 = 2.2204460493d-16*1000.d0*neqs*fnorm*abs(hcur(1))
    ! print*,yvec(15552)*sqrt(2.2204460493d-16),sig0/errwts(15552)
    
    ierr = 0
    
    yvec(1:neqs) => FN_VGetArrayPointer(sunvec_y)
    Jmat(1:neqs*lda) => FSUNBandMatrix_Data(sunmat_J)
    
    call jac_background_gas(lda*neqs, neqs, yvec, Jmat, err)
    if (len_trim(err) /= 0) then
      print*,trim(err)//". CVODE will attempt to correct the error."
      ierr = 1
    endif
    return

  end function
  
  subroutine evolve_background_atm(tstart, nq, nz, usol_start, num_t_eval, t_eval, rtol, atol, &
                                   mxsteps, solution, success, err)
    use photochem_data, only: water_sat_trop, LH2O
    use photochem_vars, only: neqs, lowerboundcond, lower_fix_mr, trop_ind, &
                              initial_dt, use_fast_jacobian, max_err_test_failures, max_order
    use photochem_wrk, only: cvode_mem
    
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, FCVodeInit, FCVodeSStolerances, &
                          FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                          FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                          FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
    use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    
    ! in/out
    real(c_double), intent(in) :: tstart
    integer, intent(in) :: nq, nz
    real(real_kind), intent(in) :: usol_start(nq,nz)
    integer, intent(in) :: num_t_eval
    real(c_double), intent(in) :: t_eval(num_t_eval)
    real(c_double), intent(in) :: rtol, atol
    integer, intent(in) :: mxsteps
    real(real_kind), intent(out) :: solution(nq,nz,num_t_eval)
    logical, intent(out) :: success
    character(len=err_len), intent(out) :: err
    
    ! local variables
    ! real(c_double) :: tstart     ! initial time
    ! real(c_double) :: rtol, atol ! relative and absolute tolerance
    ! real(c_double) :: tout       ! output time
    real(c_double) :: tcur(1)    ! current time
    integer(c_int) :: ierr       ! error flag from C functions
    ! type(c_ptr)    :: cvode_mem  ! CVODE memory
    type(N_Vector), pointer :: sunvec_y ! sundials vector
    
    ! solution vector, neq is set in the ode_mod module
    real(c_double) :: yvec(neqs)
    integer(c_long) :: neqs_long
    integer(c_long) :: mu, ml
    integer(c_long) :: mxsteps_
    type(SUNMatrix), pointer :: sunmat
    type(SUNLinearSolver), pointer :: sunlin
    
    real(real_kind) :: fH2O(trop_ind)
    integer :: i, j, k, ii
    
    ! settings
    mxsteps_ = mxsteps
    neqs_long = neqs
    tcur   = tstart
    mu = nq
    ml = nq
    
    
    ! initialize solution vector
    do j=1,nz
      do i=1,nq
        k = i + (j-1)*nq
        yvec(k) = usol_start(i,j)
      enddo
    enddo
    do i = 1,nq
      if (lowerboundcond(i) == 1) then
        yvec(i) = lower_fix_mr(i)
      endif
    enddo
    if (water_sat_trop) then
      call water_mixing_ratio(nq, nz, trop_ind, usol_start, fH2O, err)
      if (len_trim(err) /= 0) return 
      do j = 1,trop_ind
        k = LH2O + (j-1)*nq
        yvec(k) = fH2O(j)
      enddo
    endif

    ! create SUNDIALS N_Vector
    sunvec_y => FN_VMake_Serial(neqs_long, yvec)
    if (.not. associated(sunvec_y)) then
      err = "CVODE setup error."
      return
    end if
    
    ! create CVode memory
    cvode_mem = FCVodeCreate(CV_BDF)
    if (.not. c_associated(cvode_mem)) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), tstart, sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    sunmat => FSUNBandMatrix(neqs_long, mu, ml)
    sunlin => FSUNLinSol_Band(sunvec_y,sunmat)
    
    ierr = FCVodeSetLinearSolver(cvode_mem, sunlin, sunmat)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    if (use_fast_jacobian) then
      ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))
      if (ierr /= 0) then
        err = "CVODE setup error."
        return
      end if
    endif
    
    ierr = FCVodeSetMaxNumSteps(cvode_mem, mxsteps_)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetInitStep(cvode_mem, initial_dt)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxErrTestFails(cvode_mem, max_err_test_failures)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxOrd(cvode_mem, max_order)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    do ii = 1, num_t_eval
      ierr = FCVode(cvode_mem, t_eval(ii), sunvec_y, tcur, CV_NORMAL)
      if (ierr /= 0) then
        success = .false.
      else
        success = .true.
        do i = 1,nq
          if (lowerboundcond(i) == 1) then
            yvec(i) = lower_fix_mr(i)
          endif
        enddo
        do j=1,nz
          do i=1,nq  
            k = i + (j-1)*nq
            solution(i,j,ii) = yvec(k)
          enddo
        enddo
      endif
    enddo
        
    ! free memory
    call FN_VDestroy(sunvec_y)
    call FCVodeFree(cvode_mem)
    ierr = FSUNLinSolFree(sunlin)
    if (ierr /= 0) then
      err = "CVODE deallocation error"
      return
    end if
    call FSUNMatDestroy(sunmat)

  end subroutine
  
  subroutine photo_equilibrium(mxsteps, rtol, atol, success, err)
    use photochem_data, only: nq
    use photochem_vars, only: nz, usol_init, usol_out, at_photo_equilibrium
    
    integer, intent(in) :: mxsteps
    real(real_kind), intent(in) :: rtol 
    real(real_kind), intent(in) :: atol 
    logical, intent(out) :: success
    character(len=err_len), intent(out) :: err
    
    real(real_kind), pointer :: solution(:,:,:)
    
    err = ""
    
    solution(1:nq,1:nz,1:1) => usol_out
    call evolve_background_atm(0.d0, nq, nz, usol_init, 1, [1.d16],  &
                               rtol, atol, mxsteps, solution, success, err)
    if (len_trim(err) /= 0) return
    at_photo_equilibrium = success 
    
  end subroutine
  
  subroutine water_mixing_ratio(nq, nz, trop_ind, usol, fH2O, err)
    
    integer, intent(in) :: nq, nz, trop_ind
    real(real_kind), intent(in) :: usol(nq,nz)
    real(real_kind), intent(out) :: fH2O(trop_ind)
    character(len=err_len), intent(out) :: err

    real(real_kind) :: sum_usol(nz)
    real(real_kind) :: density(nz)
    real(real_kind) :: mubar(nz), pressure(nz)
    
    err = ""
    
    call prep_atm_background_gas(nq, nz, trop_ind, sum_usol, usol, &
                                 density, mubar, pressure, fH2O, err)
    if (len_trim(err) /= 0) return
    
  end subroutine
  
  subroutine diffusion_coefficients(nq, nz, dz, edd, T, den, grav, mubar, &
                                    DU, DD, DL, ADU, ADL)
    use photochem_const, only: k_boltz, N_avo
    use photochem_data, only: species_mass
    
    integer, intent(in) :: nq, nz
    real(real_kind), intent(in) :: dz(nz), edd(nz), T(nz), den(nz)
    real(real_kind), intent(in) :: grav(nz), mubar(nz)
    real(real_kind), intent(out) :: DU(nq,nz), DL(nq,nz), DD(nq,nz)
    real(real_kind), intent(out) :: ADU(nq,nz), ADL(nq,nz)
    
    real(real_kind) :: eddav_p, eddav_m, denav_p, denav_m, tav_p, tav_m
    real(real_kind) :: bx1x2_p, bx1x2_m, bx1x2_pp, bx1x2_mm, zeta_pp, zeta_mm
    
    integer :: i,j
  
    do i = 2,nz-1
      eddav_p = sqrt(edd(i)*edd(i+1))
      eddav_m = sqrt(edd(i)*edd(i-1))
      denav_p = sqrt(den(i)*den(i+1))
      denav_m = sqrt(den(i)*den(i-1))
      tav_p = sqrt(T(i)*T(i+1))
      tav_m = sqrt(T(i)*T(i-1))
      
      do j = 1,nq
        ! Equation B.4 in Catling and Kasting (2017)
        bx1x2_p = binary_diffusion_param(species_mass(j), mubar(i), tav_p)
        bx1x2_m = binary_diffusion_param(species_mass(j), mubar(i), tav_m)
  
        DU(j,i) = (eddav_p*denav_p + bx1x2_p)/(dz(i)**2.d0*den(i))
        DL(j,i) = (eddav_m*denav_m + bx1x2_m)/(dz(i)**2.d0*den(i))
        DD(j,i) = - DU(j,i) - DL(j,i)
        
        bx1x2_pp = binary_diffusion_param(species_mass(j), mubar(i+1), T(i+1))
        bx1x2_mm = binary_diffusion_param(species_mass(j), mubar(i-1), T(i-1))
        
        zeta_pp =  bx1x2_pp*((species_mass(j)*grav(i+1))/(k_boltz*T(i+1)*N_avo) &
                           - (mubar(i+1)*grav(i+1))/(k_boltz*T(i+1)*N_avo) &
                           + 0.d0) ! zeroed out thermal diffusion    
        zeta_mm =  bx1x2_mm*((species_mass(j)*grav(i-1))/(k_boltz*T(i-1)*N_avo) &
                           - (mubar(i-1)*grav(i-1))/(k_boltz*T(i-1)*N_avo) &
                           + 0.d0) ! zeroed out thermal diffusion
      
        ADU(j,i) = zeta_pp/(2.d0*dz(i)*den(i)) 
        ADL(j,i) = - zeta_mm/(2.d0*dz(i)*den(i))
      enddo      
    enddo

    ! surface layer
    eddav_p = sqrt(edd(1)*edd(2))
    denav_p = sqrt(den(1)*den(2))
    tav_p = sqrt(T(1)*T(2))
    do j = 1,nq
      bx1x2_p = binary_diffusion_param(species_mass(j), mubar(1), tav_p)
      DU(j,1) = (eddav_p*denav_p + bx1x2_p)/(dz(1)**2.d0*den(1))
      DD(j,1) = - DU(j,1)
      ! DL(j,1) = 0.d0
            
      bx1x2_pp = binary_diffusion_param(species_mass(j), mubar(2), T(2))
      zeta_pp =  bx1x2_pp*((species_mass(j)*grav(2))/(k_boltz*T(2)*N_avo) &
                         - (mubar(2)*grav(2))/(k_boltz*T(2)*N_avo) &
                         + 0.d0) ! zeroed out thermal diffusion    
      ADU(j,1) = zeta_pp/(2.d0*dz(1)*den(1)) 
      ! ADL(j,1) = 0.d0
    enddo

    ! top layer
    eddav_m = sqrt(edd(nz)*edd(nz-1))
    denav_m = sqrt(den(nz)*den(nz-1))
    tav_m = sqrt(T(nz)*T(nz-1))
    do j = 1,nq      
      bx1x2_m = binary_diffusion_param(species_mass(j), mubar(nz), tav_m)
      
      ! DU(j,nz) = 0.d0
      DL(j,nz) = (eddav_m*denav_m + bx1x2_m)/(dz(nz)**2.d0*den(nz))
      DD(j,nz) = - DL(j,nz)
            
      bx1x2_mm = binary_diffusion_param(species_mass(j), mubar(nz-1), T(nz-1))
      zeta_mm =  bx1x2_mm*((species_mass(j)*grav(nz-1))/(k_boltz*T(nz-1)*N_avo) &
                         - (mubar(nz-1)*grav(nz-1))/(k_boltz*T(nz-1)*N_avo) &
                         + 0.d0) ! zeroed out thermal diffusion    
      ADL(j,nz) = - zeta_mm/(2.d0*dz(i)*den(i))
      ! ADU(j,nz) = 0.d0
    enddo
      
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
  
  subroutine molar_weight(nq, usol_layer, sum_usol_layer, masses, background_mu, mubar_layer)
    implicit none
    integer, intent(in) :: nq
    real(real_kind), intent(in) :: usol_layer(nq)
    real(real_kind), intent(in) :: sum_usol_layer
    real(real_kind), intent(in) :: masses(nq)
    real(real_kind), intent(in) :: background_mu
    real(real_kind), intent(out) :: mubar_layer
    integer :: j
    real(real_kind) :: f_background
  
    mubar_layer = 0.d0
    do j = 1, nq
      mubar_layer = mubar_layer + usol_layer(j) * masses(j)
    enddo
    f_background = 1.d0 - sum_usol_layer
    mubar_layer = mubar_layer + f_background * background_mu
  
  end subroutine
  
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
  
  subroutine print_reaction_string(rxn)
    integer, intent(in) :: rxn
    character(len=:), allocatable :: rxstring
    call reaction_string(rxn,rxstring)
    print*,rxstring
  end subroutine
  
  subroutine reaction_string(rxn,rxstring)
    use photochem_data, only: reverse_info, nrF, nreactants, species_names, &
                              reactants_sp_inds, rxtypes, nproducts, products_sp_inds

    integer, intent(in) :: rxn
    character(len=:), allocatable, intent(out) :: rxstring
    integer j, k, i

    rxstring = ''
    if (rxn > nrF) then
      i = reverse_info(rxn)
    else
      i = rxn
    endif
    do j = 1,nreactants(rxn)-1
      k = reactants_sp_inds(j,rxn)
      rxstring = rxstring //(trim(species_names(k))//' + ')
    enddo
    
    k = reactants_sp_inds(nreactants(rxn),rxn)
    rxstring = rxstring // trim(species_names(k))//' => '
    
    if ((rxtypes(i) == 2) .or.(rxtypes(i) == 3)) then
      rxstring = rxstring(1:len(rxstring)-4) //(' + M'//' => ')
    endif
    
    do j = 1,nproducts(rxn)-1
      k = products_sp_inds(j,rxn)
      rxstring = rxstring // trim(species_names(k))//' + '
    enddo
    k = products_sp_inds(nproducts(rxn),rxn)
    rxstring = rxstring // trim(species_names(k))
    
    if ((rxtypes(i) == 2) .or.(rxtypes(i) == 3)) then
      rxstring = rxstring //' + M'
    endif
  end subroutine
  
  subroutine get_usol_init(usol)
    use photochem_vars, only: usol_init, nz
    use photochem_data, only: nq
    
    real(real_kind), intent(out) :: usol(nq,nz)
    usol = usol_init

  end subroutine

end module