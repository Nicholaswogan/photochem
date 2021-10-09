

module photochem
  use photochem_types, only: WrkBackgroundAtm
  implicit none
  
  integer, private, parameter :: real_kind = kind(1.0d0)
  integer, private, parameter :: err_len = 1024
  
contains
  
  subroutine reaction_rates(nsp, nz, nrT, temperature, density, densities, rx_rates, err)
    use photochem_data, only: rateparams, rxtypes, nreactants, &
                              nproducts, reactants_sp_inds, products_sp_inds, &
                              reverse_info, nrF, reverse, num_efficient, def_eff, &
                              eff_sp_inds, efficiencies, falloff_type, ng, npq
    use photochem_const, only: Rgas, k_boltz, smallest_real ! constants
                              
    integer, intent(in) :: nsp, nz, nrT
    real(real_kind), intent(in) :: temperature(nz)
    real(real_kind), intent(in) :: density(nz)
    real(real_kind), intent(in) :: densities(nsp+1,nz)
    real(real_kind), intent(out) :: rx_rates(nz, nrT)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k, n, l, m
    real(real_kind) :: gibbs_energy(nz,ng)
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
      call compute_gibbs_energy(nz, ng, temperature, gibbs_energy, err)
      if (len_trim(err) /= 0) return
      ! compute reverse rate
      do i = nrF+1,nrT
        
        n = reverse_info(i) ! Reaction number of the forward
        l = nreactants(n) ! number of reactants for the forward reaction
        m = nproducts(n) ! number of products for the forward reaction
        do j = 1,nz
          gibbR_forward = 0.d0
          do k = 1,l
            gibbR_forward = gibbR_forward + gibbs_energy(j,reactants_sp_inds(k,n)-npq)
          enddo
          gibbP_forward = 0.d0
          do k = 1,m
            gibbP_forward = gibbP_forward + gibbs_energy(j,products_sp_inds(k,n)-npq)
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
  
  subroutine compute_gibbs_energy(nz, ng, temperature, gibbs_energy, err)
    use photochem_data, only: thermo_data, thermo_temps, species_names
    
    integer, intent(in) :: nz, ng
    real(real_kind), intent(in) :: temperature(nz)
    real(real_kind), intent(out) :: gibbs_energy(nz,ng)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k

    err = ''
    
    do i = 1,ng
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
  
  subroutine chempl(nz, nsp, nrT, densities, rx_rates, k, xp, xl)
    use photochem_data, only: nump, numl, iprod, iloss, &
                              reactants_sp_inds, nreactants
    
    ! input
    integer, intent(in) :: nz, nsp, nrT
    real(real_kind), intent(in) :: densities(nsp+1, nz) ! molecules/cm3 of each species
    real(real_kind), intent(in) :: rx_rates(nz,nrT) ! reaction rates (various units)
    integer, intent(in) :: k ! species number
    
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
        do ii = 1,l
          iii = reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xl(j) = xl(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
  subroutine chempl_sl(nz, nsp, nrT, densities, rx_rates, k, xp, xl)
    use photochem_data, only: nump, numl, iprod, iloss, &
                              reactants_sp_inds, nreactants
    
    ! input
    integer, intent(in) :: nz, nsp, nrT
    real(real_kind), intent(in) :: densities(nsp+1, nz) ! molecules/cm3 of each species
    real(real_kind), intent(in) :: rx_rates(nz,nrT) ! reaction rates (various units)
    integer, intent(in) :: k ! species number
    
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
        do ii = 1,l
          iii = reactants_sp_inds(ii,m)
          ! We skip the short-lived species.
          if (iii /= k) then
            DD = DD * densities(iii,j)
          endif
        enddo
        xl(j) = xl(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
  subroutine dochem(neqs, nsp, np, nsl, nq, nz, trop_ind, nrT, usol, density, rx_rates, &
                    gas_sat_den, molecules_per_particle, condensation_rate, &
                    H2O_sat_mix, H2O_condensation_rate, rainout_rates, &
                    densities, xp, xl, rhs)                 
    use photochem_const, only: N_avo, pi, small_real
    use photochem_vars, only: relative_humidity_cold_trap, gas_rainout
    use photochem_data, only: ng_1, there_are_particles, particle_gas_phase_ind, &
                              particle_formation_method, fix_water_in_trop, stratospheric_cond, &
                              LH2O
    integer, intent(in) :: neqs, nsp, np, nsl, nq, nz, trop_ind, nrT
    real(real_kind), intent(in) :: usol(nq,nz), density(nz)
    real(real_kind), intent(in) :: rx_rates(nz,nrT)
    real(real_kind), intent(in) :: gas_sat_den(np,nz)
    real(real_kind), intent(in) :: molecules_per_particle(np,nz)
    real(real_kind), intent(in) :: condensation_rate(2,np)
    real(real_kind), intent(in) :: H2O_sat_mix(nz)
    real(real_kind), intent(in) :: H2O_condensation_rate(2)
    real(real_kind), intent(in) :: rainout_rates(nq, trop_ind)
    real(real_kind), intent(inout) :: densities(nsp+1,nz), xp(nz), xl(nz)
    real(real_kind), intent(out) :: rhs(neqs)
    
    real(real_kind) :: dn_gas_dt, dn_particle_dt, H2O_cold_trap, cond_rate
    
    integer :: i, ii, j, k, kk
    
    do j = 1,nz
      do k = 1,np
        densities(k,j) = max(usol(k,j)*(density(j)/molecules_per_particle(k,j)),small_real)
      enddo
      do k = ng_1,nq
        densities(k,j) = usol(k,j)*density(j)
      enddo
      densities(nsp,j) = (1.d0-sum(usol(ng_1:,j)))*density(j) ! background gas
      densities(nsp+1,j) = 1.d0 ! for hv
    enddo
    
    ! short lived
    do k = nq+1,nq+nsl
      call chempl_sl(nz, nsp, nrT, densities, rx_rates, k, xp, xl) 
      densities(k,:) = xp/xl
    enddo
    
    rhs = 0.d0
    
    ! long lived              
    do i = ng_1,nq
      call chempl(nz, nsp, nrT, densities, rx_rates, i, xp, xl)
      do j = 1,nz
        k = i + (j - 1) * nq
        rhs(k) = xp(j)/density(j) - xl(j)/density(j)
      enddo
    enddo
    
    if (gas_rainout) then
      ! rainout rates
      do j = 1,trop_ind
        do i = 1,nq
          k = i + (j - 1) * nq
          rhs(k) = rhs(k) - rainout_rates(i,j)*usol(i,j)
        enddo
      enddo
    endif
    
    ! if stratospheric water condesation is on, then condense
    ! water in the stratosphere.
    if (fix_water_in_trop) then
      if (stratospheric_cond) then
        do j = trop_ind+1,nz
          k = LH2O + (j - 1) * nq
          H2O_cold_trap = relative_humidity_cold_trap*H2O_sat_mix(j)
          if (usol(LH2O,j) >= H2O_cold_trap) then
            
            cond_rate = damp_condensation_rate(H2O_condensation_rate(1), &
                                               H2O_condensation_rate(2), &
                                               usol(LH2O,j)/H2O_cold_trap)
            rhs(k) = rhs(k) - cond_rate*(usol(LH2O,j) - H2O_cold_trap)
             
          endif
        enddo
      endif
    endif
    
    if (there_are_particles) then
      ! formation from reaction
      do i = 1,np
        call chempl(nz, nsp, nrT, densities, rx_rates, i, xp, xl)
        do j = 1,nz
          k = i + (j - 1) * nq
          rhs(k) = (xp(j) - xl(j))/density(j)
        enddo
      enddo
    
      ! particle condensation
      do j = 1,nz
        do i = 1,np
          ! if this particle forms from condensation
          if (particle_formation_method(i) == 1) then
            ii = particle_gas_phase_ind(i) ! index of gas phase
            kk = ii + (j - 1) * nq ! gas phase rhs index
            k = i + (j - 1) * nq ! particle rhs index
            
            ! if the gas phase is super-saturated
            if (densities(ii,j) >= gas_sat_den(i,j)) then
              ! compute condensation rate
              cond_rate = damp_condensation_rate(condensation_rate(1,i), &
                                                 condensation_rate(2,i), &
                                                 densities(ii,j)/gas_sat_den(i,j))
            
              ! rate the gas molecules are going to particles
              ! in molecules/cm3/s
              dn_gas_dt = - cond_rate* &
                            (densities(ii,j) - gas_sat_den(i,j))
              ! add to rhs vector, convert to change in mixing ratio/s
              rhs(kk) = rhs(kk) + dn_gas_dt/density(j)            
            
              ! rate of particle production from gas condensing
              ! in particles/cm3/s
              dn_particle_dt = - dn_gas_dt/density(j)
              ! add to rhs vector, convert to moles/cm3/s
              rhs(k) = dn_particle_dt
            else
              ! particles don't change!
              rhs(k) = 0.d0
            endif
          endif          
        enddo
      enddo
      
    endif
    
  end subroutine
  

  subroutine photorates(nz, nsp, kj, nw, np, dz, densities, xs_x_qy, &
                        w0_particles, qext_particles, gt_particles, &
                        flux, photon_scale_factor, diurnal_fac, u0, Rsfc, &
                        prates, surf_radiance,err)
    use photochem_radtran, only: two_stream
    use photochem_const, only: pi
    use photochem_data, only: photonums, reactants_sp_inds, nray, sigray, raynums
    use photochem_vars, only: particle_radius
    ! input
    integer, intent(in) :: nz, nsp, kj, nw, np
    real(real_kind), intent(in) :: dz(nz)
    real(real_kind), intent(in) :: densities(nsp+1, nz)
    real(real_kind), intent(in) :: xs_x_qy(nz,kj,nw)
    real(real_kind), intent(in) :: w0_particles(np,nz,nw)
    real(real_kind), intent(in) :: qext_particles(np,nz,nw)
    real(real_kind), intent(in) :: gt_particles(np,nz,nw)
    real(real_kind), intent(in) :: flux(nw)
    real(real_kind), intent(in) :: photon_scale_factor
    real(real_kind), intent(in) :: diurnal_fac, u0, Rsfc
    
    ! output
    real(real_kind), intent(out) :: prates(nz,kj)
    real(real_kind), intent(out) :: surf_radiance(nw)
    character(len=err_len), intent(out) :: err
    
    ! local
    real(real_kind) :: partial_prates(nz,kj)
    real(real_kind) :: tausg(nz), taua(nz), tau(nz), w0(nz), gt(nz)
    real(real_kind) :: taup(nz), tausp(nz)
    real(real_kind) :: amean(nz+1), surf_rad, flx
    real(real_kind) :: amean_grd(nz)
    real(real_kind) :: taup_1, gt_1, tausp_1(np,nz)
    integer :: l, i, j, jj, k, n, ie, ierr, ierrs
    

    ierrs = 0
    prates = 0.d0
    err = ''
    !$omp parallel private(l, i, j, jj, k, n, ie, ierr, partial_prates, &
    !$omp& taup, taup_1, tausp, tausp_1, tausg, taua, tau, w0, gt, gt_1, &
    !$omp& amean, surf_rad, &
    !$omp& amean_grd, flx)
    ierr = 0
    partial_prates = 0.d0
    !$omp do
    do l = 1,nw
      
      ! rayleigh scattering
      tausg = 0.d0
      do i = 1,nray
        j = raynums(i)
        do k = 1,nz
          n = nz+1-k
          tausg(n) = tausg(n) + sigray(i,l)*densities(j,k)*dz(k)
        enddo
      enddo
      
      ! photolysis
      taua = 0.d0
      do i = 1,kj
        jj = photonums(i)
        j = reactants_sp_inds(1,jj)
        do k = 1,nz
          n = nz+1-k
          taua(n) = taua(n) + xs_x_qy(k,i,l)*densities(j,k)*dz(k)
        enddo
      enddo
      
      ! particles
      taup = 0.d0
      tausp = 0.d0
      tausp_1 = 0.d0
      do k = 1,nz
        n = nz+1-k
        do i = 1,np
          taup_1 = qext_particles(i,k,l)*pi*particle_radius(i,j)**2*densities(i,k)*dz(k)
          taup(n) = taup(n) + taup_1
          tausp_1(i,n) = w0_particles(i,k,l)*taup_1
          tausp(n) = tausp(n) + tausp_1(i,n)
        enddo
      enddo
      gt = 0.d0
      do k = 1,nz
        n = nz+1-k
        gt_1 = 0.d0
        do i = 1,np
          gt_1 = gt_1 + gt_particles(i,k,l)*tausp_1(i,n) &
                  /(tausp(n)+tausg(n))
        enddo
        gt(n) = min(gt_1,0.999999d0)
      enddo
      
      ! sum of all contributions
      tau = tausg + taua + taup + tausp
      do i = 1,nz
        w0(i) = min(0.99999d0,(tausg(i) + tausp(i))/tau(i))
      enddo
      
      call two_stream(nz, tau, w0, gt, u0, Rsfc, amean, surf_rad, ie)
      surf_radiance(l) = surf_rad
      ierr = ierr + ie
      do i = 1, nz+1
        if (amean(i) < -1.d-5) then
          ierr = ierr + 1
        endif
        amean(i) = abs(amean(i))
        ! amean(i) = max(amean(i),0.d0)
      enddo
      
      ! convert amean to photolysis grid
      do i = 1,nz
        n = nz+1-i
        amean_grd(i) = sqrt(amean(n)*amean(n+1))        
      enddo
      
      flx = flux(l)*diurnal_fac*photon_scale_factor ! photons/cm2/s
      
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
    
    ! if openmp is used, we must round precision for reproducibility
    !$ do i=1,kj
    !$   do j=1,nz
    !$     call round(prates(j,i),-8)
    !$   enddo
    !$ enddo
    
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
                                     density, mubar, pressure, fH2O, H2O_sat_mix, err)
    use, intrinsic :: iso_c_binding, only : c_loc, c_ptr
    use cminpack2fort, only: hybrd1 ! interface to hybrd1 from cminpack
    use photochem_data, only: nll, ng_1, species_mass, back_gas_mu, fix_water_in_trop, LH2O          
    use photochem_vars, only: temperature, grav, dz, surface_pressure, &
                              use_manabe, relative_humidity
    integer, intent(in) :: nq, nz, trop_ind
    real(real_kind), intent(inout), target :: usol(nq,nz)
    real(real_kind), intent(out) :: sum_usol(nz)
    real(real_kind), intent(out) :: density(nz)
    real(real_kind), intent(out) :: mubar(nz), pressure(nz), fH2O(trop_ind), H2O_sat_mix(nz)
    character(len=err_len), intent(out) :: err
    
    real(real_kind) :: rel
    integer :: i
    
    ! hybrd1 varibles for cminpack
    integer :: info
    real(real_kind), parameter :: tol = 1.d-8
    integer :: lwa
    real(real_kind), allocatable :: fvec(:)
    real(real_kind), allocatable :: wa(:)
    type(c_ptr) :: ptr ! void c pointer for cminpack
    ! end hybrd1 varibles
    
    err = ''
    
    ! If water is fixed in the troposhere, lets set it to zero.
    ! This will be to get an initial guess for the non-linear solve.
    if (fix_water_in_trop) then
      do i = 1,trop_ind
        usol(LH2O,i) = 0.d0
      enddo
    endif 
    
    do i = 1,nz
      sum_usol(i) = sum(usol(ng_1:,i))
      if (sum_usol(i) > 1.0d0) then
        err = 'Mixing ratios sum to >1.0 at some altitude (should be <=1).' // &
              ' The atmosphere is probably in a run-away state'
        return
      endif
    enddo
    
    do i = 1,nz
      call molar_weight(nll, usol(ng_1:,i), sum_usol(i), species_mass(ng_1:), back_gas_mu, mubar(i))
    enddo
    
    call press_and_den(nz, temperature, grav, surface_pressure*1.d6, dz, &
                       mubar, pressure, density)
                       
    if (fix_water_in_trop) then
      ! Compute initial guess for fH2O
      do i = 1,trop_ind
        if (use_manabe) then
          ! manabe formula
          rel = 0.77d0*(pressure(i)/pressure(1)-0.02d0)/0.98d0
        else
          rel = relative_humidity 
        endif
        
        fH2O(i) = rel*sat_pressure_H2O(temperature(i))/pressure(i)
      enddo
      ! Here we compute self-consistent water profile.
      ! hybrd1 is nonlinear solver from cminpack. Takes negligable time.
      ptr = c_loc(usol) ! void pointer to usol. Be careful!!!
      lwa = (trop_ind*(3*trop_ind+13))/2 + 2
      allocate(fvec(trop_ind))
      allocate(wa(lwa))
      ! Call hybrd1 from cminpack (c re-write of minpack). I wrote a fortran
      ! interface to cminpack. We pass usol, the atmosphere, via a pointer.
      call hybrd1(fcn_fH2O, ptr, trop_ind, fH2O, fvec, tol, info, wa, lwa)
      if (info /= 1) then
        err = "Non-linear solve for water profile failed."
        return
      endif
      deallocate(fvec, wa)
      
      ! use the solution for tropospheric H2O to compute molar weight
      ! and pressure and density.
      do i = 1,trop_ind
        usol(LH2O,i) = fH2O(i)
      enddo
      
      do i = 1,nz
        sum_usol(i) = sum(usol(ng_1:,i))
      enddo
      
      do i = 1,nz
        call molar_weight(nll, usol(ng_1:,i), sum_usol(i), species_mass(ng_1:), back_gas_mu, mubar(i))
      enddo
      
      call press_and_den(nz, temperature, grav, surface_pressure*1.d6, dz, &
                         mubar, pressure, density)
      
      ! compute H2O saturation mixing ratio at all altitudes
      do i = 1,nz
        H2O_sat_mix(i) = sat_pressure_H2O(temperature(i))/pressure(i)
      enddo
      
    endif 
    
  end subroutine

  ! For computing self-consistent water profile. Called by minpack.
  integer function fcn_fH2O(ptr, n, x, fvec, iflag) result(res)
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr
    use photochem_data, only: nll, ng_1, species_mass, back_gas_mu, LH2O, nq
    use photochem_vars, only: temperature, grav, dz, surface_pressure, &
                              nz, use_manabe, relative_humidity    
    
    type(c_ptr) :: ptr ! void pointer. Be careful!
    integer, value :: n, iflag ! n == trop_ind
    real(real_kind), intent(in) :: x(n) ! x == fH2O
    real(real_kind), intent(out) :: fvec(n) ! fvec == residual
    
    real(real_kind), pointer :: usol(:,:)
    real(real_kind) :: sum_usol(nz) 
    real(real_kind) :: fH2O(n)
    real(real_kind) :: mubar(nz)
    real(real_kind) :: density(nz)
    real(real_kind) :: pressure(nz)
    
    integer :: i
    real(real_kind) :: rel
    
    ! dereferences the pointer. Takes the void ptr, then
    ! declares it to be pointing to a usol shaped array of doubles.
    call c_f_pointer(ptr, usol, [nq,nz])
  
    do i = 1,n
      usol(LH2O,i) = x(i)
    enddo
  
    do i = 1,nz
      sum_usol(i) = sum(usol(ng_1:,i))
    enddo
  
    do i = 1,nz
      call molar_weight(nll, usol(ng_1:,i), sum_usol(i), species_mass(ng_1:), back_gas_mu, mubar(i))
    enddo
  
    call press_and_den(nz, temperature, grav, surface_pressure*1.d6, dz, &
                       mubar, pressure, density)
  
    do i = 1,n
      if (use_manabe) then
        ! manabe formula ()
        rel = 0.77d0*(pressure(i)/pressure(1)-0.02d0)/0.98d0
      else
        rel = relative_humidity 
      endif
      fH2O(i) = rel*sat_pressure_H2O(temperature(i))/pressure(i)
    enddo  
    
    fvec = x - fH2O
    res = 0
  end function
  
  subroutine prep_all_background_gas(wrk, err)
    use photochem_rainout, only: rainout
    use photochem_const, only: pi, k_boltz, N_avo, small_real
    use photochem_data, only: ng_1, photonums, LH, LH2, &
                              back_gas_name, diff_H_escape, npq, &
                              there_are_particles, particle_sat_params, &
                              species_mass, particle_density, particle_formation_method
    use photochem_vars, only: temperature, grav, dz, edd, &
                              xs_x_qy, photon_flux, diurnal_fac, solar_zenith, &
                              surface_albedo, lowerboundcond, lower_fix_mr, &
                              photon_scale_factor, particle_radius, &
                              w0_particles, qext_particles, gt_particles, gas_rainout
    
    type(WrkBackgroundAtm), intent(inout) :: wrk
    character(len=err_len), intent(out) :: err
    
    real(real_kind) :: u0
    integer :: i, j, k
    
    err = ''
    
    do i = 1,wrk%nq
      if (lowerboundcond(i) == 1) then
        wrk%usol(i,1) = lower_fix_mr(i)
      endif
    enddo
    
    call prep_atm_background_gas(wrk%nq, wrk%nz, wrk%trop_ind, wrk%sum_usol, wrk%usol, &
                                 wrk%density, wrk%mubar, wrk%pressure, wrk%fH2O, wrk%H2O_sat_mix, err)  
    if (len_trim(err) /= 0) return  
    
    ! diffusion coefficients
    call diffusion_coefficients(wrk%nq, npq, wrk%nz, dz, edd, temperature, wrk%density, grav, wrk%mubar, &
                                particle_radius, wrk%DU, wrk%DD, wrk%DL, wrk%ADU, wrk%ADL, &
                                wrk%wfall, wrk%VH2_esc, wrk%VH_esc)
    
    ! surface scale height
    wrk%surface_scale_height = (k_boltz*temperature(1)*N_avo)/(wrk%mubar(1)*grav(1))

    ! H and H2 escape
    if (diff_H_escape) then
      if (back_gas_name /= "H2") then
        wrk%upper_veff_copy(LH2) = wrk%VH2_esc                     
      endif
      wrk%upper_veff_copy(LH) = wrk%VH_esc 
    endif
  
    if (there_are_particles) then
      do i = 1,wrk%np
        wrk%lower_vdep_copy(i) = wrk%lower_vdep_copy(i) + wrk%wfall(i,1)
        wrk%upper_veff_copy(i) = wrk%upper_veff_copy(i) + wrk%wfall(i,wrk%nz)
      enddo

      ! compute The saturation density
      do j = 1,wrk%nz
        do i = 1,wrk%np
          if (particle_formation_method(i) == 1) then
            ! problem is that H2SO4 makes aerosols with water.
            ! so pure condensation won't work I think.
            wrk%gas_sat_den(i,j) = saturation_density(temperature(j), &
                                                 particle_sat_params(1,i), &
                                                 particle_sat_params(2,i), &
                                                 particle_sat_params(3,i))
          endif
          wrk%molecules_per_particle(i,j) = (4.d0/3.d0)*pi*particle_radius(i,j)**3.d0* &
                                            particle_density(i)*(1/species_mass(i))*N_avo
        enddo
      enddo
    endif
    
    ! densities include particle densities
    ! we track mol/cm3 for particles (instead of particles/cm3) 
    ! because it makes the numbers smaller and closer to 1, which
    ! is similar to mixing ratios.
    do j = 1,wrk%nz
      do i = 1,npq
        wrk%densities(i,j) = max(wrk%usol(i,j)*(wrk%density(j)/wrk%molecules_per_particle(i,j)), small_real)
      enddo
      do i = ng_1,wrk%nq
        wrk%densities(i,j) = wrk%usol(i,j)*wrk%density(j)
      enddo
      wrk%densities(wrk%nsp,j) = (1.d0-wrk%sum_usol(j))*wrk%density(j) ! background gas
      wrk%densities(wrk%nsp+1,j) = 1.d0 ! for hv
    enddo

    call reaction_rates(wrk%nsp, wrk%nz, wrk%nrT, temperature, wrk%density, &
                        wrk%densities, wrk%rx_rates, err)
    if (len_trim(err) /= 0) return
    
    u0 = cos(solar_zenith*pi/180.d0)
    call photorates(wrk%nz, wrk%nsp, wrk%kj, wrk%nw, wrk%np, dz, wrk%densities, xs_x_qy, &
                    w0_particles, qext_particles, gt_particles, &
                    photon_flux, photon_scale_factor, diurnal_fac, u0, surface_albedo, &
                    wrk%prates, wrk%surf_radiance, err)
    if (len_trim(err) /= 0) return
    
    do i = 1,wrk%kj
      k = photonums(i)
      wrk%rx_rates(:,k) = wrk%prates(:,i) 
    enddo 
    
    ! rainout rates
    if (gas_rainout) then
      call rainout(wrk%nq, wrk%nz, wrk%trop_ind, wrk%usol, temperature, wrk%density, wrk%rainout_rates)
    endif

  end subroutine
  
  subroutine rhs_background_gas(neqs, user_data, usol_flat, rhs, err)
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real  
    use photochem_data, only: np, nq, nll, nsp, nsl, nrT, LH2O, fix_water_in_trop
    use photochem_vars, only: nz, z, dz, trop_ind, edd, condensation_rate, &
                              lowerboundcond, upperboundcond, lower_vdep, lower_flux, &
                              lower_dist_height, upper_veff, upper_flux, H2O_condensation_rate
  
    integer, intent(in) :: neqs
    type(c_ptr), intent(in) :: user_data
    real(real_kind), intent(in) :: usol_flat(neqs)
    real(real_kind), intent(out) :: rhs(neqs)
    character(len=err_len), intent(out) :: err
    
    type(WrkBackgroundAtm), pointer :: wrk
    
    real(real_kind) :: disth, ztop, ztop1    
    integer :: i, k, j, jdisth
    
    err = ''
    
    ! dereference pointer to work arrays
    call c_f_pointer(user_data, wrk)
    
    if (any(usol_flat /= usol_flat)) then
      err = 'Input mixing ratios to the rhs contains NaNs. This is typically '//&
            'related to some mixing ratios getting too negative.'
      return 
    endif
    
    ! make a copy of the mixing ratios. You can not alter the input mixing ratios
    ! this will make CVODE unstable. Guard against division by zero
    do j = 1,nz
      do i = 1,nq
        k = i + (j-1)*nq
        if (usol_flat(k) < 0.d0) then
          wrk%usol(i,j) = min(usol_flat(k),-small_real)
        else
          wrk%usol(i,j) = max(usol_flat(k),small_real)
        endif
      enddo
    enddo
    wrk%upper_veff_copy = upper_veff
    wrk%lower_vdep_copy = lower_vdep
    
    ! alters usol to deal with fixing mixing ratios, and H2O in troposphere.
    call prep_all_background_gas(wrk, err)
    if (len_trim(err) /= 0) return
    
    call dochem(neqs, nsp, np, nsl, nq, nz, trop_ind, nrT, wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, condensation_rate, &
                wrk%H2O_sat_mix, H2O_condensation_rate, wrk%rainout_rates, &
                wrk%densities, wrk%xp, wrk%xl, rhs) 
    
    ! diffusion (interior grid points)
    do j = 2,nz-1
      do i = 1,nq
        k = i + (j-1)*nq
        rhs(k) = rhs(k) + wrk%DU(i,j)*wrk%usol(i,j+1) + wrk%ADU(i,j)*wrk%usol(i,j+1) &
                        + wrk%DD(i,j)*wrk%usol(i,j) &
                        + wrk%DL(i,j)*wrk%usol(i,j-1) + wrk%ADL(i,j)*wrk%usol(i,j-1)
      enddo
    enddo
    
    ! Lower boundary
    do i = 1,nq
      if (lowerboundcond(i) == 0 .or. lowerboundcond(i) == 3) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        - wrk%DU(i,1)*wrk%usol(i,1) &
                        - wrk%lower_vdep_copy(i)*wrk%usol(i,1)/dz(1)
      elseif (lowerboundcond(i) == 1) then
        ! rhs(i) = 0.d0
        rhs(i) = 0.d0
        ! rhs(i) = - 1.d-5*atan((usol(i,1) - lower_fix_mr(i)))
      elseif (lowerboundcond(i) == 2) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        - wrk%DU(i,1)*wrk%usol(i,1) &
                        + lower_flux(i)/(wrk%density(1)*dz(1))
      ! Moses (2001) boundary condition for gas giants
      ! A deposition velocity controled by how quickly gases
      ! turbulantly mix vertically
      elseif (lowerboundcond(i) == -1) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        - wrk%DU(i,1)*wrk%usol(i,1) &
                        - (edd(1)/wrk%surface_scale_height)*wrk%usol(i,1)/dz(1)
      endif
    enddo

    ! Upper boundary
    do i = 1,nq
      k = i + (nz-1)*nq
      if (upperboundcond(i) == 0) then
        rhs(k) = rhs(k) - wrk%DL(i,nz)*wrk%usol(i,nz) &
                        + wrk%DL(i,nz)*wrk%usol(i,nz-1) + wrk%ADL(i,nz)*wrk%usol(i,nz-1) &
                        - wrk%upper_veff_copy(i)*wrk%usol(i,nz)/dz(nz)    
      elseif (upperboundcond(i) == 2) then
        rhs(k) = rhs(k) - wrk%DL(i,nz)*wrk%usol(i,nz) &
                        + wrk%DL(i,nz)*wrk%usol(i,nz-1) + wrk%ADL(i,nz)*wrk%usol(i,nz-1) &
                        - upper_flux(i)/(wrk%density(nz)*dz(nz))
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
          rhs(k) = rhs(k) + 2.d0*lower_flux(i)*(ztop1-z(j))/(wrk%density(j)*ztop**2.d0)
        enddo
      endif
    enddo 
    
    if (fix_water_in_trop) then
      do j = 1,trop_ind
        k = LH2O + (j-1)*nq
        rhs(k) = 0.d0
      enddo
    endif
    
  end subroutine
  
  subroutine jac_background_gas(lda_neqs, neqs, user_data, usol_flat, jac, err)
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real
    
    use photochem_data, only: lda, kd, ku, kl, nq, nll, nsp, np, nsl, nrT,  &
                              fix_water_in_trop, LH2O
    use photochem_vars, only: nz, dz, epsj, trop_ind, edd, condensation_rate, &
                              lowerboundcond, upperboundcond, lower_vdep, &
                              upper_veff, H2O_condensation_rate
  
    integer, intent(in) :: lda_neqs, neqs
    type(c_ptr), intent(in) :: user_data
    real(real_kind), intent(in) :: usol_flat(neqs)
    real(real_kind), intent(out), target :: jac(lda_neqs)
    character(len=err_len), intent(out) :: err
    
    real(real_kind), pointer :: djac(:,:)
    real(real_kind) :: usol_perturb(nq,nz)
    real(real_kind) :: R(nz)
    real(real_kind) :: rhs(neqs)
    real(real_kind) :: rhs_perturb(neqs)
    real(real_kind) :: densities(nsp+1,nz), xl(nz), xp(nz)
    
    type(WrkBackgroundAtm), pointer :: wrk
    ! we need these work arrays for parallel jacobian claculation.
    ! It is probably possible to use memory in "wrk", but i will ignore
    ! this for now.
    
    integer :: i, k, j, m, mm
    
    err = ''
    
    call c_f_pointer(user_data, wrk)
    
    if (any(usol_flat /= usol_flat)) then
      err = 'Input mixing ratios to the rhs contains NaNs. This is typically '//&
            'related to some mixing ratios getting too negative.'
      return 
    endif
    
    ! make a copy of the mixing ratios. You can not alter the input mixing ratios
    ! this will make CVODE unstable. Guard against division by zero
    do j = 1,nz
      do i = 1,nq
        k = i + (j-1)*nq
        if (usol_flat(k) < 0.d0) then
          wrk%usol(i,j) = min(usol_flat(k),-small_real)
        else
          wrk%usol(i,j) = max(usol_flat(k),small_real)
        endif
      enddo
    enddo
    wrk%upper_veff_copy = upper_veff
    wrk%lower_vdep_copy = lower_vdep
    
    ! pointer to jac. We reshape to make accessing more intuitive.
    djac(1:lda,1:neqs) => jac
    
    ! alters usol.
    call prep_all_background_gas(wrk, err)
    if (len_trim(err) /= 0) return
    
    ! compute chemistry contribution to jacobian using forward differences
    jac = 0.d0
    call dochem(neqs, nsp, np, nsl, nq, nz, trop_ind, nrT, wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, condensation_rate, &
                wrk%H2O_sat_mix, H2O_condensation_rate, wrk%rainout_rates, &
                densities, xp, xl, rhs) 
    !$omp parallel private(i, j, k, m, mm, usol_perturb, R, densities, xl, xp, rhs_perturb)
    usol_perturb = wrk%usol
    !$omp do
    do i = 1,nq
      do j = 1,nz
        R(j) = epsj*abs(wrk%usol(i,j))
        usol_perturb(i,j) = wrk%usol(i,j) + R(j)
      enddo

      call dochem(neqs, nsp, np, nsl, nq, nz, trop_ind, nrT, usol_perturb, wrk%density, wrk%rx_rates, &
                  wrk%gas_sat_den, wrk%molecules_per_particle, condensation_rate, &
                  wrk%H2O_sat_mix, H2O_condensation_rate, wrk%rainout_rates, &
                  densities, xp, xl, rhs_perturb) 

      do m = 1,nq
        mm = m - i + kd
        do j = 1,nz
          k = i + (j-1)*nq
          djac(mm,k) = (rhs_perturb(m + (j-1)*nq) - rhs(m + (j-1)*nq))/R(j)
        enddo
      enddo
      
      do j= 1,nz
        usol_perturb(i,j) = wrk%usol(i,j)
      enddo
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! diffusion (interior grid points)
    do j = 2,nz-1
      do i = 1,nq
        k = i + (j-1)*nq      
        djac(ku,k+nq) = wrk%DU(i,j) + wrk%ADU(i,j)
        djac(kd,k)    = djac(kd,k)    + wrk%DD(i,j)        
        djac(kl,k-nq) = wrk%DL(i,j) + wrk%ADL(i,j)
      enddo
    enddo
    
    ! Lower boundary
    do i = 1,nq
      if (lowerboundcond(i) == 0 .or. lowerboundcond(i) == 3) then
        ! rhs(i) = rhs(i) + DU(i,1)*usol(i,2) + ADU(i,1)*usol(i,2) &
        !                 - DU(i,1)*usol(i,1) &
        !                 - lower_vdep(i)*usol(i,1)/dz(1)
                        
        djac(ku,i+nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(kd,i)    = djac(kd,i)    - wrk%DU(i,1) - wrk%lower_vdep_copy(i)/dz(1)
      elseif (lowerboundcond(i) == 1) then
        ! rhs(i) = - atan((usol(i,1) - lower_fix_mr(i)))
        do m=1,nq
          mm = kd + i - m
          djac(mm,m) = 0.d0
        enddo
        djac(ku,i+nq) = 0.d0
        ! For some reason this term makes the integration
        ! much happier. I will keep it. Jacobians don't need to be perfect.
        djac(kd,i) = - wrk%DU(i,1)

        ! djac(kd,i) = 1.d0/((lower_fix_mr(i) - usol(i,1))**2.d0 + 1.d0)

      elseif (lowerboundcond(i) == 2) then
        ! rhs(i) = rhs(i) + DU(i,1)*usol(i,2) + ADU(i,1)*usol(i,2) &
        !                 - DU(i,1)*usol(i,1) &
        !                 + lower_flux(i)/(density(1)*dz(1))
                                              
        djac(ku,i+nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(kd,i)    = djac(kd,i)    - wrk%DU(i,1)
      elseif (lowerboundcond(i) == -1) then
        djac(ku,i+nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(kd,i)    = djac(kd,i)    - wrk%DU(i,1) - (edd(1)/wrk%surface_scale_height)/dz(1)
      endif
    enddo

    ! Upper boundary
    do i = 1,nq
      k = i + (nz-1)*nq
      if (upperboundcond(i) == 0) then
        ! rhs(k) = rhs(k) - DL(i,nz)*usol(i,nz) &
        !                 + DL(i,nz)*usol(i,nz-1) + ADL(i,nz)*usol(i,nz-1) &
        !                 - upper_veff(i)*usol(i,nz)/dz(nz)    
        
        djac(kd,k) = djac(kd,k) - wrk%DL(i,nz) - wrk%upper_veff_copy(i)/dz(nz) 
        djac(kl,k-nq) = wrk%DL(i,nz) + wrk%ADL(i,nz)
      elseif (upperboundcond(i) == 2) then
        djac(kd,k) = djac(kd,k) - wrk%DL(i,nz)
        djac(kl,k-nq) = wrk%DL(i,nz) + wrk%ADL(i,nz)
      endif
    enddo
    
    if (fix_water_in_trop) then
      do j = 1,trop_ind
        k = LH2O + (j-1)*nq
        do m = 1,nq
          mm = m - LH2O + kd
          djac(mm,k) = 0.d0
        enddo
        ! For some reason this term makes the integration
        ! much happier. I will keep it. Jacobians don't need to be perfect.
        djac(kd,k) = - wrk%DU(LH2O,j)
        djac(ku,k+nq) = 0.d0
        if (j /= 1) then
          djac(kl,k-nq) = 0.d0
        endif
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
    use photochem_wrk, only: atol_global
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
    integer :: k,i,j,ii
    
    ierr = 0
    
    ! get data arrays from SUNDIALS vectors
    yvec(1:neqs) => FN_VGetArrayPointer(sunvec_y)
    fvec(1:neqs) => FN_VGetArrayPointer(sunvec_f)
    
    ! fill RHS vector
    call rhs_background_gas(neqs, user_data, yvec, fvec, err)
    loc_ierr = FCVodeGetNumSteps(cvode_mem, nsteps)
    
    if (nsteps(1) /= nsteps_previous .and. verbose > 0) then
      loc_ierr = FCVodeGetCurrentStep(cvode_mem, hcur)
      
      if (verbose == 1) then
        print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,'max(dy/dt) = ',es11.5)", &
             nsteps, tn, hcur(1),maxval(abs(fvec))
             
      elseif (verbose == 2) then
        ! Find the fastest changing variable
        tmp = 0.d0
        mx = tmp
        k = 1
        do ii = 1,neqs
          tmp = abs(fvec(ii)/yvec(ii))
          if (tmp > mx .and. abs(yvec(ii)) > atol_global) then
            mx = tmp
            k = ii
          endif
        enddo
        j = k/nq
        i = k-j*nq
        
        print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,"// &
             "'dy/dt =',es12.5,3x,' y =',es12.5,3x,a8,3x,' z =',f6.2,' km')", &
             nsteps, tn, hcur(1),fvec(k),yvec(k),trim(species_names(i)),z(j+1)/1.d5
      endif
      
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
    
    ierr = 0    
    yvec(1:neqs) => FN_VGetArrayPointer(sunvec_y)
    Jmat(1:neqs*lda) => FSUNBandMatrix_Data(sunmat_J)
    call jac_background_gas(lda*neqs, neqs, user_data, yvec, Jmat, err)
    if (len_trim(err) /= 0) then
      print*,trim(err)//". CVODE will attempt to correct the error."
      ierr = 1
    endif
    return

  end function
  
  subroutine evolve_background_atm(tstart, nq, nz, usol_start, num_t_eval, t_eval, rtol, atol, &
                                   mxsteps, solution, success, err)
    use photochem_data, only: fix_water_in_trop, LH2O, nsp, nrT, kj, nw, np
    use photochem_vars, only: no_water_profile, neqs, lowerboundcond, lower_fix_mr, trop_ind, &
                              initial_dt, use_fast_jacobian, max_err_test_failures, max_order, trop_ind
    use photochem_wrk, only: cvode_mem, wrk_out
    
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, FCVodeInit, FCVodeSStolerances, &
                          FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                          FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                          FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd, &
                          FCVodeSetUserData
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
    
    type(c_ptr)    :: user_data
    type(WrkBackgroundAtm), target :: wrk
  
    err = ''
    
    ! settings
    mxsteps_ = mxsteps
    neqs_long = neqs
    tcur   = tstart
    mu = nq
    ml = nq
    
    call wrk%init(nsp, np, nq, nz, nrT, kj, nw, trop_ind)
    user_data = c_loc(wrk)
    
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
    if (fix_water_in_trop) then
      call water_mixing_ratio(nq, nz, trop_ind, usol_start, fH2O, err)
      if (len_trim(err) /= 0) return 
      do j = 1,trop_ind
        k = LH2O + (j-1)*nq
        yvec(k) = fH2O(j)
      enddo
      ! if there is no water profile, then we should extrapolate H2O
      ! above the tropopause
      if (no_water_profile) then
        do j = trop_ind+1,nz
          k = LH2O + (j-1)*nq
          yvec(k) = fH2O(trop_ind)
        enddo
      endif
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
    
    ! set user data
    ierr = FCVodeSetUserData(cvode_mem, user_data)
    if (ierr /= 0) then
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
        do j=1,nz
          do i=1,nq  
            k = i + (j-1)*nq
            solution(i,j,ii) = yvec(k)
          enddo
        enddo
        
        ! this will alter solution(:,:,ii) with proper fixed mixing ratios and H2O
        wrk%usol = solution(:,:,ii)
        call prep_all_background_gas(wrk, err)
        solution(:,:,ii) = wrk%usol
                                     
      endif
    enddo
    
    ! save the last wrk
    wrk_out = wrk
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
    use photochem_vars, only: nz, usol_init, usol_out, at_photo_equilibrium, &
                              equilibrium_time
    use photochem_wrk, only: atol_global
    integer, intent(in) :: mxsteps
    real(real_kind), intent(in) :: rtol 
    real(real_kind), intent(in) :: atol 
    logical, intent(out) :: success
    character(len=1024), intent(out) :: err ! err_len (f2py needs explicit length)
    
    real(real_kind), pointer :: solution(:,:,:)
    
    err = ""
    
    atol_global = atol
    solution(1:nq,1:nz,1:1) => usol_out
    call evolve_background_atm(0.d0, nq, nz, usol_init, 1, [equilibrium_time],  &
                               rtol, atol, mxsteps, solution, success, err)
    if (len_trim(err) /= 0) return
    at_photo_equilibrium = success 
    
  end subroutine
  
  subroutine water_mixing_ratio(nq, nz, trop_ind, usol, fH2O, err)
    
    integer, intent(in) :: nq, nz, trop_ind
    real(real_kind), intent(in) :: usol(nq,nz)
    real(real_kind), intent(out) :: fH2O(trop_ind)
    character(len=err_len), intent(out) :: err
    
    real(real_kind) :: usol_copy(nq,nz)

    real(real_kind) :: sum_usol(nz)
    real(real_kind) :: density(nz)
    real(real_kind) :: mubar(nz), pressure(nz), H2O_sat_mix(nz)
    
    err = ""
    
    usol_copy = usol
    
    call prep_atm_background_gas(nq, nz, trop_ind, sum_usol, usol_copy, &
                                 density, mubar, pressure, fH2O, H2O_sat_mix, err)
    if (len_trim(err) /= 0) return
    
  end subroutine
  
  subroutine output_surface_fluxes(nq, surface_flux, err)
    use photochem_vars, only: at_photo_equilibrium, usol_out, nz
    integer, intent(in) :: nq
    real(real_kind), intent(out) :: surface_flux(nq)
    character(len=1024), intent(out) :: err
    err = ''
    
    if (at_photo_equilibrium) then
      call compute_surface_fluxes(nq, nz, usol_out, surface_flux, err)
    else
      err = "Can not compute surface fluxes without first converging to photochemical equilibrium."
      return
    endif
    
  end subroutine
  
  subroutine compute_surface_fluxes(nq, nz, usol, surface_flux, err)
    use photochem_data, only: nsp, nrT, nll, kj, nw, nsl, fix_water_in_trop, LH2O, np
    use photochem_vars, only: trop_ind, neqs, upper_veff, lower_vdep, dz, &
                              condensation_rate, H2O_condensation_rate
                              
    integer, intent(in) :: nq, nz
    real(real_kind), intent(in) :: usol(nq,nz)
    real(real_kind), intent(out) :: surface_flux(nq)
    character(len=1024), intent(out) :: err
  
    real(real_kind) :: rhs(neqs)  
    type(WrkBackgroundAtm) :: wrk
    real(real_kind) :: diffusive_production
    real(real_kind) :: chemical_production
  
    integer :: i
  
    err = ''
    
    call wrk%init(nsp, np, nq, nz, nrT, kj, nw, trop_ind)
  
    wrk%usol = usol
    wrk%upper_veff_copy = upper_veff
    wrk%lower_vdep_copy = lower_vdep
    call prep_all_background_gas(wrk, err)
    if (len_trim(err) /= 0) return
  
    call dochem(neqs, nsp, np, nsl, nq, nz, trop_ind, nrT, wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, condensation_rate, &
                wrk%H2O_sat_mix, H2O_condensation_rate, wrk%rainout_rates, &
                wrk%densities, wrk%xp, wrk%xl, rhs) 
                          
    ! surface flux is molecules required to sustain the lower boundary
    ! chemical production + diffusion production = total change in lower cell    
    do i = 1,nq
      diffusive_production =   (wrk%DU(i,1)*usol(i,2) + wrk%ADU(i,1)*usol(i,2) &
                                - wrk%DU(i,1)*usol(i,1)) &
                                *wrk%density(1)*dz(1)
      chemical_production = rhs(i)*wrk%density(1)*dz(1)
      surface_flux(i) = -(diffusive_production + chemical_production)
      ! We don't count chemical production for water
      if (fix_water_in_trop) then
        if (i == LH2O) then
          surface_flux(i) = - diffusive_production
        endif
      endif
    enddo
    
  end subroutine
  
  subroutine diffusion_coefficients(nq, npq, nz, dz, edd, T, den, grav, mubar, &
                                    r_particles, DU, DD, DL, ADU, ADL, wfall, VH2_esc, VH_esc)
    use photochem_const, only: k_boltz, N_avo
    use photochem_data, only: species_mass, LH2, LH, back_gas_name, diff_H_escape, particle_density, ng_1
    
    integer, intent(in) :: nq, npq, nz
    real(real_kind), intent(in) :: dz(nz), edd(nz), T(nz), den(nz)
    real(real_kind), intent(in) :: grav(nz), mubar(nz)
    real(real_kind), intent(in) :: r_particles(npq,nz)
    real(real_kind), intent(out) :: DU(nq,nz), DL(nq,nz), DD(nq,nz)
    real(real_kind), intent(out) :: ADU(nq,nz), ADL(nq,nz)
    real(real_kind), intent(out) :: wfall(npq,nz)
    real(real_kind), intent(out) :: VH2_esc, VH_esc
    
    real(real_kind) :: eddav_p, eddav_m, denav_p, denav_m, tav_p, tav_m
    real(real_kind) :: bx1x2_p, bx1x2_m, bx1x2_pp, bx1x2_mm, zeta_pp, zeta_mm
    real(real_kind) :: bx1x2
    
    ! for particles
    real(real_kind) :: air_density_pp, air_density_mm
    real(real_kind) :: wfall_pp, wfall_mm
    real(real_kind) :: viscosity_pp, viscosity_mm
    
    integer :: i, j
    
    do i = 1,nz
      do j = 1,npq
        air_density_mm = (den(i)/N_avo)*mubar(i)
        viscosity_mm = dynamic_viscosity_air(T(i))
        wfall(j,i) = fall_velocity(grav(i), r_particles(j,i), particle_density(j), air_density_mm, viscosity_mm) &
                     *slip_correction_factor(r_particles(j,i), den(i))
      enddo
    enddo
  
    do i = 2,nz-1
      eddav_p = sqrt(edd(i)*edd(i+1))
      eddav_m = sqrt(edd(i)*edd(i-1))
      denav_p = sqrt(den(i)*den(i+1))
      denav_m = sqrt(den(i)*den(i-1))
      tav_p = sqrt(T(i)*T(i+1))
      tav_m = sqrt(T(i)*T(i-1))
      
      ! gases
      do j = ng_1,nq
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
      ! ! particles
      do j = 1,npq
      
        DU(j,i) = (eddav_p*denav_p)/(dz(i)**2.d0*den(i))
        DL(j,i) = (eddav_m*denav_m)/(dz(i)**2.d0*den(i))
        DD(j,i) = - DU(j,i) - DL(j,i)
      
        air_density_pp = (den(i+1)/N_avo)*mubar(i+1)
        viscosity_pp = dynamic_viscosity_air(T(i+1))
        wfall_pp = fall_velocity(grav(i+1), r_particles(j,i+1), particle_density(j), air_density_pp, viscosity_pp) &
                   *slip_correction_factor(r_particles(j,i+1), den(i+1))
      
        air_density_mm = (den(i-1)/N_avo)*mubar(i-1)
        viscosity_mm = dynamic_viscosity_air(T(i-1))
        wfall_mm = fall_velocity(grav(i-1), r_particles(j,i-1), particle_density(j), air_density_mm, viscosity_mm) &
                   *slip_correction_factor(r_particles(j,i-1), den(i-1))
      
        ADU(j,i) = wfall_pp*den(i+1)/(2.d0*dz(i)*den(i))
        ADL(j,i) = -wfall_mm*den(i-1)/(2.d0*dz(i)*den(i))
      
      enddo
    enddo

    ! surface layer
    eddav_p = sqrt(edd(1)*edd(2))
    denav_p = sqrt(den(1)*den(2))
    tav_p = sqrt(T(1)*T(2))
    do j = ng_1,nq
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
    ! particles
    do j = 1,npq
      DU(j,1) = (eddav_p*denav_p)/(dz(1)**2.d0*den(1))
      DD(j,1) = - DU(j,1)
      
      air_density_pp = (den(2)/N_avo)*mubar(2)
      viscosity_pp = dynamic_viscosity_air(T(2))
      wfall_pp = fall_velocity(grav(2), r_particles(j,2), particle_density(j), air_density_pp, viscosity_pp) &
                  *slip_correction_factor(r_particles(j,2), den(2))
      
      ADU(j,1) = wfall_pp*den(2)/(2.d0*dz(1)*den(1))
    enddo
    

    ! top layer
    eddav_m = sqrt(edd(nz)*edd(nz-1))
    denav_m = sqrt(den(nz)*den(nz-1))
    tav_m = sqrt(T(nz)*T(nz-1))
    do j = ng_1,nq      
      bx1x2_m = binary_diffusion_param(species_mass(j), mubar(nz), tav_m)
      
      ! DU(j,nz) = 0.d0
      DL(j,nz) = (eddav_m*denav_m + bx1x2_m)/(dz(nz)**2.d0*den(nz))
      DD(j,nz) = - DL(j,nz)
            
      bx1x2_mm = binary_diffusion_param(species_mass(j), mubar(nz-1), T(nz-1))
      zeta_mm =  bx1x2_mm*((species_mass(j)*grav(nz-1))/(k_boltz*T(nz-1)*N_avo) &
                         - (mubar(nz-1)*grav(nz-1))/(k_boltz*T(nz-1)*N_avo) &
                         + 0.d0) ! zeroed out thermal diffusion    
      ADL(j,nz) = - zeta_mm/(2.d0*dz(nz)*den(nz))
      ! ADU(j,nz) = 0.d0
    enddo
    ! particles
    do j = 1,npq
      DL(j,nz) = (eddav_m*denav_m)/(dz(nz)**2.d0*den(nz))
      DD(j,nz) = - DL(j,nz)
      
      air_density_mm = (den(nz-1)/N_avo)*mubar(nz-1)
      viscosity_mm = dynamic_viscosity_air(T(nz-1))
      wfall_mm = fall_velocity(grav(nz-1), r_particles(j,nz-1), particle_density(j), air_density_mm, viscosity_mm) &
                  *slip_correction_factor(r_particles(j,nz-1), den(nz-1))
    
      ADL(j,nz) = -wfall_mm*den(nz-1)/(2.d0*dz(nz)*den(nz))
    enddo

    ! H2 escape
    if (diff_H_escape) then
      if (back_gas_name /= "H2") then
        bx1x2 = binary_diffusion_param(species_mass(LH2), mubar(nz), T(nz))
        VH2_esc = bx1x2/den(nz)*(-(species_mass(LH2)*grav(nz))/(k_boltz*T(nz)*N_avo) &
                                 + (mubar(nz)*grav(nz))/(k_boltz*T(nz)*N_avo))                     
      endif
      bx1x2 = binary_diffusion_param(species_mass(LH), mubar(nz), T(nz))
      VH_esc = bx1x2/den(nz)*(-(species_mass(LH)*grav(nz))/(k_boltz*T(nz)*N_avo) &
                              + (mubar(nz)*grav(nz))/(k_boltz*T(nz)*N_avo))
    endif
      
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

end module