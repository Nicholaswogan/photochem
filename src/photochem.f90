

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
                                         * eff_den(j) ! we multiply density here!
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
        
        ! compute rate (effective density is included)
        do j = 1,nz
          rx_rates(j,i) = falloff_rate(kinf(j), Pr(j), F(j)) &
                          * eff_den(j) ! we multiply density here  
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
          ! print*,dexp(-Dg_forward/(Rgas * temperature(j)))
          Dg_forward = max(-1000000.0d0,Dg_forward)
          rx_rates(j,i) = rx_rates(j,n) * &
                          (1.d0/dexp(-Dg_forward/(Rgas * temperature(j)))) * &
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
    real(real_kind) :: enthalpy, entropy, TT
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
    entropy = coeffs(1)*dlog(TT) + coeffs(2)*TT &
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
    real(real_kind), intent(out) :: xp(nz), xl(nz) ! molecules/cm3/s
    
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
  

  subroutine photorates(nz, nsp, kj, nw, dz, densities, xs_x_qy, &
                        flux, diurnal_fac, u0, Rsfc, &
                        prates, surf_radiance,err)
    use photochem_radtran, only: two_stream
    use photochem_data, only: photonums, reactants_sp_inds, nray, sigray, raynums, wavl

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
    !$omp       & tausg, taua, tau, w0, amean, surf_rad, &
    !$omp       & amean_grd, flx)
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
        w0 = min(0.99999d0,tausg(i)/tau(i))
      enddo
      
      call two_stream(nz, tau, w0, u0, Rsfc, amean, surf_rad, ie)

      
      
      surf_radiance(l) = surf_rad
      ierr = ierr + ie
      
      ! convert amean to photolysis grid
      do i = 1,nz
        n = nz+1-i
        amean_grd(i) = dsqrt(amean(n)*amean(n+1))
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


  end subroutine
  
  
  subroutine rhs_background_gas(neqs, usol_flat, rhs, err)
    use photochem_const, only: pi
    
    use photochem_data, only: nq, nsp, nsl, nrT, kj, nw, species_mass, back_gas_mu, &
                              wavl, photonums, species_names
    use photochem_vars, only: nz, temperature, grav, dz, surface_pressure, &
                              xs_x_qy, photon_flux, diurnal_fac, solar_zenith, &
                              surface_albedo
  
    integer, intent(in) :: neqs
    real(real_kind), intent(in), target :: usol_flat(neqs)
    real(real_kind), intent(out) :: rhs(neqs)
    character(len=err_len), intent(out) :: err
    
    real(real_kind), pointer :: usol(:,:)
    
    real(real_kind) :: sum_usol(nz)
    real(real_kind) :: mubar(nz), pressure(nz), density(nz)
    real(real_kind) :: densities(nsp+1,nz)
    real(real_kind) :: rx_rates(nz,nrT)
    real(real_kind) :: prates(nz,kj), surf_radiance(nw)
    real(real_kind) :: xp(nz), xl(nz)
    
    real(real_kind) :: u0
    
    integer :: i, k, j
    
    err = ''
    ! reshape usol_flat with a pointer (no copying; same memory)
    usol(1:nq,1:nz) => usol_flat
    
    
    do i = 1,nz
      sum_usol(i) = sum(usol(:,i))
      if (sum_usol(i) > 1.0d0) then
        err = 'Mixing ratios sum to >1.0 at some altitude (should be <=1).' // &
              ' The atmosphere is probably in a run-away state.'
        return
      endif
    enddo
    
    do i = 1,nz
      call molar_weight(nq, usol(:,i), sum_usol(i), species_mass, back_gas_mu, mubar(i))
    enddo
    
    
    call press_and_den(nz, temperature, grav, surface_pressure*1.d6, dz, &
                       mubar, pressure, density)
    
    do i = 1,nq
      densities(i,:) = usol(i,:)*density
    enddo
    densities(nsp,:) = (1.d0-sum_usol)*density ! background gas
    densities(nsp+1,:) = 1.d0 ! for hv
    
    call reaction_rates(nsp, nz, nrT, temperature, density, &
                        densities, rx_rates, err)
    if (len_trim(err) /= 0) return

    u0 = dcos(solar_zenith*pi/180.d0)
    call photorates(nz, nsp, kj, nw, dz, densities, xs_x_qy, &
                    photon_flux, diurnal_fac, u0, surface_albedo, &
                    prates, surf_radiance, err)
    do i = 1,kj
      k = photonums(i)
      rx_rates(:,k) = prates(:,i) 
    enddo 

    ! short lived (need to fix)
    do i = nq+1,nq+nsl
      call chempl(nz, nsp, nrT, densities, rx_rates, i, xp, xl)
      densities(i,:) = xp/(xl/density)
    enddo

    ! long lived              
    do i = 1,nq
      call chempl(nz, nsp, nrT, densities, rx_rates, i, xp, xl)
      do j = 1,nz
        k = i + (j - 1) * nq
        rhs(k) = (xp(j) - xl(j))/density(j)
      enddo
    enddo

    ! diffusion
    ! call diffusion_coeffs()
  
  
  end subroutine
  
  ! approximate jacobian when there is a background gas.
  ! subroutine jac_backrnd_gas(neq, usol_flat, jac, err)
  ! 
  ! end subroutine
  
  ! subroutine photochem_equil_backrnd_gas
  !   use photochem_vars, only: usol_init
  !   ! the stuff that DOES NOT depend on usol
  !   ! problem dimension CAN NOT change.
  !   call stuff_before_integration(T) !etc.
  !   ! Solve the initial value problem
  !   call solve_ivp()
  !   call stuff_after_integration() ! we set up the output.
  ! end subroutine
  
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
    pressure(1) = Psurf * dexp(-((mubar(1) * grav(1))/(N_avo * k_boltz * T_temp)) * 0.5d0 * dz(1))
    density(1) = pressure(1)/(k_boltz * T(1))
    ! other layers
    do i = 2,nz
      T_temp = (T(i) + T(i-1))/2.d0
      pressure(i) = pressure(i-1) * dexp(-((mubar(i) * grav(i))/(N_avo * k_boltz * T_temp))* dz(i))
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
  
  function arrhenius_rate(A, b, Ea, T) result(k)
    real(real_kind), intent(in) :: A, b, Ea, T
    real(real_kind) :: k
    k = A * T**b * dexp(-Ea/T)
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
    
    log10Fcent = dlog10((1.d0-A)*dexp(-T/T3) + A*dexp(-T/T1))
    C = -0.4d0 - 0.67d0*log10Fcent
    N = 0.75d0 - 1.27d0*log10Fcent
    f1 = (dlog10(Pr) + C)/(N - 0.14d0*(dlog10(Pr + C)))
    F = 10.d0**((log10Fcent)/(1.d0 + f1**2.d0))
  end function

  function Troe_withT2(A, T1, T2, T3, T, Pr) result(F)
    real(real_kind), intent(in) :: A, T1, T2, T3, T, Pr
    real(real_kind) :: F
    
    real(real_kind) :: log10Fcent, f1, C, N
    
    log10Fcent = dlog10((1.d0-A)*dexp(-T/T3) + A*dexp(-T/T1) + dexp(-T2/T))
    C = -0.4d0 - 0.67d0*log10Fcent
    N = 0.75d0 - 1.27d0*log10Fcent
    f1 = (dlog10(Pr) + C)/(N - 0.14d0*(dlog10(Pr + C)))
    F = 10.d0**((log10Fcent)/(1.d0 + f1**2.d0))
  end function
  
  
    
end module