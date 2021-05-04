

module photochem
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  
contains
  
  subroutine compute_reaction_rates(temperature, density, nz, nrT, reaction_rates, err)
    use photochem_data, only: rateparams, rxtypes, nreactants, &
                              nproducts, reactants_sp_inds, products_sp_inds, &
                              reverse_info, nrF, reverse ! all protected vars
    use photochem_const, only: Rgas, k_boltz ! constants
    use photochem_wrk, only: real_nz_nsp ! pre-allocated work array
                              
    real(real_kind), intent(in) :: temperature(nz)
    real(real_kind), intent(in) :: density(nz)
    integer, intent(in) :: nz, nrT
    real(real_kind), intent(out) :: reaction_rates(nz, nrT)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k, n, l, m
    real(real_kind) :: Troe, A1, A2
    real(real_kind) :: gibbR_forward, gibbP_forward
    real(real_kind) :: Dg_forward
    err = ''
    
    do i = 1,nrF
      if (rxtypes(i) == "falloff") then
        do j = 1,nz
          Troe = 1.d0 ! Here we must compute troe.
          A1  = rateparams(1,i) * temperature(j)**rateparams(2,i) &
                * dexp(-rateparams(3,i)/temperature(j)) * density(j)
          B1 = rateparams(4,i) * temperature(j)**rateparams(5,i) &
                * dexp(-rateparams(6,i)/temperature(j))
          reaction_rates(j,i) = A1 * Troe/(A1/B1 + 1.d0)
        enddo
      else ! three-body or elementary are the same
        do j = 1,nz 
          reaction_rates(j,i) = rateparams(1,i) * temperature(j)**rateparams(2,i) &
                                * dexp(-rateparams(3,i)/temperature(j))
        enddo
      endif
    enddo
    
    if (reverse) then ! if there are reverse reactions
      ! compute gibbs energy at all altitudes
      call compute_gibbs_energy(temperature, nz, nsp, real_nz_nsp, err)
      if (len_trim(err) /= 0) return
      ! compute reverse rate
      do i = nrF+1,nrT
        n = reverse_info(i) ! Reaction number of the forward
        l = nreactants(n) ! number of reactants for the forward reaction
        m = nproducts(n) ! number of products for the forward reaction
        if (rxtypes(n) /= "elementary") then
          l = l - 1 ! "M" doesn't count as a reactant for three-body and fall-off
          m = m - 1
        endif
        do j = 1,nz
          gibbR_forward = 0.d0
          do k = 1,l
            gibbR_forward = gibbR_forward + real_nz_nsp(reactants_sp_inds(k,n),j)
          enddo
          gibbP_forward = 0.d0
          do k = 1,m
            gibbP_forward = gibbP_forward + real_nz_nsp(products_sp_inds(k,n),j)
          enddo
          Dg_forward = gibbsP_forward - gibbsRf_forward ! DG of the forward reaction (J/mol)
          ! compute the reverse rate
          reaction_rates(j,i) = reaction_rates(j,n) * &
                                (1.d0/dexp(-Dg_forward/(Rgas * temperature(j)))) * &
                                (k_boltz*temperature(j)/1.d6)**(m-l)
        enddo
      enddo
    endif

  end subroutine
  
  subroutine compute_gibbs_energy(temperature, nz, nsp, gibbs_energy, err)
    use photochem_data, only: thermo_data, thermo_temps, species_names
    
    real(real_kind), intent(in) :: temperature(nz)
    integer, intent(in) :: nz, nsp
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
  
  subroutine chempl(nz, nsp, nrT, densities, reaction_rates, k, xp, xl)
    use photochem_data, only: nump, numl, iprod, iloss, &
                              reactants_sp_inds, nreactants
    
    ! input
    integer, intent(in) :: nz, nsp, nrT
    real(real_kind), intent(in) :: densities(nz,nsp) ! molecules/cm3 of each species
    real(real_kind), intent(in) :: reaction_rates(nz,nrT) ! reaction rates (various units)
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
      l = numreactants(m) ! l is the number of reactants
      do j = 1,nz
        DD = 1.d0
        do ii = 1,l
          iii = reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xp(j) = xp(j) + reaction_rates(j,m) * DD
      enddo
    enddo
    
    nl = numl(k) ! k is a species
    ! nl is number of reactions that destroy species k
    do i=1,nl
      m = iloss(i,k) ! This will JUST be reaction number
      l = numreactants(m) ! number of reactants
      do j = 1,nz
        DD = 1.d0
        do ii = 1,l
          iii = reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xl(j) = xl(j) + reaction_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
  ! we must pass EVERYTHING into radiative transfer (no globals allowed)
  subroutine compute_photolysis_rates(nz, nsp, kj, nw, xs_x_qy, densities, wavl, flux, &
                                      usol, prates)
                                      
    ! input
    integer, intent(in) :: nz, nsp, kj, nw
    real(real_kind), intent(in) :: xs_x_qy(nz,kj,nw)
    real(real_kind), intent(in) :: densities(nz,nsp)
    real(real_kind), intent(in) :: wavl(nw+1)
    real(real_kind), intent(in) :: usol(nsp,nz)
    
    ! output
    real(real_kind), intent(out) :: prates(nz,kj)
    
    ! local
    real(real_kind) :: partial_prates(nz,kj)
    integer :: i, j, l
    
    ! initialize
    prates = 0.d0
    
    call rayleigh_prep() ! Need to make this into input files
    !$omp parallel
    !$omp private(partial_prates, flx, S)
    partial_prates = 0.d0
    !$omp do
    do l = 1,nw
      call compute_rayleigh()
      call compute_tau()
      call two_stream()
      flx = flux(L)*agl*alp
      do i=1,kj
        do j=1,nz
          partial_prates(j,i) = partial_prates(j,i) + flx*partial_prates(j,i,l)*s(j)
        enddo
      enddo
    enddo
    !$omp enddo
    !$omp critical
    prates = prates + real_nz_kj
    !$omp end critical
    !$omp end parallel

  end subroutine
  
  subroutine right_hand_side(nsp, nq, usol, neq, rhs)
  
    call compute_density(usol,T, den)
  
    call compute_reaction_rates(temperature, density, nz, nrT, reaction_rates, err)
    call compute_photolysis_rates(nz, nsp, kj, nw, xs_x_qy, densities, wavl, flux, &
                                  usol, prates)
    reaction_rates = prates
    do i = 1,nsp
      call chempl(nz, nsp, nrT, densities, reaction_rates, i, xp, xl)
      do j = 1,nz
        k = i + (j - 1) * nsp
        rhs(k) = sp(j)/den(j) - xl(j)/den(j)
      enddo
    enddo
  
    ! compute other terms
  
  end subroutine
  
  subroutine photochem_equilibrium
    use photochem_vars, only: usol_init
    ! the stuff that DOES NOT depend on usol
    ! problem dimension CAN NOT change.
    call stuff_before_integration(T) !etc.
    ! Solve the initial value problem
    call solve_ivp()
    call stuff_after_integration() ! we set up the output.
  end subroutine
  
  
  
  subroutine compute_PandD(nsp, nz, usol, temperature, gravity, Psurf, dz, &
                          mubar, pressure, density)
    use photochem_const, only: k_boltz
    
    integer, intent(in) :: nq, nz
    real(real_kind), intent(in) :: usol(nq,nz)
    real(real_kind), intent(in) :: temperature(nz), gravity
    real(real_kind), intent(in) :: Psurf, dz(nq), mubar(nq)
    
    real(real_kind), intent(out) :: pressure(nz)
    real(real_kind), intent(out) :: density(nz)
     
    real(real_kind) :: T_temp
    
    ! first layer
    ! T_temp = (T(2) - T(1))/dz(1) ! linear extrapolation
    P(1) = Psurf * dexp(((mubar(1) * gravity)/(k_boltz * T_temp)) * 0.5d0 * dz(1))
    density(1) = pressure(1)/(k_boltz * T(i))
    ! other layers
    do i = 2,nz
      T_temp = (T(i) + T(i-1))/2.d0
      pressure(i) = pressure(i-1) * dexp(((mubar(i) * gravity)/(k_boltz * T_temp)) * 0.5d0 * dz(i))
      density(i) = pressure(i)/(k_boltz * T(i))
    enddo
  
  end subroutine
  
  subroutine molar_weight_z(nq, nz, usol, masses, background_mu, mubar)
    integer, intent(in) :: nq, nz
    real(real_kind), intent(in) :: usol(nq,nz)
    real(real_kind), intent(in) :: masses(nq)
    real(real_kind), intent(in) :: background_mu
    real(real_kind), intent(out) :: mubar(nz)
    integer :: i
    
    do i = 1,nz
      call molar_weight(nq, usol(:,i), masses, background_mu, mubar(i))
    enddo
    
  end subroutine
  
  subroutine molar_weight(nq, usol_layer, masses, background_mu, mubar_layer)
    implicit none
    integer, intent(in) :: nq
    real(real_kind), intent(in) :: usol_layer(nq)
    real(real_kind), intent(in) :: masses(nq)
    real(real_kind), intent(in) :: background_mu
    real(real_kind), intent(out) :: mubar_layer
    integer :: j
    real(real_kind) :: f_background

    mubar_layer = 0.d0
    do j = 1, nq
      mubar_layer = mubar_layer + usol_layer(j) * masses(j)
    enddo
    f_background = 1 - sum(usol_layer)
    mubar_layer = mubar_layer + f_background * background_mu
    
  end subroutine
    
  
end module