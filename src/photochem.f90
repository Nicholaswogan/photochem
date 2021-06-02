

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
        w0(i) = min(0.99999d0,tausg(i)/tau(i))
      enddo
      
      call two_stream(nz, tau, w0, u0, Rsfc, amean, surf_rad, ie)
      surf_radiance(l) = surf_rad
      ierr = ierr + ie
      do i = 1, nz+1
        amean(i) = abs(amean(i))
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

  end subroutine
  
  
  subroutine rhs_background_gas(neqs, usol_flat, rhs, err)
    use photochem_const, only: pi
    
    use photochem_data, only: nq, nsp, nsl, nrT, kj, nw, species_mass, back_gas_mu, &
                              photonums, water_sat_trop
    use photochem_vars, only: nz, temperature, grav, z, dz, edd, surface_pressure, &
                              xs_x_qy, photon_flux, diurnal_fac, solar_zenith, &
                              surface_albedo, trop_ind, &
                              lowerboundcond, upperboundcond, lower_vdep, lower_flux, &
                              lower_dist_height, lower_fix_mr, upper_veff, upper_flux
  
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
    real(real_kind) :: DU(nq,nz), DD(nq,nz), DL(nq,nz), ADU(nq,nz), ADL(nq,nz)
    
    real(real_kind) :: u0, disth, ztop, ztop1    
    integer :: i, k, j, jdisth
    
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
    
    ! short lived
    do i = nq+1,nq+nsl
      call chempl(nz, nsp, nrT, densities, rx_rates, i, 2, xp, xl) 
      densities(i,:) = xp/xl
    enddo
    
    ! long lived              
    do i = 1,nq
      call chempl(nz, nsp, nrT, densities, rx_rates, i, 1, xp, xl)
      do j = 1,nz
        k = i + (j - 1) * nq
        rhs(k) = (xp(j) - xl(j))/density(j)
      enddo
    enddo
    
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
                        + DD(i,1)*usol(i,1) &
                        - lower_vdep(i)*usol(i,1)/dz(1)
      elseif (lowerboundcond(i) == 1) then
        rhs(i) = 0.d0
      else ! (lowerboundcond(i) == 2) then
        rhs(i) = rhs(i) + DU(i,1)*usol(i,2) + ADU(i,1)*usol(i,2) &
                        + DD(i,1)*usol(i,1) &
                        + lower_flux(i)/(density(1)*dz(1))
      endif
    enddo

    ! Upper boundary
    do i = 1,nq
      k = i + (nz-1)*nq
      if (upperboundcond(i) == 0) then
        rhs(k) = rhs(k) + DD(i,nz)*usol(i,nz) &
                        + DL(i,nz)*usol(i,nz-1) + ADL(i,nz)*usol(i,nz-1) &
                        - upper_veff(i)*usol(i,nz)/dz(nz)    
      elseif (upperboundcond(i) == 2) then
        rhs(k) = rhs(k) + DD(i,nz)*usol(i,nz) &
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
        do j=2,jdisth
          k = i + (j-1)*nq
          rhs(k) = rhs(k) + 2.d0*lower_flux(i)*(ztop1-z(j))/(density(j)*ztop**2.d0)
        enddo
      endif
    enddo  
    
    ! if (water_sat_trop) then ! compute it...
    !   do j=1,trop_ind
    ! 
    !   enddo
    ! endif
    
  end subroutine
  
  integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
                                result(ierr) bind(c, name='RhsFn')
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use photochem_vars, only: neqs
    
    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined dat
    
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)
    character(len=err_len) :: err
    
    ierr = 0
    
    ! get data arrays from SUNDIALS vectors
    yvec => FN_VGetArrayPointer(sunvec_y)
    fvec => FN_VGetArrayPointer(sunvec_f)
    
    ! fill RHS vector
    call rhs_background_gas(neqs, yvec, fvec, err)
    print*, tn
    if (len_trim(err) /= 0) then
      ierr = -1
    endif
    return
  end function
  
  subroutine photo_equilibrium(err)
    use photochem_data, only: nq
    use photochem_vars, only: nz, neqs, usol_init
    
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, FCVodeInit, FCVodeSStolerances, &
                          FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
    use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    
    ! in/out
    character(len=err_len), intent(out) :: err
    
    ! local variables
    real(c_double) :: tstart     ! initial time
    real(c_double) :: tend       ! final time
    real(c_double) :: rtol, atol ! relative and absolute tolerance
    real(c_double) :: dtout      ! output time interval
    real(c_double) :: tout       ! output time
    real(c_double) :: tcur(1)    ! current time
    integer(c_int) :: ierr       ! error flag from C functions
    integer(c_int) :: nout       ! number of outputs
    integer(c_int) :: outstep    ! output loop counter
    type(c_ptr)    :: cvode_mem  ! CVODE memory
    type(N_Vector), pointer :: sunvec_y ! sundials vector
    
    ! solution vector, neq is set in the ode_mod module
    real(c_double) :: yvec(neqs)
    integer(c_long) :: neqs_long
    integer(c_long) :: mu, ml
    type(SUNMatrix), pointer :: sunmat
    type(SUNLinearSolver), pointer :: sunlin
    integer :: i, j, k
    
    ! settings
    neqs_long = neqs
    tstart = 0.0d0
    tcur   = tstart
    tout   = 1.d0
    rtol = 1.d-3
    atol = 1.d-27
    mu = nq
    ml = nq
    
    ! initialize solution vector
    DO i=1,nq
      DO j=1,nz
        k = i + (j-1)*nq
        yvec(k) = usol_init(i,j)
      enddo
    enddo
    
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
    
    ierr = FCVode(cvode_mem, tout, sunvec_y, tcur, CV_NORMAL)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    endif
    
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
        bx1x2_p = 1.52d18*(1.d0/species_mass(j)+1.d0/mubar(i))**0.5d0*(tav_p**0.5)
        bx1x2_m = 1.52d18*(1.d0/species_mass(j)+1.d0/mubar(i))**0.5d0*(tav_m**0.5)
  
        DU(j,i) = (eddav_p*denav_p + bx1x2_p)/(dz(i)**2.d0*den(i))
        DL(j,i) = (eddav_m*denav_m + bx1x2_m)/(dz(i)**2.d0*den(i))
        DD(j,i) = - DU(j,i) - DU(j,i)
        
        bx1x2_pp = 1.52d18*(1.d0/species_mass(j)+1.d0/mubar(i+1))**0.5d0*(T(i+1)**0.5)
        bx1x2_mm = 1.52d18*(1.d0/species_mass(j)+1.d0/mubar(i-1))**0.5d0*(T(i-1)**0.5)
        
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
      bx1x2_p = 1.52d18*(1.d0/species_mass(j)+1.d0/mubar(1))**0.5d0*(tav_p**0.5)
      DU(j,1) = (eddav_p*denav_p + bx1x2_p)/(dz(1)**2.d0*den(1))
      DD(j,1) = - DU(j,1)
      ! DL(j,1) = 0.d0
      
      bx1x2_pp = 1.52d18*(1.d0/species_mass(j)+1.d0/mubar(2))**0.5d0*(T(2)**0.5)
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
      bx1x2_m = 1.52d18*(1.d0/species_mass(j)+1.d0/mubar(nz))**0.5d0*(tav_m**0.5)
      ! DU(j,nz) = 0.d0
      DL(j,nz) = (eddav_m*denav_m + bx1x2_m)/(dz(nz)**2.d0*den(nz))
      DD(j,nz) = - DL(j,nz)
      
      bx1x2_mm = 1.52d18*(1.d0/species_mass(j)+1.d0/mubar(nz-1))**0.5d0*(T(nz-1)**0.5)
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

end module