submodule(photochem_atmosphere) photochem_atmosphere_rhs
  implicit none
  
  ! Contains routines to compute the right-hand-side and jacobian
  ! of the system of ODEs describing photochemistry. There are two main components
  
  ! prep_all_background_gas - computes the reaction rates, photolysis rates, 
  ! diffusion coefficients, etc.
  
  ! dochem - computes the chemistry contribution to the right-hand-side
  
contains
  
  subroutine reaction_rates(self, nsp, nz, nrT, density, densities, rx_rates, err)
    use photochem_eqns, only: arrhenius_rate, Troe_noT2, Troe_withT2, falloff_rate
    use photochem_const, only: Rgas, k_boltz, smallest_real ! constants
    
    class(Atmosphere), intent(in) :: self
    integer, intent(in) :: nsp, nz, nrT
    real(real_kind), intent(in) :: density(nz)
    real(real_kind), intent(in) :: densities(nsp+1,nz)
    real(real_kind), intent(out) :: rx_rates(nz, nrT)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k, n, l, m
    real(real_kind) :: gibbs_energy(nz,self%dat%ng)
    real(real_kind) :: eff_den(nz), F(nz), k0, kinf(nz), Pr(nz)
    real(real_kind) :: gibbR_forward, gibbP_forward
    real(real_kind) :: Dg_forward
    err = ''
    
    do i = 1,self%dat%nrF
      if (self%dat%rxtypes(i) == 1) then ! elementary        
        do j = 1,nz 
          rx_rates(j,i) = arrhenius_rate(self%dat%rateparams(1,i), self%dat%rateparams(2,i), &
                                         self%dat%rateparams(3,i), self%var%temperature(j))
          
        enddo        
      elseif (self%dat%rxtypes(i) == 2) then ! three-body
        n = self%dat%num_efficient(i)
        do j = 1,self%var%nz
          eff_den(j) = density(j) * self%dat%def_eff(i)
          ! subtract the default efficiency, then add the perscribed one
          do k = 1,n ! if n is 0 then it will be skipped
            l = self%dat%eff_sp_inds(k,i)
            eff_den(j) = eff_den(j) - self%dat%def_eff(i)*densities(l,j) &
                                    + self%dat%efficiencies(k,i)*densities(l,j)
          enddo
          rx_rates(j,i) = arrhenius_rate(self%dat%rateparams(1,i), self%dat%rateparams(2,i), &
                                         self%dat%rateparams(3,i), self%var%temperature(j)) &
                                         * eff_den(j) ! we multiply density here
        enddo
        
      
      elseif (self%dat%rxtypes(i) == 3) then ! falloff
        ! compute eff_den, kinf, and Pr at all altitudes
        n = self%dat%num_efficient(i)
        do j = 1,nz
          eff_den(j) = density(j) * self%dat%def_eff(i)
          ! subtract the default efficiency, then add the perscribed one
          do k = 1,n
            l = self%dat%eff_sp_inds(k,i)
            eff_den(j) = eff_den(j) - self%dat%def_eff(i)*densities(l,j) &
                                    + self%dat%efficiencies(k,i)*densities(l,j)
          enddo
          k0  = arrhenius_rate(self%dat%rateparams(1,i), self%dat%rateparams(2,i), &
                               self%dat%rateparams(3,i), self%var%temperature(j))
          kinf(j) = arrhenius_rate(self%dat%rateparams(4,i), self%dat%rateparams(5,i), &
                                   self%dat%rateparams(6,i), self%var%temperature(j))
          kinf(j) = max(kinf(j),smallest_real)
          Pr(j) = k0*eff_den(j)/kinf(j)                        
        enddo
        
        ! compute falloff function
        if (self%dat%falloff_type(i) == 0) then ! no falloff function
          F = 1.d0
        elseif (self%dat%falloff_type(i) == 1) then ! Troe falloff without T2
          do j = 1,nz
            F(j) = Troe_noT2(self%dat%rateparams(7,i), self%dat%rateparams(8,i), &
                             self%dat%rateparams(10,i), self%var%temperature(j), Pr(j))
          enddo
        elseif (self%dat%falloff_type(i) == 2) then ! Troe falloff with T2
          do j = 1,nz
            F(j) = Troe_withT2(self%dat%rateparams(7,i), self%dat%rateparams(8,i), self%dat%rateparams(9,i), &
                               self%dat%rateparams(10,i), self%var%temperature(j), Pr(j))
          enddo
        elseif (self%dat%falloff_type(i) == 3) then ! JPL falloff function
          do j = 1,nz
            F(j) = 0.6d0**(1.d0/(1.d0 + (log10(Pr(j)))**2.d0 ))
          enddo
        endif
        
        ! compute rate
        do j = 1,nz
          rx_rates(j,i) = falloff_rate(kinf(j), Pr(j), F(j))
        enddo
    
      endif
    enddo ! end loop over forward reactions

    if (self%dat%reverse) then ! if there are reverse reactions
      ! compute gibbs energy at all altitudes
      call compute_gibbs_energy(self, self%var%nz, self%dat%ng, gibbs_energy, err)
      if (len_trim(err) /= 0) return
      ! compute reverse rate
      do i = self%dat%nrF+1,self%dat%nrT
        
        n = self%dat%reverse_info(i) ! Reaction number of the forward
        l = self%dat%nreactants(n) ! number of reactants for the forward reaction
        m = self%dat%nproducts(n) ! number of products for the forward reaction
        do j = 1,self%var%nz
          gibbR_forward = 0.d0
          do k = 1,l
            gibbR_forward = gibbR_forward + gibbs_energy(j,self%dat%reactants_sp_inds(k,n)-self%dat%npq)
          enddo
          gibbP_forward = 0.d0
          do k = 1,m
            gibbP_forward = gibbP_forward + gibbs_energy(j,self%dat%products_sp_inds(k,n)-self%dat%npq)
          enddo
          Dg_forward = gibbP_forward - gibbR_forward ! DG of the forward reaction (J/mol)
          ! compute the reverse rate
          rx_rates(j,i) = rx_rates(j,n) * &
                          (1.d0/exp(-Dg_forward/(Rgas * self%var%temperature(j)))) * &
                          (k_boltz*self%var%temperature(j)/1.d6)**(m-l)
        enddo
      enddo
    endif

  end subroutine
  
  subroutine compute_gibbs_energy(self, nz, ng, gibbs_energy, err)
    use photochem_eqns, only: gibbs_energy_shomate
    
    class(Atmosphere), target, intent(in) :: self
    integer, intent(in) :: nz, ng
    real(real_kind), intent(out) :: gibbs_energy(nz,ng)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k
    logical :: found
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var

    err = ''
    
    do i = 1,dat%ng
      do j = 1,var%nz
        found = .false.
        do k = 1,dat%thermo_data(i)%ntemps
          if (var%temperature(j) >= dat%thermo_data(i)%temps(k) .and. &
              var%temperature(j) <  dat%thermo_data(i)%temps(k+1)) then
              
            found = .true.
            call gibbs_energy_shomate(dat%thermo_data(i)%data(1:7,k), var%temperature(j), gibbs_energy(j,i))
            exit
          endif
        enddo
        if (.not. found) then
          err = 'The temperature is not within the ranges '// &
                'given for the thermodynamic data for '//trim(dat%species_names(i+dat%npq))
          return
        endif
      enddo
    enddo

  end subroutine
  
  subroutine chempl(self, nz, nsp, nrT, densities, rx_rates, k, xp, xl)
    
    ! input
    class(Atmosphere), intent(in) :: self
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
    
    np = self%dat%nump(k) ! k is a species
    ! np is number of reactions that produce species k
    do i=1,np
      m = self%dat%iprod(i,k) ! m is reaction number
      l = self%dat%nreactants(m) ! l is the number of reactants
      do j = 1,self%var%nz
        DD = 1.d0
        do ii = 1,l
          iii = self%dat%reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xp(j) = xp(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
    nl = self%dat%numl(k) ! k is a species
    ! nl is number of reactions that destroy species k
    do i=1,nl
      m = self%dat%iloss(i,k) ! This will JUST be reaction number
      l = self%dat%nreactants(m) ! number of reactants
      do j = 1,self%var%nz
        DD = 1.d0
        do ii = 1,l
          iii = self%dat%reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xl(j) = xl(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
  subroutine chempl_sl(self, nz, nsp, nrT, densities, rx_rates, k, xp, xl)
    
    ! input
    class(Atmosphere), intent(in) :: self
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
    
    np = self%dat%nump(k) ! k is a species
    ! np is number of reactions that produce species k
    do i=1,np
      m = self%dat%iprod(i,k) ! m is reaction number
      l = self%dat%nreactants(m) ! l is the number of reactants
      do j = 1,self%var%nz
        DD = 1.d0
        do ii = 1,l
          iii = self%dat%reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xp(j) = xp(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
    nl = self%dat%numl(k) ! k is a species
    ! nl is number of reactions that destroy species k
    do i=1,nl
      m = self%dat%iloss(i,k) ! This will JUST be reaction number
      l = self%dat%nreactants(m) ! number of reactants
      do j = 1,self%var%nz
        DD = 1.d0
        do ii = 1,l
          iii = self%dat%reactants_sp_inds(ii,m)
          ! We skip the short-lived species.
          if (iii /= k) then
            DD = DD * densities(iii,j)
          endif
        enddo
        xl(j) = xl(j) + rx_rates(j,m) * DD
      enddo
    enddo
    
  end subroutine
  
  module subroutine right_hand_side_chem(self, usol, rhs, err)
    class(Atmosphere), target, intent(inout) :: self
    real(real_kind), intent(in) :: usol(:,:)
    real(real_kind), intent(out) :: rhs(:)
    character(len=err_len), intent(out) :: err
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    err = ""
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz .or. size(rhs) /= var%neqs) then
      err = "Input usol or rhs to dochem_implicit has the wrong dimensions"
      return
    endif

    call self%prep_atmosphere(usol, err)
    if (len_trim(err) /= 0) return
    
    call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, &
                var%trop_ind, dat%nrT, wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%rainout_rates, &
                wrk%densities, wrk%xp, wrk%xl, rhs)
                              
  end subroutine
  
  subroutine dochem(self, neqs, nsp, np, nsl, nq, nz, trop_ind, nrT, usol, density, rx_rates, &
                    gas_sat_den, molecules_per_particle, &
                    H2O_sat_mix, rainout_rates, &
                    densities, xp, xl, rhs)                 
    use photochem_eqns, only: damp_condensation_rate
    use photochem_const, only: N_avo, pi, small_real
    
    class(Atmosphere), target, intent(in) :: self
    integer, intent(in) :: neqs, nsp, np, nsl, nq, nz, trop_ind, nrT
    real(real_kind), intent(in) :: usol(nq,nz), density(nz)
    real(real_kind), intent(in) :: rx_rates(nz,nrT)
    real(real_kind), intent(in) :: gas_sat_den(np,nz)
    real(real_kind), intent(in) :: molecules_per_particle(np,nz)
    real(real_kind), intent(in) :: H2O_sat_mix(nz)
    real(real_kind), intent(in) :: rainout_rates(nq, trop_ind)
    real(real_kind), intent(inout) :: densities(nsp+1,nz), xp(nz), xl(nz)
    real(real_kind), intent(out) :: rhs(neqs)
    
    real(real_kind) :: dn_gas_dt, dn_particle_dt, H2O_cold_trap, cond_rate
    integer :: i, ii, j, k, kk
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var
    
    do j = 1,var%nz
      do k = 1,dat%np
        densities(k,j) = max(usol(k,j)*(density(j)/molecules_per_particle(k,j)),small_real)
      enddo
      do k = dat%ng_1,dat%nq
        densities(k,j) = usol(k,j)*density(j)
      enddo
      densities(nsp,j) = (1.d0-sum(usol(dat%ng_1:,j)))*density(j) ! background gas
      densities(nsp+1,j) = 1.d0 ! for hv
    enddo
    
    ! short lived
    do k = dat%nq+1,dat%nq+dat%nsl
      call chempl_sl(self, nz, nsp, nrT, densities, rx_rates, k, xp, xl) 
      densities(k,:) = xp/xl
    enddo
    
    rhs = 0.d0
    
    ! long lived              
    do i = dat%ng_1,dat%nq
      call chempl(self, nz, nsp, nrT, densities, rx_rates, i, xp, xl)
      do j = 1,var%nz
        k = i + (j - 1) * dat%nq
        rhs(k) = xp(j)/density(j) - xl(j)/density(j)
      enddo
    enddo
    
    if (var%gas_rainout) then
      ! rainout rates
      do j = 1,var%trop_ind
        do i = 1,dat%nq
          k = i + (j - 1) * dat%nq
          rhs(k) = rhs(k) - rainout_rates(i,j)*usol(i,j)
        enddo
      enddo
    endif
    
    ! if stratospheric water condesation is on, then condense
    ! water in the stratosphere.
    if (dat%stratospheric_cond) then
      if (dat%fix_water_in_trop) then
        ! need to start above tropopause if fixing
        ! water in troposphere
        i = var%trop_ind+1
      else
        i = 1
      endif
      do j = i,var%nz
        k = dat%LH2O + (j - 1) * dat%nq
        H2O_cold_trap = var%H2O_condensation_rate(2)*H2O_sat_mix(j)
        if (usol(dat%LH2O,j) >= H2O_cold_trap) then
          
          cond_rate = damp_condensation_rate(var%H2O_condensation_rate(1), &
                                             var%H2O_condensation_rate(2), &
                                             var%H2O_condensation_rate(3), &
                                             usol(dat%LH2O,j)/H2O_sat_mix(j))
          rhs(k) = rhs(k) - cond_rate*(usol(dat%LH2O,j) - H2O_cold_trap)
           
        endif
      enddo
    endif
    
    if (dat%there_are_particles) then
      ! formation from reaction
      do i = 1,dat%np
        call chempl(self, nz, nsp, nrT, densities, rx_rates, i, xp, xl)
        do j = 1,var%nz
          k = i + (j - 1) * nq
          rhs(k) = rhs(k) + (xp(j) - xl(j))/density(j)
        enddo
      enddo
    
      ! particle condensation
      do j = 1,var%nz
        do i = 1,dat%np
          ! if this particle forms from condensation
          if (dat%particle_formation_method(i) == 1) then
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
              rhs(kk) = rhs(kk) + dn_gas_dt/density(j)            
            
              ! rate of particle production from gas condensing
              ! in particles/cm3/s
              dn_particle_dt = - dn_gas_dt/density(j)
              ! add to rhs vector, convert to moles/cm3/s
              rhs(k) = rhs(k) + dn_particle_dt
            else
              ! particles don't change!
              rhs(k) = rhs(k) + 0.d0
            endif
          endif          
        enddo
      enddo
      
    endif
    
  end subroutine
  
  subroutine photorates(self, nz, nsp, kj, nw, densities, &
                        prates, surf_radiance,err)
    use photochem_radtran, only: two_stream
    !$ use photochem_eqns, only: round
    use photochem_const, only: pi
    ! input
    class(Atmosphere), target, intent(in) :: self
    integer, intent(in) :: nz, nsp, kj, nw
    real(real_kind), intent(in) :: densities(nsp+1, nz)
    
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
    real(real_kind) :: taup_1, gt_1, tausp_1(self%dat%np,nz)
    real(real_kind) :: u0
    integer :: l, i, j, jj, k, n, ie, ierr, ierrs
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var
    
    u0 = cos(var%solar_zenith*pi/180.d0)

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
    do l = 1,dat%nw
      
      ! rayleigh scattering
      tausg = 0.d0
      do i = 1,dat%nray
        j = dat%raynums(i)
        do k = 1,var%nz
          n = var%nz+1-k
          tausg(n) = tausg(n) + dat%sigray(i,l)*densities(j,k)*var%dz(k)
        enddo
      enddo
      
      ! photolysis
      taua = 0.d0
      do i = 1,dat%kj
        jj = dat%photonums(i)
        j = dat%reactants_sp_inds(1,jj)
        do k = 1,var%nz
          n = var%nz+1-k
          taua(n) = taua(n) + var%xs_x_qy(k,i,l)*densities(j,k)*var%dz(k)
        enddo
      enddo
      
      ! particles
      taup = 0.d0
      tausp = 0.d0
      tausp_1 = 0.d0
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
      gt = 0.d0
      do k = 1,var%nz
        n = nz+1-k
        gt_1 = 0.d0
        do i = 1,dat%np    
          if (var%particle_xs(i)%ThereIsData) then
            gt_1 = gt_1 + var%particle_xs(i)%gt(k,l)*tausp_1(i,n) &
                    /(tausp(n)+tausg(n))
          endif
        enddo
        gt(n) = min(gt_1,0.999999d0)
      enddo
      
      ! sum of all contributions
      tau = tausg + taua + taup + tausp
      do i = 1,var%nz
        w0(i) = min(0.99999d0,(tausg(i) + tausp(i))/tau(i))
      enddo
      
      call two_stream(nz, tau, w0, gt, u0, var%surface_albedo, amean, surf_rad, ie)
      surf_radiance(l) = surf_rad
      ierr = ierr + ie
      do i = 1, var%nz+1
        if (amean(i) < -1.d-5) then
          ierr = ierr + 1
        endif
        amean(i) = abs(amean(i))
        ! amean(i) = max(amean(i),0.d0)
      enddo
      
      ! convert amean to photolysis grid
      do i = 1,var%nz
        n = var%nz+1-i
        amean_grd(i) = sqrt(amean(n)*amean(n+1))        
      enddo
      
      flx = var%photon_flux(l)*var%diurnal_fac*var%photon_scale_factor ! photons/cm2/s
      
      do i=1,dat%kj
        do j=1,var%nz
          partial_prates(j,i) = partial_prates(j,i) + flx*amean_grd(j)*var%xs_x_qy(j,i,l)
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
    !$ do i=1,dat%kj
    !$   do j=1,var%nz
    !$     call round(prates(j,i),-8)
    !$   enddo
    !$ enddo
    
  end subroutine
  
  subroutine prep_atm_background_gas(self, nq, nz, trop_ind, sum_usol, usol, &
                                     density, mubar, pressure, fH2O, H2O_sat_mix, err)
    use photochem_eqns, only: sat_pressure_H2O, molar_weight, press_and_den
    use, intrinsic :: iso_c_binding, only : c_loc, c_ptr
    use cminpack2fort, only: hybrd1 ! interface to hybrd1 from cminpack
    
    class(Atmosphere), target, intent(in) :: self
    integer, intent(in) :: nq, nz, trop_ind
    real(real_kind), intent(inout), target :: usol(nq,nz)
    real(real_kind), intent(out) :: sum_usol(nz)
    real(real_kind), intent(out) :: density(nz)
    real(real_kind), intent(out) :: mubar(nz), pressure(nz), fH2O(trop_ind), H2O_sat_mix(nz)
    character(len=err_len), intent(out) :: err
    
    real(real_kind) :: rel
    integer :: i
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(Atmosphere), pointer :: self_ptr
    
    ! hybrd1 varibles for cminpack
    integer :: info
    real(real_kind), parameter :: tol = 1.d-8
    integer :: lwa
    real(real_kind), allocatable :: fvec(:)
    real(real_kind), allocatable :: wa(:)
    type(c_ptr) :: ptr ! void c pointer for cminpack
    ! end hybrd1 varibles
    
    err = ''
    
    dat => self%dat
    var => self%var
    
    ! If water is fixed in the troposhere, lets set it to zero.
    ! This will be to get an initial guess for the non-linear solve.
    if (dat%fix_water_in_trop) then
      do i = 1,var%trop_ind
        usol(dat%LH2O,i) = 0.d0
      enddo
    endif 
    
    do i = 1,var%nz
      sum_usol(i) = sum(usol(dat%ng_1:,i))
      if (sum_usol(i) > 1.0d0) then
        err = 'Mixing ratios sum to >1.0 at some altitude (should be <=1).' // &
              ' The atmosphere is probably in a run-away state'
        return
      endif
    enddo
    
    do i = 1,var%nz
      call molar_weight(dat%nll, usol(dat%ng_1:,i), sum_usol(i), dat%species_mass(dat%ng_1:), dat%back_gas_mu, mubar(i))
    enddo
    
    call press_and_den(var%nz, var%temperature, var%grav, var%surface_pressure*1.d6, var%dz, &
                       mubar, pressure, density)
                       
    if (dat%fix_water_in_trop) then
      ! Compute initial guess for fH2O
      do i = 1,var%trop_ind
        if (var%use_manabe) then
          ! manabe formula
          rel = 0.77d0*(pressure(i)/pressure(1)-0.02d0)/0.98d0
        else
          rel = var%relative_humidity 
        endif
        
        fH2O(i) = rel*sat_pressure_H2O(var%temperature(i))/pressure(i)
      enddo
      ! Here we compute self-consistent water profile.
      ! hybrd1 is nonlinear solver from cminpack. Takes negligable time.
      self_ptr => self
      ptr = c_loc(self_ptr) ! void pointer to usol. Be careful!!!
      lwa = (var%trop_ind*(3*var%trop_ind+13))/2 + 2
      allocate(fvec(trop_ind))
      allocate(wa(lwa))
      ! Call hybrd1 from cminpack (c re-write of minpack). I wrote a fortran
      ! interface to cminpack. We pass usol, the atmosphere, via a pointer.
      call hybrd1(fcn_fH2O, ptr, var%trop_ind, fH2O, fvec, tol, info, wa, lwa)
      if (info /= 1) then
        err = "Non-linear solve for water profile failed."
        return
      endif
      deallocate(fvec, wa)
      
      ! use the solution for tropospheric H2O to compute molar weight
      ! and pressure and density.
      do i = 1,var%trop_ind
        usol(dat%LH2O,i) = fH2O(i)
      enddo
      
      do i = 1,var%nz
        sum_usol(i) = sum(usol(dat%ng_1:,i))
      enddo
      
      do i = 1,var%nz
        call molar_weight(dat%nll, usol(dat%ng_1:,i), sum_usol(i), &
                          dat%species_mass(dat%ng_1:), dat%back_gas_mu, mubar(i))
      enddo
      
      call press_and_den(var%nz, var%temperature, var%grav, var%surface_pressure*1.d6, var%dz, &
                         mubar, pressure, density)
      
    endif 
    
    if (dat%stratospheric_cond) then
      ! compute H2O saturation mixing ratio at all altitudes
      do i = 1,var%nz
        H2O_sat_mix(i) = sat_pressure_H2O(var%temperature(i))/pressure(i)
      enddo
    endif
    
  end subroutine

  ! For computing self-consistent water profile. Called by minpack.
  module function fcn_fH2O(ptr, n, x, fvec, iflag) result(res) bind(c)
    use photochem_eqns, only: sat_pressure_H2O, molar_weight, press_and_den
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr
  
    type(c_ptr) :: ptr ! void pointer. Be careful!
    integer, value :: n, iflag ! n == trop_ind
    real(real_kind), intent(in) :: x(n) ! x == fH2O
    real(real_kind), intent(out) :: fvec(n) ! fvec == residual
    integer :: res
  
    real(real_kind) :: fH2O(n)
    real(real_kind), pointer :: usol(:,:)
    real(real_kind), pointer :: sum_usol(:) 
    real(real_kind), pointer :: mubar(:)
    real(real_kind), pointer :: density(:)
    real(real_kind), pointer :: pressure(:)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(Atmosphere), pointer :: self
  
    integer :: i
    real(real_kind) :: rel
  
    ! dereferences the pointer. Takes the void ptr, then
    ! declares it to be pointing to a usol shaped array of doubles.
    call c_f_pointer(ptr, self)
    
    dat => self%dat
    var => self%var
    sum_usol => self%wrk%sum_usol
    mubar => self%wrk%mubar
    density => self%wrk%density
    pressure => self%wrk%pressure
    usol => self%wrk%usol
  
    do i = 1,n
      usol(dat%LH2O,i) = x(i)
    enddo
  
    do i = 1,var%nz
      sum_usol(i) = sum(usol(dat%ng_1:,i))
    enddo
  
    do i = 1,var%nz
      call molar_weight(dat%nll, usol(dat%ng_1:,i), sum_usol(i), &
                        dat%species_mass(dat%ng_1:), dat%back_gas_mu, mubar(i))
    enddo
  
    call press_and_den(var%nz, var%temperature, var%grav, var%surface_pressure*1.d6, var%dz, &
                       mubar, pressure, density)
  
    do i = 1,n
      if (var%use_manabe) then
        ! manabe formula ()
        rel = 0.77d0*(pressure(i)/pressure(1)-0.02d0)/0.98d0
      else
        rel = var%relative_humidity 
      endif
      fH2O(i) = rel*sat_pressure_H2O(var%temperature(i))/pressure(i)
    enddo  
  
    fvec = x - fH2O
    res = 0
  end function
  
  subroutine diffusion_coefficients(self, nq, npq, nz, den, mubar, &
                                    DU, DD, DL, ADU, ADL, wfall, VH2_esc, VH_esc)
    use photochem_eqns, only: dynamic_viscosity_air, fall_velocity, slip_correction_factor, &
                              binary_diffusion_param
    use photochem_const, only: k_boltz, N_avo
    
    class(Atmosphere), target, intent(in) :: self
    integer, intent(in) :: nq, npq, nz
    real(real_kind), intent(in) :: den(nz)
    real(real_kind), intent(in) :: mubar(nz)
    real(real_kind), intent(out) :: DU(nq,nz), DL(nq,nz), DD(nq,nz)
    real(real_kind), intent(out) :: ADU(nq,nz), ADL(nq,nz)
    real(real_kind), intent(out) :: wfall(npq,nz)
    real(real_kind), intent(out) :: VH2_esc, VH_esc
    
    real(real_kind) :: eddav_p, eddav_m, denav_p, denav_m, tav_p, tav_m
    real(real_kind) :: bx1x2_p, bx1x2_m, zeta_p, zeta_m
    real(real_kind) :: grav_p, grav_m, mubar_p, mubar_m
    real(real_kind) :: bx1x2
    
    ! for particles
    real(real_kind) :: air_density_p, air_density_m
    real(real_kind) :: wfall_p, wfall_m
    real(real_kind) :: viscosity_p, viscosity_m
    real(real_kind) :: radius_p, radius_m
    
    integer :: i, j
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var
    
    ! eddy diffusion. particles and gases
    ! middle
    do i = 2,var%nz-1
      eddav_p = sqrt(var%edd(i)*var%edd(i+1))
      eddav_m = sqrt(var%edd(i)*var%edd(i-1))
      denav_p = sqrt(den(i)*den(i+1))
      denav_m = sqrt(den(i)*den(i-1))
      do j = 1,dat%nq
        DU(j,i) = (eddav_p*denav_p)/(den(i)*var%dz(i)**2.d0)
        DL(j,i) = (eddav_m*denav_m)/(den(i)*var%dz(i)**2.d0)
        DD(j,i) = - DU(j,i) - DL(j,i)
      enddo
    enddo
    ! top and bottom
    eddav_p = sqrt(var%edd(1)*var%edd(2))
    eddav_m = sqrt(var%edd(var%nz)*var%edd(var%nz-1))
    denav_p = sqrt(den(1)*den(2))
    denav_m = sqrt(den(var%nz)*den(var%nz-1))
    do j = 1,dat%nq
      DU(j,1) = (eddav_p*denav_p)/(den(1)*var%dz(1)**2.d0)
      DL(j,var%nz) = (eddav_m*denav_m)/(den(var%nz)*var%dz(var%nz)**2.d0)
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
        
        DU(j,i) = DU(j,i) + bx1x2_p/(var%dz(i)**2.d0*den(i))
        DL(j,i) = DL(j,i) + bx1x2_m/(var%dz(i)**2.d0*den(i))
        DD(j,i) = - DU(j,i) - DL(j,i)
        
        zeta_p =  bx1x2_p*((dat%species_mass(j)*grav_p)/(k_boltz*tav_p*N_avo) &
                           - (mubar_p*grav_p)/(k_boltz*tav_p*N_avo) &
                           + 0.d0) ! zeroed out thermal diffusion   
        zeta_m =  bx1x2_m*((dat%species_mass(j)*grav_m)/(k_boltz*tav_m*N_avo) &
                          - (mubar_m*grav_m)/(k_boltz*tav_m*N_avo) &
                          + 0.d0) ! zeroed out thermal diffusion
                          
        ADU(j,i) = zeta_p/(2.d0*var%dz(i)*den(i)) 
        ADL(j,i) = - zeta_m/(2.d0*var%dz(i)*den(i))
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
      DU(j,i) = DU(j,i) + bx1x2_p/(var%dz(i)**2.d0*den(i))
      zeta_p =  bx1x2_p*((dat%species_mass(j)*grav_p)/(k_boltz*tav_p*N_avo) &
                         - (mubar_p*grav_p)/(k_boltz*tav_p*N_avo) &
                         + 0.d0) ! zeroed out thermal diffusion   
      ADU(j,i) = zeta_p/(2.d0*var%dz(i)*den(i)) 
    enddo
    ! upper boundary
    i = var%nz
    do j = dat%ng_1, dat%nq
      bx1x2_m = binary_diffusion_param(dat%species_mass(j), mubar_m, tav_m)
      
      DL(j,i) = DL(j,i) + bx1x2_m/(var%dz(i)**2.d0*den(i))
      
      zeta_m =  bx1x2_m*((dat%species_mass(j)*grav_m)/(k_boltz*tav_m*N_avo) &
                        - (mubar_m*grav_m)/(k_boltz*tav_m*N_avo) &
                        + 0.d0) ! zeroed out thermal diffusion
                        
      ADL(j,i) = - zeta_m/(2.d0*var%dz(i)*den(i))
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
      
        ADU(j,i) = wfall_p*denav_p/(2.d0*var%dz(i)*den(i))
        ADL(j,i) = - wfall_m*denav_m/(2.d0*var%dz(i)*den(i))
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
    
      ADU(j,i) = wfall_p*denav_p/(2.d0*var%dz(i)*den(i))
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
    
      ADL(j,i) = - wfall_m*denav_m/(2.d0*var%dz(i)*den(i))
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
  
  subroutine rainout(self, nq, nz, trop_ind, usol, den, rainout_rates)
    use photochem_const, only: k_boltz, N_avo, small_real
    use photochem_eqns, only: henrys_law
    
    class(Atmosphere), target, intent(in) :: self
    integer, intent(in) :: nq, nz, trop_ind
    real(real_kind), intent(in) :: usol(nq, nz)
    real(real_kind), intent(in) :: den(nz)
    real(real_kind), intent(out) :: rainout_rates(nq, trop_ind)
    
    integer :: i, j
    
    real(real_kind) :: wH2O(trop_ind)
    real(real_kind) :: slope, intercept
    real(real_kind) :: denav_p, eddav_p, denav_m, eddav_m
    real(real_kind) :: scale_factor
    real(real_kind) :: k_bar, Q_i, H_coeff
    
    real(real_kind), parameter :: earth_rainfall_rate = 1.1d17 ! molecules/cm2/s
    real(real_kind), parameter :: fz = 0.05d0 ! fraction of the time it rains
    real(real_kind), parameter :: gamma = 4.d5 ! average time of storm cycles (s)
    real(real_kind), parameter :: LLL = 1.d0 ! g H2O/m3 of clouds
    real(real_kind), parameter :: C1 = 1.d-6 !dynes/bar
    real(real_kind), parameter :: C2 = 1.d-9 ![m3/cm3][L H2O/g H2O]
    real(real_kind), parameter :: MH2O = 18.d0 ! g H2O/mol H2O
    real(real_kind), parameter :: rho_H2O = 1000.d0 ! g H2O/L H2O
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var
    
    !!!!!!! calculate raining rate !!!!!!!
    ! middle of atmosphere
    wH2O = 0.d0
    do i = 2,var%trop_ind-1
      denav_p = sqrt(den(i+1)*den(i))
      eddav_p = sqrt(var%edd(i+1)*var%edd(i))
      denav_m = sqrt(den(i-1)*den(i))
      eddav_m = sqrt(var%edd(i-1)*var%edd(i))
      wH2O(i) = (eddav_p*denav_p/var%dz(i)**2.d0) * usol(dat%LH2O,i+1) &
              - (eddav_p*denav_p/var%dz(i)**2.d0 + eddav_m*denav_m/var%dz(i)**2.d0) * usol(dat%lH2O,i) &
              + (eddav_m*denav_m/var%dz(i)**2.d0) * usol(dat%lH2O,i-1)
      if (wH2O(i) < 0.d0) then
        wH2O(i) = 1.d-20
      endif
    enddo
    ! lets just linearly extrapolate wH2O to bottom and top grid cell
    !!! lower boundary !!!
    slope = (wH2O(3) - wH2O(2))/(var%dz(2))
    intercept = wH2O(2) - slope*var%z(2)
    wH2O(1) = slope*var%z(1) + intercept
    if (wH2O(1) < 0.d0) then
      wH2O(1) = 1.d-20
    endif
    !!! upper boundary !!!
    slope = (wH2O(var%trop_ind-1) - wH2O(var%trop_ind-2))/(var%dz(var%trop_ind-1))
    intercept = wH2O(var%trop_ind-1) - slope*var%z(var%trop_ind-1)
    wH2O(var%trop_ind) = slope*var%z(var%trop_ind) + intercept
    if (wH2O(var%trop_ind) < 0.d0) then
      wH2O(var%trop_ind) = 1.d-20
    endif
    ! Earth's globally averaged rainfall is 1.1e+17 molecules/cm2/s (Giorgi+1985)
    ! here we rescale the raining rate so that the total integrated value is
    ! the same as Earth's.
    scale_factor = earth_rainfall_rate/sum(wH2O*var%dz(1))
    wH2O = wH2O*scale_factor
    !!!!!!! end calculate raining rate !!!!!!!
    
    !!!!!!! dissolve gas in the rain !!!!!!!!!
    do j = 1,var%trop_ind
      do i = 1,dat%nq
        H_coeff = henrys_law(max(var%temperature(j),273.15d0),dat%henry_data(1,i),dat%henry_data(2,i))*(1.d5)
        H_coeff = max(H_coeff, small_real)
        k_bar = (C1*k_boltz*var%temperature(j)*H_coeff/ &
                (1.d0+C1*C2*N_avo*LLL*k_boltz*var%temperature(j)*H_coeff)) &
                * (WH2O(j)*MH2O/rho_H2O) 
        Q_i = (1.d0-fz) + (fz/(gamma*k_bar))*(1.d0 - exp(-k_bar*gamma))
        rainout_rates(i,j) = (1.d0/(gamma*Q_i)) * (1.d0 - exp(-k_bar*gamma))
      enddo
    enddo
    !!!!!!! end dissolve gas in the rain !!!!!!!!!
    
  end subroutine
  
  module subroutine prep_all_background_gas(self, usol_in, err)
    use photochem_eqns, only: saturation_density
    use photochem_const, only: pi, k_boltz, N_avo, small_real
  
    class(Atmosphere), target, intent(inout) :: self
    real(real_kind), intent(in) :: usol_in(:,:)
    character(len=err_len), intent(out) :: err
  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    integer :: i, j, k
    real(real_kind) :: P_H2SO4
  
    err = ''
  
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    ! make a copy of the mixing ratios. You can not alter the input mixing ratios
    ! this will make CVODE unstable. Guard against division by zero
    do j = 1,var%nz
      do i = 1,dat%nq
        if (usol_in(i,j) < 0.d0) then
          wrk%usol(i,j) = min(usol_in(i,j),-small_real)
        else
          wrk%usol(i,j) = max(usol_in(i,j), small_real)
        endif
      enddo
    enddo

    wrk%upper_veff_copy = var%upper_veff
    wrk%lower_vdep_copy = var%lower_vdep
  
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == 1) then
        wrk%usol(i,1) = var%lower_fix_mr(i)
      endif
    enddo
  
    call prep_atm_background_gas(self, dat%nq, var%nz, var%trop_ind, wrk%sum_usol, wrk%usol, &
                                 wrk%density, wrk%mubar, wrk%pressure, wrk%fH2O, wrk%H2O_sat_mix, err)  
    if (len_trim(err) /= 0) return  
  
    ! diffusion coefficients
    call diffusion_coefficients(self, dat%nq, dat%npq, var%nz, wrk%density, wrk%mubar, &
                                wrk%DU, wrk%DD, wrk%DL, wrk%ADU, wrk%ADL, &
                                wrk%wfall, wrk%VH2_esc, wrk%VH_esc)
  
    ! surface scale height
    wrk%surface_scale_height = (k_boltz*var%temperature(1)*N_avo)/(wrk%mubar(1)*var%grav(1))
  
    ! H and H2 escape
    if (dat%diff_H_escape) then
      if (dat%back_gas_name /= "H2") then
        wrk%upper_veff_copy(dat%LH2) = wrk%VH2_esc                     
      endif
      wrk%upper_veff_copy(dat%LH) = wrk%VH_esc 
    endif
  
    if (dat%there_are_particles) then
      do i = 1,dat%np
        wrk%lower_vdep_copy(i) = wrk%lower_vdep_copy(i) + wrk%wfall(i,1)
        ! wrk%upper_veff_copy(i) = wrk%upper_veff_copy(i) + wrk%wfall(i,var%nz)
      enddo
  
      ! compute The saturation density
      do j = 1,var%nz
        do i = 1,dat%np
          if (dat%particle_formation_method(i) == 1) then
            if (dat%particle_sat_type(i) == 1) then ! arrhenius
              wrk%gas_sat_den(i,j) = saturation_density(var%temperature(j), &
                                                   dat%particle_sat_params(1,i), &
                                                   dat%particle_sat_params(2,i), &
                                                   dat%particle_sat_params(3,i))
            elseif (dat%particle_sat_type(i) == 2) then ! interpolate to H2SO4 data
              call dat%H2SO4_sat%evaluate(var%temperature(j), &
                                          log10(wrk%usol(dat%LH2O,j)*wrk%pressure(j)/1.d6),P_H2SO4)
              P_H2SO4 = 10.d0**(P_H2SO4)
              wrk%gas_sat_den(i,j) = (P_H2SO4*1.d6)/(k_boltz*var%temperature(j))
            endif
          endif
          wrk%molecules_per_particle(i,j) = (4.d0/3.d0)*pi*var%particle_radius(i,j)**3.d0* &
                                            dat%particle_density(i)*(1/dat%species_mass(i))*N_avo
        enddo
      enddo
    endif
  
    ! densities include particle densities
    ! we track mol/cm3 for particles (instead of particles/cm3) 
    ! because it makes the numbers smaller and closer to 1, which
    ! is similar to mixing ratios.
    do j = 1,var%nz
      do i = 1,dat%npq
        wrk%densities(i,j) = max(wrk%usol(i,j)*(wrk%density(j)/wrk%molecules_per_particle(i,j)), small_real)
      enddo
      do i = dat%ng_1,dat%nq
        wrk%densities(i,j) = wrk%usol(i,j)*wrk%density(j)
      enddo
      wrk%densities(dat%nsp,j) = (1.d0-wrk%sum_usol(j))*wrk%density(j) ! background gas
      wrk%densities(dat%nsp+1,j) = 1.d0 ! for hv
    enddo
    
    call reaction_rates(self, dat%nsp, var%nz, dat%nrT, wrk%density, &
                        wrk%densities, wrk%rx_rates, err)
    if (len_trim(err) /= 0) return
  
    call photorates(self, var%nz, dat%nsp, dat%kj, dat%nw, wrk%densities, &
                    wrk%prates, wrk%surf_radiance, err)
    if (len_trim(err) /= 0) return
  
    do i = 1,dat%kj
      k = dat%photonums(i)
      wrk%rx_rates(:,k) = wrk%prates(:,i) 
    enddo 
  
    ! rainout rates
    if (var%gas_rainout) then
      call rainout(self, dat%nq, var%nz, var%trop_ind, wrk%usol, wrk%density, wrk%rainout_rates)
    endif
  
  end subroutine
  
  module subroutine rhs_background_gas(self, neqs, usol_flat, rhs, err)
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real  
    
    class(Atmosphere), target, intent(inout) :: self
    integer, intent(in) :: neqs
    real(real_kind), target, intent(in) :: usol_flat(neqs)
    real(real_kind), intent(out) :: rhs(neqs)
    character(len=err_len), intent(out) :: err
    
    real(real_kind) :: disth, ztop, ztop1    
    integer :: i, k, j, jdisth
    
    real(real_kind), pointer :: usol_in(:,:)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    err = ''
    
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
    call prep_all_background_gas(self, usol_in, err)
    if (len_trim(err) /= 0) return
    
    call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, &
                var%trop_ind, dat%nrT, wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%rainout_rates, &
                wrk%densities, wrk%xp, wrk%xl, rhs) 
    
    ! diffusion (interior grid points)
    do j = 2,var%nz-1
      do i = 1,dat%nq
        k = i + (j-1)*dat%nq
        rhs(k) = rhs(k) + wrk%DU(i,j)*wrk%usol(i,j+1) + wrk%ADU(i,j)*wrk%usol(i,j+1) &
                        + wrk%DD(i,j)*wrk%usol(i,j) + (wrk%ADU(i,j) + wrk%ADL(i,j))*wrk%usol(i,j) &
                        + wrk%DL(i,j)*wrk%usol(i,j-1) + wrk%ADL(i,j)*wrk%usol(i,j-1)
      enddo
    enddo
    
    ! Lower boundary
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == 0 .or. var%lowerboundcond(i) == 3) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        - wrk%DU(i,1)*wrk%usol(i,1) + wrk%ADU(i,1)*wrk%usol(i,1) &
                        - wrk%lower_vdep_copy(i)*wrk%usol(i,1)/var%dz(1)
      elseif (var%lowerboundcond(i) == 1) then
        rhs(i) = 0.d0
      elseif (var%lowerboundcond(i) == 2) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        - wrk%DU(i,1)*wrk%usol(i,1) + wrk%ADU(i,1)*wrk%usol(i,1) &
                        + var%lower_flux(i)/(wrk%density(1)*var%dz(1))
      ! Moses (2001) boundary condition for gas giants
      ! A deposition velocity controled by how quickly gases
      ! turbulantly mix vertically
      elseif (var%lowerboundcond(i) == -1) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        - wrk%DU(i,1)*wrk%usol(i,1) + wrk%ADU(i,1)*wrk%usol(i,1) &
                        - (var%edd(1)/wrk%surface_scale_height)*wrk%usol(i,1)/var%dz(1)
      endif
    enddo

    ! Upper boundary
    do i = 1,dat%nq
      k = i + (var%nz-1)*dat%nq
      if (var%upperboundcond(i) == 0) then
        rhs(k) = rhs(k) - wrk%DL(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz) &
                        + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1) &
                        - wrk%upper_veff_copy(i)*wrk%usol(i,var%nz)/var%dz(var%nz)    
      elseif (var%upperboundcond(i) == 2) then
        rhs(k) = rhs(k) - wrk%DL(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz) &
                        + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1) &
                        - var%upper_flux(i)/(wrk%density(var%nz)*var%dz(var%nz))
      endif
    enddo
    
    ! Distributed (volcanic) sources
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == 3) then
        disth = var%lower_dist_height(i)*1.d5        
        jdisth = minloc(var%Z,1, var%Z >= disth) - 1
        jdisth = max(jdisth,2)
        ztop = var%z(jdisth)-var%z(1)
        ztop1 = var%z(jdisth) + 0.5d0*var%dz(jdisth)
        do j = 2,jdisth
          k = i + (j-1)*dat%nq
          rhs(k) = rhs(k) + 2.d0*var%lower_flux(i)*(ztop1-var%z(j))/(wrk%density(j)*ztop**2.d0)
        enddo
      endif
    enddo 
    
    if (dat%fix_water_in_trop) then
      do j = 1,var%trop_ind
        k = dat%LH2O + (j-1)*dat%nq
        rhs(k) = 0.d0
      enddo
    endif
    
  end subroutine
  
  module subroutine jac_background_gas(self, lda_neqs, neqs, usol_flat, jac, err)
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real
    
    class(Atmosphere), target, intent(inout) :: self
    integer, intent(in) :: lda_neqs, neqs
    real(real_kind), target, intent(in) :: usol_flat(neqs)
    real(real_kind), intent(out), target :: jac(lda_neqs)
    character(len=err_len), intent(out) :: err
    
    real(real_kind), pointer :: usol_in(:,:)
    real(real_kind), pointer :: djac(:,:)
    real(real_kind) :: usol_perturb(self%dat%nq,self%var%nz)
    real(real_kind) :: R(self%var%nz)
    real(real_kind) :: rhs(self%var%neqs)
    real(real_kind) :: rhs_perturb(self%var%neqs)
    real(real_kind) :: densities(self%dat%nsp+1,self%var%nz), xl(self%var%nz), xp(self%var%nz)
    ! we need these work arrays for parallel jacobian claculation.
    ! It is probably possible to use memory in "wrk", but i will ignore
    ! this for now.
  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    integer :: i, k, j, m, mm
  
    err = ''
    
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
  
    call prep_all_background_gas(self, usol_in, err)
    if (len_trim(err) /= 0) return
  
    ! compute chemistry contribution to jacobian using forward differences
    jac = 0.d0
    call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, var%trop_ind, dat%nrT, &
                wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%rainout_rates, &
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
                  wrk%H2O_sat_mix, wrk%rainout_rates, &
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
  
    ! diffusion (interior grid points)
    do j = 2,var%nz-1
      do i = 1,dat%nq
        k = i + (j-1)*dat%nq      
        djac(dat%ku,k+dat%nq) = wrk%DU(i,j) + wrk%ADU(i,j)
        djac(dat%kd,k) = djac(dat%kd,k) + wrk%DD(i,j) + (wrk%ADU(i,j) + wrk%ADL(i,j))     
        djac(dat%kl,k-dat%nq) = wrk%DL(i,j) + wrk%ADL(i,j)
      enddo
    enddo
  
    ! Lower boundary
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == 0 .or. var%lowerboundcond(i) == 3) then

        djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(dat%kd,i) = djac(dat%kd,i) - wrk%DU(i,1) + wrk%ADU(i,1) - wrk%lower_vdep_copy(i)/var%dz(1)
      elseif (var%lowerboundcond(i) == 1) then

        do m=1,dat%nq
          mm = dat%kd + i - m
          djac(mm,m) = 0.d0
        enddo
        djac(dat%ku,i+dat%nq) = 0.d0
        ! For some reason this term makes the integration
        ! much happier. I will keep it. Jacobians don't need to be perfect.
        djac(dat%kd,i) = - wrk%DU(i,1)
  
      elseif (var%lowerboundcond(i) == 2) then
        djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(dat%kd,i) = djac(dat%kd,i) - wrk%DU(i,1) + wrk%ADU(i,1)
      elseif (var%lowerboundcond(i) == -1) then
        djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(dat%kd,i) = djac(dat%kd,i) - wrk%DU(i,1) + wrk%ADU(i,1) - &
                         (var%edd(1)/wrk%surface_scale_height)/var%dz(1)
      endif
    enddo
  
    ! Upper boundary
    do i = 1,dat%nq
      k = i + (var%nz-1)*dat%nq
      if (var%upperboundcond(i) == 0) then
        ! rhs(k) = rhs(k) - DL(i,nz)*usol(i,nz) &
        !                 + DL(i,nz)*usol(i,nz-1) + ADL(i,nz)*usol(i,nz-1) &
        !                 - upper_veff(i)*usol(i,nz)/dz(nz)    
  
        djac(dat%kd,k) = djac(dat%kd,k) - wrk%DL(i,var%nz) + wrk%ADL(i,var%nz) &
                        - wrk%upper_veff_copy(i)/var%dz(var%nz) 
        djac(dat%kl,k-dat%nq) = wrk%DL(i,var%nz) + wrk%ADL(i,var%nz)
      elseif (var%upperboundcond(i) == 2) then
        djac(dat%kd,k) = djac(dat%kd,k) - wrk%DL(i,var%nz) + wrk%ADL(i,var%nz)
        djac(dat%kl,k-dat%nq) = wrk%DL(i,var%nz) + wrk%ADL(i,var%nz)
      endif
    enddo
  
    if (dat%fix_water_in_trop) then
      do j = 1,var%trop_ind
        k = dat%LH2O + (j-1)*dat%nq
        do m = 1,dat%nq
          mm = m - dat%LH2O + dat%kd
          djac(mm,k) = 0.d0
        enddo
        ! For some reason this term makes the integration
        ! much happier. I will keep it. Jacobians don't need to be perfect.
        djac(dat%kd,k) = - wrk%DU(dat%LH2O,j)
        djac(dat%ku,k+dat%nq) = 0.d0
        if (j /= 1) then
          djac(dat%kl,k-dat%nq) = 0.d0
        endif
      enddo
    endif
  
  end subroutine
  
  subroutine chempl_t(self, nz, nsp, nrT, np, nl, densities, rx_rates, k, xpT, xlT)
  
    ! input
    class(Atmosphere), intent(in) :: self
    integer, intent(in) :: nz, nsp, nrT
    integer, intent(in) :: np, nl
    real(real_kind), intent(in) :: densities(nsp+1, nz) ! molecules/cm3 of each species
    real(real_kind), intent(in) :: rx_rates(nz,nrT) ! reaction rates (various units)
    integer, intent(in) :: k ! species number
  
    ! output
    real(real_kind), intent(out) :: xpT(nz,np), xlT(nz,nl) ! molecules/cm3/s. if loss_start_ind = 2, then xl is in units of 1/s
  
    ! local
    real(real_kind) :: DD
    integer :: i, ii, iii, m, l, j
    xpT = 0.d0
    xlT = 0.d0
  
    ! np = self%dat%nump(k) ! k is a species
    ! np is number of reactions that produce species k
    do i=1,np
      m = self%dat%iprod(i,k) ! m is reaction number
      l = self%dat%nreactants(m) ! l is the number of reactants
      do j = 1,self%var%nz
        DD = 1.d0
        do ii = 1,l
          iii = self%dat%reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xpT(j,i) = rx_rates(j,m) * DD
      enddo
    enddo
  
    ! nl = self%dat%numl(k) ! k is a species
    ! nl is number of reactions that destroy species k
    do i=1,nl
      m = self%dat%iloss(i,k) ! This will JUST be reaction number
      l = self%dat%nreactants(m) ! number of reactants
      do j = 1,self%var%nz
        DD = 1.d0
        do ii = 1,l
          iii = self%dat%reactants_sp_inds(ii,m)
          DD = DD * densities(iii,j)
        enddo
        xlT(j,i) = rx_rates(j,m) * DD
      enddo
    enddo
  
  end subroutine
  
  module subroutine production_and_loss(self, species, usol, pl, err)     
    use sorting, only: argsort            
    use photochem_types, only: ProductionLoss
    use photochem_const, only: small_real
  
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: species
    real(real_kind), intent(in) :: usol(:,:)
    type(ProductionLoss), intent(out) :: pl
    character(len=err_len), intent(out) :: err
  
    real(real_kind) :: xl(self%var%nz), xp(self%var%nz)
    integer, allocatable :: prod_inds(:), loss_inds(:)
    integer :: ind(1), sp_ind
    integer :: i, j, k, np, nl
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
  
    dat => self%dat
    var => self%var
    wrk => self%wrk
  
    err = ""
    
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = "Input usol to production_and_loss has the wrong dimensions"
      return
    endif
  
    ind = findloc(dat%species_names(1:dat%nsp),species)
    sp_ind = ind(1)
    if (sp_ind == 0) then
      err = "Species "//trim(species)//" is not in the list of species."
      return
    endif
    
    call self%prep_atmosphere(usol, err)
    if (len_trim(err) /= 0) return
  
    np = self%dat%nump(sp_ind)
    nl = self%dat%numl(sp_ind)
    allocate(pl%production(var%nz,np))
    allocate(pl%loss(var%nz,nl))
    allocate(pl%integrated_production(np), pl%integrated_loss(nl))
    allocate(pl%loss_rx(nl),pl%production_rx(np))
    allocate(prod_inds(np), loss_inds(nl))
  
    do j = 1,var%nz
      do k = 1,dat%np
        wrk%densities(k,j) = max(wrk%usol(k,j)* &
                            (wrk%density(j)/wrk%molecules_per_particle(k,j)),small_real)
      enddo
      do k = dat%ng_1,dat%nq
        wrk%densities(k,j) = wrk%usol(k,j)*wrk%density(j)
      enddo
      wrk%densities(dat%nsp,j) = (1.d0-sum(wrk%usol(dat%ng_1:,j)))*wrk%density(j) ! background gas
      wrk%densities(dat%nsp+1,j) = 1.d0 ! for hv
    enddo
  
    if (sp_ind <= dat%nq) then ! long lived or particle
      do k = dat%nq+1,dat%nq+dat%nsl
        call chempl_sl(self, var%nz, dat%nsp, dat%nrT, wrk%densities, wrk%rx_rates, &
                       k, xp, xl) 
        wrk%densities(k,:) = xp/xl
      enddo
    endif
    
    call chempl_t(self, var%nz, dat%nsp, dat%nrT, np, nl, &
                  wrk%densities, wrk%rx_rates, sp_ind, pl%production, pl%loss)
  
    do i = 1,np
      pl%integrated_production(i) = sum(pl%production(:,i)*var%dz)
      k = dat%iprod(i,sp_ind) ! reaction number
      pl%production_rx(i) = dat%reaction_equations(k)
    enddo
    do i = 1,nl
      pl%integrated_loss(i) = sum(pl%loss(:,i)*var%dz)
      k = dat%iloss(i,sp_ind) ! reaction number
      pl%loss_rx(i) = dat%reaction_equations(k)
    enddo
    
    ! sort 
    prod_inds = argsort(pl%integrated_production)
    loss_inds = argsort(pl%integrated_loss)
    prod_inds = prod_inds(np:1:-1)
    loss_inds = loss_inds(nl:1:-1)
  
    pl%integrated_production = pl%integrated_production(prod_inds)
    pl%integrated_loss = pl%integrated_loss(loss_inds)
    
    pl%production = pl%production(:,prod_inds)
    pl%loss = pl%loss(:,loss_inds)
    
    pl%production_rx = pl%production_rx(prod_inds)
    pl%loss_rx = pl%loss_rx(loss_inds)
    
  end subroutine
  
end submodule