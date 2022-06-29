submodule(photochem_atmosphere) photochem_atmosphere_rhs
  implicit none
  
  ! Contains routines to compute the right-hand-side and jacobian
  ! of the system of ODEs describing photochemistry. There are two main components
  
  ! prep_all_background_gas - computes the reaction rates, photolysis rates, 
  ! diffusion coefficients, etc.
  
  ! dochem - computes the chemistry contribution to the right-hand-side
  
contains
  
  subroutine dochem(self, neqs, nsp, np, nsl, nq, nz, trop_ind, nrT, usol, density, rx_rates, &
                    gas_sat_den, molecules_per_particle, &
                    H2O_sat_mix, rainout_rates, &
                    densities, xp, xl, rhs)                 
    use photochem_enum, only: CondensingParticle
    use photochem_common, only: chempl, chempl_sl
    use photochem_eqns, only: damp_condensation_rate
    use photochem_const, only: N_avo, pi, small_real
    
    class(Atmosphere), target, intent(in) :: self
    integer, intent(in) :: neqs, nsp, np, nsl, nq, nz, trop_ind, nrT
    real(dp), intent(in) :: usol(nq,nz), density(nz)
    real(dp), intent(in) :: rx_rates(nz,nrT)
    real(dp), intent(in) :: gas_sat_den(np,nz)
    real(dp), intent(in) :: molecules_per_particle(np,nz)
    real(dp), intent(in) :: H2O_sat_mix(nz)
    real(dp), intent(in) :: rainout_rates(nq, trop_ind)
    real(dp), intent(inout) :: densities(nsp+1,nz), xp(nz), xl(nz)
    real(dp), intent(out) :: rhs(neqs)
    
    real(dp) :: dn_gas_dt, dn_particle_dt, H2O_cold_trap, cond_rate
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
      densities(nsp,j) = (1.0_dp-sum(usol(dat%ng_1:,j)))*density(j) ! background gas
      densities(nsp+1,j) = 1.0_dp ! for hv
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
        rhs(k) = xp(j)/density(j) - xl(j)/density(j)
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
    
    ! if stratospheric water condesation is on, then condense
    ! water in the stratosphere.
    if (dat%water_cond) then
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
        call chempl(self%dat, self%var, densities, rx_rates, i, xp, xl)
        do j = 1,var%nz
          k = i + (j - 1) * nq
          rhs(k) = rhs(k) + (xp(j) - xl(j))/density(j)
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
              rhs(kk) = rhs(kk) + dn_gas_dt/density(j)            
            
              ! rate of particle production from gas condensing
              ! in particles/cm3/s
              dn_particle_dt = - dn_gas_dt/density(j)
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
  
  subroutine prep_atm_background_gas(self, usol, sum_usol, &
                                     density, mubar, pressure, fH2O, H2O_sat_mix, err)
    use photochem_eqns, only: sat_pressure_H2O, molar_weight, press_and_den
    use, intrinsic :: iso_c_binding, only : c_loc, c_ptr
    use cminpack_wrapper, only: hybrd1 ! interface to hybrd1 from cminpack
    
    type(Atmosphere), target, intent(inout) :: self
    real(dp), intent(inout) :: usol(:,:)
    real(dp), intent(out) :: sum_usol(:)
    real(dp), intent(out) :: density(:)
    real(dp), intent(out) :: mubar(:), pressure(:), fH2O(:), H2O_sat_mix(:)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: rel
    integer :: i
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(Atmosphere), pointer :: self_ptr
    
    ! hybrd1 varibles for cminpack
    integer :: info
    real(dp), parameter :: tol = 1.d-8
    integer :: lwa
    real(dp), allocatable :: fvec(:)
    real(dp), allocatable :: wa(:)
    type(c_ptr) :: ptr ! void c pointer for cminpack
    ! end hybrd1 varibles
    
    dat => self%dat
    var => self%var
    
    ! If water is fixed in the troposhere, lets set it to zero.
    ! This will be to get an initial guess for the non-linear solve.
    if (dat%fix_water_in_trop) then
      do i = 1,var%trop_ind
        usol(dat%LH2O,i) = 0.0_dp
      enddo
    endif 
    
    do i = 1,var%nz
      sum_usol(i) = sum(usol(dat%ng_1:,i))
      if (sum_usol(i) > 1.0e0_dp) then
        err = 'Mixing ratios sum to >1.0 at some altitude (should be <=1).' // &
              ' The atmosphere is probably in a run-away state'
        return
      endif
    enddo
    
    do i = 1,var%nz
      call molar_weight(dat%nll, usol(dat%ng_1:,i), sum_usol(i), dat%species_mass(dat%ng_1:), dat%back_gas_mu, mubar(i))
    enddo
    
    call press_and_den(var%nz, var%temperature, var%grav, var%surface_pressure*1.e6_dp, var%dz, &
                       mubar, pressure, density)
                       
    if (dat%fix_water_in_trop) then
      ! Compute initial guess for fH2O
      do i = 1,var%trop_ind
        if (var%use_manabe) then
          ! manabe formula
          rel = 0.77e0_dp*(pressure(i)/pressure(1)-0.02e0_dp)/0.98e0_dp
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
      allocate(fvec(var%trop_ind))
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
      
      call press_and_den(var%nz, var%temperature, var%grav, var%surface_pressure*1.e6_dp, var%dz, &
                         mubar, pressure, density)
      
    endif 
    
    if (dat%water_cond) then
      ! compute H2O saturation mixing ratio at all altitudes
      do i = 1,var%nz
        H2O_sat_mix(i) = sat_pressure_H2O(var%temperature(i))/pressure(i)
      enddo
    endif
    
  end subroutine
  
  module function fcn_fH2O(ptr, n, x, fvec, iflag) result(res) bind(c)
    use photochem_eqns, only: sat_pressure_H2O, molar_weight, press_and_den
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_double, c_int
  
    type(c_ptr) :: ptr ! void pointer. Be careful!
    integer(c_int), value :: n, iflag ! n == trop_ind
    real(c_double), intent(in) :: x(n) ! x == fH2O
    real(c_double), intent(out) :: fvec(n) ! fvec == residual
    integer(c_int) :: res
  
    real(dp) :: fH2O(n)
    real(dp), pointer :: usol(:,:)
    real(dp), pointer :: sum_usol(:) 
    real(dp), pointer :: mubar(:)
    real(dp), pointer :: density(:)
    real(dp), pointer :: pressure(:)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(Atmosphere), pointer :: self
  
    integer :: i
    real(dp) :: rel
  
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
  
    call press_and_den(var%nz, var%temperature, var%grav, var%surface_pressure*1.e6_dp, var%dz, &
                       mubar, pressure, density)
  
    do i = 1,n
      if (var%use_manabe) then
        ! manabe formula ()
        rel = 0.77e0_dp*(pressure(i)/pressure(1)-0.02e0_dp)/0.98e0_dp
      else
        rel = var%relative_humidity 
      endif
      fH2O(i) = rel*sat_pressure_H2O(var%temperature(i))/pressure(i)
    enddo  
  
    fvec = x - fH2O
    res = 0
  end function
  
  module subroutine prep_all_background_gas(self, usol_in, err)
    use photochem_enum, only: CondensingParticle
    use photochem_enum, only: ArrheniusSaturation, H2SO4Saturation
    use photochem_enum, only: MixingRatioBC
    use photochem_common, only: reaction_rates, diffusion_coefficients, rainout, photorates
    use photochem_eqns, only: saturation_density
    use photochem_const, only: pi, k_boltz, N_avo, small_real
  
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    character(:), allocatable, intent(out) :: err
  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    integer :: i, j, k
    real(dp) :: P_H2SO4
  
  
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    ! make a copy of the mixing ratios. You can not alter the input mixing ratios
    ! this will make CVODE unstable. Guard against division by zero
    do j = 1,var%nz
      do i = 1,dat%nq
        if (usol_in(i,j) < 0.0_dp) then
          wrk%usol(i,j) = min(usol_in(i,j),-small_real)
        else
          wrk%usol(i,j) = max(usol_in(i,j), small_real)
        endif
      enddo
    enddo

    wrk%upper_veff_copy = var%upper_veff
    wrk%lower_vdep_copy = var%lower_vdep
  
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == MixingRatioBC) then
        wrk%usol(i,1) = var%lower_fix_mr(i)
      endif
    enddo
  
    call prep_atm_background_gas(self, wrk%usol, wrk%sum_usol, &
                                 wrk%density, wrk%mubar, wrk%pressure, wrk%fH2O, wrk%H2O_sat_mix, err)  
    if (allocated(err)) return  
  
    ! diffusion coefficients
    call diffusion_coefficients(dat, var, wrk%density, wrk%mubar, &
                                wrk%DU, wrk%DD, wrk%DL, wrk%ADU, wrk%ADL, wrk%ADD, &
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
          if (dat%particle_formation_method(i) == CondensingParticle) then
            if (dat%particle_sat_type(i) == ArrheniusSaturation) then ! arrhenius
              wrk%gas_sat_den(i,j) = saturation_density(var%temperature(j), &
                                                   dat%particle_sat_params(1,i), &
                                                   dat%particle_sat_params(2,i), &
                                                   dat%particle_sat_params(3,i))
            elseif (dat%particle_sat_type(i) == H2SO4Saturation) then ! interpolate to H2SO4 data
              call dat%H2SO4_sat%evaluate(var%temperature(j), &
                                          log10(wrk%usol(dat%LH2O,j)*wrk%pressure(j)/1.e6_dp),P_H2SO4)
              P_H2SO4 = 10.0_dp**(P_H2SO4)
              wrk%gas_sat_den(i,j) = (P_H2SO4*1.e6_dp)/(k_boltz*var%temperature(j))
            endif
          endif
          wrk%molecules_per_particle(i,j) = (4.0_dp/3.0_dp)*pi*var%particle_radius(i,j)**3.0_dp* &
                                            dat%particle_density(i)*(1/dat%species_mass(i))*N_avo
        enddo
      enddo
    endif
  
    ! densities include particle densities
    do j = 1,var%nz
      do i = 1,dat%npq
        wrk%densities(i,j) = max(wrk%usol(i,j)*(wrk%density(j)/wrk%molecules_per_particle(i,j)), small_real)
      enddo
      do i = dat%ng_1,dat%nq
        wrk%densities(i,j) = wrk%usol(i,j)*wrk%density(j)
      enddo
      wrk%densities(dat%nsp,j) = (1.0_dp-wrk%sum_usol(j))*wrk%density(j) ! background gas
      wrk%densities(dat%nsp+1,j) = 1.0_dp ! for hv
    enddo
    
    call reaction_rates(self%dat, self%var, wrk%density, wrk%densities, wrk%rx_rates)

    call photorates(dat, var, wrk%densities, &
                    wrk%prates, wrk%surf_radiance, wrk%amean_grd, wrk%optical_depth, err)
    if (allocated(err)) return
  
    do i = 1,dat%kj
      k = dat%photonums(i)
      wrk%rx_rates(:,k) = wrk%prates(:,i) 
    enddo 
  
    ! rainout rates
    if (dat%gas_rainout) then
      call rainout(self%dat, self%var, &
                   wrk%usol, wrk%density, wrk%rainout_rates)
    endif
  
  end subroutine
  
  module subroutine rhs_background_gas(self, neqs, usol_flat, rhs, err)
    use photochem_enum, only: MosesBC, VelocityBC, MixingRatioBC, FluxBC, VelocityDistributedFluxBC
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real  
    
    class(Atmosphere), target, intent(inout) :: self
    integer, intent(in) :: neqs
    real(dp), target, intent(in) :: usol_flat(neqs)
    real(dp), intent(out) :: rhs(neqs)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: disth, ztop, ztop1    
    integer :: i, k, j, jdisth
    
    real(dp), pointer :: usol_in(:,:)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    
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
    if (allocated(err)) return
    
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
      elseif (var%lowerboundcond(i) == MixingRatioBC) then
        rhs(i) = 0.0_dp
      elseif (var%lowerboundcond(i) == FluxBC) then
        rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                        + wrk%DD(i,1)*wrk%usol(i,1) + wrk%ADD(i,1)*wrk%usol(i,1) &
                        + var%lower_flux(i)/(wrk%density(1)*var%dz(1))
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
                        - var%upper_flux(i)/(wrk%density(var%nz)*var%dz(var%nz))
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
          rhs(k) = rhs(k) + 2.0_dp*var%lower_flux(i)*(ztop1-var%z(j))/(wrk%density(j)*ztop**2.0_dp)
        enddo
      endif
    enddo 
    
    if (dat%fix_water_in_trop) then
      do j = 1,var%trop_ind
        k = dat%LH2O + (j-1)*dat%nq
        rhs(k) = 0.0_dp
      enddo
    endif
    
  end subroutine
  
  module subroutine jac_background_gas(self, lda_neqs, neqs, usol_flat, jac, err)
    use photochem_enum, only: MosesBC, VelocityBC, MixingRatioBC, FluxBC, VelocityDistributedFluxBC
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real
    
    class(Atmosphere), target, intent(inout) :: self
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
    real(dp) :: densities(self%dat%nsp+1,self%var%nz), xl(self%var%nz), xp(self%var%nz)
    ! we need these work arrays for parallel jacobian claculation.
    ! It is probably possible to use memory in "wrk", but i will ignore
    ! this for now.
  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

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
  
    call prep_all_background_gas(self, usol_in, err)
    if (allocated(err)) return
  
    ! compute chemistry contribution to jacobian using forward differences
    jac = 0.0_dp
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
        djac(dat%kd,k) = djac(dat%kd,k) + wrk%DD(i,j) + wrk%ADD(i,j)     
        djac(dat%kl,k-dat%nq) = wrk%DL(i,j) + wrk%ADL(i,j)
      enddo
    enddo
  
    ! Lower boundary
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityBC .or. var%lowerboundcond(i) == VelocityDistributedFluxBC) then

        djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
        djac(dat%kd,i) = djac(dat%kd,i) + wrk%DD(i,1) + wrk%ADD(i,1) - wrk%lower_vdep_copy(i)/var%dz(1)
      elseif (var%lowerboundcond(i) == MixingRatioBC) then

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
        ! rhs(k) = rhs(k) - DL(i,nz)*usol(i,nz) &
        !                 + DL(i,nz)*usol(i,nz-1) + ADL(i,nz)*usol(i,nz-1) &
        !                 - upper_veff(i)*usol(i,nz)/dz(nz)    
  
        djac(dat%kd,k) = djac(dat%kd,k) + wrk%DD(i,var%nz) + wrk%ADD(i,var%nz) &
                        - wrk%upper_veff_copy(i)/var%dz(var%nz) 
        djac(dat%kl,k-dat%nq) = wrk%DL(i,var%nz) + wrk%ADL(i,var%nz)
      elseif (var%upperboundcond(i) == FluxBC) then
        djac(dat%kd,k) = djac(dat%kd,k) + wrk%DD(i,var%nz) + wrk%ADD(i,var%nz)
        djac(dat%kl,k-dat%nq) = wrk%DL(i,var%nz) + wrk%ADL(i,var%nz)
      endif
    enddo
  
    if (dat%fix_water_in_trop) then
      do j = 1,var%trop_ind
        k = dat%LH2O + (j-1)*dat%nq
        do m = 1,dat%nq
          mm = m - dat%LH2O + dat%kd
          djac(mm,k) = 0.0_dp
        enddo
        ! For some reason this term makes the integration
        ! much happier. I will keep it. Jacobians don't need to be perfect.
        djac(dat%kd,k) = - wrk%DU(dat%LH2O,j)
        djac(dat%ku,k+dat%nq) = 0.0_dp
        if (j /= 1) then
          djac(dat%kl,k-dat%nq) = 0.0_dp
        endif
      enddo
    endif
  
  end subroutine
  
  module subroutine right_hand_side_chem(self, usol, rhs, err)
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(out) :: rhs(:)
    character(:), allocatable, intent(out) :: err
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz .or. size(rhs) /= var%neqs) then
      err = "Input usol or rhs to dochem_implicit has the wrong dimensions"
      return
    endif

    call self%prep_atmosphere(usol, err)
    if (allocated(err)) return
    
    call dochem(self, var%neqs, dat%nsp, dat%np, dat%nsl, dat%nq, var%nz, &
                var%trop_ind, dat%nrT, wrk%usol, wrk%density, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%H2O_sat_mix, wrk%rainout_rates, &
                wrk%densities, wrk%xp, wrk%xl, rhs)
                              
  end subroutine
  
  module subroutine production_and_loss(self, species, usol, pl, err)     
    use futils, only: argsort            
    use photochem_common, only: chempl_sl, chempl_t
    use photochem_types, only: ProductionLoss
    use photochem_const, only: small_real
  
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: species
    real(dp), intent(in) :: usol(:,:)
    type(ProductionLoss), intent(out) :: pl
    character(:), allocatable, intent(out) :: err
  
    real(dp) :: xl(self%var%nz), xp(self%var%nz)
    integer, allocatable :: prod_inds(:), loss_inds(:)
    integer :: ind(1), sp_ind
    integer :: i, j, k, np, nl, nlT
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
  
    dat => self%dat
    var => self%var
    wrk => self%wrk
  
    
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
      do k = 1,dat%np
        wrk%densities(k,j) = max(wrk%usol(k,j)* &
                            (wrk%density(j)/wrk%molecules_per_particle(k,j)),small_real)
      enddo
      do k = dat%ng_1,dat%nq
        wrk%densities(k,j) = wrk%usol(k,j)*wrk%density(j)
      enddo
      wrk%densities(dat%nsp,j) = (1.0_dp-sum(wrk%usol(dat%ng_1:,j)))*wrk%density(j) ! background gas
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
          wrk%rainout_rates(sp_ind,1:var%trop_ind)*wrk%densities(sp_ind,1:var%trop_ind)
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