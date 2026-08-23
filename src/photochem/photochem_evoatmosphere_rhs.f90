#:set TYPES = ['real(dp)', 'type(dual)']
#:set NAMES = ['real', 'dual']
#:set TYPES_NAMES = list(zip(TYPES, NAMES))

submodule(photochem_evoatmosphere) photochem_evoatmosphere_rhs
  use differentia, only: dual, assignment(=), operator(+), operator(-), &
                         operator(*), operator(/), operator(**), &
                         operator(>), operator(<=), max, atan
  implicit none
  
  interface dochem
    module procedure :: dochem_real, dochem_dual
  end interface

  interface prepare_chemistry_state
    module procedure :: prepare_chemistry_state_real, prepare_chemistry_state_dual
  end interface

  interface set_gas_reaction_tendencies
    module procedure :: set_gas_reaction_tendencies_real, set_gas_reaction_tendencies_dual
  end interface

  interface add_rainout_tendencies
    module procedure :: add_rainout_tendencies_real, add_rainout_tendencies_dual
  end interface

  interface add_particle_reaction_tendencies
    module procedure :: add_particle_reaction_tendencies_real, add_particle_reaction_tendencies_dual
  end interface

  interface add_condensation_tendencies
    module procedure :: add_condensation_tendencies_real, add_condensation_tendencies_dual
  end interface

  interface chempl
    module procedure :: chempl_real, chempl_dual
  end interface

  interface chempl_sl
    module procedure :: chempl_sl_real, chempl_sl_dual
  end interface

  interface damp_condensation_rate
    module procedure :: damp_condensation_rate_real, damp_condensation_rate_dual
  end interface

contains

  #:for TYPE1, NAME in TYPES_NAMES
  subroutine dochem_${NAME}$(self, usol, rx_rates, &
                    gas_sat_den, molecules_per_particle, &
                    rainout_rates, scale_height, wfall, &
                    density, mix, densities, xp, xl, rhs)                 
    class(EvoAtmosphere), target, intent(in) :: self
    ${TYPE1}$, intent(in) :: usol(:,:)
    real(dp), intent(in) :: rx_rates(:,:)
    real(dp), intent(in) :: gas_sat_den(:,:)
    real(dp), intent(in) :: molecules_per_particle(:,:)
    real(dp), intent(in) :: rainout_rates(:,:), scale_height(:), wfall(:,:)
    real(dp), intent(in) :: density(:)
    ${TYPE1}$, intent(inout) :: mix(:,:)
    ${TYPE1}$, intent(inout) :: densities(:,:), xp(:), xl(:)
    ${TYPE1}$, intent(inout) :: rhs(:) ! neqs

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var

    ! Preserve this process order when adding the individual tendencies.
    ! Maintenance rule: every new composition-dependent tendency added here
    ! must have a matching derivative in analytical_chemistry_jacobian and
    ! analytical-versus-autodiff coverage in test_jacobian. An unsupported
    ! analytical term must be rejected rather than silently omitted; callers
    ! can explicitly select the autodiff method until support is implemented.
    call prepare_chemistry_state(dat, var, usol, rx_rates, &
                                 molecules_per_particle, density, mix, &
                                 densities, xp, xl)
    rhs = 0.0_dp
    call set_gas_reaction_tendencies(dat, var, densities, rx_rates, &
                                     xp, xl, rhs)
    call add_rainout_tendencies(dat, var, usol, rainout_rates, rhs)
    if (dat%there_are_particles) then
      call add_particle_reaction_tendencies(dat, var, densities, rx_rates, &
                                            xp, xl, rhs)
      call add_condensation_tendencies(dat, var, usol, gas_sat_den, &
                                       scale_height, wfall, rhs)
    endif

  end subroutine

  ! Prepare number densities, including diagnostic short-lived species.
  subroutine prepare_chemistry_state_${NAME}$(dat, var, usol, rx_rates, &
                                              molecules_per_particle, density, &
                                              mix, densities, xp, xl)
    use photochem_const, only: small_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: usol(:,:)
    real(dp), intent(in) :: rx_rates(:,:), molecules_per_particle(:,:)
    real(dp), intent(in) :: density(:)
    ${TYPE1}$, intent(inout) :: mix(:,:), densities(:,:), xp(:), xl(:)

    integer :: i, j, k

    do j = 1,var%nz
      mix(:,j) = usol(:,j)/density(j)
    enddo

    do j = 1,var%nz
      do i = 1,dat%npq
        densities(i,j) = max(usol(i,j)* &
                             (1.0_dp/molecules_per_particle(i,j)), &
                             small_real)
      enddo
      do i = dat%ng_1,dat%nq
        densities(i,j) = usol(i,j)
      enddo
      densities(dat%nsp+1,j) = 1.0_dp ! for hv
    enddo

    do k = dat%nq+1,dat%nq+dat%nsl
      call chempl_sl(dat, var, densities, rx_rates, k, xp, xl)
      densities(k,:) = xp/xl
    enddo
  end subroutine

  ! Set the reaction tendencies for long-lived gas species.
  subroutine set_gas_reaction_tendencies_${NAME}$(dat, var, densities, &
                                                  rx_rates, xp, xl, rhs)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: densities(:,:)
    real(dp), intent(in) :: rx_rates(:,:)
    ${TYPE1}$, intent(inout) :: xp(:), xl(:), rhs(:)

    integer :: i, j, k

    do i = dat%ng_1,dat%nq
      call chempl(dat, var, densities, rx_rates, i, xp, xl)
      do j = 1,var%nz
        k = i + (j - 1) * dat%nq
        rhs(k) = xp(j) - xl(j)
      enddo
    enddo
  end subroutine

  ! Add tropospheric gas rainout losses.
  subroutine add_rainout_tendencies_${NAME}$(dat, var, usol, &
                                             rainout_rates, rhs)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: usol(:,:)
    real(dp), intent(in) :: rainout_rates(:,:)
    ${TYPE1}$, intent(inout) :: rhs(:)

    integer :: i, j, k

    if (dat%gas_rainout) then
      do j = 1,var%trop_ind
        do i = 1,dat%nq
          k = i + (j - 1) * dat%nq
          rhs(k) = rhs(k) - rainout_rates(i,j)*usol(i,j)
        enddo
      enddo
    endif
  end subroutine

  ! Add particle production and loss from chemical reactions.
  subroutine add_particle_reaction_tendencies_${NAME}$(dat, var, &
                                                       densities, rx_rates, &
                                                       xp, xl, rhs)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: densities(:,:)
    real(dp), intent(in) :: rx_rates(:,:)
    ${TYPE1}$, intent(inout) :: xp(:), xl(:), rhs(:)

    integer :: i, j, k

    do i = 1,dat%np
      call chempl(dat, var, densities, rx_rates, i, xp, xl)
      do j = 1,var%nz
        k = i + (j - 1) * dat%nq
        rhs(k) = rhs(k) + (xp(j) - xl(j))
      enddo
    enddo
  end subroutine

  ! Add gas-particle transfer from condensation and evaporation.
  subroutine add_condensation_tendencies_${NAME}$(dat, var, usol, &
                                                  gas_sat_den, scale_height, &
                                                  wfall, rhs)
    use photochem_enum, only: CondensingParticle
    use photochem_const, only: small_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: usol(:,:)
    real(dp), intent(in) :: gas_sat_den(:,:), scale_height(:), wfall(:,:)
    ${TYPE1}$, intent(inout) :: rhs(:)

    real(dp) :: cond_rate0
    ${TYPE1}$ :: rh, dn_gas_dt, dn_particle_dt, cond_rate
    integer :: i, ii, j, k, kk

    #:if NAME == 'dual'
    rh = dual(size(usol(1,1)%der))
    dn_gas_dt = dual(size(usol(1,1)%der))
    dn_particle_dt = dual(size(usol(1,1)%der))
    cond_rate = dual(size(usol(1,1)%der))
    #:endif

    do j = 1,var%nz
      do i = 1,dat%np
        if (dat%particle_formation_method(i) == CondensingParticle) then
          ii = dat%particle_gas_phase_ind(i)
          kk = ii + (j - 1) * dat%nq
          k = i + (j - 1) * dat%nq

          rh = max(usol(ii,j)/gas_sat_den(i,j),small_real)
          if (rh > var%cond_params(i)%RHc) then
            cond_rate0 = var%cond_params(i)%k_cond* &
                         (var%edd(j)/scale_height(j)**2.0_dp)
            cond_rate = damp_condensation_rate( &
                cond_rate0, var%cond_params(i)%RHc, &
                (1.0_dp + var%cond_params(i)%smooth_factor)* &
                    var%cond_params(i)%RHc, rh)

            dn_gas_dt = -cond_rate*usol(ii,j)
            rhs(kk) = rhs(kk) + dn_gas_dt
            dn_particle_dt = -dn_gas_dt
            rhs(k) = rhs(k) + dn_particle_dt

          elseif (rh <= var%cond_params(i)%RHc .and. var%evaporation) then
            cond_rate0 = var%cond_params(i)%k_evap* &
                         (wfall(i,j)/scale_height(j))
            cond_rate = damp_condensation_rate( &
                cond_rate0, 1.0_dp/var%cond_params(i)%RHc, &
                (1.0_dp + var%cond_params(i)%smooth_factor)/ &
                    var%cond_params(i)%RHc, 1.0_dp/rh)

            dn_gas_dt = cond_rate*usol(i,j)
            rhs(kk) = rhs(kk) + dn_gas_dt
            dn_particle_dt = -dn_gas_dt
            rhs(k) = rhs(k) + dn_particle_dt
          endif
        endif
      enddo
    enddo
  end subroutine

  ! Sum chemical production and loss rates for one species.
  subroutine chempl_${NAME}$(dat, var, densities, rx_rates, k, xp, xl)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: densities(:,:)
    real(dp), intent(in) :: rx_rates(:,:)
    integer, intent(in) :: k
    ${TYPE1}$, intent(inout) :: xp(:), xl(:)

    ${TYPE1}$ :: DD
    integer :: i, ii, iii, m, l, j

    #:if NAME == 'dual'
    DD = dual(size(densities(1,1)%der))
    #:endif
    xp = 0.0_dp
    xl = 0.0_dp

    do i = 1,dat%pl(k)%nump
      m = dat%pl(k)%iprod(i)
      l = dat%rx(m)%nreact
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD*densities(iii,j)
        enddo
        xp(j) = xp(j) + rx_rates(j,m)*DD
      enddo
    enddo

    do i = 1,dat%pl(k)%numl
      m = dat%pl(k)%iloss(i)
      l = dat%rx(m)%nreact
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD*densities(iii,j)
        enddo
        xl(j) = xl(j) + rx_rates(j,m)*DD
      enddo
    enddo
  end subroutine

  ! Compute production and pseudo-first-order loss for a short-lived species.
  subroutine chempl_sl_${NAME}$(dat, var, densities, rx_rates, k, xp, xl)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    ${TYPE1}$, intent(in) :: densities(:,:)
    real(dp), intent(in) :: rx_rates(:,:)
    integer, intent(in) :: k
    ${TYPE1}$, intent(inout) :: xp(:), xl(:)

    ${TYPE1}$ :: DD
    integer :: i, ii, iii, m, l, j

    #:if NAME == 'dual'
    DD = dual(size(densities(1,1)%der))
    #:endif
    xp = 0.0_dp
    xl = 0.0_dp

    do i = 1,dat%pl(k)%nump
      m = dat%pl(k)%iprod(i)
      l = dat%rx(m)%nreact
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD*densities(iii,j)
        enddo
        xp(j) = xp(j) + rx_rates(j,m)*DD
      enddo
    enddo

    do i = 1,dat%pl(k)%numl
      m = dat%pl(k)%iloss(i)
      l = dat%rx(m)%nreact
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          if (iii /= k) then
            DD = DD*densities(iii,j)
          endif
        enddo
        xl(j) = xl(j) + rx_rates(j,m)*DD
      enddo
    enddo
  end subroutine

  ! Smooth a condensation or evaporation coefficient across saturation.
  pure function damp_condensation_rate_${NAME}$(A, rhc, rh0, rh) result(k)
    use photochem_const, only: pi
    real(dp), intent(in) :: A, rhc, rh0
    ${TYPE1}$, intent(in) :: rh
    ${TYPE1}$ :: k

    k = A*(2.0_dp/pi)*atan((rh - rhc)/(rh0 - rhc))
  end function

  #:endfor

  ! Compute reaction-resolved chemical production and loss for one species.
  pure subroutine chempl_t(dat, var, densities, rx_rates, k, xpT, xlT)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: densities(:,:)
    real(dp), intent(in) :: rx_rates(:,:)
    integer, intent(in) :: k
    real(dp), intent(out) :: xpT(:,:), xlT(:,:)

    real(dp) :: DD
    integer :: i, ii, iii, m, l, j

    xpT = 0.0_dp
    xlT = 0.0_dp

    do i = 1,dat%pl(k)%nump
      m = dat%pl(k)%iprod(i)
      l = dat%rx(m)%nreact
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD*densities(iii,j)
        enddo
        xpT(j,i) = rx_rates(j,m)*DD
      enddo
    enddo

    do i = 1,dat%pl(k)%numl
      m = dat%pl(k)%iloss(i)
      l = dat%rx(m)%nreact
      do j = 1,var%nz
        DD = 1.0_dp
        do ii = 1,l
          iii = dat%rx(m)%react_sp_inds(ii)
          DD = DD*densities(iii,j)
        enddo
        xlT(j,i) = rx_rates(j,m)*DD
      enddo
    enddo
  end subroutine

  subroutine diffusion_coefficients_evo(dat, var, den, mubar, &
                                    DU, DD, DL, ADU, ADL, ADD, wfall, VH2_esc, VH_esc)
    use photochem_eqns, only: dynamic_viscosity_air, fall_velocity, slip_correction_factor, &
                              default_binary_diffusion_param
    use photochem_const, only: k_boltz, N_avo
    use photochem_enum, only: DiffusionLimHydrogenEscape
    use photochem_vars, only: binary_diffusion_fcn

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
    real(dp) :: gamma_i_gas_D(dat%nll,var%nz)
    real(dp) :: grav_av, mubar_av
    real(dp) :: bx1x2

    ! for particles
    real(dp) :: air_density_pp, air_density
    real(dp) :: wfall_pp, wfall_i
    real(dp) :: viscosity_pp, viscosity
    real(dp) :: FF2, FF1

    ! molecular diffusion function
    procedure(binary_diffusion_fcn), pointer :: binary_diffusion_param
    
    integer :: j, i, k

    if (associated(var%custom_binary_diffusion_fcn)) then
      binary_diffusion_param => var%custom_binary_diffusion_fcn
    else
      binary_diffusion_param => default_binary_diffusion_param
    endif

    if (var%upwind_molec_diff) then; block
    real(dp) :: dTdz_tmp, scale_height_i_tmp, b1x2_tmp, grav_tmp, mubar_tmp
    ! If using upwind scheme for molecular diffusion, then compute the needed
    ! advective velocities.
    do j = 1,var%nz-1
      dTdz_tmp = (var%temperature(j+1) - var%temperature(j))/var%dz(j)
      grav_tmp = var%grav(j)
      mubar_tmp = mubar(j)
      do i = dat%ng_1,dat%nq
        k = i - dat%np
        scale_height_i_tmp = (N_avo*k_boltz*var%temperature(j))/(grav_tmp*dat%species_mass(i))
        b1x2_tmp = binary_diffusion_param(dat%species_mass(i), mubar_tmp, var%temperature(j))
        gamma_i_gas_D(k,j) = (b1x2_tmp/den(j))*(1.0_dp/scale_height_i_tmp + (1.0_dp/var%temperature(j))*dTdz_tmp)
      enddo
    enddo
    j = var%nz
    dTdz_tmp = (var%temperature(var%nz) - var%temperature(var%nz-1))/var%dz(var%nz)
    grav_tmp = var%grav(j)
    mubar_tmp = mubar(j)
    do i = dat%ng_1,dat%nq
      k = i - dat%np
      scale_height_i_tmp = (N_avo*k_boltz*var%temperature(j))/(grav_tmp*dat%species_mass(i))
      b1x2_tmp = binary_diffusion_param(dat%species_mass(i), mubar_tmp, var%temperature(j))
      gamma_i_gas_D(k,j) = (b1x2_tmp/den(j))*(1.0_dp/scale_height_i_tmp + (1.0_dp/var%temperature(j))*dTdz_tmp)     
    enddo
    endblock; endif

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

        gamma_i_gas_av(k,j) = eddav(j)*(1.0_dp/scale_height_av(j) + (1.0_dp/tav(j))*dTdz(j))

        ! If we don't use upwind scheme for molecular diffusion, then we add the centered
        ! molecular diffusion terms here.
        if (.not.var%upwind_molec_diff) then
          gamma_i_gas_av(k,j) = gamma_i_gas_av(k,j) + &
            (b1x2av(k,j)/denav(j))*(1.0_dp/scale_height_i_av + (1.0_dp/tav(j))*dTdz(j))
        endif

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

        if (var%upwind_molec_diff) then
          ADU(i,j) = ADU(i,j) + gamma_i_gas_D(k,j+1)/var%dz(j)
          ADL(i,j) = ADL(i,j) + 0.0_dp
          ADD(i,j) = ADD(i,j) + (- gamma_i_gas_D(k,j)/var%dz(j))
        endif
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

      if (var%upwind_molec_diff) then
        ADU(i,j) = ADU(i,j) + gamma_i_gas_D(k,j+1)/var%dz(j)
        ADD(i,j) = ADD(i,j) + 0.0_dp
      endif
    enddo
    ! upper boundary
    j = var%nz
    do i = dat%ng_1,dat%nq
      k = i - dat%np
      DL(i,j) = (eddav(j-1) + (b1x2av(k,j-1)/denav(j-1)))/(var%dz(j)**2.0_dp)
      DD(i,j) = - DL(i,j)

      ADL(i,j) = - gamma_i_gas_av(k,j-1)/(2.0_dp*var%dz(j))
      ADD(i,j) = ADL(i,j)

      if (var%upwind_molec_diff) then
        ADL(i,j) = ADL(i,j) + 0.0_dp
        ADD(i,j) = ADD(i,j) + (- gamma_i_gas_D(k,j)/var%dz(j))
      endif
    enddo
    
    ! ! particles (eddy diffusion)
    ! ! middle
    ! do j = 2,var%nz-1
    !   do i = 1,dat%np
    !     ! diffusion
    !     DU(i,j) = eddav(j)/(var%dz(j)**2.0_dp)
    !     DL(i,j) = eddav(j-1)/(var%dz(j)**2.0_dp)
    !     DD(i,j) = - DU(i,j) - DL(i,j)

    !     ! advection
    !     ADU(i,j) = gamma_i_part_av(i,j)/(2.0_dp*var%dz(j))
    !     ADL(i,j) = - gamma_i_part_av(i,j-1)/(2.0_dp*var%dz(j))
    !     ADD(i,j) = ADU(i,j) + ADL(i,j)
    !   enddo
    ! enddo
    ! ! lower boundary
    ! j = 1
    ! do i = 1,dat%np
    !   DU(i,j) = eddav(j)/(var%dz(j)**2.0_dp)
    !   DD(i,j) = - DU(i,j)

    !   ADU(i,j) = gamma_i_part_av(i,j)/(2.0_dp*var%dz(j))
    !   ADD(i,j) = ADU(i,j)
    ! enddo
    ! ! upper boundary
    ! j = var%nz
    ! do i = 1,dat%np
    !   DL(i,j) = eddav(j-1)/(var%dz(j)**2.0_dp)
    !   DD(i,j) = - DL(i,j)

    !   ADL(i,j) = - gamma_i_part_av(i,j-1)/(2.0_dp*var%dz(j))
    !   ADD(i,j) = ADL(i,j)
    ! enddo

    ! We do not include eddy diffusion for particles
    do j = 1,var%nz
      do i = 1,dat%np
        DU(i,j) = 0.0_dp
        DL(i,j) = 0.0_dp
        DD(i,j) = 0.0_dp
        ADU(i,j) = 0.0_dp
        ADL(i,j) = 0.0_dp
        ADD(i,j) = 0.0_dp
      enddo
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

    ! H2 escape
    if (dat%H_escape_type == DiffusionLimHydrogenEscape) then
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

  module subroutine apply_lower_boundary_conditions(self, temperature, usol_bottom, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_const, only: k_boltz
    use photochem_enum, only: DensityBC, PressureBC
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: temperature
    real(dp), intent(inout) :: usol_bottom(:)
    character(:), allocatable, intent(out) :: err

    real(dp) :: Psat
    integer :: i, particle_ind
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    if (size(usol_bottom) /= dat%nq) then
      err = 'The bottom-layer state has the wrong dimensions'
      return
    endif
    if (.not. ieee_is_finite(temperature) .or. temperature <= 0.0_dp) then
      err = 'The bottom-layer temperature is not finite and positive'
      return
    endif

    do i = 1,dat%nq
      if (var%lowerboundcond(i) == DensityBC) then
        usol_bottom(i) = var%lower_fix_den(i)
      elseif (var%lowerboundcond(i) == PressureBC) then
        Psat = huge(1.0_dp)
        if (dat%gas_particle_ind(i) /= 0) then
          particle_ind = dat%gas_particle_ind(i)
          Psat = dat%particle_sat(particle_ind)%sat_pressure(temperature)* &
                 var%cond_params(particle_ind)%RHc
        endif
        usol_bottom(i) = min(var%lower_fix_press(i), Psat)/(k_boltz*temperature)
      endif
    enddo

  end subroutine

  module subroutine clip_usol(usol_in, usol_out)
    use photochem_const, only: small_real
    real(dp), intent(in) :: usol_in(:,:)
    real(dp), intent(out) :: usol_out(:,:)

    usol_out = usol_in
    where (usol_out >= 0.0_dp)
      usol_out = max(usol_out, small_real)
    elsewhere
      usol_out = min(usol_out, -small_real)
    endwhere

  end subroutine

  module subroutine prepare_atmosphere_structure(self, usol_in, usol, &
                                                  molecules_per_particle, pressure, density, mix, mubar, &
                                                  pressure_hydro, density_hydro, apply_persistent_profile, err)
    use photochem_eqns, only: press_and_den
    use photochem_evoatmosphere_chemistry, only: molec_per_particle
    use photochem_const, only: N_avo, k_boltz
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    real(dp), intent(out) :: usol(:,:)
    real(dp), intent(out) :: molecules_per_particle(:,:)
    real(dp), intent(out) :: pressure(:), density(:), mix(:,:), mubar(:)
    real(dp), intent(out) :: pressure_hydro(:), density_hydro(:)
    logical, optional, intent(in) :: apply_persistent_profile
    character(:), allocatable, intent(out) :: err

    logical :: apply_profile
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    integer :: j

    dat => self%dat
    var => self%var

    apply_profile = .true.
    if (present(apply_persistent_profile)) apply_profile = apply_persistent_profile

    ! A persistent pressure-based profile depends on the trial composition.
    ! Apply it before boundary conditions, hydrostatics, transport, chemistry,
    ! and saturation quantities are prepared.
    if (apply_profile) then
      call apply_press_temp_edd_profile(self, usol_in, err)
      if (allocated(err)) return
    endif

    call clip_usol(usol_in, usol)

    call self%apply_lower_boundary_conditions(var%temperature(1), usol(:,1), err)
    if (allocated(err)) return

    !!! molecules/particle
    if (dat%there_are_particles) then
      call molec_per_particle(dat, var, molecules_per_particle)
    endif

    !!! pressure, density and mean molcular weight
    do j = 1,var%nz
      density(j) = sum(usol(dat%ng_1:,j))
      mix(:,j) = usol(:,j)/density(j) ! mixing ratios
      mubar(j) = sum(dat%species_mass(dat%ng_1:dat%nq)*mix(dat%ng_1:,j))
    enddo
    
    ! surface pressure by adding up all the mass in the atmosphere (bars)
    self%wrk%surface_pressure = sum(density(:)*mubar(:)*var%grav(:)*var%dz(:))/N_avo/1.0e6_dp
    call press_and_den(var%temperature, var%grav, self%wrk%surface_pressure*1.0e6_dp, var%dz, &
                       mubar, pressure_hydro, density_hydro)
    pressure(:) = density(:)*k_boltz*var%temperature(:)

  end subroutine

  module subroutine prep_atmosphere_unchecked(self, usol_in, apply_persistent_profile, err)

    use photochem_evoatmosphere_chemistry, only: reaction_rates, rainout, photorates
    use photochem_evoatmosphere_chemistry, only: gas_saturation_density
    use photochem_const, only: pi, N_avo, small_real, k_boltz
    use photochem_enum, only: DiffusionLimHydrogenEscape

    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    logical, optional, intent(in) :: apply_persistent_profile
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    integer :: i, j, k

    dat => self%dat
    var => self%var
    wrk => self%wrk

    call prepare_atmosphere_structure(self, usol_in, wrk%usol, &
                                      wrk%molecules_per_particle, wrk%pressure, wrk%density, &
                                      wrk%mix, wrk%mubar, wrk%pressure_hydro, wrk%density_hydro, &
                                      apply_persistent_profile, err)
    if (allocated(err)) return

    !!! diffusion and advection coefficients
    call diffusion_coefficients_evo(dat, var, wrk%density, wrk%mubar, &
    wrk%DU, wrk%DD, wrk%DL, wrk%ADU, wrk%ADL, wrk%ADD, wrk%wfall, wrk%VH2_esc, wrk%VH_esc)
    
    wrk%scale_height = (k_boltz*var%temperature(:)*N_avo)/(wrk%mubar(:)*var%grav(:))

    !!! H and H2 escape
    wrk%upper_veff_copy = var%upper_veff
    wrk%lower_vdep_copy = var%lower_vdep
    if (dat%H_escape_type == DiffusionLimHydrogenEscape) then
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

      call gas_saturation_density(dat, var, wrk%gas_sat_den)
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
    call reaction_rates(self%dat, self%var, wrk%pressure, wrk%density, wrk%densities, wrk%rx_rates)
    
    ! Update the photon_flux if the function is associated.
    ! we use time wrk%tn, which MUST be updated.
    if (associated(var%photon_flux_fcn)) then
      call var%photon_flux_fcn(wrk%tn, dat%nw, var%photon_flux)
    endif
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

  module subroutine prep_atmosphere(self, usol_in, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    character(:), allocatable, intent(out) :: err

    call self%require_atmosphere_initialized('prep_atmosphere', err)
    if (allocated(err)) return

    call prep_atmosphere_unchecked(self, usol_in, err=err)
    if (allocated(err)) return

  end subroutine

  module subroutine right_hand_side(self, neqs, tn, usol_flat, rhs, err)
    use photochem_enum, only: MosesBC, VelocityBC, DensityBC, PressureBC, FluxBC, VelocityDistributedFluxBC
    use photochem_enum, only: ZahnleHydrogenEscape
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real  
    
    class(EvoAtmosphere), target, intent(inout) :: self
    integer, intent(in) :: neqs
    real(dp), intent(in) :: tn
    real(dp), target, intent(in) :: usol_flat(neqs)
    real(dp), intent(out) :: rhs(neqs)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: disth, ztop, ztop1, boundary_correction
    integer :: i, k, j, jdisth
    
    real(dp), pointer :: usol_in(:,:)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    call self%require_atmosphere_initialized('right_hand_side', err)
    if (allocated(err)) return

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

    ! time
    wrk%tn = tn
    
    ! fills self%wrk with data
    call prep_atmosphere_unchecked(self, usol_in, err=err)
    if (allocated(err)) return

    call dochem(self, wrk%usol, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                wrk%density, wrk%mix, wrk%densities, wrk%xp, wrk%xl, rhs)  

    ! Extra functions specifying production or destruction
    do i = 1,dat%nq
      if (associated(var%rate_fcns(i)%fcn)) then
        call var%rate_fcns(i)%fcn(tn, var%nz, wrk%xp) ! using wrk%xp space.
        do j = 1,var%nz
          k = i + (j-1)*dat%nq
          rhs(k) = rhs(k) + wrk%xp(j) ! (molecules/cm^3/s)
        enddo
      endif
    enddo
    ! zahnle hydrogen escape
    if (dat%H_escape_type == ZahnleHydrogenEscape) then
      ! for Zahnle hydrogen escape, we pull H2 out of 
      ! the bottom grid cell of the model.
      rhs(dat%LH2) = rhs(dat%LH2) &
      - dat%H_escape_coeff*wrk%mix(dat%LH2,1)/var%dz(1)
    endif

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
      select case (var%lowerboundcond(i))
      case (DensityBC, PressureBC)
        ! Fixed lower boundaries replace the differential equation.
        rhs(i) = 0.0_dp
        cycle
      case (VelocityBC, VelocityDistributedFluxBC)
        boundary_correction = -wrk%lower_vdep_copy(i)*wrk%usol(i,1)/var%dz(1)
      case (FluxBC)
        boundary_correction = var%lower_flux(i)/var%dz(1)
      ! Moses (2001) boundary condition for gas giants
      ! A deposition velocity controled by how quickly gases
      ! turbulantly mix vertically
      case (MosesBC)
        boundary_correction = -(var%edd(1)/wrk%scale_height(1))*wrk%usol(i,1)/var%dz(1)
      case default
        err = 'Invalid lower boundary condition type'
        return
      end select

      rhs(i) = rhs(i) + wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                      + wrk%DD(i,1)*wrk%usol(i,1) + wrk%ADD(i,1)*wrk%usol(i,1) &
                      + boundary_correction
    enddo

    ! Upper boundary
    do i = 1,dat%nq
      k = i + (var%nz-1)*dat%nq
      select case (var%upperboundcond(i))
      case (VelocityBC)
        boundary_correction = -wrk%upper_veff_copy(i)*wrk%usol(i,var%nz)/var%dz(var%nz)
      case (FluxBC)
        boundary_correction = -var%upper_flux(i)/var%dz(var%nz)
      case default
        err = 'Invalid upper boundary condition type'
        return
      end select

      rhs(k) = rhs(k) + wrk%DD(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADD(i,var%nz)*wrk%usol(i,var%nz) &
                      + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1) &
                      + boundary_correction
    enddo

    ! Distributed (volcanic) sources
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityDistributedFluxBC) then
        disth = var%lower_dist_height(i)*1.e5_dp        
        if (disth < var%z(1) - 0.5_dp*var%dz(1)) then
        ! If the height is below the model domain, then we will put all flux into
        ! lowest layer.
        rhs(i) = rhs(i) + var%lower_flux(i)/var%dz(1)
        else
        ! If the height is within the model domain, then we will distribute the flux
        ! throught the model.
        jdisth = minloc(var%Z,1, var%Z >= disth) - 1
        jdisth = max(jdisth,2)
        ztop = var%z(jdisth)-var%z(1)
        ztop1 = var%z(jdisth) + 0.5e0_dp*var%dz(jdisth)
        do j = 2,jdisth
          k = i + (j-1)*dat%nq
          rhs(k) = rhs(k) + 2.0_dp*var%lower_flux(i)*(ztop1-var%z(j))/(ztop**2.0_dp)
        enddo
        endif
      endif
    enddo 

  end subroutine

  subroutine autodiff_chemistry_jacobian(self, usol, rhs, djac, err)
    use differentia, only: jacobian, dual, BlockDiagonalJacobian, initialize_dual_array
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), target, contiguous, intent(in) :: usol(:,:)
    real(dp), intent(out) :: rhs(:), djac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), pointer :: usol_flat(:)
    type(dual), allocatable :: mix(:,:), densities(:,:), xp(:), xl(:)
    integer :: blocksize

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    dat => self%dat
    var => self%var
    wrk => self%wrk

    usol_flat(1:var%neqs) => usol(:,:) ! reshape input
    blocksize = dat%nq
    ! allocate some work memory
    allocate(mix(dat%nq,var%nz))
    call initialize_dual_array(mix, blocksize)
    allocate(densities(size(wrk%densities,1),size(wrk%densities,2)))
    call initialize_dual_array(densities, blocksize)
    allocate(xp(size(wrk%xp)))
    call initialize_dual_array(xp, blocksize)
    allocate(xl(size(wrk%xl)))
    call initialize_dual_array(xl, blocksize)

    call jacobian(fcn, usol_flat, rhs, djac, jt=BlockDiagonalJacobian, blocksize=blocksize, err=err)
    if (allocated(err)) return

  contains
    subroutine fcn(x_, f_)
      type(dual), target, intent(in) :: x_(:)
      type(dual), target, intent(out) :: f_(:)
      type(dual), pointer :: usol_(:,:)
      usol_(1:dat%nq,1:var%nz) => x_(:)
      call initialize_dual_array(f_, blocksize)
      call dochem(self, usol_, wrk%rx_rates, &
                  wrk%gas_sat_den, wrk%molecules_per_particle, &
                  wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                  wrk%density, mix, densities, xp, xl, f_) 
    end subroutine
  end subroutine

  subroutine analytical_chemistry_jacobian(self, usol, djac, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), target, contiguous, intent(in) :: usol(:,:)
    real(dp), intent(out) :: djac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: density_derivatives(:,:)
    real(dp), allocatable :: reaction_derivative(:)
    real(dp), allocatable :: production_derivative(:), loss_derivative(:)
    integer :: j

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    dat => self%dat
    var => self%var
    wrk => self%wrk

    allocate(density_derivatives(dat%nsp+1,dat%nq))
    allocate(reaction_derivative(dat%nq))
    allocate(production_derivative(dat%nq))
    allocate(loss_derivative(dat%nq))
    djac = 0.0_dp

    ! Preserve the chemistry process order used by dochem without retaining a
    ! three-dimensional density-derivative workspace for every model layer.
    do j = 1,var%nz
      call prepare_chemistry_jacobian_layer( &
          dat, usol, wrk%rx_rates, wrk%molecules_per_particle, &
          wrk%densities, j, density_derivatives, reaction_derivative, &
          production_derivative, loss_derivative)
      call add_reaction_jacobian_layer( &
          dat, wrk%rx_rates, wrk%densities, density_derivatives, j, &
          reaction_derivative, djac)
      call add_rainout_jacobian_layer(dat, var, wrk%rainout_rates, j, djac)
      if (dat%there_are_particles) then
        call add_condensation_jacobian_layer( &
            dat, var, usol, wrk%gas_sat_den, wrk%scale_height, wrk%wfall, &
            j, djac)
      endif
    enddo

  end subroutine

  ! Prepare density derivatives, including diagnostic short-lived species,
  ! for one model layer.
  subroutine prepare_chemistry_jacobian_layer( &
      dat, usol, rx_rates, molecules_per_particle, densities, layer, &
      density_derivatives, reaction_derivative, production_derivative, &
      loss_derivative)
    use photochem_const, only: small_real
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: usol(:,:), rx_rates(:,:)
    real(dp), intent(in) :: molecules_per_particle(:,:)
    real(dp), intent(inout) :: densities(:,:)
    integer, intent(in) :: layer
    real(dp), intent(out) :: density_derivatives(:,:)
    real(dp), intent(out) :: reaction_derivative(:)
    real(dp), intent(out) :: production_derivative(:)
    real(dp), intent(out) :: loss_derivative(:)

    real(dp) :: reaction_rate, production, loss
    integer :: i, k, m, species

    density_derivatives = 0.0_dp

    ! Gas chemistry uses evolved gas densities directly. Evolved particle
    ! variables represent molecules/cm^3, while reaction rates use
    ! particles/cm^3. Match the lower clamp in prepare_chemistry_state.
    do i = 1,dat%npq
      if (usol(i,layer)/molecules_per_particle(i,layer) > small_real) then
        density_derivatives(i,i) = &
            1.0_dp/molecules_per_particle(i,layer)
      endif
    enddo
    do i = dat%ng_1,dat%nq
      density_derivatives(i,i) = 1.0_dp
    enddo

    ! Short-lived species are eliminated algebraically as n = P/L. Their
    ! derivatives therefore follow dn = (dP - n*dL)/L. Model validation
    ! forbids dependencies between short-lived species, so no implicit solve
    ! is required.
    do species = dat%nq+1,dat%nq+dat%nsl
      production = 0.0_dp
      loss = 0.0_dp
      production_derivative = 0.0_dp
      loss_derivative = 0.0_dp

      do k = 1,dat%pl(species)%nump
        m = dat%pl(species)%iprod(k)
        call reaction_rate_and_derivative( &
            dat, rx_rates, densities, density_derivatives, m, layer, 0, &
            reaction_rate, reaction_derivative)
        production = production + reaction_rate
        production_derivative = production_derivative + reaction_derivative
      enddo
      do k = 1,dat%pl(species)%numl
        m = dat%pl(species)%iloss(k)
        call reaction_rate_and_derivative( &
            dat, rx_rates, densities, density_derivatives, m, layer, &
            species, reaction_rate, reaction_derivative)
        loss = loss + reaction_rate
        loss_derivative = loss_derivative + reaction_derivative
      enddo

      densities(species,layer) = production/loss
      density_derivatives(species,:) = &
          (production_derivative - &
           densities(species,layer)*loss_derivative)/loss
    enddo

  end subroutine

  ! Add mass-action reaction derivatives for one model layer.
  subroutine add_reaction_jacobian_layer( &
      dat, rx_rates, densities, density_derivatives, layer, &
      reaction_derivative, djac)
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: rx_rates(:,:), densities(:,:)
    real(dp), intent(in) :: density_derivatives(:,:)
    integer, intent(in) :: layer
    real(dp), intent(out) :: reaction_derivative(:)
    real(dp), intent(inout) :: djac(:,:)

    real(dp) :: reaction_rate
    integer :: k, m, n, species

    n = dat%nq*(layer-1) + 1
    do m = 1,dat%nrT
      call reaction_rate_and_derivative( &
          dat, rx_rates, densities, density_derivatives, m, layer, 0, &
          reaction_rate, reaction_derivative)

      ! Reaction metadata repeats species indices for stoichiometric
      ! coefficients greater than one, so scatter once per occurrence.
      do k = 1,dat%rx(m)%nprod
        species = dat%rx(m)%prod_sp_inds(k)
        if (species <= dat%nq) then
          djac(species,n:n+dat%nq-1) = &
              djac(species,n:n+dat%nq-1) + reaction_derivative
        endif
      enddo
      do k = 1,dat%rx(m)%nreact
        species = dat%rx(m)%react_sp_inds(k)
        if (species <= dat%nq) then
          djac(species,n:n+dat%nq-1) = &
              djac(species,n:n+dat%nq-1) - reaction_derivative
        endif
      enddo
    enddo

  end subroutine

  ! Add fixed-coefficient gas-rainout derivatives for one model layer.
  subroutine add_rainout_jacobian_layer(dat, var, rainout_rates, layer, djac)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: rainout_rates(:,:)
    integer, intent(in) :: layer
    real(dp), intent(inout) :: djac(:,:)

    integer :: i, n

    if (.not.dat%gas_rainout .or. layer > var%trop_ind) return

    ! The tendency is -k_rain*n, so its only derivative is -k_rain on the
    ! corresponding species diagonal.
    n = dat%nq*(layer-1) + 1
    do i = 1,dat%nq
      djac(i,n+i-1) = djac(i,n+i-1) - rainout_rates(i,layer)
    enddo

  end subroutine

  ! Add condensation and evaporation derivatives for one model layer.
  subroutine add_condensation_jacobian_layer( &
      dat, var, usol, gas_sat_den, scale_height, wfall, layer, djac)
    use photochem_const, only: small_real
    use photochem_enum, only: CondensingParticle
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    real(dp), intent(in) :: usol(:,:), gas_sat_den(:,:)
    real(dp), intent(in) :: scale_height(:), wfall(:,:)
    integer, intent(in) :: layer
    real(dp), intent(inout) :: djac(:,:)

    real(dp) :: rh, rh_derivative, cond_rate0, cond_rate
    real(dp) :: rate_derivative, transfer_derivative
    integer :: i, ii, n

    n = dat%nq*(layer-1) + 1
    do i = 1,dat%np
      if (dat%particle_formation_method(i) /= CondensingParticle) cycle

      ii = dat%particle_gas_phase_ind(i)
      rh = max(usol(ii,layer)/gas_sat_den(i,layer),small_real)
      rh_derivative = 0.0_dp
      if (usol(ii,layer)/gas_sat_den(i,layer) > small_real) then
        rh_derivative = 1.0_dp/gas_sat_den(i,layer)
      endif

      if (rh > var%cond_params(i)%RHc) then
        ! Condensation C = k_cond(RH)*n_gas. Apply -dC to the gas equation and
        ! +dC to the corresponding particle equation.
        cond_rate0 = var%cond_params(i)%k_cond* &
                     (var%edd(layer)/scale_height(layer)**2.0_dp)
        cond_rate = damp_condensation_rate( &
            cond_rate0, var%cond_params(i)%RHc, &
            (1.0_dp + var%cond_params(i)%smooth_factor)* &
                var%cond_params(i)%RHc, rh)
        rate_derivative = damp_condensation_rate_derivative( &
            cond_rate0, var%cond_params(i)%RHc, &
            (1.0_dp + var%cond_params(i)%smooth_factor)* &
                var%cond_params(i)%RHc, rh)*rh_derivative
        transfer_derivative = cond_rate + usol(ii,layer)*rate_derivative

        djac(ii,n+ii-1) = djac(ii,n+ii-1) - transfer_derivative
        djac(i,n+ii-1) = djac(i,n+ii-1) + transfer_derivative

      elseif (var%evaporation) then
        ! Evaporation E = k_evap(1/RH)*n_particle. Apply +dE to gas and -dE
        ! to particles.
        cond_rate0 = var%cond_params(i)%k_evap* &
                     (wfall(i,layer)/scale_height(layer))
        cond_rate = damp_condensation_rate( &
            cond_rate0, 1.0_dp/var%cond_params(i)%RHc, &
            (1.0_dp + var%cond_params(i)%smooth_factor)/ &
                var%cond_params(i)%RHc, 1.0_dp/rh)
        rate_derivative = damp_inverse_rh_rate_derivative( &
            cond_rate0, var%cond_params(i)%RHc, &
            var%cond_params(i)%smooth_factor, rh)*rh_derivative
        transfer_derivative = usol(i,layer)*rate_derivative

        djac(ii,n+ii-1) = djac(ii,n+ii-1) + transfer_derivative
        djac(i,n+ii-1) = djac(i,n+ii-1) - transfer_derivative
        djac(ii,n+i-1) = djac(ii,n+i-1) + cond_rate
        djac(i,n+i-1) = djac(i,n+i-1) - cond_rate
      endif
    enddo

  end subroutine

  ! Evaluate one mass-action rate and its evolved-state derivative.
  subroutine reaction_rate_and_derivative( &
      dat, rx_rates, densities, density_derivatives, reaction, layer, &
      excluded_species, rate, derivative)
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: rx_rates(:,:), densities(:,:)
    real(dp), intent(in) :: density_derivatives(:,:)
    integer, intent(in) :: reaction, layer, excluded_species
    real(dp), intent(out) :: rate, derivative(:)

    real(dp) :: rate_without_reactant
    integer :: reactant, other_reactant, reactant_species

    ! Optionally omit an algebraically eliminated species to obtain its
    ! pseudo-first-order loss coefficient.
    rate = rx_rates(layer,reaction)
    do reactant = 1,dat%rx(reaction)%nreact
      reactant_species = dat%rx(reaction)%react_sp_inds(reactant)
      if (reactant_species /= excluded_species) then
        rate = rate*densities(reactant_species,layer)
      endif
    enddo

    ! Differentiate with a leave-one-reactant-out product. This avoids
    ! division by density and remains valid for repeated or zero reactants.
    derivative = 0.0_dp
    do reactant = 1,dat%rx(reaction)%nreact
      reactant_species = dat%rx(reaction)%react_sp_inds(reactant)
      if (reactant_species == excluded_species) cycle

      rate_without_reactant = rx_rates(layer,reaction)
      do other_reactant = 1,dat%rx(reaction)%nreact
        if (other_reactant /= reactant .and. &
            dat%rx(reaction)%react_sp_inds(other_reactant) /= &
                excluded_species) then
          rate_without_reactant = rate_without_reactant* &
              densities( &
                  dat%rx(reaction)%react_sp_inds(other_reactant),layer)
        endif
      enddo
      derivative = derivative + rate_without_reactant* &
                   density_derivatives(reactant_species,:)
    enddo

  end subroutine

  pure function damp_condensation_rate_derivative(A, rhc, rh0, rh) &
      result(derivative)
    use photochem_const, only: pi
    real(dp), intent(in) :: A, rhc, rh0, rh
    real(dp) :: derivative, argument

    ! Derivative with respect to RH of the arctangent smoothing function.
    argument = (rh-rhc)/(rh0-rhc)
    derivative = A*(2.0_dp/pi)/((rh0-rhc)*(1.0_dp+argument**2))

  end function

  pure function damp_inverse_rh_rate_derivative(A, rhc, smooth_factor, &
                                                 rh) result(derivative)
    use photochem_const, only: pi
    real(dp), intent(in) :: A, rhc, smooth_factor, rh
    real(dp) :: derivative, inverse_rhc, width

    ! Derivative with respect to RH of the evaporation smoothing function
    ! evaluated at 1/RH. This form avoids explicitly forming 1/RH**2.
    inverse_rhc = 1.0_dp/rhc
    width = smooth_factor/rhc
    derivative = -A*(2.0_dp/pi)*width/ &
                 ((width*rh)**2 + (1.0_dp-inverse_rhc*rh)**2)

  end function

  module subroutine jacobian(self, lda_neqs, neqs, usol_flat, jac, err)
    use photochem_enum, only: MosesBC, VelocityBC, DensityBC, PressureBC, FluxBC, VelocityDistributedFluxBC
    use photochem_enum, only: ZahnleHydrogenEscape
    use photochem_enum, only: AnalyticalJacobian, AutodiffJacobian, &
                              FiniteDifferenceJacobian
    use iso_c_binding, only: c_ptr, c_f_pointer
    use photochem_const, only: pi, small_real
    
    class(EvoAtmosphere), target, intent(inout) :: self
    integer, intent(in) :: lda_neqs, neqs
    real(dp), target, intent(in) :: usol_flat(neqs)
    real(dp), intent(out), target :: jac(lda_neqs)
    character(:), allocatable, intent(out) :: err
    
    real(dp), pointer :: usol_in(:,:)
    real(dp), pointer :: djac(:,:)
    real(dp) :: rhs(self%var%neqs), boundary_correction

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    integer :: i, k, j, m, mm
    logical :: insert_chemistry_jacobian
  
    call self%require_atmosphere_initialized('jacobian', err)
    if (allocated(err)) return

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
  
    call prep_atmosphere_unchecked(self, usol_in, err=err)
    if (allocated(err)) return
  
    jac = 0.0_dp
    insert_chemistry_jacobian = .false.
  
    select case (var%jacobian_method)
    case (FiniteDifferenceJacobian); block
    real(dp) :: usol_perturb(dat%nq,var%nz)
    real(dp) :: R(var%nz)
    real(dp) :: rhs_perturb(var%neqs)
    real(dp) :: mix(dat%nq,var%nz)
    real(dp) :: densities(dat%nsp+1,self%var%nz), xl(var%nz), xp(var%nz)
  
    ! Finite differenced Jacobian

    ! compute chemistry contribution to jacobian using forward differences
    call dochem(self, wrk%usol, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                wrk%density, mix, densities, xp, xl, rhs) 

    !$omp parallel private(i, j, k, m, mm, usol_perturb, R, mix, densities, xl, xp, rhs_perturb)
    usol_perturb = wrk%usol
    !$omp do
    do i = 1,dat%nq
      do j = 1,var%nz
        R(j) = var%epsj*abs(wrk%usol(i,j))
        usol_perturb(i,j) = wrk%usol(i,j) + R(j)
      enddo
      
      call dochem(self, usol_perturb, wrk%rx_rates, &
                  wrk%gas_sat_den, wrk%molecules_per_particle, &
                  wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                  wrk%density, mix, densities, xp, xl, rhs_perturb) 
  
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
    
    endblock

    case (AutodiffJacobian)

    call autodiff_chemistry_jacobian(self, wrk%usol, rhs, wrk%djac_chem, err)
    if (allocated(err)) return
    insert_chemistry_jacobian = .true.

    case (AnalyticalJacobian)

    call analytical_chemistry_jacobian(self, wrk%usol, wrk%djac_chem, err)
    if (allocated(err)) return

    insert_chemistry_jacobian = .true.

    case default
      err = 'Invalid chemistry Jacobian method'
      return
    end select

    if (insert_chemistry_jacobian) then
      do mm = 1,dat%nq
        do m = 1,var%nz
          j = mm + dat%nq*(m-1)
          do i = 1,dat%nq
            djac(i + 2*dat%nq - (mm-1),j) = wrk%djac_chem(i,j)
          enddo
        enddo
      enddo
    endif

    ! zahnle hydrogen escape
    if (dat%H_escape_type == ZahnleHydrogenEscape) then

      djac(dat%kd,dat%LH2) = djac(dat%kd,dat%LH2) &
      - (dat%H_escape_coeff/var%dz(1)) &
        *((wrk%density(1) - wrk%usol(dat%LH2,1))/wrk%density(1)**2.0_dp)

      do m = dat%ng_1,dat%nq 
        if (m /= dat%LH2) then
          mm = dat%kd + dat%LH2 - m
          djac(mm,m) = djac(mm,m) &
          - (dat%H_escape_coeff/var%dz(1)) &
            *(-wrk%usol(dat%LH2,1)/wrk%density(1)**2.0_dp)
        endif
      enddo
    endif
  
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
      select case (var%lowerboundcond(i))
      case (DensityBC, PressureBC)
        ! Fixed lower boundaries replace the differential equation.
        do m=1,dat%nq
          mm = dat%kd + i - m
          djac(mm,m) = 0.0_dp
        enddo
        djac(dat%ku,i+dat%nq) = 0.0_dp
        ! For some reason this term makes the integration
        ! much happier. I will keep it. Jacobians don't need to be perfect.
        djac(dat%kd,i) = - wrk%DU(i,1)
        cycle
      case (VelocityBC, VelocityDistributedFluxBC)
        boundary_correction = -wrk%lower_vdep_copy(i)/var%dz(1)
      case (FluxBC)
        boundary_correction = 0.0_dp
      case (MosesBC)
        boundary_correction = -(var%edd(1)/wrk%scale_height(1))/var%dz(1)
      case default
        err = 'Invalid lower boundary condition type'
        return
      end select

      djac(dat%ku,i+dat%nq) = wrk%DU(i,1) + wrk%ADU(i,1)
      djac(dat%kd,i) = djac(dat%kd,i) + wrk%DD(i,1) + wrk%ADD(i,1) + &
                       boundary_correction
    enddo
  
    ! Upper boundary
    do i = 1,dat%nq
      k = i + (var%nz-1)*dat%nq
      select case (var%upperboundcond(i))
      case (VelocityBC)
        boundary_correction = -wrk%upper_veff_copy(i)/var%dz(var%nz)
      case (FluxBC)
        boundary_correction = 0.0_dp
      case default
        err = 'Invalid upper boundary condition type'
        return
      end select

      djac(dat%kd,k) = djac(dat%kd,k) + wrk%DD(i,var%nz) + wrk%ADD(i,var%nz) + &
                       boundary_correction
      djac(dat%kl,k-dat%nq) = wrk%DL(i,var%nz) + wrk%ADL(i,var%nz)
    enddo
  
  end subroutine

  module subroutine chemistry_right_hand_side(self, usol, rhs, err)
    use photochem_enum, only: ZahnleHydrogenEscape
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(out) :: rhs(:)
    character(:), allocatable, intent(out) :: err
    
    integer :: i, j, k
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    call self%require_atmosphere_initialized('chemistry_right_hand_side', err)
    if (allocated(err)) return

    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz .or. size(rhs) /= var%neqs) then
      err = "Input usol or rhs has the wrong dimensions"
      return
    endif

    call self%prep_atmosphere(usol, err)
    if (allocated(err)) return
    
    call dochem(self, wrk%usol, wrk%rx_rates, &
                wrk%gas_sat_den, wrk%molecules_per_particle, &
                wrk%rainout_rates, wrk%scale_height, wrk%wfall, &
                wrk%density, wrk%mix, wrk%densities, wrk%xp, wrk%xl, rhs) 

    ! We additionally have to add the below. I treat them as "chemistry"
  
    ! Extra functions specifying production or destruction
    do i = 1,dat%nq
      if (associated(var%rate_fcns(i)%fcn)) then
        call var%rate_fcns(i)%fcn(wrk%tn, var%nz, wrk%xp) ! using wrk%xp space.
        do j = 1,var%nz
          k = i + (j-1)*dat%nq
          rhs(k) = rhs(k) + wrk%xp(j) ! (molecules/cm^3/s)
        enddo
      endif
    enddo
    ! Zahnle hydrogen escape
    if (dat%H_escape_type == ZahnleHydrogenEscape) then
      ! for Zahnle hydrogen escape, we pull H2 out of 
      ! the bottom grid cell of the model.
      rhs(dat%LH2) = rhs(dat%LH2) &
      - dat%H_escape_coeff*wrk%mix(dat%LH2,1)/var%dz(1)
    endif
                              
  end subroutine

  module subroutine production_and_loss(self, species, usol, pl, err)     
    use futils, only: argsort            
    use photochem_const, only: small_real
  
    class(EvoAtmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: species
    real(dp), intent(in) :: usol(:,:)
    type(ProductionLoss), intent(out) :: pl
    character(:), allocatable, intent(out) :: err
  
    real(dp) :: xl(self%var%nz), xp(self%var%nz)
    integer, allocatable :: prod_inds(:), loss_inds(:)
    integer :: sp_ind
    integer :: i, j, k, np, nl, nlT
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
  
    call self%require_atmosphere_initialized('production_and_loss', err)
    if (allocated(err)) return

    dat => self%dat
    var => self%var
    wrk => self%wrk
  
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = "Input usol to production_and_loss has the wrong dimensions"
      return
    endif
  
    sp_ind = findloc(dat%species_names(1:dat%nsp),species,1)
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
      do i = 1,dat%npq
        wrk%densities(i,j) = max(wrk%usol(i,j)*(1.0_dp/wrk%molecules_per_particle(i,j)), small_real)
      enddo
      do i = dat%ng_1,dat%nq
        wrk%densities(i,j) = wrk%usol(i,j)
      enddo
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
          wrk%rainout_rates(sp_ind,1:var%trop_ind)*wrk%usol(sp_ind,1:var%trop_ind)
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
