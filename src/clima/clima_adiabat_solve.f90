
submodule(clima_adiabat) clima_adiabat_solve
  implicit none

contains

  module subroutine AdiabatClimate_make_profile_rc(self, P_i_surf, T_in, err)
    use clima_const, only: k_boltz, N_avo
    use clima_adiabat_rc, only: make_profile_rc
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: T_in(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: P_e(:), z_e(:), f_i_e(:,:), lapse_rate_e(:)
    logical, allocatable :: super_saturated_e(:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    allocate(P_e(2*self%nz+1),z_e(2*self%nz+1),f_i_e(2*self%nz+1,self%sp%ng),lapse_rate_e(2*self%nz+1))
    allocate(super_saturated_e(2*self%nz+1))
    allocate(density(self%nz))

    if (size(P_i_surf) /= self%sp%ng) then
      err = "P_i_surf has the wrong dimension"
      return
    endif
    if (size(T_in) /= self%nz+1) then
      err = "T_in has the wrong dimension"
      return
    endif

    self%T_surf = T_in(1)
    self%T = T_in(2:)
    call make_profile_rc(self%T_surf, self%T, P_i_surf, &
                         self%sp_custom, self%mix_custom, &
                         self%convecting_with_below, &
                         self%sp, self%nz, self%planet_mass, &
                         self%planet_radius, self%P_top, self%RH, &
                         self%rtol, self%atol, &
                         self%ocean_fcns, self%ocean_args_p, &
                         P_e, z_e, f_i_e, lapse_rate_e, super_saturated_e, &
                         self%N_surface, self%N_ocean, &
                         err)
    if (allocated(err)) return

    z_e(2*self%nz+1) = z_e(2*self%nz) + (z_e(2*self%nz) - z_e(2*self%nz-1))

    self%f_i_surf = f_i_e(1,:)
    self%P_surf = P_e(1)
    do i = 1,self%nz
      self%P(i) = P_e(2*i)
      do j =1,self%sp%ng
        self%f_i(i,j) = f_i_e(2*i,j)
      enddo
    enddo

    call self%compute_altitude(err)
    if (allocated(err)) return

    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    call self%interpolate_particles(self%P, err)
    if (allocated(err)) return

    do i = 1,self%sp%ng
      ! mol/cm^2 in atmosphere
      self%N_atmos(i) = sum(density*self%f_i(:,i)*self%dz)/N_avo
    enddo

    ! lapse rates
    self%lapse_rate_intended(1) = lapse_rate_e(1)
    do i = 2,self%nz
      self%lapse_rate_intended(i) = lapse_rate_e(2*i-2)
    enddo

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo

    do i = 1,self%nz
      self%super_saturated(i) = super_saturated_e(2*i)
    enddo

  end subroutine

  !> initializes custom mixing ratios
  subroutine intialize_custom_inputs(self, sp_custom, P_custom, mix_custom, err)
    class(AdiabatClimate), intent(inout) :: self
    character(*), optional, intent(in) :: sp_custom(:)
    real(dp), optional, intent(in) :: P_custom(:)
    real(dp), optional, intent(in) :: mix_custom(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: mix_custom_copy(:,:)
    real(dp), allocatable :: log10P_interp(:), log10mix_interp(:)
    integer :: i, ind, istat

    ! Ensure all are input if just one is input
    if (present(sp_custom) .or. present(P_custom) .or. present(mix_custom)) then
      if (.not.present(sp_custom)) then
        err = '`sp_custom` must be an input if `P_custom` and `mix_custom` are inputs'
        return
      endif
      if (.not.present(P_custom)) then
        err = '`P_custom` must be an input if `sp_custom` and `mix_custom` are inputs'
        return
      endif
      if (.not.present(mix_custom)) then
        err = '`mix_custom` must be an input if `P_custom` and `sp_custom` are inputs'
        return
      endif
    endif

    self%sp_custom = .false.
    if (.not.present(sp_custom)) return ! no custom inputs

    ! Check shapes
    if (size(sp_custom) /= size(mix_custom,2)) then
      err = '`sp_custom` and `mix_custom` have incompatible shapes'
      return
    endif
    if (size(P_custom) /= size(mix_custom,1)) then
      err = '`P_custom` and `mix_custom` have incompatible shapes'
      return
    endif
    if (any(mix_custom < 0.0_dp)) then
      err = '`mix_custom` can not have negative values'
      return
    endif
    if (any(P_custom <= 0.0_dp)) then
      err = '`P_custom` must be >= 0 for all values'
      return
    endif

    ! Normalize mixing ratios so they sum to 1
    mix_custom_copy = mix_custom
    do i = 1,size(mix_custom,1)
      mix_custom_copy(i,:) = mix_custom_copy(i,:)/sum(mix_custom_copy(i,:))
    enddo

    ! Pressures. Ensure constant extrapolation
    log10P_interp = [[huge(1.0_dp),P_custom], tiny(1.0_dp)]
    log10P_interp = log10(log10P_interp(size(log10P_interp):1:-1))
    allocate(log10mix_interp(size(log10P_interp)))

    do i = 1,size(sp_custom)
      ind = findloc(self%species_names, sp_custom(i), 1)
      if (ind == 0) then
        err = 'Custom species "'//trim(sp_custom(i))//'" is not in the list of species'
        return
      endif

      log10mix_interp(2:size(log10mix_interp)-1) = mix_custom_copy(:,i)
      log10mix_interp(1) = mix_custom_copy(1,i)
      log10mix_interp(size(log10mix_interp)) = mix_custom_copy(size(mix_custom_copy,1),i)
      log10mix_interp = log10(log10mix_interp(size(log10mix_interp):1:-1))

      self%sp_custom(ind) = .true.
      call self%mix_custom(ind)%initialize(log10P_interp, log10mix_interp, istat)
      if (istat /= 0) then
        err = 'Interpolation initialization for custom mixing ratios failed'
        return
      endif
    enddo

  end subroutine

  module function AdiabatClimate_RCE(self, P_i_surf, T_surf_guess, T_guess, convecting_with_below, &
    sp_custom, P_custom, mix_custom, err) result(converged)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_surf_guess
    real(dp), intent(in) :: T_guess(:)
    logical, optional, intent(in) :: convecting_with_below(:)
    character(*), optional, intent(in) :: sp_custom(:)
    real(dp), optional, intent(in) :: P_custom(:)
    real(dp), optional, intent(in) :: mix_custom(:,:)
    character(:), allocatable, intent(out) :: err
    logical :: converged

    logical, allocatable :: convecting_with_below_save(:,:)
    real(dp), allocatable :: T_in(:), x_init(:), dFdt(:), x_stage(:), x_sol(:), f_sol(:)
    integer :: i, j, k, mode_update
    logical :: solver_ok, perform_solve
    logical :: mask_changed

    converged = .false.

    if (.not.self%double_radiative_grid) then
      err = 'AdiabatClimate must be initialized with "double_radiative_grid" set to True '// &
            'in order to call RCE.'
      return
    endif
    if (size(T_guess) /= self%nz) then
      err = "T_guess has the wrong dimension"
      return
    endif
    if (self%max_rc_iters < 1) then
      return
    endif

    call intialize_custom_inputs(self, sp_custom, P_custom, mix_custom, err)
    if (allocated(err)) return
    
    allocate(convecting_with_below_save(self%nz,0))
    allocate(T_in(self%nz+1))

    T_in(1) = T_surf_guess
    T_in(2:) = T_guess(:)
    self%T_surf = T_surf_guess
    self%T(:) = T_guess(:)
    self%prevent_overconvection_lock = 0

    ! setup the convecting zone
    if (present(convecting_with_below)) then
      call AdiabatClimate_set_convecting_zones(self, convecting_with_below, err)
      if (allocated(err)) return
    else
      self%convecting_with_below = .false.
      call AdiabatClimate_update_convecting_zones(self, P_i_surf, T_in, mode=1, err=err)
      if (allocated(err)) return
    endif

    perform_solve = .true.
    mode_update = 1
    if (self%max_rc_iters_convection <= 1) then
      mode_update = 2
    endif

    do i = 1,self%max_rc_iters
      j = i

      if (self%verbose) then
        print"(1x,'Iteration =',i3,', Mode =',i3)", i, mode_update
      endif

      if (perform_solve) then

      if (allocated(x_init)) deallocate(x_init)
      allocate(x_init(size(self%inds_Tx)))
      x_init(1) = self%T_surf
      do k = 2,size(self%inds_Tx)
        x_init(k) = self%T(self%inds_Tx(k)-1)
      enddo
      if (allocated(dFdt)) deallocate(dFdt)
      allocate(dFdt(size(self%inds_Tx)))
      if (allocated(x_stage)) deallocate(x_stage)
      allocate(x_stage(size(self%inds_Tx)))
      if (allocated(x_sol)) deallocate(x_sol)
      allocate(x_sol(size(self%inds_Tx)))
      if (allocated(f_sol)) deallocate(f_sol)
      allocate(f_sol(size(self%inds_Tx)))

      select case (self%rce_solve_strategy)
      case (RCE_SOLVE_HYBRJ_ONLY)
        call run_hybrj(self, P_i_surf, x_init, x_sol, f_sol, dFdt, solver_ok)
        if (.not. solver_ok) then
          err = 'hybrj root solve failed in RCE (HYBRJ_ONLY).'
          return
        endif

      case (RCE_SOLVE_PTC_THEN_HYBRJ)
        call run_ptc(self, P_i_surf, x_init, x_stage, f_sol, dFdt, solver_ok)
        if (solver_ok) then
          ! If converged, then we are good
          x_sol = x_stage
        else
          ! If not converged, we tighten with hybrd
          call run_hybrj(self, P_i_surf, x_stage, x_sol, f_sol, dFdt, solver_ok)
        endif
        if (.not. solver_ok) then
          err = 'hybrj root solve failed in RCE (PTC_THEN_HYBRJ).'
          return
        endif

      case (RCE_SOLVE_HYBRJ_THEN_PTC_THEN_HYBRJ)
        call run_hybrj(self, P_i_surf, x_init, x_stage, f_sol, dFdt, solver_ok)
        if (solver_ok) then
          x_sol = x_stage
        else
          call run_ptc(self, P_i_surf, x_init, x_stage, f_sol, dFdt, solver_ok)
          if (solver_ok) then
            ! If converged, then we are good
            x_sol = x_stage
          else
            ! If not converged, we tighten with hybrd
            call run_hybrj(self, P_i_surf, x_stage, x_sol, f_sol, dFdt, solver_ok)
          endif
        endif
        if (.not. solver_ok) then
          err = 'hybrj root solve failed in RCE (HYBRJ_THEN_PTC_THEN_HYBRJ).'
          return
        endif

      case default
        err = 'Invalid rce_solve_strategy.'
        return
      end select

      ! Update all variables to the current root
      call AdiabatClimate_objective(self, P_i_surf, x_sol, dFdt, f_sol, err)
      if (allocated(err)) return

      endif
      perform_solve = .true.

      ! Save the current convective zones
      convecting_with_below_save = reshape(convecting_with_below_save,shape=[self%nz,i],pad=self%convecting_with_below)

      call AdiabatClimate_update_convecting_zones(self, P_i_surf, [self%T_surf, self%T], mode_update, err)
      if (allocated(err)) return
      mask_changed = .not. all(convecting_with_below_save(:,i) .eqv. self%convecting_with_below)

      ! Change modes and check convergence
      select case (mode_update)
      case (1)
        if (.not.mask_changed) then
          if (self%require_mode2) then
            ! We must pass through mode 2
            mode_update = 2
            perform_solve = .false.
            cycle
          endif
          if (self%prevent_overconvection) then
            ! We must skip to mode 3
            mode_update = 3
            perform_solve = .false.
            cycle
          endif
          ! If we get here, then we are converged
          converged = .true.
          exit
        else
          ! Mask is still changing
          if (i >= self%max_rc_iters_convection - 1) then
            ! We skip to mode 2
            mode_update = 2
            cycle
          endif
        endif
      case (2)
        if (.not.mask_changed) then
          if (self%prevent_overconvection) then
            ! We move on to mode 3
            mode_update = 3
            perform_solve = .false.
            cycle
          endif
          ! Otherwise we are converged
          converged = .true.
          exit
        endif
      case (3)
        if (.not.mask_changed) then
          converged = .true.
          exit
        endif
      end select

    enddo

    if (converged .and. self%verbose) then
      print'(1x,A)','CONVERGED'
    endif

    ! Return all information to what is was prior to checking for the root
    call AdiabatClimate_set_convecting_zones(self, convecting_with_below_save(:,j), err)
    if (allocated(err)) return
    call AdiabatClimate_objective(self, P_i_surf, x_sol, dFdt, f_sol, err)
    if (allocated(err)) return

  end function

  subroutine run_hybrj(self, P_i_surf, x_seed, x_out, f_out, dFdt, ok)
    use clima_useful, only: MinpackHybrj
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: x_seed(:)
    real(dp), intent(out) :: x_out(:), f_out(:)
    real(dp), intent(inout) :: dFdt(:)
    logical, intent(out) :: ok

    type(MinpackHybrj) :: mv
    integer :: k
    real(dp) :: perturbation
    logical :: minpack_custom_converged, have_base
    real(dp), allocatable :: dTdt_base(:), x_base(:), x_conv(:), f_conv(:)
    character(:), allocatable :: err

    allocate(dTdt_base(size(self%inds_Tx)))
    allocate(x_base(size(self%inds_Tx)))
    allocate(x_conv(size(self%inds_Tx)))
    allocate(f_conv(size(self%inds_Tx)))

    mv = MinpackHybrj(fcn, size(self%inds_Tx))
    mv%xtol = 1.0e-12_dp
    mv%nprint = 1
    mv%maxfev = 100

    k = 0
    do
      if (mod(k,2) == 0) then
        perturbation = real(k,dp)*1.0_dp
      else
        perturbation = -real(k,dp)*1.0_dp
      endif

      if (self%verbose .and. k > 0) then
        print'(3x,"Perturbation = ",f7.1)',perturbation
      endif

      minpack_custom_converged = .false.
      have_base = .false.
      x_conv(:) = 0.0_dp
      f_conv(:) = 0.0_dp
      mv%x = x_seed + perturbation
      call mv%hybrj()
      if (minpack_custom_converged) then
        x_out = x_conv
        f_out = f_conv
        ok = .true.
        return
      endif

      if (k > 2) then
        ok = .false.
        return
      endif

      k = k + 1
    enddo
    
  contains

    subroutine fcn(n_, x_, fvec_, fjac_, ldfjac_, iflag_)
      implicit none
      integer, intent(in) :: n_
      real(dp), dimension(n_), intent(in) :: x_
      integer, intent(in) :: ldfjac_
      real(dp), dimension(n_), intent(inout) :: fvec_
      real(dp), dimension(ldfjac_, n_), intent(inout) :: fjac_
      integer, intent(inout) :: iflag_
      real(dp) :: max_flux_imbalance_wm2_, max_flux_ratio_
      logical :: compute_objective_

      if (iflag_ == 1) then
        ! Compute right-hand-side
        call AdiabatClimate_objective(self, P_i_surf, x_, dFdt, fvec_, err)
        if (allocated(err)) then
          deallocate(err)
          iflag_ = -1
          return
        endif
        dTdt_base(:) = fvec_
        x_base(:) = x_
        have_base = .true.
        if (custom_flux_converged(self, dFdt)) then
          minpack_custom_converged = .true.
          x_conv(:) = x_
          f_conv(:) = fvec_
          iflag_ = -77
          return
        endif
      elseif (iflag_ == 2) then
        ! Compute jacobian
        compute_objective_ = .not. have_base
        if (.not. compute_objective_) then
          compute_objective_ = any(x_ /= x_base)
        endif
        if (compute_objective_) then
          call AdiabatClimate_objective(self, P_i_surf, x_, dFdt, dTdt_base, err)
          if (allocated(err)) then
            deallocate(err)
            iflag_ = -1
            return
          endif
          x_base(:) = x_
          have_base = .true.
        endif
        call AdiabatClimate_jacobian_from_base(self, P_i_surf, x_, dFdt, dTdt_base, fjac_, err)
        if (allocated(err)) then
          deallocate(err)
          iflag_ = -1
          return
        endif
        ! Base state was consumed to build this Jacobian; force refresh on next Jacobian call.
        have_base = .false.
      endif

      if (iflag_ == 0 .and. self%verbose) then
        call get_flux_metrics(self, dFdt, max_flux_imbalance_wm2_, max_flux_ratio_)
        print"(3x,'step =',i4,3x,'njev =',i10,3x,'max|F| = ',es9.2,3x,'max|F/F0| = ',es9.2,3x,'max(T) = ',f7.1,3x,'min(T) = ',f7.1)", &
              mv%nfev, mv%njev, max_flux_imbalance_wm2_, &
              max_flux_ratio_, maxval(x_), minval(x_)
      endif

    end subroutine

  end subroutine

  subroutine run_ptc(self, P_i_surf, x_seed, x_out, f_out, dFdt, ok)
    use clima_ptc, only: PTCSolver, PTC_JAC_DENSE, PTC_CONVERGED_USER
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: x_seed(:)
    real(dp), intent(out) :: x_out(:)
    real(dp), intent(out) :: f_out(:)
    real(dp), intent(inout) :: dFdt(:)
    logical, intent(out) :: ok

    type(PTCSolver) :: solver
    real(dp), allocatable :: dTdt_base(:), x_base(:)
    logical :: have_base
    character(:), allocatable :: err

    allocate(dTdt_base(size(self%inds_Tx)))
    allocate(x_base(size(self%inds_Tx)))
    have_base = .false.

    call solver%initialize(x_seed, f_ptc, jac_ptc, PTC_JAC_DENSE, dt_increment=self%dt_increment, max_steps=300)
    if (self%verbose) call solver%set_progress(progress_ptc)
    call solver%set_custom_convergence(convergence_ptc)
    call solver%solve()
    x_out = solver%x
    f_out = solver%fvec
    ok = (solver%reason == PTC_CONVERGED_USER)

  contains
    
    subroutine f_ptc(solver_, u_, udot_, ierr_)
      use clima_ptc, only: PTCSolver, wp
      class(PTCSolver), intent(in) :: solver_
      real(wp), intent(in) :: u_(:)
      real(wp), intent(out) :: udot_(:)
      integer, intent(out) :: ierr_

      ierr_ = 0

      call AdiabatClimate_objective(self, P_i_surf, u_, dFdt, udot_, err)
      if (allocated(err)) then
        deallocate(err)
        ierr_ = 1
        return
      endif
      dTdt_base(:) = udot_
      x_base(:) = u_
      have_base = .true.

    end subroutine

    subroutine jac_ptc(solver_, u_, jac_, ierr_)
      use clima_ptc, only: PTCSolver, wp
      class(PTCSolver), intent(in) :: solver_
      real(wp), intent(in) :: u_(:)
      real(wp), intent(out) :: jac_(:, :)
      integer, intent(out) :: ierr_

      logical :: compute_objective_

      ierr_ = 0
      if (allocated(err)) then
        deallocate(err)
        ierr_ = 1
        return
      endif

      compute_objective_ = .not. have_base
      if (.not. compute_objective_) then
        compute_objective_ = any(u_ /= x_base)
      endif
      if (compute_objective_) then
        call AdiabatClimate_objective(self, P_i_surf, u_, dFdt, dTdt_base, err)
        if (allocated(err)) then
          deallocate(err)
          ierr_ = 1
          return
        endif
        x_base(:) = u_
        have_base = .true.
      endif

      call AdiabatClimate_jacobian_from_base(self, P_i_surf, u_, dFdt, dTdt_base, jac_, err)
      if (allocated(err)) then
        deallocate(err)
        ierr_ = 1
        return
      endif
      ! Base state was consumed to build this Jacobian; force refresh on next Jacobian call.
      have_base = .false.

    end subroutine

    subroutine convergence_ptc(solver_, converged_, ierr_)
      use clima_ptc, only: PTCSolver
      class(PTCSolver), intent(in) :: solver_
      logical, intent(out) :: converged_
      integer, intent(out) :: ierr_
      ierr_ = 0
      converged_ = custom_flux_converged(self, dFdt)
    end subroutine

    subroutine progress_ptc(solver_)
      use clima_ptc, only: PTCSolver
      class(PTCSolver), intent(in) :: solver_
      real(dp) :: max_flux_imbalance_wm2_, max_flux_ratio_
      call get_flux_metrics(self, dFdt, max_flux_imbalance_wm2_, max_flux_ratio_)
      print"(3x,'step =',i4,3x,'dt   =',es10.3,3x,'max|F| = ',es9.2,3x,'max|F/F0| = ',es9.2,3x,'max(T) = ',f7.1,3x,'min(T) = ',f7.1)", &
            solver_%steps, solver_%dt, &
            max_flux_imbalance_wm2_, max_flux_ratio_, &
            maxval(solver_%x), minval(solver_%x)
    end subroutine

  end subroutine

  subroutine get_flux_metrics(self, dFdt, max_flux_imbalance_wm2, max_flux_ratio)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: dFdt(:)
    real(dp), intent(out) :: max_flux_imbalance_wm2, max_flux_ratio
    real(dp) :: characteristic_flux_wm2

    ! Characteristic flux scale of the atmosphere in W/m^2.
    characteristic_flux_wm2 = abs(self%rad%bolometric_flux()/4.0_dp + self%surface_heat_flow*1.0e-3_dp)
    characteristic_flux_wm2 = max(characteristic_flux_wm2, 1.0e-6_dp)

    ! Maximum flux imbalance in W/m^2.
    max_flux_imbalance_wm2 = maxval(abs(dFdt))*1.0e-3_dp
    max_flux_ratio = max_flux_imbalance_wm2/characteristic_flux_wm2

  end subroutine

  function custom_flux_converged(self, dFdt) result(converged_)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: dFdt(:)
    logical :: converged_
    real(dp) :: max_flux_imbalance_wm2, max_flux_ratio

    ! Converged if energy is conserved.
    call get_flux_metrics(self, dFdt, max_flux_imbalance_wm2, max_flux_ratio)
    converged_ = max_flux_ratio < self%xtol_rc

  end function

  subroutine AdiabatClimate_objective(self, P_i_surf, x, dFdt, dTdt, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: dFdt(:), dTdt(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: T_in(:)
    integer :: i

    allocate(T_in(self%nz+1))
    do i = 1,size(self%inds_Tx)
      T_in(self%inds_Tx(i)) = x(i)
    enddo

    ! Sets, self%T_surf, self%T, self%N_surface, self%N_ocean, self%P_surf, self%P,
    ! self%z, self%dz, self%f_i, self%densities, self%N_atmos, self%lapse_rate_intended
    ! self%lapse_rate
    call self%make_profile_rc(P_i_surf, T_in, err)
    if (allocated(err)) return

    T_in(1) = self%T_surf
    T_in(2:) = self%T

    ! resets self%T_surf, self%T, self%densities, self%lapse_rate
    ! also does radiative transfer and computes res
    call AdiabatClimate_objective_(self, P_i_surf, T_in, .true., .true., dFdt, dTdt, err)
    if (allocated(err)) return

  end subroutine

  subroutine AdiabatClimate_objective_(self, P_i_surf, T_in, compute_solar, compute_opacity, dFdt, dTdt, err)
    use clima_const, only: k_boltz
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    logical, intent(in) :: compute_solar, compute_opacity
    real(dp), intent(out) :: dFdt(:), dTdt(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: f_total(:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    ! work storage
    allocate(f_total(self%nz+1))
    allocate(density(self%nz))

    ! Stuff set in make_profile_rc that does not change
    ! - self%N_surface, self%N_ocean, self%P_surf, self%P, self%z, 
    ! - self%dz, self%f_i, self%N_atmos, self%lapse_rate_intended
    ! Stuff that does change
    ! - self%T_surf, self%T, self%densities, self%lapse_rate
    self%T_surf = T_in(1)
    self%T = T_in(2:)

    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo

    ! radiative transfer
    call self%copy_atm_to_radiative_grid()
    call self%rad%radiate(self%T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
      self%pdensities_r, self%pradii_r, compute_solar=compute_solar, compute_opacity=compute_opacity, err=err)
    if (allocated(err)) return

    ! Considers a tidally locked dayside climate. We only compute this for
    ! compute_solar == True.
    if (self%tidally_locked_dayside .and. compute_solar) then; block
      real(dp) :: tau_LW, k_term, f_term, rad_enhancement
      call self%heat_redistribution_parameters(tau_LW, k_term, f_term, err)
      if (allocated(err)) return

      rad_enhancement = 4.0_dp*f_term
      call self%rad%apply_radiation_enhancement(rad_enhancement)
    endblock; endif

    do i = 1,self%nz+1
      f_total(i) = self%rad%f_total(2*i-1)
    enddo
    f_total(1) = f_total(1) + self%surface_heat_flow

    call AdiabatClimate_residuals_with_convection(self, f_total, self%lapse_rate, self%lapse_rate_intended, dFdt, dTdt, err)
    if (allocated(err)) return

  end subroutine

  subroutine AdiabatClimate_jacobian(self, P_i_surf, x, jac, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: jac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: dFdt(:), dTdt(:)

    ! Check inputs
    if (size(jac,1) /= size(self%inds_Tx) .or. size(jac,2) /= size(self%inds_Tx)) then
      err = "jac has the wrong dimension"
      return
    endif

    ! allocate work
    allocate(dFdt(size(self%inds_Tx)), dTdt(size(self%inds_Tx)))

    ! First evaluate res at T.
    call AdiabatClimate_objective(self, P_i_surf, x, dFdt, dTdt, err)
    if (allocated(err)) return

    call AdiabatClimate_jacobian_from_base(self, P_i_surf, x, dFdt, dTdt, jac, err)
    if (allocated(err)) return

  end subroutine

  subroutine AdiabatClimate_jacobian_from_base(self, P_i_surf, x, dFdt, dTdt, jac, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: dFdt(:), dTdt(:)
    real(dp), intent(out) :: jac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: dFdt_perturb(:), dTdt_perturb(:), T_perturb(:), T_in(:)
    real(dp) :: deltaT
    integer :: i, ind

    ! Check inputs
    if (size(jac,1) /= size(self%inds_Tx) .or. size(jac,2) /= size(self%inds_Tx)) then
      err = "jac has the wrong dimension"
      return
    endif
    if (size(dFdt) /= size(self%inds_Tx) .or. size(dTdt) /= size(self%inds_Tx)) then
      err = "Base residual vectors have wrong dimension in AdiabatClimate_jacobian_from_base"
      return
    endif

    ! allocate work
    allocate(dFdt_perturb(size(self%inds_Tx)), dTdt_perturb(size(self%inds_Tx)))
    allocate(T_in(self%nz+1), T_perturb(self%nz+1))

    T_in(1) = self%T_surf
    T_in(2:) = self%T

    T_perturb(:) = T_in(:)
    do i = 1,size(x)

      ! Perturb temperature
      deltaT = self%epsj*abs(x(i))
      T_perturb(self%inds_Tx(i)) = T_in(self%inds_Tx(i)) + deltaT

      ! We perturb the temp also in a convecting region
      ind = findloc(self%ind_conv_lower_x, i, 1)
      if (ind /= 0) then
        T_perturb(self%ind_conv_lower(ind):self%ind_conv_upper(ind)) = &
        T_in(self%ind_conv_lower(ind):self%ind_conv_upper(ind)) + deltaT
      endif

      call AdiabatClimate_objective_(self, P_i_surf, T_perturb, self%compute_solar_in_jac, &
                                    .false., dFdt_perturb, dTdt_perturb, err)
      if (allocated(err)) return

      ! Compute jacobian
      jac(:,i) = (dTdt_perturb(:) - dTdt(:))/deltaT

      ! unperturb T
      T_perturb(:) = T_in(:)
    enddo

  end subroutine

  subroutine AdiabatClimate_set_convecting_zones(self, convecting_with_below, err)
    class(AdiabatClimate), intent(inout) :: self
    !> Length nz. States whether the layer below is convecting
    !> with the current layer. `convecting(1)` determines whether
    !> the first atmospheric layer is convecting with the ground.
    logical, intent(in) :: convecting_with_below(:)
    character(:), allocatable, intent(out) :: err

    integer :: i,j,k,ind

    if (size(convecting_with_below) /= self%nz) then
      err = 'Input "convecting_with_below" has the wrong dimension'
      return
    endif

    self%convecting_with_below = convecting_with_below

    self%n_convecting_zones = 0
    if (allocated(self%ind_conv_lower)) deallocate(self%ind_conv_lower)
    if (allocated(self%ind_conv_upper)) deallocate(self%ind_conv_upper)
    allocate(self%ind_conv_lower(0),self%ind_conv_upper(0))
    i = 1
    do while(i <= self%nz)
      if (convecting_with_below(i)) then
        ! identified a convecting zone
        self%n_convecting_zones = self%n_convecting_zones + 1
        ! Save the lower index
        self%ind_conv_lower = [self%ind_conv_lower, i]
        ! Search for the upper bound of convecting zone
        do j = i,self%nz
          if (convecting_with_below(j)) then
            k = j + 1
          else
            ! convecting zone stops
            exit
          endif
        enddo
        self%ind_conv_upper = [self%ind_conv_upper, k]
        i = j
        cycle
      endif
      i = i + 1
    enddo

    if (allocated(self%inds_Tx)) deallocate(self%inds_Tx)
    allocate(self%inds_Tx(1))
    self%inds_Tx(1) = 1
    do i = 1,size(convecting_with_below)
      if (convecting_with_below(i)) then
        ! nothing
      else
        self%inds_Tx = [self%inds_Tx, i+1]
      endif
    enddo

    if (allocated(self%ind_conv_lower_x)) deallocate(self%ind_conv_lower_x)
    allocate(self%ind_conv_lower_x(self%n_convecting_zones))
    do i = 1,size(self%ind_conv_lower)
      ind = findloc(self%inds_Tx, self%ind_conv_lower(i), 1)
      if (ind == 0) then
        err = 'Problem setting a convective zone'
        return
      endif
      self%ind_conv_lower_x(i) = ind
    enddo

  end subroutine

  !> Given P_surf and T profile, this routine determines which layers are convecting
  !> and which are purely radiative. To accomplish this, the routine does a small
  !> Newton step for a purely radiative atmosphere. In a layer, if the step tends towards
  !> a T profile stable to convection, then a layer is radiative. If instead, the step
  !> tends to a T profile unstable to convection, then the layer is convective. The
  !> routine specifically updates self%convecting_with_below, self%n_convecting_zones, 
  !> self%ind_conv_lower and self%ind_conv_upper. It also updates all atmosphere varibles.
  subroutine AdiabatClimate_update_convecting_zones(self, P_i_surf, T_in, mode, err)
    use clima_useful, only: linear_solve
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    integer, intent(in) :: mode
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: dFdt(:), dTdt(:), dTdt_dT(:,:), deltaT(:), T_perturb(:), x_in(:)
    real(dp), allocatable :: lapse_rate_perturb(:), difference(:)
    logical, allocatable :: convecting_with_below_save(:)
    logical, allocatable :: convecting_with_below_candidate(:)
    integer :: i, ierr
    integer :: bt
    real(dp) :: alpha
    character(:), allocatable :: err_trial
    logical :: got_perturb
    real(dp) :: thresh_on, thresh_off
    integer :: n_flip_on, n_flip_off, n_zones_prev
    integer :: l, r

    if (all(mode /= [1, 2, 3])) then
      err = 'Invalid mode in AdiabatClimate_update_convecting_zones: expected 1, 2, or 3.'
      return
    endif

    ! work storage
    allocate(deltaT(size(T_in)),T_perturb(size(T_in)))
    allocate(lapse_rate_perturb(self%nz),difference(self%nz))
    allocate(convecting_with_below_candidate(self%nz))

    convecting_with_below_save = self%convecting_with_below
    if (mode /= 3) then
      self%convecting_with_below = .false.
      call AdiabatClimate_set_convecting_zones(self, self%convecting_with_below, err)
      if (allocated(err)) return
    endif

    allocate(x_in(size(self%inds_Tx)))
    do i = 1,size(self%inds_Tx)
      x_in(i) = T_in(self%inds_Tx(i))
    enddo
    allocate(dFdt(size(self%inds_Tx)), dTdt(size(self%inds_Tx)), dTdt_dT(size(self%inds_Tx),size(self%inds_Tx)))

    call AdiabatClimate_objective(self, P_i_surf, x_in, dFdt, dTdt, err)
    if (allocated(err)) return

    if (mode == 1) then

      call AdiabatClimate_jacobian_from_base(self, P_i_surf, x_in, dFdt, dTdt, dTdt_dT, err)
      if (allocated(err)) return

      deltaT = -dTdt
      call linear_solve(dTdt_dT, deltaT, ierr)
      if (ierr /= 0) then
        err = 'Linear solved failed in "update_convecting_zones"'
        return
      endif

      ! Trial Newton step (used only for convective classification). Safeguard the step
      ! size so `make_profile_rc` does not see invalid temperature profiles.
      alpha = min(max(0.0_dp, self%convective_newton_step_size), 1.0_dp)
      got_perturb = .false.
      do bt = 1,20
        T_perturb = deltaT*alpha + T_in

        ! Ensure there are no invalid temperatures
        if (minval(T_perturb) < 1.0_dp) then    
          alpha = 0.5_dp*alpha
          cycle
        endif

        ! Try to get a perturbed T profile.
        call self%make_profile_rc(P_i_surf, T_perturb, err_trial)
        if (.not. allocated(err_trial)) then
          lapse_rate_perturb = self%lapse_rate
          got_perturb = .true.
          exit
        else
          deallocate(err_trial)
          alpha = 0.5_dp*alpha
        endif
        
        if (alpha < 1.0e-8_dp) exit
      enddo

      if (.not. got_perturb) then
        err = 'Failed to update convecting zones.'
        return
      endif

      ! Re-update all variables at T_in, including self%lapse_rate_intended
      call AdiabatClimate_objective(self, P_i_surf, x_in, dFdt, dTdt, err)
      if (allocated(err)) return

      difference = lapse_rate_perturb - self%lapse_rate_intended

      do i = 1,self%nz
        thresh_on = max(self%convective_hysteresis_min, &
                         self%convective_hysteresis_frac_on*abs(self%lapse_rate_intended(i)))
        thresh_off = max(self%convective_hysteresis_min, &
                          self%convective_hysteresis_frac_off*abs(self%lapse_rate_intended(i)))
        if (convecting_with_below_save(i)) then
          if (difference(i) < -thresh_off) then
            self%convecting_with_below(i) = .false.
          else
            self%convecting_with_below(i) = .true.
          endif
        else
          if (difference(i) > thresh_on) then
            self%convecting_with_below(i) = .true.
          else
            self%convecting_with_below(i) = .false.
          endif
        endif
      enddo

      ! Save hysteresis result as candidate mask for boundary limiting / nucleation.
      convecting_with_below_candidate = self%convecting_with_below

      ! Apply limiter to control mask changes (boundary motion and nucleation).
      call AdiabatClimate_apply_convective_mask_limiter(self, convecting_with_below_save, &
          convecting_with_below_candidate, difference, .false.)

    elseif (mode == 2) then

      difference = self%lapse_rate - self%lapse_rate_intended

      self%convecting_with_below = convecting_with_below_save
      do i = 1,self%nz
        if (.not.self%convecting_with_below(i)) then
          thresh_on = max(self%convective_hysteresis_min, &
                           self%convective_hysteresis_frac_on*abs(self%lapse_rate_intended(i)))
          if (difference(i) > thresh_on) then
            self%convecting_with_below(i) = .true.
          endif
        endif
      enddo

      ! Save hysteresis result as candidate mask for boundary limiting / nucleation.
      convecting_with_below_candidate = self%convecting_with_below

      ! Apply limiter to control mask changes (boundary motion and nucleation).
      call AdiabatClimate_apply_convective_mask_limiter(self, convecting_with_below_save, &
          convecting_with_below_candidate, difference, .true.)

    elseif (mode == 3) then

      difference = self%lapse_rate - self%lapse_rate_intended
      where (self%prevent_overconvection_lock > 0)
        self%prevent_overconvection_lock = self%prevent_overconvection_lock - 1
      end where

      i = 1
      do while (i <= self%nz)
        if (self%convecting_with_below(i)) then
          ! We are in a convecting zone
          l = i ! The bottom of a convecting zone
          do while (i <= self%nz)
            if (.not. self%convecting_with_below(i)) exit
            i = i + 1
          enddo
          r = i - 1 ! The top of a convecting zone
          if (r >= self%nz) exit ! guard against bounds error

          thresh_on = max(self%convective_hysteresis_min, &
                       self%convective_hysteresis_frac_on*abs(self%lapse_rate_intended(r+1)))
          thresh_off = max(self%convective_hysteresis_min, &
                        self%convective_hysteresis_frac_off*abs(self%lapse_rate_intended(r+1)))

          if (difference(r+1) > thresh_on) then
            ! Hard invariant for mode 3:
            ! do not allow a superadiabatic first radiative layer above a convective top.
            ! Promote the first radiative link above the current top.
            self%convecting_with_below(r+1) = .true.
            self%prevent_overconvection_lock(r+1) = 2
          elseif (self%lapse_rate(r+1) < -thresh_off) then
            if (self%prevent_overconvection_lock(r) == 0) then
              ! If inversion just above convecting zone, then we shrink the top
              ! of the convecting zone
              self%convecting_with_below(r) = .false.
            endif
          endif
        else
          i = i + 1
        endif
      enddo

    endif

    call AdiabatClimate_set_convecting_zones(self, self%convecting_with_below, err)
    if (allocated(err)) return

    if (self%verbose) then
      n_zones_prev = 0
      i = 1
      do while (i <= self%nz)
        if (convecting_with_below_save(i)) then
          n_zones_prev = n_zones_prev + 1
          do while (i <= self%nz)
            if (.not. convecting_with_below_save(i)) exit
            i = i + 1
          enddo
        else
          i = i + 1
        endif
      enddo
      n_flip_on = count(.not.convecting_with_below_save .and. self%convecting_with_below)
      n_flip_off = count(convecting_with_below_save .and. .not.self%convecting_with_below)
      print"(1x,'Conv mask: +',i0,'  -',i0,'  zones ',i0,'->',i0)", &
            n_flip_on, n_flip_off, n_zones_prev, self%n_convecting_zones
    endif

  end subroutine

  !> Apply boundary-motion and nucleation limits to the convective mask.
  !> If convective_max_boundary_shift < 0, the candidate mask is adopted directly.
  !> If convective_max_boundary_shift == 0, the previous mask is retained.
  !> When no_convection_to_radiation is true, only radiative->convective growth is allowed.
  subroutine AdiabatClimate_apply_convective_mask_limiter(self, convecting_with_below_save, &
      convecting_with_below_candidate, difference, no_convection_to_radiation)
    class(AdiabatClimate), intent(inout) :: self
    logical, intent(in) :: convecting_with_below_save(:)
    logical, intent(in) :: convecting_with_below_candidate(:)
    real(dp), intent(in) :: difference(:)
    logical, intent(in) :: no_convection_to_radiation
    integer :: i, l, r, length, shift

    if (self%convective_max_boundary_shift < 0) then
      ! No limiter: adopt candidate mask directly.
      self%convecting_with_below = convecting_with_below_candidate
      return
    endif

    ! Start from previous mask and only allow limited boundary motion.
    self%convecting_with_below = convecting_with_below_save
    shift = self%convective_max_boundary_shift

    if (shift == 0) return

    i = 1
    do while (i <= self%nz)
      if (convecting_with_below_save(i)) then
        ! existing convective zone in previous mask
        l = i ! lower edge of convective zone
        do while (i <= self%nz)
          if (.not. convecting_with_below_save(i)) exit
          i = i + 1
        enddo
        r = i - 1 ! upper edge of convective zone

        ! candidate expansion limited by shift
        if (convecting_with_below_candidate(l)) then
          ! lower boundary can move down by <= shift
          if (l - shift >= 1) then
            if (any(convecting_with_below_candidate(l-shift:l-1))) then
              self%convecting_with_below(l-shift:l-1) = .true.
            endif
          endif
        endif
        if (convecting_with_below_candidate(r)) then
          ! upper boundary can move up by <= shift
          if (r + shift <= self%nz) then
            if (any(convecting_with_below_candidate(r+1:r+shift))) then
              self%convecting_with_below(r+1:r+shift) = .true.
            endif
          endif
        endif

        ! allow contraction inside the zone based on candidate (only within shift)
        if (.not.no_convection_to_radiation) then
          if (shift < (r-l+1)) then
            if (.not.any(convecting_with_below_candidate(l:l+shift-1))) then
              self%convecting_with_below(l:l+shift-1) = .false.
            endif
            if (.not.any(convecting_with_below_candidate(r-shift+1:r))) then
              self%convecting_with_below(r-shift+1:r) = .false.
            endif
          endif
        endif
      else
        i = i + 1
      endif
    enddo

    ! Allow new convective islands if instability is strong enough.
    i = 1
    do while (i <= self%nz)
      if (.not.convecting_with_below_save(i) .and. convecting_with_below_candidate(i)) then
        l = i
        do while (i <= self%nz)
          if (convecting_with_below_candidate(i) .and. .not.convecting_with_below_save(i)) then
            i = i + 1
          else
            exit
          endif
        enddo
        r = i - 1
        length = r - l + 1

        ! require sufficiently strong instability within the candidate island
        if (maxval(difference(l:r)) > max(self%convective_hysteresis_min, &
            self%convective_hysteresis_frac_on*maxval(abs(self%lapse_rate_intended(l:r))))) then
          ! cap nucleation size to 2*shift per iteration
          self%convecting_with_below(l:min(r, l+2*shift-1)) = .true.
        endif
      else
        i = i + 1
      endif
    enddo

  end subroutine

  subroutine AdiabatClimate_residuals_with_convection(self, f_total, lapse_rate, lapse_rate_intended, dFdt, dTdt, err)
    use clima_eqns, only: heat_capacity_eval
    use clima_const, only: k_boltz, N_avo
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: f_total(:) !! fluxes at the edges of layers
    real(dp), intent(in) :: lapse_rate(:)
    real(dp), intent(in) :: lapse_rate_intended(:)
    real(dp), intent(out) :: dFdt(:), dTdt(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: fluxes(:), mubar(:), cp(:), rho(:), density(:)
    real(dp) :: f_lower, f_upper, cp_tmp, c_eff, c_surface
    integer :: i, j, ind_lower, ind_upper, ind_zone, k_lower, k_upper
    logical :: found

    if (size(dFdt) /= size(self%inds_Tx) .or. size(dTdt) /= size(self%inds_Tx)) then
      err = 'dFdt/dTdt has the wrong shape in residuals_with_convection'
      return
    endif

    ! work storage
    allocate(fluxes(self%nz+1), mubar(self%nz), cp(self%nz), rho(self%nz), density(self%nz))

    ! Radiative energy going into each layer (ergs/(cm^2*s))
    fluxes(1) = f_total(1)
    do i = 2,self%nz+1
      fluxes(i) = (f_total(i) - f_total(i-1))
    enddo

    ! mean molecular weight (g/mol)
    do j = 1,self%nz
      mubar(j) = 0.0_dp
      do i = 1,self%sp%ng
        mubar(j) = mubar(j) + self%f_i(j,i)*self%sp%g(i)%mass
      enddo
    enddo

    ! molecules/cm^3
    density = self%P/(k_boltz*self%T)
    ! [molecules/cm^3]*[mol/molecules]*[g/mol] = [g/cm^3]
    rho = density*(1.0_dp/N_avo)*mubar 

    ! Heat capacity in erg/(g*K)
    do j = 1,self%nz
      cp(j) = 0.0_dp
      do i = 1,self%sp%ng
        call heat_capacity_eval(self%sp%g(i)%thermo, self%T(j), found, cp_tmp) ! J/(mol*K)
        if (.not. found) then
          err = "not found"
          return
        endif
        ! J/(mol*K)
        cp(j) = cp(j) + cp_tmp*self%f_i(j,i) ! J/(mol*K)
      enddo
      ! J/(mol*K) * (mol/kg) = J/(kg*K)
      cp(j) = cp(j)*(1.0_dp/(mubar(j)*1.0e-3_dp))
      ! J/(kg*K) * (erg/J) * (kg/g) = erg/(g*K)
      cp(j) = cp(j)*1.0e4_dp
    enddo

    ! Default radiative residual in [erg/(cm^2*s)] for each active temperature DOF.
    do i = 1,size(self%inds_Tx)
      dFdt(i) = fluxes(self%inds_Tx(i))
    enddo

    ! Change residual to account for convection
    do i = 1,self%n_convecting_zones

      ind_lower = self%ind_conv_lower(i)
      if (ind_lower == 1) then
        f_lower = 0.0_dp
      else
        f_lower = f_total(ind_lower-1)
      endif

      ind_upper = self%ind_conv_upper(i)
      if (ind_lower == 1) then
        f_upper = f_total(ind_upper) + self%surface_heat_flow
      else
        f_upper = f_total(ind_upper)
      endif

      ! Net radiative energy going into the convective zone [erg/(cm^2*s)].
      dFdt(self%ind_conv_lower_x(i)) = f_upper - f_lower

    enddo

    ! Convert residual from [erg/(cm^2*s)] to [K/s] using the effective
    ! areal heat capacity associated with each active temperature DOF.
    c_surface = rho(1)*cp(1)*self%dz(1)
    do i = 1,size(self%inds_Tx)
      ind_zone = findloc(self%ind_conv_lower_x, i, 1)
      if (ind_zone > 0) then
        ind_lower = self%ind_conv_lower(ind_zone)
        ind_upper = self%ind_conv_upper(ind_zone)

        k_lower = max(1, ind_lower-1)
        k_upper = ind_upper - 1
        c_eff = sum(rho(k_lower:k_upper)*cp(k_lower:k_upper)*self%dz(k_lower:k_upper))

        ! Include a surface thermal mass term for surface-connected zones.
        if (ind_lower == 1) c_eff = c_eff + c_surface
      else
        if (self%inds_Tx(i) == 1) then
          ! Surface DOF: assume same rho, cp, dz as the first atmospheric layer.
          c_eff = c_surface
        else
          j = self%inds_Tx(i) - 1
          c_eff = rho(j)*cp(j)*self%dz(j)
        endif
      endif

      dTdt(i) = dFdt(i)/max(c_eff, tiny(1.0_dp))
    enddo

  end subroutine

end submodule
