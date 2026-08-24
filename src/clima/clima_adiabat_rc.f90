module clima_adiabat_rc
  use clima_const, only: dp
  use clima_types, only: Species
  use dop853_module, only: dop853_class
  use futils, only: brent_class
  use clima_eqns, only: ocean_solubility_fcn
  use clima_adiabat_general, only: OceanFunction, DrySpeciesType, CondensingSpeciesType, CustomDrySpeciesType
  use linear_interpolation_module, only: linear_interp_1d
  implicit none

  type :: AdiabatRCProfileData
    ! Input data
    real(dp), pointer :: T_surf
    type(linear_interp_1d), pointer :: mix_custom(:)
    logical, pointer :: convecting_with_below(:)
    type(Species), pointer :: sp
    integer, pointer :: nz
    real(dp), pointer :: planet_mass
    real(dp), pointer :: planet_radius
    real(dp), pointer :: P_top
    real(dp), pointer :: RH(:)
    real(dp), pointer :: rtol
    real(dp), pointer :: atol
    ! Ouput
    real(dp), pointer :: P(:)
    real(dp), pointer :: z(:)
    real(dp), pointer :: f_i(:,:)
    real(dp), pointer :: lapse_rate(:)
    logical, pointer :: super_saturated(:)

    ! For interpolating the input T profile.
    type(linear_interp_1d) :: T
    real(dp), allocatable :: P_interp(:)
    real(dp), allocatable :: T_interp(:)
    real(dp), allocatable :: T_layers(:) ! size nz+1 (includes ground)

    real(dp) :: P_surf !! surface pressure from inputs

    !> indicates whether species is dry or condensing (length ng)
    integer, allocatable :: sp_type(:)
    integer :: stopping_reason !! Reason why we stop integration
    logical :: in_convecting_region
    !> Index of root
    integer :: ind_root

    !> Work space for root finding. All length (ng)
    real(dp), allocatable :: gout(:)
    real(dp), allocatable :: gout_old(:)
    real(dp), allocatable :: gout_tmp(:)
    logical, allocatable :: root_found(:)
    real(dp), allocatable :: P_roots(:)
    !> When we find a root, and exit integration, these
    !> guys will give use the root
    real(dp) :: P_root  = huge(1.0_dp)
    real(dp), allocatable :: u_root(:)
    real(dp) :: epsj = 1.0e-8

    !> Keeps track of how the dry atmosphere is partitioned
    !> between different species
    real(dp), allocatable :: f_i_dry(:)

    !> index for saving z, mixing ratios, and lapse rate
    integer :: j

    ! work variables. All dimension (ng)
    real(dp), allocatable :: P_i_cur(:)
    real(dp), allocatable :: f_i_cur(:)
    real(dp), allocatable :: cp_i_cur(:)
    real(dp), allocatable :: L_i_cur(:)

    !> This helps us propogate error messages outside of the integration
    character(:), allocatable :: err

  end type

  ! stopping_reason
  enum, bind(c)
    enumerator :: ReachedPtop, ReachedRoot
  end enum

  type, extends(dop853_class) :: dop853_rc
    type(AdiabatRCProfileData), pointer :: d => NULL()
  end type

contains

  ! inputs: P_i_surf, T_surf, T in each grid cell, P_top, planet_mass, planet_radius, nz, sp, ocean_fcns
  !         indexes where convection is occuring

  ! outputs: Pressure, mixing ratios, z, lapse rate everywhere

  subroutine make_profile_rc(T_surf, T, P_i_surf, &
                             sp_custom, mix_custom, &
                             convecting_with_below, &
                             sp, nz, planet_mass, &
                             planet_radius, P_top, RH, &
                             rtol, atol, &
                             ocean_fcns, args_p, &
                             P, z, f_i, lapse_rate, super_saturated, &
                             N_surface, N_ocean, &
                             err)
    use futils, only: linspace
    use clima_eqns, only: gravity, ocean_solubility_fcn
    use iso_c_binding, only: c_ptr

    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(inout) :: T(:) !! (nz)
    real(dp), intent(in) :: P_i_surf(:) !! (ng) dynes/cm2
    logical, intent(in) :: sp_custom(:)
    type(linear_interp_1d), target, intent(inout) :: mix_custom(:)
    logical, target, intent(in) :: convecting_with_below(:)

    type(Species), target, intent(inout) :: sp
    integer, target, intent(in) :: nz
    real(dp), target, intent(in) :: planet_mass, planet_radius
    real(dp), target, intent(in) :: P_top, RH(:)
    real(dp), target, intent(in) :: rtol, atol
    type(OceanFunction), intent(in) :: ocean_fcns(:)
    type(c_ptr), value, intent(in) :: args_p

    real(dp), target, intent(out) :: P(:), z(:)
    real(dp), target, intent(out) :: f_i(:,:)
    real(dp), target, intent(out) :: lapse_rate(:)
    logical, target, intent(out) :: super_saturated(:)
    real(dp), target, intent(out) :: N_surface(:), N_ocean(:,:)
    character(:), allocatable, intent(out) :: err

    type(AdiabatRCProfileData) :: d
    integer :: i, j
    real(dp) :: P_sat, grav, P_surface_inventory, f_tmp, P_custom

    ! check inputs
    if (size(T) /= nz) then
      err = 'make_profile: Input "T" has the wrong shape'
      return
    endif
    if (any(P_i_surf < 0.0_dp)) then
      err = 'make_profile: Surface pressures (input "P_i_surf") must be positive'
      return
    endif
    if (size(P_i_surf) /= sp%ng) then
      err = 'make_profile: Input "P_i_surf" has the wrong shape'
      return
    endif
    if (size(convecting_with_below) /= nz) then
      err = 'make_profile: Input "convecting_with_below" has the wrong shape'
      return
    endif
    if (size(RH) /= sp%ng) then
      err = 'make_profile: Input "RH" has the wrong dimension.'
      return
    endif
    if (size(P) /= 2*nz+1) then
      err = 'make_profile: Input "P" has the wrong shape'
      return
    endif
    if (size(z) /= 2*nz+1) then
      err = 'make_profile: Input "z" has the wrong shape'
      return
    endif
    if (size(f_i, 1) /= 2*nz+1 .or. size(f_i, 2) /= sp%ng) then
      err = 'make_profile: Input "f_i" has the wrong shape'
      return
    endif
    if (size(lapse_rate) /= 2*nz+1) then
      err = 'make_profile: Input "lapse_rate" has the wrong shape'
      return
    endif
    if (size(super_saturated) /= 2*nz+1) then
      err = 'make_profile: Input "lapse_rate" has the wrong shape'
      return
    endif
    if (size(N_surface) /= sp%ng) then
      err = 'make_profile: Input "N_surface" has the wrong dimension.'
      return
    endif
    if (size(N_ocean,1) /= sp%ng .or. size(N_ocean,2) /= sp%ng) then
      err = 'make_profile: Input "N_ocean" has the wrong dimension.'
      return
    endif

    ! associate
    ! inputs
    d%T_surf => T_surf
    d%mix_custom => mix_custom
    d%convecting_with_below => convecting_with_below
    d%sp => sp
    d%nz => nz
    d%planet_mass => planet_mass
    d%planet_radius => planet_radius
    d%P_top => P_top
    d%RH => RH
    d%rtol => rtol
    d%atol => atol
    ! outputs
    d%P => P
    d%z => z
    d%f_i => f_i
    d%lapse_rate => lapse_rate
    d%super_saturated => super_saturated

    ! allocate work memory
    allocate(d%sp_type(sp%ng))
    allocate(d%gout(sp%ng))
    allocate(d%gout_old(sp%ng))
    allocate(d%gout_tmp(sp%ng))
    allocate(d%P_roots(sp%ng))
    allocate(d%root_found(sp%ng))
    allocate(d%f_i_dry(sp%ng))
    allocate(d%P_i_cur(sp%ng))
    allocate(d%f_i_cur(sp%ng))
    allocate(d%cp_i_cur(sp%ng))
    allocate(d%L_i_cur(sp%ng))

    ! gravity at the surface of planet
    grav = gravity(planet_radius, planet_mass, 0.0_dp)

    d%P_surf = 0.0_dp
    P_custom = 0.0_dp
    do i = 1,d%sp%ng

      if (sp_custom(i)) then
        ! The total custom pressure.
        P_custom = P_custom + P_i_surf(i)
        cycle
      endif

      ! compute saturation vapor pressures
      P_sat = huge(1.0_dp)
      if (allocated(d%sp%g(i)%sat)) then
        P_sat = d%RH(i)*d%sp%g(i)%sat%sat_pressure(T_surf)
      endif

      ! determine if species are condensing, or not
      if (P_i_surf(i) > P_sat) then
        d%P_i_cur(i) = P_sat
        ! the pressure of the ocean is everything not in the atmosphere
        P_surface_inventory = P_i_surf(i) - P_sat 
        ! The surface mol/cm^2 can then be computed
        ! from the "surface pressure"
        N_surface(i) = P_surface_inventory/(sp%g(i)%mass*grav)
        d%sp_type(i) = CondensingSpeciesType
      else
        d%P_i_cur(i) = P_i_surf(i)
        N_surface(i) = 0.0_dp ! no surface reservoir if not condensing at the surface
        d%sp_type(i) = DrySpeciesType
      endif

      d%P_surf = d%P_surf + d%P_i_cur(i)
    enddo

    d%P_surf = d%P_surf + P_custom
    if (P_top > d%P_surf) then
      err = 'make_profile: "P_top" is bigger than the surface pressure'
      return
    endif
    do i = 1,d%sp%ng
      if (.not.sp_custom(i)) cycle

      call d%mix_custom(i)%evaluate(log10(d%P_surf), f_tmp)
      d%P_i_cur(i) = P_custom*10.0_dp**f_tmp
      N_surface(i) = 0.0_dp
      d%sp_type(i) = CustomDrySpeciesType
    enddo

    d%f_i_cur(:) = d%P_i_cur(:)/d%P_surf ! mixing ratios
    call update_f_i_dry(d, d%P_surf, d%f_i_cur)

    ! compute gases dissolved in oceans. This allows for multiple ocean.
    do j = 1,sp%ng
      if (associated(ocean_fcns(j)%fcn)) then; block
        real(dp), allocatable :: m_i_cur(:)

        ! Compute mol/kg of each gas dissolved in the ocean
        allocate(m_i_cur(sp%ng))
        call ocean_fcns(j)%fcn(T_surf, sp%ng, d%P_i_cur/1.0e6_dp, m_i_cur, args_p)
        
        ! Compute mol/cm^2 of each gas dissolved in ocean
        N_ocean(j,j) = 0.0_dp ! ocean can not dissolve into itself
        do i = 1,sp%ng
          if (i /= j) then
            N_ocean(i,j) = m_i_cur(i)*N_surface(j)*(sp%g(j)%mass/1.0e3_dp)
          endif
        enddo
      endblock; else
        ! Nothing dissolved in ocean
        N_ocean(:,j) = 0.0_dp
      endif
    enddo

    ! Make P profile
    call linspace(log10(d%P_surf),log10(P_top),P)
    P(:) = 10.0_dp**P(:)
    P(1) = d%P_surf
    P(2*nz+1) = P_top
    if (minval(P(1:2*nz) - P(2:2*nz+1)) <= 10.0_dp*spacing(P(1))) then
      err = 'make_profile: pressure grid spacing too fine (P_top too close to P_surf)'
      return
    endif

    ! If convecting with below, we don't know the temp.
    ! We set it to negative as a place holder
    do i = 1,nz
      if (d%convecting_with_below(i)) then
        T(i) = -1.0_dp
      endif
    enddo

    ! Copy temps to T_layers
    allocate(d%T_layers(nz+1))
    d%T_layers(1) = T_surf
    d%T_layers(2:size(d%T_layers)) = T(:)

    allocate(d%T_interp(nz+1),d%P_interp(nz+1))

    call integrate(d, err)
    T(:) = d%T_layers(2:)
    if (allocated(err)) return

  end subroutine

  subroutine initialize_interpolator(d, err)
    type(AdiabatRCProfileData), target, intent(inout) :: d
    character(:), allocatable, intent(out) :: err

    integer :: i, ierr

    d%T_interp(1) = d%T_layers(1)
    d%P_interp(1) = log10(d%P_surf)
    do i = 1,d%nz
      d%T_interp(i+1) = d%T_layers(i+1)
      d%P_interp(i+1) = log10(d%P(2*i))
    enddo
    d%T_interp = d%T_interp(size(d%T_interp):1:-1)
    d%P_interp = d%P_interp(size(d%P_interp):1:-1)
    call d%T%initialize(d%P_interp, d%T_interp, ierr)
    if (ierr /= 0) then
      err = 'Failed to initialize interpolator in "make_profile_rc"'
      return
    endif

  end subroutine

  subroutine integrate(d, err)
    type(AdiabatRCProfileData), target, intent(inout) :: d
    character(:), allocatable, intent(out) :: err

    type(dop853_rc) :: dop
    logical :: status_ok, convecting_with_below
    integer :: idid, i, j, k, n, kk, j_prev
    real(dp) :: Pn, T_root, Ptop
    real(dp), allocatable :: u(:)
    integer, allocatable :: icomp(:)
    character(6) :: tmp_char

    d%j = 2
    ! Surface values
    d%lapse_rate(1) = general_adiabat_lapse_rate(d, d%T_surf, d%P_surf)
    d%f_i(1,:) = d%f_i_cur(:)
    d%z(1) = 0.0_dp

    do
      ! i is the layer of the atmosphere we are starting from.
      ! i == 1 is the ground.
      if (d%j == 2) then
        i = 1
      else
        i = (d%j-1)/2 + 1
      endif

      ! Check if layer above is convective
      convecting_with_below = d%convecting_with_below(i)

      if (convecting_with_below) then
        d%in_convecting_region = .true.
        n = 2
        icomp = [1, 2]
        u = [d%z(d%j-1), d%T_layers(i)]
        Pn = d%P(d%j-1)

        ! Find the pressure to integrate to
        do j = i,d%nz
          if (d%convecting_with_below(j)) then
            k = j*2
          else
            exit
          endif
        enddo

        Ptop = d%P(k)
        
      else
        d%in_convecting_region = .false.
        n = 1
        icomp = [1]
        u = [d%z(d%j-1)]
        Pn = d%P(d%j-1)
        ! Interpolator for T
        call initialize_interpolator(d, err)
        if (allocated(err)) return

        ! Find the pressure to integrate to
        do j = i,d%nz
          if (.not.d%convecting_with_below(j)) then
            k = j*2
          else
            exit
          endif
        enddo

        Ptop = d%P(k)

      endif

      ! Some logic to identify cold traps at the surface.
      if (d%j == 2 .and. .not.d%in_convecting_region) then
        call d%T%evaluate(log10(Pn), T_root)
        call root_fcn(d, Pn, T_root, d%gout)
        do i = 1,d%sp%ng
          if (d%sp_type(i) == CondensingSpeciesType .and. d%gout(i) <= 0.0_dp) then
            ! If d%gout(i) is <= 0, this means a moist species is increasing with decreasing pressure
            ! so, we must switch it to a dry species.
            d%sp_type(i) = DrySpeciesType
          endif
        enddo
        call update_f_i_dry(d, Pn, d%f_i(1,:))
      endif

      j_prev = d%j
      call dop%initialize(fcn=right_hand_side_dop, solout=solout_dop, n=n, &
                        iprint=0, icomp=icomp, status_ok=status_ok)
      dop%d => d
      d%stopping_reason = ReachedPtop
      if (.not. status_ok) then
        err = 'dop853 initialization failed'
        return
      endif
      
      kk = 1
      do
        call dop%integrate(Pn, u, Ptop, [d%rtol], [d%atol], iout=2, idid=idid)
        if (allocated(d%err)) then
          err = d%err
          return
        endif
        if (idid < 0) then
          write(tmp_char,'(i6)') idid
          err = 'dop853 integration failed: '//trim(tmp_char)
          return
        endif 

        if (d%stopping_reason == ReachedPtop) then
          exit
        elseif (d%stopping_reason == ReachedRoot) then; block
          real(dp) :: f_dry
          logical :: super_saturated
          logical :: switched_other
          ! A root was hit. We restart integration
          Pn = d%P_root ! root pressure
          u = d%u_root ! root altitude

          ! Get T at the root
          if (d%in_convecting_region) then
            T_root = dop%contd8(2, Pn)
          else
            call d%T%evaluate(log10(Pn), T_root)
          endif
          
          call mixing_ratios(d, Pn, T_root, d%f_i_cur, f_dry, super_saturated) ! get mixing ratios at the root

          if (d%sp_type(d%ind_root) == DrySpeciesType) then
            ! A gas has reached saturation.
            d%sp_type(d%ind_root) = CondensingSpeciesType
          elseif (d%sp_type(d%ind_root) == CondensingSpeciesType) then
            ! A gas has reached a cold trap
            d%sp_type(d%ind_root) = DrySpeciesType
          endif
          call update_f_i_dry(d, Pn, d%f_i_cur)

          ! Check for other roots detected in this step that may be missed after switching.
          call root_fcn(d, Pn, T_root, d%gout)
          switched_other = .false.
          do i = 1,d%sp%ng
            if (.not.d%root_found(i)) cycle
            if (i == d%ind_root) cycle

            if (d%sp_type(i) == CondensingSpeciesType) then
              ! If mixing ratio would increase with decreasing pressure, force dry.
              if (d%gout(i) < 0.0_dp) then
                d%sp_type(i) = DrySpeciesType
                switched_other = .true.
              endif
            elseif (d%sp_type(i) == DrySpeciesType) then
              ! If supersaturated, force condensing.
              if (d%gout(i) > 0.0_dp) then
                d%sp_type(i) = CondensingSpeciesType
                switched_other = .true.
              endif
            endif
          enddo
          if (switched_other) then
            call update_f_i_dry(d, Pn, d%f_i_cur)
          endif

          d%root_found(:) = .false.
          d%P_roots(:) = 0.0_dp
          d%stopping_reason = ReachedPtop

          ! Ensure we make progress away from the root. It is possible (especially with noisy
          ! prescribed radiative T profiles) to repeatedly detect roots extremely close together,
          ! which can lead to an excessive number of restarts. Nudge the restart pressure toward
          ! Ptop by a few ulps and recompute the state from the previous segment's dense output.
          if (Pn > Ptop) then
            Pn = Pn - 3.0_dp*spacing(Pn)
            if (d%j <= size(d%P)) then
              Pn = max(Pn, Ptop + 3.0_dp*spacing(Ptop), d%P(d%j) + 3.0_dp*spacing(d%P(d%j)))
            else
              Pn = max(Pn, Ptop + 3.0_dp*spacing(Ptop))
            endif
            if (d%in_convecting_region) then
              u = [dop%contd8(1, Pn), dop%contd8(2, Pn)]
            else
              u = [dop%contd8(1, Pn)]
            endif
          endif

        endblock; endif

        if (kk > 100000) then
          err = 'Terrible integration error in clima_adiabat_rc. Report this bug!'
          return
        endif

        kk = kk + 1
      enddo

      if (d%j == j_prev) then
        err = 'Integration stalled in clima_adiabat_rc: no progress in j'
        return
      endif

      if (k > 2*d%nz) then
        print*,'Terrible out-of-bounds error in clima_adiabat_rc. Report this bug!'
        stop 1
      endif

      if (k == 2*d%nz) exit

    enddo

  end subroutine

  subroutine solout_dop(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout
    
    type(AdiabatRCProfileData), pointer :: d
    real(dp) :: z_cur, P_cur, T_cur, P_old
    real(dp) :: T_i, PP, f_dry, P_tol
    logical :: super_saturated
    integer :: i
    
    P_old = xold
    P_cur = x
    z_cur = y(1)

    select type (self)
    class is (dop853_rc)
      d => self%d
    end select

    if (d%in_convecting_region) then
      T_cur = y(2)
    else
      call d%T%evaluate(log10(P_cur), T_cur)
    endif

    if (allocated(d%err)) then
      ! if there is an error (from RHS function)
      ! we return immediately
      irtrn = -10
      return
    endif

    call root_fcn(d, P_cur, T_cur, d%gout)
    if (nr == 1) then
      d%gout_old(:) = d%gout(:)
    endif
    do i = 1,size(d%gout_old)
      if ((d%gout_old(i) < 0 .and. d%gout(i) > 0) .or. &
          (d%gout_old(i) > 0 .and. d%gout(i) < 0)) then
        d%root_found(i) = .true. 
        call find_root(self, d, P_old, P_cur, i, d%P_roots(i))
        if (allocated(d%err)) then
          irtrn = -10
          return
        endif
      else
        d%root_found(i) = .false. 
        d%P_roots(i) = 0.0_dp
      endif
    enddo

    if (any(d%root_found)) then
      ! lets use the highest pressure root.
      ! this makes sure we don't skip an interesting root
      irtrn = -1
      d%ind_root = maxloc(d%P_roots,1) ! highest pressure root
      d%P_root = d%P_roots(d%ind_root)
      if (d%in_convecting_region) then
        d%u_root = [self%contd8(1, d%P_root), self%contd8(2, d%P_root)]
      else
        d%u_root = [self%contd8(1, d%P_root)]
      endif
      d%stopping_reason = ReachedRoot
    endif

    d%gout_old(:) = d%gout(:) ! save for next step

    ! save the results
    if (d%j <= size(d%P)) then
      if (irtrn == -1) then
        PP = d%P_root
      else
        PP = P_cur
      endif

      ! Tolerance for P comparisons: use a relative scale plus a machine-precision guard.
      P_tol = max(1.0e-10_dp*max(abs(P_old), abs(PP), abs(d%P(d%j))), &
                  10.0_dp*spacing(max(P_old, PP)))
      if (d%P(d%j) <= P_old .and. d%P(d%j) >= PP - P_tol) then
        do while (d%P(d%j) >= PP - P_tol)
          d%z(d%j) = self%contd8(1, d%P(d%j))
          if (d%in_convecting_region) then
            T_i = self%contd8(2, d%P(d%j))
            if (mod(d%j,2) == 0) then
              d%T_layers(d%j/2+1) = T_i
            endif
          else
            call d%T%evaluate(log10(d%P(d%j)), T_i)
          endif
          call mixing_ratios(d, d%P(d%j), T_i, d%f_i(d%j,:), f_dry, super_saturated)
          d%lapse_rate(d%j) = general_adiabat_lapse_rate(d, T_i, d%P(d%j))
          d%super_saturated(d%j) = super_saturated
          d%j = d%j + 1
          if (d%j > size(d%P)) exit
        enddo
      endif

    endif

  end subroutine

  subroutine find_root(dop, d, P_old, P_cur, ind, P_root)
    use futils, only: brent_class
    type(dop853_class), intent(inout) :: dop
    type(AdiabatRCProfileData), intent(inout) :: d
    real(dp), intent(in) :: P_old
    real(dp), intent(in) :: P_cur
    integer, intent(in) :: ind
    real(dp), intent(out) :: P_root

    real(dp), parameter :: tol = 1.0e-8_dp
    real(dp) :: fzero
    integer :: info

    type(brent_class) :: b

    call b%set_function(fcn)
    call b%find_zero(P_old, P_cur, tol, P_root, fzero, info)
    if (info /= 0) then
      d%err = 'brent failed in "find_root"'
      return
    endif

  contains
    function fcn(me_, x_) result(f_)
      class(brent_class), intent(inout) :: me_
      real(dp), intent(in) :: x_
      real(dp) :: f_
      real(dp) :: T_
      if (d%in_convecting_region) then
        T_ = dop%contd8(2, x_)
      else
        call d%T%evaluate(log10(x_), T_)
      endif
      call root_fcn(d, x_, T_, d%gout_tmp)
      f_ = d%gout_tmp(ind)
    end function
  end subroutine

  subroutine root_fcn(d, P, T, gout)
    type(AdiabatRCProfileData), intent(inout) :: d
    real(dp), intent(in) :: P, T
    real(dp), intent(out) :: gout(:) 

    real(dp) :: f_dry, P_sat, dTdlog10P
    logical :: super_saturated
    integer :: i

    call mixing_ratios(d, P, T, d%f_i_cur, f_dry, super_saturated)
    d%P_i_cur = d%f_i_cur*P

    if (any(d%sp_type == CondensingSpeciesType) .and. .not.d%in_convecting_region) then
      call d%T%evaluate_derivative(log10(P), dTdlog10P)
    endif

    gout(:) = 1.0_dp

    do i = 1,d%sp%ng
      ! compute saturation vapor pressures
      P_sat = huge(1.0_dp)
      if (allocated(d%sp%g(i)%sat)) then
        P_sat = d%RH(i)*d%sp%g(i)%sat%sat_pressure(T)
      endif

      if (d%sp_type(i) == CondensingSpeciesType .and. .not.d%in_convecting_region) then; block
        real(dp) :: dPi_dT, dTdP, dPi_dP, dfi_dP, dlog10fi_dP

        ! Derivative of SVP curve for species i
        dPi_dT = d%RH(i)*d%sp%g(i)%sat%sat_pressure_derivative(T)

        ! Derivative of temperature profile
        dTdP = dTdlog10P*(1.0_dp/(P*log(10.0_dp)))

        ! Get dP_i/dP
        dPi_dP = dPi_dT*dTdP
        
        ! Convert to df_i/dP, where f_i is mixing ratio of condensible
        dfi_dP = (1.0_dp/P)*(dPi_dP) - P_sat/P**2.0_dp

        ! Convert to the derivative of the log of the mixing ratio
        dlog10fi_dP = dfi_dP*(1/(d%f_i_cur(i)*log(10.0_dp)))

        ! The root occurs when the mixing ratio begins to increase in concentration
        gout(i) = dlog10fi_dP - 1.0e-8_dp

      end block; elseif (d%sp_type(i) == DrySpeciesType) then
        ! Dry species can become condensing species if
        ! they reach saturation.
        gout(i) = d%P_i_cur(i)/P_sat - (1.0_dp + 1.0e-8_dp)
      endif

    enddo 

  end subroutine

  subroutine right_hand_side_dop(self, P, u, du)
    use clima_eqns, only: gravity, heat_capacity_eval
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: P
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du

    select type (self)
    class is (dop853_rc)
      call right_hand_side(self%d, P, u, du)
    end select

  end subroutine

  subroutine update_f_i_dry(d, P, f_i_layer)
    type(AdiabatRCProfileData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(in) :: f_i_layer(:)

    integer :: i
    real(dp) :: P_dry

    d%P_i_cur(:) = f_i_layer(:)*P
    P_dry = 0.0_dp
    do i = 1,d%sp%ng
      if (d%sp_type(i) == DrySpeciesType .or. d%sp_type(i) == CustomDrySpeciesType) then
        P_dry = P_dry + d%P_i_cur(i)
      endif
    enddo
    d%f_i_dry(:) = d%P_i_cur(:)/P_dry

  end subroutine

  subroutine mixing_ratios(d, P, T, f_i_layer, f_dry, super_saturated)
    type(AdiabatRCProfileData), intent(inout) :: d
    real(dp), intent(in) :: P, T
    real(dp), intent(out) :: f_i_layer(:), f_dry
    logical, intent(out) :: super_saturated

    real(dp) :: f_moist, f_dry_tot, f_custom_tot, f_tmp, f_custom(d%sp%ng)
    integer :: i

    super_saturated = .false.

    ! moist mixing ratios
    f_moist = 0.0_dp
    do i = 1,d%sp%ng
      if (d%sp_type(i) == CondensingSpeciesType) then
        f_i_layer(i) = min(d%RH(i)*d%sp%g(i)%sat%sat_pressure(T)/P,1.0_dp)
        f_moist = f_moist + f_i_layer(i)
      endif
    enddo

    ! fraction of the atmosphere that is dry
    f_dry = max(1.0_dp - f_moist, 1.0e-40_dp)
    ! mixing ratios of dry species
    f_dry_tot = 0.0_dp ! Total DrySpeciesType dry mixing ratio
    f_custom_tot = 0.0_dp ! Total mixing ratio from interpolation
    f_custom = 0.0_dp ! Save all the mixing ratios from interpolation
    do i = 1,d%sp%ng
      if (d%sp_type(i) == DrySpeciesType) then
        f_dry_tot = f_dry_tot + d%f_i_dry(i)
        f_i_layer(i) = f_dry*d%f_i_dry(i)
      elseif (d%sp_type(i) == CustomDrySpeciesType) then
        call d%mix_custom(i)%evaluate(log10(P), f_tmp)
        f_tmp = 10.0_dp**f_tmp
        f_custom_tot = f_custom_tot + f_tmp
        f_custom(i) = f_tmp
      endif
    enddo
    f_custom = f_custom/max(f_custom_tot, tiny(1.0_dp)) ! Normalize to 1
    do i = 1,d%sp%ng
      if (d%sp_type(i) == CustomDrySpeciesType) then
        ! (1.0_dp - f_dry_tot) is the fraction of the dry atmosphere yet not filled.
        ! So we will distribution custom interpolated gases into this. Then to get true mixing
        ! ratio, multiply by f_dry
        f_i_layer(i) = f_custom(i)*(1.0_dp - f_dry_tot)*f_dry
      endif
    enddo

  end subroutine

  function general_adiabat_lapse_rate(d, T, P) result(dlnT_dlnP)
    use clima_eqns, only: heat_capacity_eval
    use clima_eqns_water, only: Rgas_cgs => Rgas
    type(AdiabatRCProfileData), intent(inout) :: d
    real(dp), intent(in) :: P, T
    real(dp) :: dlnT_dlnP

    real(dp) :: f_dry, cp_dry
    real(dp) :: cp, L
    real(dp) :: first_sumation, second_sumation, Beta_i
    logical :: found, super_saturated
    integer :: i

    real(dp), parameter :: Rgas_si = Rgas_cgs/1.0e7_dp ! ideal gas constant in SI units (J/(mol*K))

    call mixing_ratios(d, P, T, d%f_i_cur, f_dry, super_saturated)

    ! heat capacity and latent heat
    cp_dry = tiny(0.0_dp)
    do i = 1,d%sp%ng
      if (d%sp_type(i) == CondensingSpeciesType) then
        L = d%sp%g(i)%sat%latent_heat(T) ! erg/g
        L = L*d%sp%g(i)%mass*1.0e-7_dp ! convert to J/mol
        d%L_i_cur(i) = L
      endif
      
      ! heat capacity in J/(mol*K)
      call heat_capacity_eval(d%sp%g(i)%thermo, T, found, cp)
      if (.not. found) then
        d%err = "Failed to compute heat capacity"
        return
      endif
      d%cp_i_cur(i) = cp
      if (d%sp_type(i) == DrySpeciesType .or. d%sp_type(i) == CustomDrySpeciesType) then
        cp_dry = cp_dry + d%f_i_dry(i)*cp ! J/(mol*K)
      endif
    enddo

    ! Two sums in Equation 1 of Graham et al. (2021)
    first_sumation = 0.0_dp
    second_sumation = 0.0_dp
    do i = 1,d%sp%ng
      if (d%sp_type(i) == CondensingSpeciesType) then
        Beta_i = d%L_i_cur(i)/(Rgas_si*T)
        first_sumation = first_sumation + &
          d%f_i_cur(i)*(d%cp_i_cur(i) - Rgas_si*Beta_i + Rgas_si*Beta_i**2.0_dp)

        second_sumation = second_sumation + &
          Beta_i*d%f_i_cur(i)
      endif
    enddo

    ! Equation 1 in Graham et al. (2021), except simplified to assume no condensate present.
    ! Units are SI
    dlnT_dlnP = 1.0_dp/(f_dry*((cp_dry*f_dry + first_sumation)/(Rgas_si*(f_dry + second_sumation))) + second_sumation)

  end function

  subroutine right_hand_side(d, P, u, du)
    use clima_eqns, only: heat_capacity_eval, gravity
    use clima_eqns_water, only: Rgas_cgs => Rgas
    type(AdiabatRCProfileData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(in) :: u(:)
    real(dp), intent(out) :: du(:)

    real(dp) :: z
    real(dp) :: T, f_dry, mubar, grav, dz_dP
    logical :: super_saturated
    integer :: i

    ! unpack u
    z = u(1)

    if (d%in_convecting_region) then
      T = u(2)
    else
      call d%T%evaluate(log10(P),T)
    endif
    call mixing_ratios(d, P, T, d%f_i_cur, f_dry, super_saturated)

    mubar = 0.0_dp
    do i = 1,d%sp%ng
      mubar = mubar + d%f_i_cur(i)*d%sp%g(i)%mass
    enddo

    grav = gravity(d%planet_radius, d%planet_mass, z)
    dz_dP = -(Rgas_cgs*T)/(grav*P*mubar)

    du(1) = dz_dP

    if (d%in_convecting_region) then; block
      real(dp) :: dlnT_dlnP, dT_dP
      dlnT_dlnP = general_adiabat_lapse_rate(d, T, P)
      dT_dP = dlnT_dlnP*(T/P)
      du(2) = dT_dP
    endblock; endif

  end subroutine

end module
