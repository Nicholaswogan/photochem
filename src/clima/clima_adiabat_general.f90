module clima_adiabat_general
  use clima_const, only: dp
  use clima_types, only: Species
  use dop853_module, only: dop853_class
  use futils, only: brent_class
  use clima_eqns, only: ocean_solubility_fcn
  implicit none
  private

  public :: make_profile
  public :: OceanFunction
  public :: DrySpeciesType, CondensingSpeciesType, CustomDrySpeciesType

  type :: OceanFunction
    procedure(ocean_solubility_fcn), nopass, pointer :: fcn => null()
  end type

  type :: AdiabatProfileData
    !> Input data
    type(Species), pointer :: sp
    real(dp), pointer :: planet_mass
    real(dp), pointer :: planet_radius
    real(dp), pointer :: P_top
    real(dp), pointer :: T_trop
    real(dp), pointer :: RH(:)
    real(dp), pointer :: T_surf
    real(dp), pointer :: rtol
    real(dp), pointer :: atol
    !> Ouput
    real(dp), pointer :: P(:)
    real(dp), pointer :: z(:)
    real(dp), pointer :: T(:)
    real(dp), pointer :: f_i(:,:)
    real(dp), pointer :: P_trop

    real(dp) :: P_surf !! surface pressure from inputs

    !> indicates whether species is dry or condensing (length ng)
    integer, allocatable :: sp_type(:)
    integer :: stopping_reason !! Reason why we stop integration (see enum below)
    !> Index of gas that reached saturation if stopping_reason == ReachedGasSaturation
    integer :: ind_gas_reached_saturation

    !> Work space for root finding. All length (ng+1)
    real(dp), allocatable :: gout(:)
    real(dp), allocatable :: gout_old(:)
    real(dp), allocatable :: gout_tmp(:)
    logical, allocatable :: root_found(:)
    real(dp), allocatable :: P_roots(:)
    !> When we find a root, and exit integration, these
    !> guys will give use the root
    real(dp) :: P_root 
    real(dp) :: u_root(2)

    !> Keeps track of how the dry atmosphere is partitioned
    !> between different species
    real(dp), allocatable :: f_i_dry(:)

    !> index for saving T, z and mixing ratios
    integer :: j

    !> work variables. All dimension (ng)
    real(dp), allocatable :: P_i_cur(:)
    real(dp), allocatable :: f_i_cur(:)
    real(dp), allocatable :: cp_i_cur(:)
    real(dp), allocatable :: L_i_cur(:)
    
    !> This helps us propogate error messages outside of the integration
    character(:), allocatable :: err
    
  end type

  ! stopping_reason
  enum, bind(c)
    enumerator :: ReachedPtop, ReachedTropopause, ReachedGasSaturation
  end enum

  ! sp_type
  enum, bind(c)
    enumerator :: DrySpeciesType, CondensingSpeciesType, CustomDrySpeciesType
  end enum

  type, extends(dop853_class) :: dop853_custom
    type(AdiabatProfileData), pointer :: d => NULL()
  end type

contains

  subroutine make_profile(T_surf, P_i_surf, &
                          sp, nz, planet_mass, &
                          planet_radius, P_top, T_trop, RH, &
                          rtol, atol, &
                          ocean_fcns, args_p, &
                          P, z, T, f_i, P_trop, &
                          N_surface, N_ocean, &
                          err)
    use futils, only: linspace
    use clima_eqns, only: gravity, ocean_solubility_fcn
    use iso_c_binding, only: c_ptr

    real(dp), target, intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_i_surf(:) !! (ng) dynes/cm2

    type(Species), target, intent(inout) :: sp
    integer, intent(in) :: nz
    real(dp), target, intent(in) :: planet_mass, planet_radius
    real(dp), target, intent(in) :: P_top, T_trop, RH(:)
    real(dp), target, intent(in) :: rtol, atol
    type(OceanFunction), intent(in) :: ocean_fcns(:)
    type(c_ptr), value, intent(in) :: args_p

    real(dp), target, intent(out) :: P(:), z(:), T(:) ! (ng)
    real(dp), target, intent(out) :: f_i(:,:) ! (nz,ng)
    real(dp), target, intent(out) :: P_trop !! Tropopause pressure (dynes/cm^2). 
                                            !! Value is negative if no tropopause is found.
    real(dp), target, intent(out) :: N_surface(:), N_ocean(:,:)
    character(:), allocatable, intent(out) :: err

    type(AdiabatProfileData) :: d
    integer :: i, j
    real(dp) :: P_sat, grav, P_surface_inventory

    ! check inputs
    if (T_surf < T_trop) then
      err = 'make_profile: Input "T_surf" is less than input "T_trop"'
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
    if (T_trop < 0.0_dp) then
      err = 'make_profile: Input "T_trop" is less than 0'
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
    if (size(T) /= 2*nz+1) then
      err = 'make_profile: Input "T" has the wrong shape'
      return
    endif
    if (size(f_i, 1) /= 2*nz+1 .or. size(f_i, 2) /= sp%ng) then
      err = 'make_profile: Input "f_i" has the wrong shape'
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
    d%sp => sp
    d%planet_mass => planet_mass
    d%planet_radius => planet_radius
    d%P_top => P_top
    d%T_trop => T_trop
    d%RH => RH
    d%T_surf => T_surf
    d%rtol => rtol
    d%atol => atol
    ! outputs
    d%P => P
    d%z => z
    d%T => T
    d%f_i => f_i
    d%P_trop => P_trop
    ! allocate
    allocate(d%sp_type(sp%ng))
    allocate(d%gout(sp%ng+1))
    allocate(d%gout_old(sp%ng+1))
    allocate(d%gout_tmp(sp%ng+1))
    allocate(d%P_roots(sp%ng+1))
    allocate(d%root_found(sp%ng+1))
    allocate(d%f_i_dry(sp%ng))
    allocate(d%P_i_cur(sp%ng))
    allocate(d%f_i_cur(sp%ng))
    allocate(d%cp_i_cur(sp%ng))
    allocate(d%L_i_cur(sp%ng))

    ! gravity at the surface of planet
    grav = gravity(planet_radius, planet_mass, 0.0_dp)

    do i = 1,d%sp%ng

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
    enddo

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

    d%P_surf = sum(d%P_i_cur) ! total pressure
    if (P_top > d%P_surf) then
      err = 'make_profile: "P_top" is bigger than the surface pressure'
      return
    endif
    d%f_i_cur(:) = d%P_i_cur(:)/d%P_surf ! mixing ratios
    call update_f_i_dry(d, d%P_surf, d%f_i_cur)

    ! Make P profile
    call linspace(log10(d%P_surf),log10(P_top),P)
    P(:) = 10.0_dp**P(:)
    P(1) = d%P_surf
    P(2*nz+1) = P_top

    ! integrate
    call integrate(d, err)
    if (allocated(err)) return

    if (any(z < 0.0_dp)) then
      err = '"make_profile" yielded negative altitudes. '// &
            'This may be caused by the lack of a hydrostatic solution to the entered atmosphere.'
      return
    endif

  end subroutine

  subroutine integrate(d, err)
    type(AdiabatProfileData), target, intent(inout) :: d
    character(:), allocatable, intent(out) :: err

    type(dop853_custom) :: dop
    logical :: status_ok
    integer :: i, idid
    real(dp) :: Pn, u(2)
    character(6) :: tmp_char

    call dop%initialize(fcn=right_hand_side_dop, solout=solout_dop, n=2, &
                        iprint=0, icomp=[1,2], status_ok=status_ok)
    dop%d => d
    d%j = 2
    d%f_i(1,:) = d%f_i_cur(:)
    d%T(1) = d%T_surf
    d%z(1) = 0.0_dp
    d%stopping_reason = ReachedPtop
    if (.not. status_ok) then
      err = 'dop853 initialization failed'
      return
    endif

    ! intial conditions
    u = [d%T_surf, 0.0_dp]
    Pn = d%P_surf
    do
      call dop%integrate(Pn, u, d%P_top, [d%rtol], [d%atol], iout=2, idid=idid)
      if (allocated(d%err)) then
        err = d%err
        return
      endif
      if (idid < 0) then
        write(tmp_char,'(i6)') idid
        err = 'dop853 integration failed: '//trim(tmp_char)
        return
      endif 

      if (d%stopping_reason == ReachedPtop .or. &
          d%stopping_reason == ReachedTropopause) then
        exit
      elseif (d%stopping_reason == ReachedGasSaturation) then; block
        real(dp) :: f_dry

        Pn = d%P_root
        u = d%u_root
        call mixing_ratios(d, Pn, u(1), d%f_i_cur, f_dry)
        d%sp_type(d%ind_gas_reached_saturation) = CondensingSpeciesType
        call update_f_i_dry(d, Pn, d%f_i_cur)
        d%stopping_reason = ReachedPtop

      endblock; endif

    enddo

    if (d%stopping_reason == ReachedPtop) then
      d%P_trop = -1.0_dp ! no tropopause identified
    elseif (d%stopping_reason == ReachedTropopause) then; block
      real(dp) :: mubar, f_dry

      call mixing_ratios(d, d%P_root, d%u_root(1), d%f_i_cur, f_dry)

      ! Propopause pressure
      d%P_trop = d%P_root

      ! Values above the tropopause
      mubar = 0.0_dp
      do i = 1,d%sp%ng
        mubar = mubar + d%f_i_cur(i)*d%sp%g(i)%mass
      enddo
      do i = d%j,size(d%z)
        d%T(i) = d%T_trop
        d%f_i(i,:) = d%f_i_cur(:)
        d%z(i) = altitude_vs_P(d%P(i), d%T_trop, mubar, d%P_root, &
                 d%u_root(2), d%planet_mass, d%planet_radius) 
      enddo

    endblock; endif

  end subroutine

  subroutine solout_dop(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout
    
    type(AdiabatProfileData), pointer :: d
    real(dp) :: T_cur, z_cur, P_cur, P_old
    real(dp) :: PP, f_dry
    integer :: i
    
    P_old = xold
    P_cur = x
    T_cur = y(1)
    z_cur = y(2)

    select type (self)
    class is (dop853_custom)
      d => self%d
    end select

    if (allocated(d%err)) then
      ! if there is an error (from RHS function)
      ! we return immediately
      irtrn = -10
      return
    endif

    ! Look for roots
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

    if (any(d%root_found)) then; block
      integer :: ind_root
      ! lets use the highest pressure root.
      ! this makes sure we don't skip an interesting root
      irtrn = -1
      ind_root = maxloc(d%P_roots,1)
      d%P_root = d%P_roots(ind_root)
      d%u_root(1) = self%contd8(1, d%P_root)
      d%u_root(2) = self%contd8(2, d%P_root)

      if (ind_root == 1) then
        ! Tropopause
        d%stopping_reason = ReachedTropopause
      elseif (ind_root >= 2 .and. ind_root <= d%sp%ng+1) then
        ! Some species is supersaturated
        d%stopping_reason = ReachedGasSaturation
        d%ind_gas_reached_saturation = ind_root - 1
      endif
      
    endblock; endif

    d%gout_old(:) = d%gout(:) ! save for next step

    ! save the results
    if (d%j <= size(d%P)) then
      if (irtrn == -1) then
        PP = d%P_root
      else
        PP = P_cur
      endif

      if (d%P(d%j) <= P_old .and. d%P(d%j) >= PP) then
        do while (d%P(d%j) >= PP)
          d%T(d%j) = self%contd8(1, d%P(d%j))
          d%z(d%j) = self%contd8(2, d%P(d%j))
          call mixing_ratios(d, d%P(d%j), d%T(d%j), d%f_i(d%j,:), f_dry)
          d%j = d%j + 1
          if (d%j > size(d%P)) exit
        enddo
      endif

    endif

  end subroutine

  subroutine find_root(dop, d, P_old, P_cur, ind, P_root)
    use futils, only: brent_class
    type(dop853_class), intent(inout) :: dop
    type(AdiabatProfileData), intent(inout) :: d
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
      call root_fcn(d, x_, dop%contd8(1, x_), d%gout_tmp)
      f_ = d%gout_tmp(ind)
    end function
  end subroutine

  subroutine root_fcn(d, P, T, gout)
    type(AdiabatProfileData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(in) :: T
    real(dp), intent(out) :: gout(:) 

    real(dp) :: f_dry, P_sat
    integer :: i

    call mixing_ratios(d, P, T, d%f_i_cur, f_dry)
    d%P_i_cur = d%f_i_cur*P

    do i = 1,d%sp%ng
      ! compute saturation vapor pressures
      P_sat = huge(1.0_dp)
      if (allocated(d%sp%g(i)%sat)) then
        P_sat = d%RH(i)*d%sp%g(i)%sat%sat_pressure(T)
      endif

      if (d%sp_type(i) == CondensingSpeciesType) then
        gout(i+1) = 1.0
      elseif (d%sp_type(i) == DrySpeciesType) then
        gout(i+1) = P_sat - d%P_i_cur(i)
      endif

    enddo

    ! also stop at tropopause
    gout(1) = T - d%T_trop    

  end subroutine

  subroutine right_hand_side_dop(self, P, u, du)
    use clima_eqns, only: gravity, heat_capacity_eval
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: P
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du

    select type (self)
    class is (dop853_custom)
      call right_hand_side(self%d, P, u, du)
    end select

  end subroutine

  subroutine update_f_i_dry(d, P, f_i_layer)
    type(AdiabatProfileData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(in) :: f_i_layer(:)

    integer :: i
    real(dp) :: P_dry

    d%P_i_cur(:) = f_i_layer(:)*P
    P_dry = 0.0_dp
    do i = 1,d%sp%ng
      if (d%sp_type(i) == DrySpeciesType) then
        P_dry = P_dry + d%P_i_cur(i)
      endif
    enddo
    d%f_i_dry(:) = d%P_i_cur(:)/P_dry

  end subroutine

  subroutine mixing_ratios(d, P, T, f_i_layer, f_dry)
    type(AdiabatProfileData), intent(inout) :: d
    real(dp), intent(in) :: P, T
    real(dp), intent(out) :: f_i_layer(:), f_dry

    real(dp) :: f_moist
    integer :: i

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
    do i = 1,d%sp%ng
      if (d%sp_type(i) == DrySpeciesType) then
        f_i_layer(i) = f_dry*d%f_i_dry(i)
      endif
    enddo

  end subroutine

  subroutine right_hand_side(d, P, u, du)
    use clima_eqns, only: heat_capacity_eval, gravity
    use clima_eqns_water, only: Rgas_cgs => Rgas
    type(AdiabatProfileData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(in) :: u(:)
    real(dp), intent(out) :: du(:)

    real(dp) :: T, z

    real(dp) :: f_dry, cp_dry
    real(dp) :: cp, mubar, L
    real(dp) :: first_sumation, second_sumation, Beta_i
    real(dp) :: grav
    real(dp) :: dlnT_dlnP, dT_dP, dz_dP
    logical :: found
    integer :: i

    real(dp), parameter :: Rgas_si = Rgas_cgs/1.0e7_dp ! ideal gas constant in SI units (J/(mol*K))

    ! unpack u
    T = u(1)
    z = u(2)

    call mixing_ratios(d, P, T, d%f_i_cur, f_dry)

    mubar = 0.0_dp
    do i = 1,d%sp%ng
      mubar = mubar + d%f_i_cur(i)*d%sp%g(i)%mass
    enddo

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
      if (d%sp_type(i) == DrySpeciesType) then
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
    ! Convert to dT/dP. Here we introduce CGS units
    ! P = [dynes/cm2]
    ! T = [K]
    dT_dP = dlnT_dlnP*(T/P)

    ! rate of change of altitude
    grav = gravity(d%planet_radius, d%planet_mass, z)
    dz_dP = -(Rgas_cgs*T)/(grav*P*mubar)

    du(:) = [dT_dP, dz_dP]

  end subroutine

  !> Analytical function for the alitude vs height for constant T and mubar.
  pure function altitude_vs_P(P, T, mubar, P0, z0, planet_mass, planet_radius) result(z)
    use clima_const, only: k_boltz, N_avo
    real(dp), intent(in) :: P
    real(dp), intent(in) :: T, mubar, P0, z0, planet_mass, planet_radius
    real(dp) :: z
    real(dp), parameter :: G_grav_cgs = 6.67e-8_dp
    
    z = ((N_avo*k_boltz*T)/(G_grav_cgs*planet_mass*mubar)*log(P/P0) &
        + 1/(planet_radius + z0))**(-1.0_dp) - planet_radius

  end function

end module