module clima_adiabat_dry
  use clima_const, only: dp
  use clima_types, only: Species
  use dop853_module, only: dop853_class
  use futils, only: brent_class
  use clima_eqns, only: ocean_solubility_fcn
  use clima_adiabat_general, only: OceanFunction, DrySpeciesType, CondensingSpeciesType
  use linear_interpolation_module, only: linear_interp_1d
  implicit none

  type :: AdiabatDryProfileData
    ! Input data
    real(dp), pointer :: P_in(:)
    real(dp), pointer :: T_in(:)
    real(dp), pointer :: f_i_in(:,:)
    type(Species), pointer :: sp
    integer, pointer :: nz
    real(dp), pointer :: planet_mass
    real(dp), pointer :: planet_radius
    real(dp), pointer :: P_top
    real(dp), pointer :: rtol
    real(dp), pointer :: atol
    ! Ouput
    real(dp), pointer :: P(:)
    real(dp), pointer :: z(:)
    real(dp), pointer :: T(:)
    real(dp), pointer :: f_i(:,:)
    real(dp), pointer :: lapse_rate(:)

    real(dp) :: P_surf !! surface pressure from inputs

    type(linear_interp_1d) :: T_interp !! Interpolates input T
    type(linear_interp_1d), allocatable :: f_i_interp(:) !! Interpolates mixing ratios (ng)

    !> index for saving T, z and mixing ratios
    integer :: j

    !> work variables. All dimension (ng)
    real(dp), allocatable :: f_i_cur(:)

    !> This helps us propogate error messages outside of the integration
    character(:), allocatable :: err

  end type

  type, extends(dop853_class) :: dop853_custom
    type(AdiabatDryProfileData), pointer :: d => NULL()
  end type

contains

  subroutine make_profile_dry(P_in, T_in, f_i_in, &
                              sp, nz, &
                              planet_mass, planet_radius, P_top, &
                              rtol, atol, &
                              P, z, T, f_i, lapse_rate, &
                              err)
    use futils, only: linspace
    real(dp), target, intent(in) :: P_in(:)
    real(dp), target, intent(in) :: T_in(:)
    real(dp), intent(in) :: f_i_in(:,:)
    type(Species), target, intent(inout) :: sp
    integer, target, intent(in) :: nz
    real(dp), target, intent(in) :: planet_mass, planet_radius
    real(dp), target, intent(in) :: P_top
    real(dp), target, intent(in) :: rtol, atol
    real(dp), target, intent(out) :: P(:), z(:), T(:)
    real(dp), target, intent(out) :: f_i(:,:)
    real(dp), target, intent(out) :: lapse_rate(:)
    character(:), allocatable, intent(out) :: err

    real(dp), target, allocatable :: tmp_arr(:), f_i_in_copy(:,:)
    real(dp), allocatable :: log10P_in_interp(:)
    type(AdiabatDryProfileData) :: d
    integer :: i, ierr

    ! Check inputs
    if (any(T_in < 0.0_dp)) then
      err = '`T` can not have negative elements'
      return
    endif
    if (any(P_in < 0.0_dp)) then
      err = '`P` can not have negative elements'
      return
    endif
    if (P_in(1) < P_top) then
      err = 'The first element of `P` must be greater than `P_top`'
      return
    endif
    if (size(T_in) <= 1) then
      err = '`P` must have a size > 1.'
      return
    endif
    if (size(T_in) /= size(P_in)) then
      err = '`T` and `P` must have the same length'
      return
    endif
    do i = 1,size(P_in)-1
      if (P_in(i+1) >= P_in(i)) then
        err = '`P` must be strictly decreasing'
        return
      endif
    enddo
    if (any(f_i_in < 0.0_dp)) then
      err = '`f_i` can not have negative elements'
      return
    endif
    if (size(f_i_in,1) /= size(P_in)) then
      err = 'The first dimension of `f_i` must match the size of `P`'
      return
    endif
    if (size(f_i_in,2) /= sp%ng) then
      err = 'The second dimension of `f_i` must match the the number of species'
      return
    endif

    ! First, lets renormalize the mixing ratios so they sum to 1
    f_i_in_copy = f_i_in
    do i = 1,size(f_i_in,1)
      f_i_in_copy(i,:) = f_i_in_copy(i,:)/sum(f_i_in_copy(i,:))
    enddo

    ! Associate inputs
    d%P_in => P_in
    d%T_in => T_in
    d%f_i_in => f_i_in_copy
    d%sp => sp
    d%nz => nz
    d%planet_mass => planet_mass
    d%planet_radius => planet_radius
    d%P_top => P_top
    d%rtol => rtol
    d%atol => atol
    ! outputs
    d%P => P
    d%z => z
    d%T => T
    d%f_i => f_i
    d%lapse_rate => lapse_rate

    ! allocate
    allocate(d%f_i_cur(sp%ng))

    ! Make P profile
    d%P_surf = P_in(1)
    call linspace(log10(d%P_surf),log10(P_top),P)
    P(:) = 10.0_dp**P(:)
    P(1) = d%P_surf
    P(2*nz+1) = P_top

    ! Build the interpolators
    allocate(log10P_in_interp(size(P_in)))
    log10P_in_interp = log10(P_in)
    log10P_in_interp = log10P_in_interp(size(log10P_in_interp):1:-1)

    ! Temperature
    allocate(tmp_arr(size(T_in)))
    tmp_arr = T_in
    tmp_arr = tmp_arr(size(tmp_arr):1:-1)
    call d%T_interp%initialize(log10P_in_interp, tmp_arr, ierr)
    if (ierr /= 0) then
      err = 'Failed to initialize interpolator in "make_profile_dry"'
      return
    endif

    ! Mixing ratios
    allocate(d%f_i_interp(sp%ng))
    do i = 1,sp%ng
      tmp_arr = log10(f_i_in_copy(:,i))
      tmp_arr = tmp_arr(size(tmp_arr):1:-1)
      call d%f_i_interp(i)%initialize(log10P_in_interp, tmp_arr, ierr)
      if (ierr /= 0) then
        err = 'Failed to initialize interpolator in "make_profile_dry"'
        return
      endif
    enddo

    call integrate(d, err)
    if (allocated(err)) return

  end subroutine

  subroutine integrate(d, err)
    type(AdiabatDryProfileData), target, intent(inout) :: d
    character(:), allocatable, intent(out) :: err

    type(dop853_custom) :: dop
    logical :: status_ok
    integer :: idid
    real(dp) :: Pn, u(1)

    call dop%initialize(fcn=right_hand_side_dop, solout=solout_dop, n=1, &
                        iprint=0, icomp=[1], status_ok=status_ok)
    dop%d => d
    d%j = 2
    call mixing_ratios(d, d%P_surf, d%f_i(1,:))
    call d%T_interp%evaluate(log10(d%P_surf),d%T(1))
    d%z(1) = 0.0_dp
    call dry_adiabat(d, d%P_surf, d%lapse_rate(1), err)
    if (allocated(err)) return
    if (.not. status_ok) then
      err = 'dop853 initialization failed'
      return
    endif

    u = [0.0_dp]
    Pn = d%P_surf
    call dop%integrate(Pn, u, d%P_top, [d%rtol], [d%atol], iout=2, idid=idid)
    if (allocated(d%err)) then
      err = d%err
      return
    endif

  end subroutine

  subroutine solout_dop(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout
    
    type(AdiabatDryProfileData), pointer :: d
    real(dp) :: z_cur, P_cur, P_old
    
    P_old = xold
    P_cur = x
    z_cur = y(1)

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

    ! save the results
    if (d%j <= size(d%P)) then
      if (d%P(d%j) <= P_old .and. d%P(d%j) >= P_cur) then
        do while (d%P(d%j) >= P_cur)
          d%z(d%j) = self%contd8(1, d%P(d%j))
          call d%T_interp%evaluate(log10(d%P(d%j)), d%T(d%j))
          call mixing_ratios(d, d%P(d%j), d%f_i(d%j,:))
          call dry_adiabat(d, d%P(d%j), d%lapse_rate(d%j), d%err)
          if (allocated(d%err)) return
          d%j = d%j + 1
          if (d%j > size(d%P)) exit
        enddo
      endif
    endif

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

  subroutine dry_adiabat(d, P, lapse_rate, err)
    use clima_eqns, only: heat_capacity_eval
    use clima_eqns_water, only: Rgas_cgs => Rgas
    type(AdiabatDryProfileData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(out) :: lapse_rate
    character(:), allocatable, intent(out) :: err

    real(dp) :: T, cp_tmp, cp
    logical :: found
    integer :: i

    real(dp), parameter :: Rgas_si = Rgas_cgs/1.0e7_dp ! ideal gas constant in SI units (J/(mol*K))

    ! Get temperature
    call d%T_interp%evaluate(log10(P), T)

    ! Get mixing ratios
    call mixing_ratios(d, P, d%f_i_cur)

    cp = tiny(0.0_dp)
    do i = 1,d%sp%ng
      call heat_capacity_eval(d%sp%g(i)%thermo, T, found, cp_tmp) ! J/(mol*K)
      if (.not. found) then
        err = "Failed to compute heat capacity"
        return
      endif
      ! J/(mol*K)
      cp = cp + cp_tmp*d%f_i_cur(i) ! J/(mol*K)
    enddo

    lapse_rate = Rgas_si/cp

  end subroutine

  subroutine mixing_ratios(d, P, f_i_layer)
    type(AdiabatDryProfileData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(out) :: f_i_layer(:)

    integer :: i

    do i = 1,d%sp%ng
      call d%f_i_interp(i)%evaluate(log10(P),f_i_layer(i))
      f_i_layer(i) = 10.0_dp**f_i_layer(i)
    enddo

  end subroutine

  subroutine right_hand_side(d, P, u, du)
    use clima_eqns, only: gravity
    use clima_eqns_water, only: Rgas_cgs => Rgas
    type(AdiabatDryProfileData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(in) :: u(:)
    real(dp), intent(out) :: du(:)

    real(dp) :: z
    real(dp) :: T, mubar, grav, dz_dP
    integer :: i

    ! unpack u
    z = u(1)

    ! Get temperature
    call d%T_interp%evaluate(log10(P), T)

    ! Get mixing ratios
    call mixing_ratios(d, P, d%f_i_cur)

    mubar = 0.0_dp
    do i = 1,d%sp%ng
      mubar = mubar + d%f_i_cur(i)*d%sp%g(i)%mass
    enddo

    ! rate of change of altitude
    grav = gravity(d%planet_radius, d%planet_mass, z)
    dz_dP = -(Rgas_cgs*T)/(grav*P*mubar)

    du(:) = [dz_dP]

  end subroutine

end module