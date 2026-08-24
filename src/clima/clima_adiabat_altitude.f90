submodule(clima_adiabat) clima_adiabat_altitude
  use clima_useful, only: solve_ivp_dop853
  implicit none

  type :: AltitudeIntegrationData
    type(linear_interp_1d) :: T_interp
    type(linear_interp_1d) :: mubar_interp
    real(dp), pointer :: P_e(:) => NULL()
    real(dp), pointer :: z_e(:) => NULL()
    real(dp) :: planet_mass
    real(dp) :: planet_radius
  end type

contains

  module subroutine AdiabatClimate_compute_altitude(self, err)
    use clima_eqns, only: gravity
    class(AdiabatClimate), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    type(AltitudeIntegrationData), target :: d
    integer :: i, j, ierr, k
    integer :: j_start_top, j_start_bot
    real(dp) :: z_shift, z_ref_for_radius
    real(dp) :: sentinel
    logical :: is_on_k, is_on_kp1
    real(dp), allocatable :: T_grid(:), mubar_grid(:), P_grid(:)
    real(dp), allocatable :: P_interp(:), T_interp(:), mubar_interp(:)
    real(dp), allocatable, target :: P_e(:), z_e(:)
    character(256) :: msg

    if (self%nz < 2) then
      err = 'compute_altitude requires nz >= 2'
      return
    endif
    if (self%P_surf <= 0.0_dp .or. self%P_top <= 0.0_dp) then
      err = 'compute_altitude requires positive P_surf and P_top'
      return
    endif

    allocate(P_e(2*self%nz+1), z_e(2*self%nz+1))
    allocate(T_grid(self%nz+1), mubar_grid(self%nz+1), P_grid(self%nz+1))
    allocate(P_interp(self%nz+1), T_interp(self%nz+1), mubar_interp(self%nz+1))

    P_e(1) = self%P_surf
    do i = 1,self%nz
      P_e(2*i) = self%P(i)
      if (i < self%nz) P_e(2*i+1) = sqrt(self%P(i)*self%P(i+1))
    enddo
    P_e(2*self%nz+1) = self%P_top

    if (minval(P_e(1:2*self%nz) - P_e(2:2*self%nz+1)) <= 10.0_dp*spacing(P_e(1))) then
      err = 'compute_altitude: pressure grid is not strictly decreasing.'
      return
    endif

    P_grid(1) = self%P_surf
    T_grid(1) = self%T_surf
    mubar_grid(1) = 0.0_dp
    do j = 1,self%sp%ng
      mubar_grid(1) = mubar_grid(1) + self%f_i_surf(j)*self%sp%g(j)%mass
    enddo
    do i = 1,self%nz
      P_grid(i+1) = self%P(i)
      T_grid(i+1) = self%T(i)
      mubar_grid(i+1) = 0.0_dp
      do j = 1,self%sp%ng
        mubar_grid(i+1) = mubar_grid(i+1) + self%f_i(i,j)*self%sp%g(j)%mass
      enddo
    enddo

    do i = 1,self%nz+1
      P_interp(i) = log10(P_grid(self%nz+2-i))
      T_interp(i) = T_grid(self%nz+2-i)
      mubar_interp(i) = mubar_grid(self%nz+2-i)
    enddo

    call d%T_interp%initialize(P_interp, T_interp, ierr)
    if (ierr /= 0) then
      err = 'compute_altitude: failed to initialize temperature interpolator.'
      return
    endif
    call d%mubar_interp%initialize(P_interp, mubar_interp, ierr)
    if (ierr /= 0) then
      err = 'compute_altitude: failed to initialize mubar interpolator.'
      return
    endif

    d%P_e => P_e
    d%z_e => z_e
    d%planet_mass = self%planet_mass
    d%planet_radius = self%planet_radius

    sentinel = 1.0e300_dp
    z_e(:) = sentinel

    if (self%reference_pressure > 0.0_dp) then
      if (self%reference_pressure < self%P_top .or. self%reference_pressure > self%P_surf) then
        write(msg,'(a,1pe11.4,a,1pe11.4,a,1pe11.4,a)') &
          'compute_altitude: reference_pressure=', self%reference_pressure, &
          ' outside model domain [', self%P_top, ', ', self%P_surf, ']'
        err = trim(msg)
        return
      endif
    endif

    if (self%reference_pressure <= 0.0_dp) then
      call integrate_altitude_profile(d, self%rtol, self%atol, err)
      if (allocated(err)) return
      z_ref_for_radius = 0.0_dp
    else
      do k = 1,size(P_e)-1
        if (self%reference_pressure <= P_e(k) .and. self%reference_pressure >= P_e(k+1)) then
          exit
        endif
      enddo
      if (k > size(P_e)-1) then
        write(msg,'(a,1pe11.4,a,1pe11.4,a,1pe11.4,a)') &
          'compute_altitude: could not bracket reference pressure ', self%reference_pressure, &
          ' in [', P_e(size(P_e)), ', ', P_e(1), ']'
        err = trim(msg)
        return
      endif

      is_on_k = self%reference_pressure == P_e(k)
      is_on_kp1 = self%reference_pressure == P_e(k+1)

      if (is_on_k) then
        z_e(k) = 0.0_dp
        j_start_top = k + 1
        j_start_bot = k - 1
      elseif (is_on_kp1) then
        z_e(k+1) = 0.0_dp
        j_start_top = k + 2
        j_start_bot = k
      else
        j_start_top = k + 1
        j_start_bot = k
      endif

      ! Integrate from reference pressure upward (toward lower pressures/top).
      if (j_start_top <= size(P_e)-1) then
        call integrate_altitude_segment(d, self%rtol, self%atol, &
                                        self%reference_pressure, 0.0_dp, &
                                        j_start_top, size(P_e)-1, +1, err)
        if (allocated(err)) return
      endif

      ! Integrate from reference pressure downward (toward higher pressures/surface).
      if (j_start_bot >= 1) then
        call integrate_altitude_segment(d, self%rtol, self%atol, &
                                        self%reference_pressure, 0.0_dp, &
                                        j_start_bot, 1, -1, err)
        if (allocated(err)) return
      endif

      do i = 1,size(P_e)-1
        if (z_e(i) == sentinel) then
          err = 'compute_altitude: failed to populate altitude profile from reference-pressure integration.'
          return
        endif
      enddo

      ! Keep stored altitude surface-anchored for compatibility.
      z_shift = z_e(1)
      z_e(1:size(P_e)-1) = z_e(1:size(P_e)-1) - z_shift
      z_e(size(P_e)) = z_e(size(P_e)-1) + (z_e(size(P_e)-1) - z_e(size(P_e)-2))
      z_ref_for_radius = -z_shift
    endif

    do i = 1,self%nz
      self%z(i) = z_e(2*i)
      self%dz(i) = z_e(2*i+1) - z_e(2*i-1)
      self%gravity(i) = gravity(self%planet_radius, self%planet_mass, self%z(i) - z_ref_for_radius)
    enddo
    self%gravity_surf = gravity(self%planet_radius, self%planet_mass, -z_ref_for_radius)

  end subroutine

  subroutine integrate_altitude_profile(d, rtol, atol, err)
    type(AltitudeIntegrationData), intent(inout) :: d
    real(dp), intent(in) :: rtol, atol
    character(:), allocatable, intent(out) :: err

    d%z_e(1) = 0.0_dp

    call integrate_altitude_segment(d, rtol, atol, &
                                    d%P_e(1), 0.0_dp, &
                                    2, size(d%P_e)-1, +1, err)
    if (allocated(err)) return
    d%z_e(size(d%z_e)) = d%z_e(size(d%z_e)-1) + (d%z_e(size(d%z_e)-1) - d%z_e(size(d%z_e)-2))

  end subroutine

  subroutine integrate_altitude_segment(d, rtol, atol, &
                                        P_start, z_start, &
                                        j_start, j_end, j_step, err)
    type(AltitudeIntegrationData), target, intent(inout) :: d
    real(dp), intent(in) :: rtol, atol
    real(dp), intent(in) :: P_start, z_start
    integer, intent(in) :: j_start, j_end, j_step
    character(:), allocatable, intent(out) :: err

    integer :: idid, m, i, idx
    real(dp), allocatable :: x_eval(:), y_eval(:,:)

    m = abs(j_end - j_start) + 1
    allocate(x_eval(m), y_eval(1,m))

    do i = 1,m
      idx = j_start + (i-1)*j_step
      x_eval(i) = d%P_e(idx)
    enddo

    call solve_ivp_dop853(rhs_altitude_segment, P_start, [z_start], x_eval, rtol, atol, y_eval, idid, err)
    if (allocated(err)) return
    if (idid < 0) then
      err = 'compute_altitude: altitude segment integration failed.'
      return
    endif

    do i = 1,m
      idx = j_start + (i-1)*j_step
      d%z_e(idx) = y_eval(1,i)
    enddo

  contains

    subroutine rhs_altitude_segment(P, u, du)
      use clima_eqns, only: gravity
      use clima_eqns_water, only: Rgas_cgs => Rgas
      real(dp), intent(in) :: P
      real(dp), intent(in) :: u(:)
      real(dp), intent(out) :: du(:)

      real(dp) :: T, mubar, grav
      call d%T_interp%evaluate(log10(P), T)
      call d%mubar_interp%evaluate(log10(P), mubar)
      grav = gravity(d%planet_radius, d%planet_mass, u(1))
      du(1) = -(Rgas_cgs*T)/(grav*P*mubar)
    end subroutine

  end subroutine

end submodule
