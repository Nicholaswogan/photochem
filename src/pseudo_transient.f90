module pseudo_transient
  !! Standalone pseudo-transient continuation (PTC) solver.
  !! Uses implicit Euler with Newton + backtracking line search.
  !! Supports dense or LAPACK banded Jacobians supplied by the caller.
  use iso_fortran_env, only: dp => real64, int32
  use iso_c_binding, only: c_ptr, c_null_ptr
  implicit none
  private

  public :: PTCOptions, PTCSolver

  type :: PTCOptions
    real(dp) :: rtol = 1.0e-6_dp
    real(dp) :: atol = 1.0e-20_dp
    real(dp) :: dt_init = 1.0_dp
    real(dp) :: dt_min  = epsilon(1.0_dp)
    real(dp) :: dt_max  = huge(1.0_dp)
    real(dp) :: dt_growth = 2.0_dp
    real(dp) :: dt_shrink = 0.5_dp
    real(dp) :: ls_c = 1.0e-4_dp      !! Armijo parameter
    integer  :: max_newton = 8
    integer  :: max_backtracks = 10
    integer  :: max_steps = 10000
    logical  :: use_band = .false.
    integer  :: kl = 0
    integer  :: ku = 0
    integer  :: verbosity = 0
  end type PTCOptions

  abstract interface
    subroutine ptc_rhs(y, f, udata, ierr)
      import :: dp, c_ptr
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: f(:)
      type(c_ptr), value :: udata
      integer, intent(out) :: ierr
    end subroutine ptc_rhs
  end interface

  abstract interface
    subroutine ptc_jac(y, J, udata, ierr)
      import :: dp, c_ptr
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: J(:,:)
      type(c_ptr), value :: udata
      integer, intent(out) :: ierr
    end subroutine ptc_jac
  end interface

  ! External LAPACK declarations
  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: dp, int32
      integer(int32), intent(in) :: n, nrhs, lda, ldb
      integer(int32), intent(out) :: ipiv(*)
      integer(int32), intent(out) :: info
      real(dp), intent(inout) :: a(lda,*), b(ldb,*)
    end subroutine dgesv

    subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
      import :: dp, int32
      integer(int32), intent(in) :: m, n, kl, ku, ldab
      integer(int32), intent(out) :: ipiv(*)
      integer(int32), intent(out) :: info
      real(dp), intent(inout) :: ab(ldab,*)
    end subroutine dgbtrf

    subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
      import :: dp, int32
      character(len=*), intent(in) :: trans
      integer(int32), intent(in) :: n, kl, ku, nrhs, ldab, ldb
      integer(int32), intent(in) :: ipiv(*)
      integer(int32), intent(out) :: info
      real(dp), intent(inout) :: ab(ldab,*), b(ldb,*)
    end subroutine dgbtrs
  end interface

  type :: PTCSolver
    type(PTCOptions) :: opt
    procedure(ptc_rhs), pointer :: rhs => null()
    procedure(ptc_jac), pointer :: jac => null()
    type(c_ptr) :: udata = c_null_ptr
    logical :: use_band = .false.
    logical :: use_atol_vec = .false.
    integer :: n = 0
    integer :: ldab = 0
    integer(int32) :: n_i32 = 0, ldab_i32 = 0
    real(dp) :: dt = 0.0_dp
    real(dp) :: res_norm = huge(1.0_dp)
    logical :: initialized = .false.
    real(dp), allocatable :: y(:), f(:), f_trial(:), y_trial(:), delta(:)
    real(dp), allocatable :: atol_vec(:), wt(:)
    real(dp), allocatable :: Jd(:,:), Ab(:,:)
    integer(int32), allocatable :: ipiv(:)
  contains
    procedure :: initialize
    procedure :: step
    procedure :: integrate
  end type PTCSolver

contains

  subroutine initialize(self, y0, rhs, jac, opts, udata, atol_vec, err_code)
    class(PTCSolver), intent(inout) :: self
    real(dp), intent(in) :: y0(:)
    procedure(ptc_rhs) :: rhs
    procedure(ptc_jac) :: jac
    type(PTCOptions), intent(in), optional :: opts
    type(c_ptr), value, intent(in), optional :: udata
    real(dp), intent(in), optional :: atol_vec(:)
    integer, intent(out) :: err_code

    err_code = 0
    self%opt = PTCOptions()
    if (present(opts)) self%opt = opts
    self%use_band = self%opt%use_band
    if (present(udata)) self%udata = udata
    self%rhs => rhs
    self%jac => jac

    self%n = size(y0)
    self%n_i32 = int(self%n, int32)
    self%ldab = 2*self%opt%kl + self%opt%ku + 1
    self%ldab_i32 = int(self%ldab, int32)

    allocate(self%y(self%n), self%f(self%n), self%f_trial(self%n), &
             self%y_trial(self%n), self%delta(self%n))
    allocate(self%wt(self%n))
    if (self%use_band) then
      allocate(self%Ab(self%ldab, self%n), self%ipiv(self%n))
    else
      allocate(self%Jd(self%n, self%n), self%ipiv(self%n))
    endif

    self%use_atol_vec = .false.
    if (present(atol_vec)) then
      if (size(atol_vec) /= self%n) then
        err_code = -12
        return
      endif
      allocate(self%atol_vec(self%n))
      self%atol_vec = atol_vec
      self%use_atol_vec = .true.
    endif

    self%y = y0
    self%dt = self%opt%dt_init
    call self%rhs(self%y, self%f, self%udata, err_code)
    if (err_code /= 0) return
    if (self%use_atol_vec) then
      call build_weights(self%y, self%opt%rtol, self%opt%atol, self%use_atol_vec, atol_vec=self%atol_vec, w=self%wt)
    else
      call build_weights(self%y, self%opt%rtol, self%opt%atol, self%use_atol_vec, w=self%wt)
    endif
    self%res_norm = residual_norm(self%f, self%wt)
    self%initialized = .true.
  end subroutine initialize

  subroutine step(self, converged, err_code)
    !! Perform one outer PTC step (with internal Newton/line search).
    class(PTCSolver), intent(inout) :: self
    logical, intent(out) :: converged
    integer, intent(out) :: err_code
    integer :: it, bt
    real(dp) :: res_new
    real(dp) :: dt_floor
    integer(int32) :: info

    converged = .false.
    err_code = 0
    if (.not. self%initialized) then
      err_code = -100
      return
    endif

    if (self%res_norm <= 1.0_dp) then
      converged = .true.
      return
    endif

    dt_floor = max(self%opt%dt_min, epsilon(1.0_dp)*max(1.0_dp, maxval(abs(self%y))))

    retry_step: do
      do it = 1, self%opt%max_newton
        if (self%use_band) then
          self%Ab = 0.0_dp
          call self%jac(self%y, self%Ab, self%udata, err_code)
          if (err_code > 0) then
            self%dt = max(dt_floor, self%dt*self%opt%dt_shrink)
            if (self%dt <= dt_floor) then
              err_code = -1
              return
            endif
            err_code = 0
            cycle retry_step
          elseif (err_code /= 0) then
            return
          endif
          call assemble_band_system(self%n, self%opt%kl, self%opt%ku, self%ldab, self%Ab, self%dt)
          self%delta = -self%f
          call dgbtrf(self%n_i32, self%n_i32, self%opt%kl, self%opt%ku, self%Ab, self%ldab_i32, self%ipiv, info)
          if (info /= 0) then
            err_code = 1000 + info
            exit
          endif
          call dgbtrs('N', self%n_i32, self%opt%kl, self%opt%ku, 1, self%Ab, self%ldab_i32, self%ipiv, self%delta, self%n_i32, info)
          if (info /= 0) then
            err_code = 1100 + info
            exit
          endif
        else
          call self%jac(self%y, self%Jd, self%udata, err_code)
          if (err_code > 0) then
            self%dt = max(dt_floor, self%dt*self%opt%dt_shrink)
            if (self%dt <= dt_floor) then
              err_code = -1
              return
            endif
            err_code = 0
            cycle retry_step
          elseif (err_code /= 0) then
            return
          endif
          call assemble_dense_system(self%n, self%Jd, self%dt)
          self%delta = -self%f
          call dgesv(self%n_i32, 1, self%Jd, self%n_i32, self%ipiv, self%delta, self%n_i32, info)
          if (info /= 0) then
            err_code = 2000 + info
            exit
          endif
        endif

        res_new = self%res_norm
        self%y_trial = self%y
        do bt = 0, self%opt%max_backtracks
          self%y_trial = self%y + self%delta
          call self%rhs(self%y_trial, self%f_trial, self%udata, err_code)
          if (err_code > 0) then
            self%dt = max(dt_floor, self%dt*self%opt%dt_shrink)
            if (self%dt <= dt_floor) then
              err_code = -1
              return
            endif
            err_code = 0
            cycle retry_step
          elseif (err_code /= 0) then
            return
          endif
          if (self%use_atol_vec) then
            call build_weights(self%y_trial, self%opt%rtol, self%opt%atol, self%use_atol_vec, atol_vec=self%atol_vec, w=self%wt)
          else
            call build_weights(self%y_trial, self%opt%rtol, self%opt%atol, self%use_atol_vec, w=self%wt)
          endif
          res_new = residual_norm(self%f_trial, self%wt)
          if (res_new <= (1.0_dp - self%opt%ls_c)*self%res_norm) exit
          self%delta = 0.5_dp*self%delta
        enddo

        self%y = self%y_trial
        self%f = self%f_trial
        self%res_norm = res_new

        if (self%res_norm <= 1.0_dp) then
          converged = .true.
          exit
        endif
      enddo

      exit retry_step
    enddo

    if (converged) then
      self%dt = min(self%opt%dt_max, self%dt*self%opt%dt_growth)
    else
      self%dt = max(dt_floor, self%dt*self%opt%dt_shrink)
      if (self%dt <= dt_floor) err_code = -1
    endif
  end subroutine step

  subroutine integrate(self, converged, err_code, max_steps)
    class(PTCSolver), intent(inout) :: self
    logical, intent(out) :: converged
    integer, intent(out) :: err_code
    integer, intent(in), optional :: max_steps
    integer :: nsteps, msteps
    logical :: loc_conv
    integer :: loc_err

    msteps = self%opt%max_steps
    if (present(max_steps)) msteps = max_steps
    converged = .false.
    err_code = 0

    do nsteps = 1, msteps
      call self%step(loc_conv, loc_err)
      if (loc_err /= 0) then
        err_code = loc_err
        return
      endif
      if (loc_conv) then
        converged = .true.
        exit
      endif
    enddo
    if (.not. converged .and. err_code == 0) err_code = -2
  end subroutine integrate

  subroutine build_weights(y, rtol, atol, use_atol_vec, atol_vec, w)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: rtol, atol
    logical, intent(in) :: use_atol_vec
    real(dp), intent(in), optional :: atol_vec(:)
    real(dp), intent(out) :: w(:)
    if (use_atol_vec) then
      if (present(atol_vec)) then
        w = atol_vec + rtol*abs(y)
      else
        w = atol + rtol*abs(y)
      endif
    else
      w = atol + rtol*abs(y)
    endif
  end subroutine build_weights

  function residual_norm(r, w) result(val)
    real(dp), intent(in) :: r(:), w(:)
    real(dp) :: val
    integer :: ii, nr
    nr = size(r)
    val = 0.0_dp
    do ii = 1, nr
      val = val + (r(ii)/w(ii))**2
    enddo
    val = sqrt(val/real(nr,dp))
  end function residual_norm

  subroutine assemble_dense_system(nr, J, dtloc)
    integer, intent(in) :: nr
    real(dp), intent(inout) :: J(nr,nr)
    real(dp), intent(in) :: dtloc
    integer :: ii
    J = -dtloc*J
    do ii = 1, nr
      J(ii,ii) = J(ii,ii) + 1.0_dp
    enddo
  end subroutine assemble_dense_system

  subroutine assemble_band_system(nr, kl, ku, ld, Abnd, dtloc)
    integer, intent(in) :: nr, kl, ku, ld
    real(dp), intent(inout) :: Abnd(ld, nr)
    real(dp), intent(in) :: dtloc
    integer :: jj
    Abnd = -dtloc*Abnd
    do jj = 1, nr
      Abnd(kl+ku+1, jj) = Abnd(kl+ku+1, jj) + 1.0_dp
    enddo
  end subroutine assemble_band_system

end module pseudo_transient
