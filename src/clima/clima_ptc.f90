module clima_ptc
  use clima_const, only: wp => dp
  implicit none
  private

  integer, parameter, public :: PTC_JAC_DENSE = 1
  integer, parameter, public :: PTC_JAC_BAND  = 2

  integer, parameter, public :: PTC_REASON_NONE                  = 0
  integer, parameter, public :: PTC_CONVERGED_PSEUDO_FATOL       = 1
  integer, parameter, public :: PTC_CONVERGED_PSEUDO_FRTOL       = 2
  integer, parameter, public :: PTC_CONVERGED_USER               = 3
  integer, parameter, public :: PTC_DIVERGED_STEP_REJECTED       = -1
  integer, parameter, public :: PTC_DIVERGED_CALLBACK_FATAL      = -2
  integer, parameter, public :: PTC_DIVERGED_NOT_INITIALIZED     = -3
  integer, parameter, public :: PTC_DIVERGED_INVALID_INPUT       = -4
  integer, parameter, public :: PTC_DIVERGED_MAX_STEPS           = -5
  integer, parameter, public :: PTC_DIVERGED_STAGNATION          = -6

  public :: PTCSolver, wp

  type :: PTCSolver
    integer :: neq = 0  !! Number of unknowns in the state vector.
    integer :: jacobian_type = 0  !! Jacobian mode (`PTC_JAC_DENSE` or `PTC_JAC_BAND`).
    integer :: kl = 0  !! Number of sub-diagonals for banded Jacobians.
    integer :: ku = 0  !! Number of super-diagonals for banded Jacobians.
    integer :: ldab = 0  !! Leading dimension for LAPACK banded matrix storage.

    procedure(rhs_fcn), pointer, nopass :: f => null()  !! User residual callback computing `f(x)`.
    procedure(jac_fcn), pointer, nopass :: jac => null()  !! User Jacobian callback (dense or compact banded layout).
    procedure(verify_step_fcn), pointer, nopass :: verify => null()  !! Optional step verification callback.
    procedure(timestep_fcn), pointer, nopass :: compute_dt => null()  !! Optional timestep update callback.
    procedure(convergence_fcn), pointer, nopass :: custom_convergence => null()  !! Optional custom convergence callback.
    procedure(progress_fcn), pointer, nopass :: progress => null()  !! Optional progress callback (called at step 0 and after accepted steps).

    real(wp) :: dt = 0.0_wp  !! Current pseudo-time step size.
    real(wp) :: dt_initial = 0.0_wp  !! Initial pseudo-time step size.
    real(wp) :: dt_increment = 1.1_wp  !! Growth factor used by default timestep adaptation.
    real(wp) :: dt_max = 0.0_wp  !! Maximum pseudo-time step (`<=0` disables cap).
    logical :: increment_dt_from_initial_dt = .false.  !! If true, adapt from initial `(dt, fnorm)` pair.

    real(wp) :: fatol = 1.0e-50_wp  !! Absolute convergence tolerance on `||f(x)||_2`.
    real(wp) :: frtol = 1.0e-12_wp  !! Relative convergence tolerance on `||f(x)||_2 / ||f(x0)||_2`.

    real(wp) :: fnorm = -1.0_wp  !! Current residual norm `||f(x)||_2`.
    real(wp) :: fnorm_initial = -1.0_wp  !! Residual norm at first accepted step.
    real(wp) :: fnorm_previous = -1.0_wp  !! Residual norm from previous accepted step.
    logical :: residual_valid = .false.  !! True when `fvec/fnorm` correspond to current `x`.
    real(wp) :: fnorm_old = -1.0_wp  !! Cached residual norm for rollback state `x_old`.
    logical :: residual_old_valid = .false.  !! True when `fvec_old/fnorm_old` correspond to `x_old`.
    logical :: jac_valid = .false.  !! True when `jac_mat` corresponds to current `x`.
    logical :: jac_old_valid = .false.  !! True when `jac_mat_old` corresponds to `x_old`.

    integer :: steps = 0  !! Number of accepted pseudo-steps.
    integer :: rejects_total = 0  !! Total number of rejected step attempts.
    integer :: max_reject = 10  !! Maximum rejects allowed per `step()` before failure.
    integer :: max_steps = 10000  !! Maximum accepted steps allowed in `solve()`.
    integer :: reason = PTC_REASON_NONE  !! Solver state/reason code.
    logical :: enable_stagnation_check = .true.  !! Enable generic residual-norm stagnation detection.
    integer :: stagnation_warmup_steps = 10  !! Number of accepted steps before stagnation checks begin.
    integer :: stagnation_window = 150  !! Accepted non-improving steps allowed before stagnation failure.
    real(wp) :: stagnation_rel_improve_tol = 1.0e-3_wp  !! Relative improvement threshold counted as progress.
    integer :: stagnation_count = 0  !! Counter of consecutive accepted non-improving steps.
    real(wp) :: fnorm_best = huge(1.0_wp)  !! Best residual norm seen so far.
    logical :: progress_step0_emitted = .false.  !! Internal guard to emit initial progress only once.

    logical :: initialized = .false.  !! True once arrays and callbacks are configured.

    real(wp), allocatable :: x(:)  !! Current solution iterate.
    real(wp), allocatable :: x_old(:)  !! Backup iterate for rollback on rejected steps.
    real(wp), allocatable :: fvec(:)  !! Residual workspace.
    real(wp), allocatable :: fvec_old(:)  !! Cached residual vector for rollback state `x_old`.
    real(wp), allocatable :: step_vec(:)  !! Linear correction vector workspace.
    real(wp), allocatable :: rhs_mat(:, :)  !! Right-hand-side workspace passed to LAPACK solvers.

    real(wp), allocatable :: jac_mat(:, :)  !! Jacobian workspace (dense or compact banded).
    real(wp), allocatable :: jac_mat_old(:, :)  !! Cached Jacobian for rollback state `x_old`.
    real(wp), allocatable :: a_dense(:, :)  !! Dense system matrix workspace for `(I/dt - J)`.

    ! Banded Jacobian compact storage (LAPACK standard):
    ! jac_mat(ku+1+i-j, j) = J(i,j), for max(1,j-ku) <= i <= min(n,j+kl)
    real(wp), allocatable :: a_band(:, :)  !! Banded system matrix workspace for LAPACK `dgbsv`.

    integer, allocatable :: ipiv(:)  !! Pivot indices returned by LAPACK factorizations.
  contains
    procedure :: initialize => PTCSolver_initialize
    procedure :: step => PTCSolver_step
    procedure :: solve => PTCSolver_solve
    procedure :: check_convergence => PTCSolver_check_convergence
    procedure :: update_stagnation => PTCSolver_update_stagnation

    procedure :: set_verify_timestep => PTCSolver_set_verify_timestep
    procedure :: set_compute_timestep => PTCSolver_set_compute_timestep
    procedure :: set_custom_convergence => PTCSolver_set_custom_convergence
    procedure :: set_progress => PTCSolver_set_progress
  end type PTCSolver

    abstract interface
    subroutine rhs_fcn(solver, u, udot, ierr)
      import :: PTCSolver, wp
      implicit none
      class(PTCSolver), intent(in) :: solver  !! Solver object (`self`) passed for callback-side inspection.
      real(wp), intent(in) :: u(:)  !! Input state vector.
      real(wp), intent(out) :: udot(:)  !! Residual/right-hand-side vector `f(u)`.
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).
    end subroutine rhs_fcn

    subroutine jac_fcn(solver, u, jac, ierr)
      import :: PTCSolver, wp
      implicit none
      class(PTCSolver), intent(in) :: solver  !! Solver object (`self`) passed for callback-side inspection.
      real(wp), intent(in) :: u(:)  !! Input state vector.
      real(wp), intent(out) :: jac(:, :)  !! Jacobian in configured dense or compact-banded layout.
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).
    end subroutine jac_fcn

    !> Step verification callback timing/data contract.
    !! Called after a candidate state update (`solver%x` is candidate `x_{n+1}`).
    !! Residual (`solver%fvec/solver%fnorm`) is not guaranteed up-to-date for this `x`.
    !! Jacobian (`solver%jac_mat`) is not guaranteed up-to-date for this `x`.
    subroutine verify_step_fcn(solver, accept, reject_dt, ierr)
      import :: PTCSolver, wp
      implicit none
      class(PTCSolver), intent(in) :: solver  !! Solver object (`self`) passed for callback-side inspection.
      logical, intent(out) :: accept  !! Set true to accept candidate step.
      real(wp), intent(out) :: reject_dt  !! Suggested retry timestep if rejecting (ignored when accepting).
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).
    end subroutine verify_step_fcn

    !> Custom timestep callback timing/data contract.
    !! Called after step acceptance and residual recomputation at current `solver%x`.
    !! Residual (`solver%fvec/solver%fnorm`) is up-to-date for this `x`.
    !! Jacobian (`solver%jac_mat`) is not guaranteed up-to-date for this `x`.
    subroutine timestep_fcn(solver, new_dt, ierr)
      import :: PTCSolver, wp
      implicit none
      class(PTCSolver), intent(in) :: solver  !! Solver object (`self`) passed for callback-side inspection.
      real(wp), intent(out) :: new_dt  !! Computed timestep for the next accepted step.
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).
    end subroutine timestep_fcn

    !> Custom convergence callback timing/data contract.
    !! Called after step acceptance and residual recomputation at current `solver%x`.
    !! Residual (`solver%fvec/solver%fnorm`) is up-to-date for this `x`.
    !! Jacobian (`solver%jac_mat`) is not guaranteed up-to-date for this `x`.
    subroutine convergence_fcn(solver, converged, ierr)
      import PTCSolver
      implicit none
      class(PTCSolver), intent(in) :: solver  !! Solver object (`self`) passed for callback-side inspection.
      logical, intent(out) :: converged  !! Set true to terminate as converged.
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero fatal failure).
    end subroutine convergence_fcn

    !> Progress callback timing/data contract.
    !! Called once at start of `solve()` (step 0, with residual current),
    !! and once after each accepted step.
    !! For a given reported step, this callback is invoked after `rhs_fcn`
    !! and before `jac_fcn`.
    subroutine progress_fcn(solver)
      import PTCSolver
      implicit none
      class(PTCSolver), intent(in) :: solver
    end subroutine progress_fcn

  end interface

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: wp
      implicit none
      integer, intent(in) :: n  !! System size.
      integer, intent(in) :: nrhs  !! Number of right-hand sides.
      integer, intent(in) :: lda  !! Leading dimension of `a`.
      integer, intent(in) :: ldb  !! Leading dimension of `b`.
      integer, intent(out) :: ipiv(*)  !! LAPACK pivot indices.
      real(wp), intent(inout) :: a(lda, *)  !! Coefficient matrix (overwritten by LU factors).
      real(wp), intent(inout) :: b(ldb, *)  !! Right-hand side(s), overwritten by solution(s).
      integer, intent(out) :: info  !! LAPACK status code.
    end subroutine dgesv

    subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
      import :: wp
      implicit none
      integer, intent(in) :: n  !! System size.
      integer, intent(in) :: kl  !! Number of sub-diagonals.
      integer, intent(in) :: ku  !! Number of super-diagonals.
      integer, intent(in) :: nrhs  !! Number of right-hand sides.
      integer, intent(in) :: ldab  !! Leading dimension of `ab`.
      integer, intent(in) :: ldb  !! Leading dimension of `b`.
      integer, intent(out) :: ipiv(*)  !! LAPACK pivot indices.
      real(wp), intent(inout) :: ab(ldab, *)  !! Banded coefficient matrix (overwritten by LU factors).
      real(wp), intent(inout) :: b(ldb, *)  !! Right-hand side(s), overwritten by solution(s).
      integer, intent(out) :: info  !! LAPACK status code.
    end subroutine dgbsv
  end interface

contains

  !> Initialize solver state, allocate work arrays, and register user callbacks.
  !!
  !! Configures dense or banded Jacobian storage, sets PETSc-like defaults,
  !! and optionally applies user-provided tolerances and stepping controls.
  subroutine PTCSolver_initialize(self, x0, f, jac, jacobian_type, dt0, dt0_guess_fac, kl, ku, fatol, frtol, dt_increment, dt_max, increment_dt_from_initial_dt, max_reject, max_steps, &
      enable_stagnation_check, stagnation_warmup_steps, stagnation_window, stagnation_rel_improve_tol)
    class(PTCSolver), intent(inout) :: self  !! Solver object to initialize.
    real(wp), intent(in) :: x0(:)  !! Initial state guess.
    procedure(rhs_fcn) :: f  !! User residual callback.
    procedure(jac_fcn) :: jac  !! User Jacobian callback.
    integer, intent(in) :: jacobian_type  !! Jacobian mode (`PTC_JAC_DENSE` or `PTC_JAC_BAND`).
    real(wp), intent(in), optional :: dt0  !! Initial pseudo-time step. If omitted, estimated from Jacobian diagonal.
    real(wp), intent(in), optional :: dt0_guess_fac  !! Scale factor for auto `dt0` guess; used only when `dt0` is not provided.
    integer, intent(in), optional :: kl  !! Number of sub-diagonals for banded Jacobian mode.
    integer, intent(in), optional :: ku  !! Number of super-diagonals for banded Jacobian mode.
    integer, intent(in), optional :: max_reject  !! Maximum rejections allowed per `step()` call.
    integer, intent(in), optional :: max_steps  !! Maximum accepted steps allowed in `solve()`.
    logical, intent(in), optional :: enable_stagnation_check  !! Enable generic residual stagnation detection.
    integer, intent(in), optional :: stagnation_warmup_steps  !! Accepted-step warmup before enabling stagnation checks.
    integer, intent(in), optional :: stagnation_window  !! Consecutive non-improving accepted steps allowed.
    real(wp), intent(in), optional :: stagnation_rel_improve_tol  !! Relative improvement threshold considered progress.
    real(wp), intent(in), optional :: fatol  !! Absolute residual-norm convergence tolerance.
    real(wp), intent(in), optional :: frtol  !! Relative residual-norm convergence tolerance.
    real(wp), intent(in), optional :: dt_increment  !! Default timestep growth factor.
    real(wp), intent(in), optional :: dt_max  !! Maximum allowed timestep (non-positive means no cap).
    logical, intent(in), optional :: increment_dt_from_initial_dt  !! Optional switch for initial-reference dt adaptation.
    integer :: i, ierr
    real(wp) :: maxdiag, dt0_guess_fac_

    call reset_storage(self)

    self%reason = PTC_REASON_NONE
    if (size(x0) <= 0) then
      self%reason = PTC_DIVERGED_INVALID_INPUT
      return
    end if
    if (present(dt0)) then
      if (dt0 <= 0.0_wp) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
    end if
    if (present(dt0_guess_fac)) then
      if (dt0_guess_fac <= 0.0_wp) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
    end if
    if (present(dt_increment)) then
      if (dt_increment <= 0.0_wp) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
    end if
    if (present(stagnation_warmup_steps)) then
      if (stagnation_warmup_steps < 0) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
    end if
    if (present(stagnation_window)) then
      if (stagnation_window <= 0) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
    end if
    if (present(stagnation_rel_improve_tol)) then
      if (stagnation_rel_improve_tol < 0.0_wp) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
    end if

    self%neq = size(x0)
    self%dt = 0.0_wp
    self%dt_initial = 0.0_wp
    self%steps = 0
    self%rejects_total = 0

    self%f => f

    if (present(dt_increment)) self%dt_increment = dt_increment
    if (present(dt_max)) self%dt_max = dt_max
    if (present(fatol)) self%fatol = fatol
    if (present(frtol)) self%frtol = frtol
    if (present(increment_dt_from_initial_dt)) self%increment_dt_from_initial_dt = increment_dt_from_initial_dt
    if (present(max_reject)) self%max_reject = max_reject
    if (present(max_steps)) self%max_steps = max_steps
    if (present(enable_stagnation_check)) self%enable_stagnation_check = enable_stagnation_check
    if (present(stagnation_warmup_steps)) self%stagnation_warmup_steps = stagnation_warmup_steps
    if (present(stagnation_window)) self%stagnation_window = stagnation_window
    if (present(stagnation_rel_improve_tol)) self%stagnation_rel_improve_tol = stagnation_rel_improve_tol

    self%fnorm = -1.0_wp
    self%fnorm_initial = -1.0_wp
    self%fnorm_previous = -1.0_wp
    self%residual_valid = .false.
    self%fnorm_old = -1.0_wp
    self%residual_old_valid = .false.
    self%jac_valid = .false.
    self%jac_old_valid = .false.
    self%stagnation_count = 0
    self%fnorm_best = huge(1.0_wp)

    allocate(self%x(self%neq), self%x_old(self%neq), self%fvec(self%neq), self%fvec_old(self%neq), self%step_vec(self%neq), self%rhs_mat(self%neq, 1), self%ipiv(self%neq))
    self%x = x0

    self%jacobian_type = jacobian_type
    self%jac => jac

    select case (self%jacobian_type)
    case (PTC_JAC_DENSE)
      allocate(self%jac_mat(self%neq, self%neq), self%jac_mat_old(self%neq, self%neq), self%a_dense(self%neq, self%neq))

    case (PTC_JAC_BAND)
      if (.not. present(kl) .or. .not. present(ku)) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
      if (kl < 0 .or. ku < 0) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
      self%kl = kl
      self%ku = ku
      self%ldab = 2 * self%kl + self%ku + 1
      allocate(self%jac_mat(self%kl + self%ku + 1, self%neq), self%jac_mat_old(self%kl + self%ku + 1, self%neq), self%a_band(self%ldab, self%neq))

    case default
      self%reason = PTC_DIVERGED_INVALID_INPUT
      return
    end select

    if (present(dt0)) then
      self%dt = dt0
      self%dt_initial = dt0
    else
      call self%jac(self, self%x, self%jac_mat, ierr)
      if (ierr /= 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        return
      end if

      maxdiag = 0.0_wp
      select case (self%jacobian_type)
      case (PTC_JAC_DENSE)
        do i = 1, self%neq
          maxdiag = max(maxdiag, abs(self%jac_mat(i, i)))
        end do
      case (PTC_JAC_BAND)
        do i = 1, self%neq
          maxdiag = max(maxdiag, abs(self%jac_mat(self%ku + 1, i)))
        end do
      end select

      dt0_guess_fac_ = 0.1_wp
      if (present(dt0_guess_fac)) dt0_guess_fac_ = dt0_guess_fac

      self%dt = min(dt0_guess_fac_ / max(maxdiag, tiny(1.0_wp)), 1.0e12_wp)
      self%dt_initial = self%dt
      self%jac_valid = .true.
    end if

    self%initialized = .true.
  end subroutine PTCSolver_initialize

  !> Register a callback to verify/possibly reject each candidate pseudo-step.
  subroutine PTCSolver_set_verify_timestep(self, verify)
    class(PTCSolver), intent(inout) :: self  !! Solver object to update.
    procedure(verify_step_fcn) :: verify  !! User callback for accept/reject decisions.

    self%verify => verify
  end subroutine PTCSolver_set_verify_timestep

  !> Register a callback that overrides default pseudo-time-step adaptation.
  subroutine PTCSolver_set_compute_timestep(self, compute_dt)
    class(PTCSolver), intent(inout) :: self  !! Solver object to update.
    procedure(timestep_fcn) :: compute_dt  !! User callback for computing next timestep.

    self%compute_dt => compute_dt
  end subroutine PTCSolver_set_compute_timestep

  !> Register a callback that overrides built-in convergence checks.
  !!
  !! When set, this callback is the only convergence criterion used.
  subroutine PTCSolver_set_custom_convergence(self, convergence)
    class(PTCSolver), intent(inout) :: self  !! Solver object to update.
    procedure(convergence_fcn) :: convergence  !! User callback deciding convergence.

    self%custom_convergence => convergence
  end subroutine PTCSolver_set_custom_convergence

  !> Register a callback for progress reporting.
  subroutine PTCSolver_set_progress(self, progress)
    class(PTCSolver), intent(inout) :: self  !! Solver object to update.
    procedure(progress_fcn) :: progress  !! User callback for per-step reporting.

    self%progress => progress
  end subroutine PTCSolver_set_progress

  !> Advance the solver by one accepted pseudo-step (with internal retries).
  !!
  !! Performs linearized PTC update, optional verify callback, default/custom
  !! timestep adaptation, and convergence checks.
  subroutine PTCSolver_step(self)
    class(PTCSolver), intent(inout) :: self  !! Solver object advanced by one accepted pseudo-step.

    integer :: ierr, rejections
    real(wp) :: next_dt, reject_dt
    logical :: accept

    if (.not. self%initialized) then
      self%reason = PTC_DIVERGED_NOT_INITIALIZED
      return
    end if

    if (self%reason /= PTC_REASON_NONE) return

    if (self%steps == 0) self%dt_initial = self%dt

    call PTCSolver_check_step0_convergence(self)
    if (self%reason /= PTC_REASON_NONE) return

    rejections = 0

    do
      self%x_old = self%x
      if (self%residual_valid) then
        self%fvec_old = self%fvec
        self%fnorm_old = self%fnorm
        self%residual_old_valid = .true.
      else
        self%residual_old_valid = .false.
      end if
      if (self%jac_valid) then
        self%jac_mat_old = self%jac_mat
        self%jac_old_valid = .true.
      else
        self%jac_old_valid = .false.
      end if
      call PTCSolver_take_newton_update(self, ierr)
      if (ierr < 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        self%x = self%x_old
        return
      end if

      if (ierr > 0) then
        reject_dt = max(0.5_wp * self%dt, tiny(1.0_wp))
        call PTCSolver_reject_step(self, rejections, reject_dt)
        if (self%reason /= PTC_REASON_NONE) return
        cycle
      end if

      accept = .true.
      reject_dt = self%dt
      if (associated(self%verify)) then
        call self%verify(self, accept, reject_dt, ierr)
        if (ierr < 0) then
          self%reason = PTC_DIVERGED_CALLBACK_FATAL
          self%x = self%x_old
          return
        end if
        if (ierr > 0) accept = .false.
      end if

      if (.not. accept) then
        call PTCSolver_reject_step(self, rejections, max(reject_dt, tiny(1.0_wp)))
        if (self%reason /= PTC_REASON_NONE) return
        cycle
      end if

      call PTCSolver_compute_residual(self, self%x, self%fvec, self%fnorm, ierr)
      if (ierr < 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        self%x = self%x_old
        self%residual_valid = .false.
        return
      end if
      if (ierr > 0) then
        reject_dt = max(0.5_wp * self%dt, tiny(1.0_wp))
        call PTCSolver_reject_step(self, rejections, reject_dt)
        if (self%reason /= PTC_REASON_NONE) return
        cycle
      end if
      self%residual_valid = .true.

      if (self%fnorm_initial < 0.0_wp) then
        self%fnorm_initial = self%fnorm
        self%fnorm_previous = self%fnorm
      end if

      call PTCSolver_compute_next_dt(self, next_dt, ierr)
      if (ierr < 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        self%x = self%x_old
        return
      end if
      if (ierr > 0) then
        reject_dt = max(0.5_wp * self%dt, tiny(1.0_wp))
        call PTCSolver_reject_step(self, rejections, reject_dt)
        if (self%reason /= PTC_REASON_NONE) return
        cycle
      end if

      self%dt = next_dt
      self%fnorm_previous = self%fnorm
      self%steps = self%steps + 1
      call self%update_stagnation()
      if (associated(self%progress)) call self%progress(self)

      call PTCSolver_check_convergence(self)
      return
    end do
  end subroutine PTCSolver_step

  !> For step 0, ensure residual/progress/convergence are current before any timestep.
  subroutine PTCSolver_check_step0_convergence(self)
    class(PTCSolver), intent(inout) :: self
    integer :: ierr

    if (self%steps /= 0) return

    if (.not. self%residual_valid) then
      call PTCSolver_compute_residual(self, self%x, self%fvec, self%fnorm, ierr)
      if (ierr < 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        self%residual_valid = .false.
        return
      end if
      if (ierr > 0) then
        self%residual_valid = .false.
        return
      end if
      self%residual_valid = .true.
    end if

    if (self%fnorm_initial < 0.0_wp) then
      self%fnorm_initial = self%fnorm
      self%fnorm_previous = self%fnorm
    end if

    if (.not. self%progress_step0_emitted) then
      if (associated(self%progress)) call self%progress(self)
      self%progress_step0_emitted = .true.
    end if

    call PTCSolver_check_convergence(self)
  end subroutine PTCSolver_check_step0_convergence

  !> Repeatedly call `step()` until convergence or a terminal failure reason.
  subroutine PTCSolver_solve(self)
    class(PTCSolver), intent(inout) :: self  !! Solver object advanced until termination.

    if (.not. self%initialized) then
      self%reason = PTC_DIVERGED_NOT_INITIALIZED
      return
    end if

    do while (self%reason == PTC_REASON_NONE)
      if (self%steps >= self%max_steps) then
        self%reason = PTC_DIVERGED_MAX_STEPS
        exit
      end if
      call self%step()
    end do
  end subroutine PTCSolver_solve

  !> Evaluate stopping criteria.
  !!
  !! If `custom_convergence` is set, it fully overrides the built-in
  !! `fatol/frtol` checks.
  subroutine PTCSolver_check_convergence(self)
    class(PTCSolver), intent(inout) :: self  !! Solver object whose residual norms are tested.
    logical :: converged
    integer :: ierr

    if (self%enable_stagnation_check) then
      if (self%steps >= self%stagnation_warmup_steps) then
        if (self%stagnation_count >= self%stagnation_window) then
          self%reason = PTC_DIVERGED_STAGNATION
          return
        end if
      end if
    end if

    if (associated(self%custom_convergence)) then
      call self%custom_convergence(self, converged, ierr)
      if (ierr /= 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        return
      end if

      if (converged) then
        self%reason = PTC_CONVERGED_USER
      else
        self%reason = PTC_REASON_NONE
      end if
      return
    end if

    if (self%fnorm < self%fatol) then
      self%reason = PTC_CONVERGED_PSEUDO_FATOL
      return
    end if

    if (self%fnorm_initial > 0.0_wp) then
      if ((self%fnorm / self%fnorm_initial) < self%frtol) then
        self%reason = PTC_CONVERGED_PSEUDO_FRTOL
        return
      end if
    end if

    self%reason = PTC_REASON_NONE
  end subroutine PTCSolver_check_convergence

  !> Update generic stagnation monitor using residual-norm progress.
  subroutine PTCSolver_update_stagnation(self)
    class(PTCSolver), intent(inout) :: self
    real(wp) :: improve_factor

    if (.not. self%enable_stagnation_check) return
    if (self%steps < self%stagnation_warmup_steps) return
    if (self%fnorm < 0.0_wp) return

    if (self%fnorm_best == huge(1.0_wp)) then
      self%fnorm_best = self%fnorm
      self%stagnation_count = 0
      return
    end if

    improve_factor = 1.0_wp - max(0.0_wp, self%stagnation_rel_improve_tol)
    if (self%fnorm < self%fnorm_best*improve_factor) then
      self%fnorm_best = self%fnorm
      self%stagnation_count = 0
    else
      self%stagnation_count = self%stagnation_count + 1
    end if
  end subroutine PTCSolver_update_stagnation

  !> Compute residual vector and its 2-norm at a given state.
  subroutine PTCSolver_compute_residual(self, x, fvec, fnorm, ierr)
    class(PTCSolver), intent(inout) :: self  !! Solver object providing residual callback.
    real(wp), intent(in) :: x(:)  !! State at which to evaluate residual.
    real(wp), intent(out) :: fvec(:)  !! Residual vector `f(x)`.
    real(wp), intent(out) :: fnorm  !! Euclidean norm of residual vector.
    integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).

    call self%f(self, x, fvec, ierr)
    if (ierr /= 0) then
      fnorm = -1.0_wp
      return
    end if

    fnorm = norm2(fvec)
  end subroutine PTCSolver_compute_residual

  !> Perform one linearized pseudo-transient correction solve and state update.
  !!
  !! Solves `(I/dt - J) s = f(x)` using dense (`dgesv`) or banded (`dgbsv`)
  !! LAPACK routines, then updates `x <- x + s`.
  subroutine PTCSolver_take_newton_update(self, ierr)
    class(PTCSolver), intent(inout) :: self  !! Solver object updated in place.
    integer, intent(out) :: ierr  !! Status (`0` success, positive rejectable failure, negative fatal failure).

    integer :: info, i, j, row_compact, row_solve
    real(wp) :: inv_dt

    ierr = 0
    if (self%dt <= 0.0_wp) then
      ierr = 1
      return
    end if

    if (.not. self%residual_valid) then
      call PTCSolver_compute_residual(self, self%x, self%fvec, self%fnorm, ierr)
      if (ierr /= 0) then
        self%residual_valid = .false.
        return
      end if
      self%residual_valid = .true.
    end if
    if (self%steps == 0 .and. .not. self%progress_step0_emitted) then
      if (associated(self%progress)) call self%progress(self)
      self%progress_step0_emitted = .true.
    end if

    inv_dt = 1.0_wp / self%dt

    if (.not. self%jac_valid) then
      call self%jac(self, self%x, self%jac_mat, ierr)
      if (ierr /= 0) return
      self%jac_valid = .true.
    end if

    select case (self%jacobian_type)
    case (PTC_JAC_DENSE)
      self%a_dense = -self%jac_mat
      do i = 1, self%neq
        self%a_dense(i, i) = self%a_dense(i, i) + inv_dt
      end do

      self%rhs_mat(:, 1) = self%fvec
      call dgesv(self%neq, 1, self%a_dense, self%neq, self%ipiv, self%rhs_mat, self%neq, info)
      if (info /= 0) then
        ierr = 1
        return
      end if

      self%step_vec = self%rhs_mat(:, 1)
      self%x = self%x + self%step_vec
      self%residual_valid = .false.
      self%jac_valid = .false.

    case (PTC_JAC_BAND)
      self%a_band = 0.0_wp
      do j = 1, self%neq
        do i = max(1, j - self%ku), min(self%neq, j + self%kl)
          row_compact = self%ku + 1 + i - j
          row_solve = self%kl + row_compact
          self%a_band(row_solve, j) = -self%jac_mat(row_compact, j)
        end do
        self%a_band(self%kl + self%ku + 1, j) = self%a_band(self%kl + self%ku + 1, j) + inv_dt
      end do

      self%rhs_mat(:, 1) = self%fvec
      call dgbsv(self%neq, self%kl, self%ku, 1, self%a_band, self%ldab, self%ipiv, self%rhs_mat, self%neq, info)
      if (info /= 0) then
        ierr = 1
        return
      end if

      self%step_vec = self%rhs_mat(:, 1)
      self%x = self%x + self%step_vec
      self%residual_valid = .false.
      self%jac_valid = .false.

    case default
      ierr = -1
    end select
  end subroutine PTCSolver_take_newton_update

  !> Compute the next pseudo-time step.
  !!
  !! Uses either user callback `compute_dt(self, new_dt, ierr)` or the
  !! default TSPSEUDO-style residual-ratio formula with optional cap `dt_max`.
  subroutine PTCSolver_compute_next_dt(self, next_dt, ierr)
    class(PTCSolver), intent(inout) :: self  !! Solver object providing residual history and controls.
    real(wp), intent(out) :: next_dt  !! Computed timestep for next accepted step.
    integer, intent(out) :: ierr  !! Status (`0` success, nonzero failure).

    ierr = 0

    if (associated(self%compute_dt)) then
      call self%compute_dt(self, next_dt, ierr)
      if (ierr /= 0) return
    else
      if (self%fnorm == 0.0_wp) then
        next_dt = 1.0e12_wp * self%dt_increment * self%dt
      else if (self%increment_dt_from_initial_dt) then
        next_dt = self%dt_increment * self%dt_initial * self%fnorm_initial / self%fnorm
      else
        next_dt = self%dt_increment * self%dt * self%fnorm_previous / self%fnorm
      end if
      if (self%dt_max > 0.0_wp) next_dt = min(next_dt, self%dt_max)
    end if

    if (next_dt <= 0.0_wp) then
      ierr = 1
      return
    end if
  end subroutine PTCSolver_compute_next_dt

  !> Reject current attempt, rollback state, update `dt`, and count rejection.
  subroutine PTCSolver_reject_step(self, rejections, new_dt)
    class(PTCSolver), intent(inout) :: self  !! Solver object to rollback and update.
    integer, intent(inout) :: rejections  !! Rejection counter for current `step()` call.
    real(wp), intent(in) :: new_dt  !! Timestep to use after rejection.

    self%x = self%x_old
    self%dt = max(new_dt, tiny(1.0_wp))
    if (self%residual_old_valid) then
      self%fvec = self%fvec_old
      self%fnorm = self%fnorm_old
      self%residual_valid = .true.
    else
      self%residual_valid = .false.
    end if
    if (self%jac_old_valid) then
      self%jac_mat = self%jac_mat_old
      self%jac_valid = .true.
    else
      self%jac_valid = .false.
    end if
    self%rejects_total = self%rejects_total + 1

    rejections = rejections + 1
    if (rejections > self%max_reject) then
      self%reason = PTC_DIVERGED_STEP_REJECTED
    end if
  end subroutine PTCSolver_reject_step

  !> Deallocate solver work arrays and clear callback pointers.
  subroutine reset_storage(self)
    class(PTCSolver), intent(inout) :: self  !! Solver object whose allocations/pointers are cleared.

    if (allocated(self%x)) deallocate(self%x)
    if (allocated(self%x_old)) deallocate(self%x_old)
    if (allocated(self%fvec)) deallocate(self%fvec)
    if (allocated(self%fvec_old)) deallocate(self%fvec_old)
    if (allocated(self%step_vec)) deallocate(self%step_vec)
    if (allocated(self%rhs_mat)) deallocate(self%rhs_mat)
    if (allocated(self%jac_mat)) deallocate(self%jac_mat)
    if (allocated(self%jac_mat_old)) deallocate(self%jac_mat_old)
    if (allocated(self%a_dense)) deallocate(self%a_dense)
    if (allocated(self%a_band)) deallocate(self%a_band)
    if (allocated(self%ipiv)) deallocate(self%ipiv)

    self%initialized = .false.
    self%residual_valid = .false.
    self%residual_old_valid = .false.
    self%jac_valid = .false.
    self%jac_old_valid = .false.
    self%f => null()
    self%jac => null()
    self%verify => null()
    self%compute_dt => null()
    self%custom_convergence => null()
    self%progress => null()
    self%progress_step0_emitted = .false.
  end subroutine reset_storage

end module clima_ptc
