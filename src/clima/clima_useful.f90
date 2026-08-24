module clima_useful
  use clima_const, only: dp
  use minpack_module, only: fcn_hybrj
  use dop853_module, only: dop853_class
  use h5fortran, only: hdf5_file
  implicit none
  private

  public :: hdf5_file_closer
  public :: MinpackHybrd1Vars, MinpackHybrj
  public :: linear_solve
  public :: solve_ivp_dop853

  abstract interface
    subroutine dop853_rhs_fcn(x, y, dy)
      import :: dp
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: dy(:)
    end subroutine
  end interface

  type, extends(dop853_class) :: dop853_solveivp
    procedure(dop853_rhs_fcn), pointer, nopass :: rhs => null()
    integer :: j = 1
    integer :: neq = 0
    integer :: neval = 0
    real(dp), pointer :: x_eval(:) => null()
    real(dp), pointer :: y_eval(:,:) => null()
    real(dp) :: xtol = 0.0_dp
  end type

  ! Helps close hdf5 files
  type :: hdf5_file_closer
    type(hdf5_file), pointer :: h => null()
  contains
    final :: hdf5_file_closer_final
  end type

  type :: MinpackHybrd1Vars
    integer :: n
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: fvec(:)
    real(dp) :: tol = 1.0e-8_dp
    integer :: info
    real(dp), allocatable :: wa(:)
    integer :: lwa
  end type
  interface MinpackHybrd1Vars
    procedure :: create_MinpackHybrd1Vars
  end interface

  type :: MinpackHybrj
    procedure(fcn_hybrj), nopass, pointer :: fcn => null()
    integer :: n
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: fvec(:)
    real(dp), allocatable :: fjac(:,:)
    integer :: Ldfjac
    real(dp) :: xtol
    integer :: maxfev
    real(dp), allocatable :: diag(:)
    integer :: mode
    real(dp) :: factor
    integer :: nprint
    integer :: info
    integer :: nfev
    integer :: njev
    real(dp), allocatable :: r(:)
    integer :: Lr
    real(dp), allocatable :: qtf(:)
    real(dp), allocatable :: wa1(:), wa2(:), wa3(:), wa4(:)
    integer :: lwa
  contains
    procedure :: hybrj => MinpackHybrj_hybrj
    procedure :: code_to_message => MinpackHybrj_code_to_message
  end type
  interface MinpackHybrj
    procedure :: create_MinpackHybrj
  end interface

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      real(8), intent(inout) :: a(lda,n)
      integer, intent(in) :: lda
      integer, intent(in) :: ipiv(n)
      real(8), intent(inout) :: b(ldb,nrhs)
      integer, intent(in) :: ldb
      integer, intent(out) :: info
    end subroutine
  end interface

contains

  subroutine solve_ivp_dop853(rhs, x_start, y_start, x_eval, rtol, atol, y_eval, idid, err)
    procedure(dop853_rhs_fcn) :: rhs
    real(dp), intent(in) :: x_start
    real(dp), intent(in) :: y_start(:)
    real(dp), intent(in), target :: x_eval(:)
    real(dp), intent(in) :: rtol, atol
    real(dp), intent(out), target :: y_eval(:,:)
    integer, intent(out) :: idid
    character(:), allocatable, intent(out) :: err

    type(dop853_solveivp) :: dop
    logical :: status_ok
    integer :: i, n, m
    integer, allocatable :: icomp(:)
    real(dp) :: x0, x_end
    real(dp), allocatable :: y0(:)

    n = size(y_start)
    m = size(x_eval)

    if (size(y_eval,1) /= n .or. size(y_eval,2) /= m) then
      err = 'solve_ivp_dop853: y_eval has wrong shape.'
      return
    endif
    if (m < 1) then
      err = 'solve_ivp_dop853: x_eval must have size >= 1.'
      return
    endif

    x_end = x_eval(m)
    if (x_end >= x_start) then
      do i = 2,m
        if (x_eval(i) < x_eval(i-1)) then
          err = 'solve_ivp_dop853: x_eval must be nondecreasing.'
          return
        endif
      enddo
      if (x_eval(1) < x_start) then
        err = 'solve_ivp_dop853: x_eval(1) must be >= x_start for forward integration.'
        return
      endif
    else
      do i = 2,m
        if (x_eval(i) > x_eval(i-1)) then
          err = 'solve_ivp_dop853: x_eval must be nonincreasing.'
          return
        endif
      enddo
      if (x_eval(1) > x_start) then
        err = 'solve_ivp_dop853: x_eval(1) must be <= x_start for backward integration.'
        return
      endif
    endif

    dop%xtol = max(1.0e-12_dp*max(1.0_dp, abs(x_start), maxval(abs(x_eval))), &
                   10.0_dp*spacing(max(1.0_dp, abs(x_start), maxval(abs(x_eval)))))

    y_eval = 0.0_dp
    dop%j = 1
    do while (dop%j <= m .and. abs(x_eval(dop%j)-x_start) <= dop%xtol)
      y_eval(:,dop%j) = y_start
      dop%j = dop%j + 1
    enddo
    if (dop%j > m) then
      idid = 1
      return
    endif

    allocate(icomp(n))
    allocate(y0(n))
    y0 = y_start
    do i = 1,n
      icomp(i) = i
    enddo

    call dop%initialize(fcn=rhs_dop853_solveivp, solout=solout_dop853_solveivp, n=n, &
                        iprint=0, icomp=icomp, status_ok=status_ok)
    if (.not.status_ok) then
      err = 'solve_ivp_dop853: failed to initialize dop853.'
      return
    endif

    dop%rhs => rhs
    dop%neq = n
    dop%neval = m
    dop%x_eval => x_eval
    dop%y_eval => y_eval

    x0 = x_start
    call dop%integrate(x0, y0, x_end, [rtol], [atol], iout=2, idid=idid)
    if (idid < 0) then
      err = 'solve_ivp_dop853: dop853 integration failed.'
      return
    endif

    if (dop%j <= m) then
      err = 'solve_ivp_dop853: integration finished before sampling all output points.'
      idid = -1
      return
    endif

  end subroutine

  subroutine rhs_dop853_solveivp(self, x, y, dy)
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: x
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: dy(:)

    select type (self)
    class is (dop853_solveivp)
      call self%rhs(x, y, dy)
    end select
  end subroutine

  subroutine solout_dop853_solveivp(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in) :: nr
    real(dp),intent(in) :: xold
    real(dp),intent(in) :: x
    real(dp),dimension(:),intent(in) :: y
    integer,intent(inout) :: irtrn
    real(dp),intent(out) :: xout

    real(dp) :: x_lo, x_hi, xe
    integer :: i

    xout = x
    select type (self)
    class is (dop853_solveivp)
      x_lo = min(xold, x) - self%xtol
      x_hi = max(xold, x) + self%xtol
      do while (self%j <= self%neval)
        xe = self%x_eval(self%j)
        if (xe < x_lo .or. xe > x_hi) exit
        self%y_eval(:,self%j) = [(self%contd8(i, xe), i=1,self%neq)]
        self%j = self%j + 1
      enddo
    end select
  end subroutine

  subroutine hdf5_file_closer_final(self)
    type(hdf5_file_closer), intent(inout) :: self
    if (associated(self%h)) then
      if (self%h%is_open) call self%h%close()
    endif
  end subroutine

  function create_MinpackHybrd1Vars(n, tol) result(v)
    integer, intent(in) :: n
    real(dp), optional, intent(in) :: tol
    type(MinpackHybrd1Vars) :: v
    v%n = n
    allocate(v%x(n))
    allocate(v%fvec(n))
    if (present(tol)) then
      v%tol = tol
    endif
    v%lwa = (n*(3*n+13))/2 + 1
    allocate(v%wa(v%lwa))
  end function

  function create_MinpackHybrj(fcn, n) result(v)
    procedure(fcn_hybrj) :: fcn
    integer, intent(in) :: n
    type(MinpackHybrj) :: v

    v%fcn => fcn
    v%n = n
    allocate(v%x(v%n))
    v%x = 1.0_dp
    allocate(v%fvec(v%n))
    allocate(v%fjac(v%n,v%n))
    v%ldfjac = v%n
    v%xtol = 1.0e-8_dp ! not default in hybrdj1
    v%maxfev = 100*(v%n + 1)
    allocate(v%diag(v%n))
    v%diag(:) = 1.0_dp
    v%mode = 2
    v%factor = 100.0_dp
    v%nprint = 0
    v%lr = (v%n*(v%n+1))/2
    allocate(v%r(v%lr))
    allocate(v%qtf(v%n))
    allocate(v%wa1(v%n), v%wa2(v%n), v%wa3(v%n), v%wa4(v%n))

  end function

  subroutine MinpackHybrj_hybrj(self)
    use minpack_module, only: hybrj
    class(MinpackHybrj), intent(inout) :: self

    if (.not.associated(self%fcn)) return

    call hybrj(self%fcn, self%n, self%x, self%Fvec, self%Fjac, &
               self%Ldfjac, self%Xtol, self%Maxfev, self%Diag, self%Mode, &
               self%Factor, self%Nprint, self%Info, self%Nfev, self%Njev, &
               self%r, self%Lr, self%Qtf, self%Wa1, self%Wa2, &
               self%Wa3, self%Wa4)
    
  end subroutine

  function MinpackHybrj_code_to_message(self, info) result(message)
    class(MinpackHybrj), intent(inout) :: self
    integer, intent(in) :: info
    character(:), allocatable :: message

    if (info < 0) then
      message = 'user terminated execution.'
    elseif (info == 0) then
      message = 'improper input parameters.'
    elseif (info == 1) then
      message = 'relative error between two consecutive iterates is '// &
      'at most xtol.'
    elseif (info == 2) then
      message = 'number of calls to fcn with iflag = 1 has reached maxfev.'
    elseif (info == 3) then
      message = 'xtol is too small. no further improvement in the '// &
      'approximate solution x is possible.'
    elseif (info == 4) then
      message = 'iteration is not making good progress, as measured '// &
      'by the improvement from the last five jacobian evaluations.'
    elseif (info == 5) then
      message = 'iteration is not making good progress, as measured '// &
      'by the improvement from the last ten iterations.'
    else
      message = 'unknown message'
    endif

  end function

  subroutine linear_solve(A, b, info)
    real(dp), intent(inout) :: a(:,:)
    real(dp), target, intent(inout) :: b(:)
    integer, intent(out) :: info

    integer :: n, nrhs, lda
    integer :: ipiv(size(b))
    integer :: ldb

    n = size(b)
    if (size(a,1) /= n .or. size(a,1) /= n) then
      info = -1
      return
    endif
    lda = n
    ldb = n
    nrhs = 1

    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
  end subroutine

end module
