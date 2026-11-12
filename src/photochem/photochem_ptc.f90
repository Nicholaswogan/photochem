! Experimental steady-state solver; pseudo-time is not a physical trajectory.
module photochem_ptc
  use photochem_const, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: PTCSolver

  abstract interface
    subroutine ptc_rhs(context, y, f, scale, err)
      import dp
      class(*), intent(inout) :: context
      real(dp), target, intent(in) :: y(:)
      real(dp), intent(out) :: f(:), scale(:)
      character(:), allocatable, intent(out) :: err
    end subroutine
    subroutine ptc_jac(context, y, jac, err)
      import dp
      class(*), intent(inout) :: context
      real(dp), target, intent(in) :: y(:)
      real(dp), target, intent(out) :: jac(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine
  end interface

  interface
    subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
      import dp
      integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
      integer, intent(out) :: ipiv(*), info
      real(dp), intent(inout) :: ab(ldab,*), b(ldb,*)
    end subroutine
  end interface

  !> Banded linearly implicit Euler PTC with change or residual adaptation.
  !! The caller supplies LAPACK extended-band Jacobians (diagonal kl+ku+1).
  !! Callback scales have state units and must be finite and strictly positive.
  type :: PTCSolver
    real(dp), allocatable :: y(:), f(:), scale(:), trial(:), ft(:), st(:)
    real(dp), allocatable :: jac(:,:), matrix(:,:), correction(:,:)
    integer, allocatable :: pivots(:)
    integer :: kl=0, ku=0, steps=0, rejects=0, rhs_calls=0, jac_calls=0, solves=0
    integer :: limiting_index=0, max_retries=100, clipped_components=0
    real(dp) :: maximum_clipped_scaled=0
    real(dp) :: time=0, dt=1.0e-6_dp, last_dt=0, change=0, residual=0
    real(dp) :: dt_max=1.0e30_dp, change_limit=0.25_dp
    logical :: residual_policy=.false., nonnegative=.false., initialized=.false.
  contains
    procedure :: initialize
    procedure :: step
  end type
contains
  subroutine initialize(self, y, kl, ku, initial_dt, err)
    class(PTCSolver), intent(inout) :: self
    real(dp), intent(in) :: y(:), initial_dt
    integer, intent(in) :: kl, ku
    character(:), allocatable, intent(out) :: err
    integer :: n, lda
    self%initialized=.false.
    if (size(y)<1.or.kl<0.or.ku<0.or.kl>=size(y).or.ku>=size(y).or. &
        .not.ieee_is_finite(initial_dt).or.initial_dt<=0.or..not.all(ieee_is_finite(y))) then
      err='Invalid PTC initial state, bandwidth, or timestep'; return
    endif
    if (allocated(self%y)) deallocate(self%y,self%f,self%scale,self%trial,self%ft,self%st, &
                                      self%jac,self%matrix,self%correction,self%pivots)
    n=size(y); lda=2*kl+ku+1
    allocate(self%y(n),self%f(n),self%scale(n),self%trial(n),self%ft(n),self%st(n), &
             self%jac(lda,n),self%matrix(lda,n),self%correction(n,1),self%pivots(n))
    self%y=y; self%kl=kl; self%ku=ku; self%dt=initial_dt
    self%time=0; self%steps=0; self%rejects=0
    self%last_dt=0; self%change=0; self%residual=0; self%limiting_index=0
    self%clipped_components=0; self%maximum_clipped_scaled=0
    self%rhs_calls=0; self%jac_calls=0; self%solves=0
    self%initialized=.true.
  end subroutine

  subroutine step(self, context, rhs, jacobian, err)
    class(PTCSolver), intent(inout) :: self
    class(*), intent(inout) :: context
    procedure(ptc_rhs) :: rhs
    procedure(ptc_jac) :: jacobian
    character(:), allocatable, intent(out) :: err
    character(:), allocatable :: trial_err, restore_err
    real(dp) :: h, factor, old_norm, new_norm
    integer :: attempt, i, info, n
    logical :: valid
    if (.not.self%initialized) then
      err='PTC is not initialized'; return
    endif
    if (.not.ieee_is_finite(self%dt).or.self%dt<=0.or. &
        .not.ieee_is_finite(self%dt_max).or.self%dt_max<=0.or. &
        .not.ieee_is_finite(self%change_limit).or.self%change_limit<=0) then
      err='Invalid PTC timestep controls'; return
    endif
    n=size(self%y)
    call rhs(context,self%y,self%f,self%scale,err); self%rhs_calls=self%rhs_calls+1
    if (allocated(err)) return
    if (.not.all(ieee_is_finite(self%f)).or..not.all(ieee_is_finite(self%scale)).or.any(self%scale<=0)) then
      err='Invalid PTC residual or scale'; return
    endif
    call jacobian(context,self%y,self%jac,err); self%jac_calls=self%jac_calls+1
    if (allocated(err)) return
    if (.not.all(ieee_is_finite(self%jac))) then
      err='Nonfinite PTC Jacobian'; return
    endif
    old_norm=maxval(abs(self%f)/self%scale)
    if (.not.ieee_is_finite(old_norm)) then
      err="Nonfinite scaled PTC residual"; return
    endif
    do attempt=1,self%max_retries
      h=min(self%dt,self%dt_max)
      if (self%time+h<=self%time.or.h<=tiny(h)) exit
      self%matrix=-self%jac
      do i=1,n
        self%matrix(self%kl+self%ku+1,i)=self%matrix(self%kl+self%ku+1,i)+1.0_dp/h
      enddo
      self%correction(:,1)=self%f
      call dgbsv(n,self%kl,self%ku,1,self%matrix,size(self%matrix,1),self%pivots,self%correction,n,info)
      self%solves=self%solves+1
      valid=info==0.and.all(ieee_is_finite(self%correction))
      if (valid) then
        self%limiting_index=maxloc(abs(self%correction(:,1))/self%scale,1)
        self%change=maxval(abs(self%correction(:,1))/self%scale)
        self%trial=self%y+self%correction(:,1)
        ! Trace negative values are tolerated only below 1e-4 of the scale.
        ! Only these roundoff-sized negative trials may be projected to zero.
        valid=all(ieee_is_finite(self%trial)).and.self%change<=self%change_limit
        if (self%nonnegative) valid=valid.and.all(self%trial>=-1e-4_dp*self%scale)
      endif
      if (valid) then
        if (self%nonnegative) self%trial=max(self%trial,0.0_dp)
        call rhs(context,self%trial,self%ft,self%st,trial_err); self%rhs_calls=self%rhs_calls+1
        valid=.not.allocated(trial_err)
        if (valid) valid=all(ieee_is_finite(self%ft)).and.all(ieee_is_finite(self%st)).and.all(self%st>0)
      endif
      if (valid) then
        ! Compare residuals with the same frozen scales, not moving weights.
        new_norm=maxval(abs(self%ft)/self%scale)
        valid=ieee_is_finite(new_norm)
      endif
      if (valid) then
        if (self%residual_policy) then
          if (new_norm<=old_norm/5.0_dp) then
            factor=5.0_dp
          else
            factor=max(0.5_dp,old_norm/new_norm)
          endif
        else
          factor=1.0_dp
          if(self%change<0.15_dp)factor=1.1_dp
          if(self%change<0.07_dp)factor=1.2_dp
          if(self%change<0.01_dp)factor=1.4_dp
          if(self%change<0.008_dp)factor=2.0_dp
          if(self%change<0.004_dp)factor=3.0_dp
          if(self%change<0.001_dp)factor=4.0_dp
          if(self%change<0.0005_dp)factor=5.0_dp
        endif
        if (self%nonnegative) then
          self%clipped_components=self%clipped_components+count(self%y+self%correction(:,1)<0)
          self%maximum_clipped_scaled=max(self%maximum_clipped_scaled, &
            maxval(max(-(self%y+self%correction(:,1)),0.0_dp)/self%scale))
        endif
        self%y=self%trial; self%f=self%ft; self%residual=new_norm
        self%time=self%time+h; self%last_dt=h; self%dt=min(self%dt_max,h*factor)
        self%steps=self%steps+1
        return
      endif
      if (allocated(trial_err)) deallocate(trial_err)
      ! Reprepare committed atmospheric state before retrying or returning.
      call rhs(context,self%y,self%f,self%st,restore_err); self%rhs_calls=self%rhs_calls+1
      if (allocated(restore_err)) then
        err='PTC rollback failed: '//restore_err; return
      endif
      self%rejects=self%rejects+1; self%dt=0.7_dp*h
    enddo
    err='PTC retries exhausted or pseudo-time cannot advance'
  end subroutine
end module
