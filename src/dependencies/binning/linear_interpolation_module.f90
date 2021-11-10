!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!  Multidimensional linear interpolation/extrapolation.
!
!  Uses repeated linear interpolation to evaluate
!  functions \(f(x), f(x,y), f(x,y,z), f(x,y,z,q), f(x,y,z,q,r), f(x,y,z,q,r,s) \)
!  which have been tabulated at the nodes of an n-dimensional rectangular grid.
!  If any coordinate \( (x_i, y_i, ...) \) lies outside the range of the corresponding
!  variable, then extrapolation is performed using the two nearest points.

    module linear_interpolation_module

    use iso_fortran_env,    only: wp => real64   ! working precision

    implicit none

    private

    real(wp),parameter,private :: zero = 0.0_wp  !! numeric constant
    real(wp),parameter,private :: one  = 1.0_wp  !! numeric constant

    type,public,abstract :: linear_interp_class
        !! Base class for the linear interpolation types
        private
        logical :: initialized = .false. !! if the class was properly initialized
    contains
        private
        procedure(destroy_func),deferred,public :: destroy  !! destructor
        procedure :: check_inputs
    end type linear_interp_class

    abstract interface
        pure elemental subroutine destroy_func(me)  !! interface for bspline destructor routines
        import :: linear_interp_class
        implicit none
        class(linear_interp_class),intent(inout) :: me
        end subroutine destroy_func
    end interface

    type,extends(linear_interp_class),public :: linear_interp_2d
        !! Class for 2d linear interpolation.
        private
        real(wp),dimension(:,:),allocatable :: f
        real(wp),dimension(:),allocatable :: x
        real(wp),dimension(:),allocatable :: y
        integer :: ilox = 1
        integer :: iloy = 1
        contains
        private
        procedure,public :: initialize => initialize_2d
        procedure,public :: evaluate   => interp_2d
        procedure,public :: destroy    => destroy_2d
        final :: finalize_2d
    end type linear_interp_2d
    
    type,extends(linear_interp_2d),public :: nearest_interp_2d
        !! Class for 2d nearest neighbor interpolation.
        contains
        procedure,public :: evaluate => nearest_2d
    end type nearest_interp_2d

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for a [[linear_interp_2d]] type.

    pure elemental subroutine finalize_2d(me)

    implicit none

    type(linear_interp_2d),intent(inout) :: me
    call me%destroy()

    end subroutine finalize_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for a [[linear_interp_2d]] class.

    pure elemental subroutine destroy_2d(me)

    implicit none

    class(linear_interp_2d),intent(inout) :: me

    if (allocated(me%f)) deallocate(me%f)
    if (allocated(me%x)) deallocate(me%x)
    if (allocated(me%y)) deallocate(me%y)
    me%ilox = 1
    me%iloy = 1
    me%initialized = .false.

    end subroutine destroy_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[linear_interp_2d]] class.

    pure subroutine initialize_2d(me,x,y,f,istat)

    implicit none

    class(linear_interp_2d),intent(inout) :: me
    real(wp),dimension(:),intent(in)      :: x
    real(wp),dimension(:),intent(in)      :: y
    real(wp),dimension(:,:),intent(in)    :: f
    integer,intent(out)                   :: istat  !! `0`   : no problems,
                                                    !! `1`   : `x` is not strictly increasing,
                                                    !! `2`   : `y` is not strictly increasing,
                                                    !! `10`  : `x` is not equal to size(f,1),
                                                    !! `20`  : `y` is not equal to size(f,2),
                                                    !! `100` : cannot use linear interpolation for only one point.

    call me%destroy()

    istat = 0

    if (istat==0 .and. size(x)/=size(f,1)) istat = 10
    if (istat==0 .and. size(y)/=size(f,2)) istat = 20

    if (istat==0) then
        call me%check_inputs(x=x,y=y,ierr=istat)
        if (istat==0) then
            allocate(me%f(size(x),size(y))); me%f = f
            allocate(me%x(size(x))); me%x = x
            allocate(me%y(size(y))); me%y = y
            me%initialized = .true.
        end if
    end if

    end subroutine initialize_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  2D linear interpolation routine.

    pure subroutine interp_2d(me,x,y,f,istat)

    implicit none

    class(linear_interp_2d),intent(inout) :: me
    real(wp),intent(in)                   :: x
    real(wp),intent(in)                   :: y
    real(wp),intent(out)                  :: f     !! Interpolated \( f(x,y) \)
    integer,intent(out),optional          :: istat !! `0`  : no problems,
                                                   !! `-1` : class has not been initialized

    integer,dimension(2) :: ix, iy
    real(wp) :: p1, p2
    real(wp) :: q1, q2
    integer :: mflag
    real(wp) :: fx1, fx2

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag)

        q1 = (x-me%x(ix(1)))/(me%x(ix(2))-me%x(ix(1)))
        q2 = (y-me%y(iy(1)))/(me%y(iy(2))-me%y(iy(1)))
        p1 = one-q1
        p2 = one-q2

        fx1 = p1*me%f(ix(1),iy(1)) + q1*me%f(ix(2),iy(1))
        fx2 = p1*me%f(ix(1),iy(2)) + q1*me%f(ix(2),iy(2))

        f = p2*fx1 + q2*fx2
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine interp_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the indices in `xt` that bound `x`, to use for interpolation.
!  If outside the range, then the indices are returned that can
!  be used for extrapolation.
!  Precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   iright=2,    mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   iright=i+1,  mflag=0
!         if   xt(n) <= x           then ileft=n-1, iright=n,    mflag=1
!```
!
!### History
!
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.
!  * Jacob Williams, 2/17/2016 : additional refactoring (eliminated GOTOs).
!  * Jacob Williams, 2/22/2016 : modified bspline-fortran `dintrv` routine for
!    linear interpolation/extrapolation use.
!  * Jacob Williams, 10/9/2019 : added optional `inearest` output.

    pure subroutine dintrv(xt,x,ilo,ileft,iright,mflag,inearest)

    implicit none

    real(wp),dimension(:),intent(in) :: xt       !! a knot or break point vector
    real(wp),intent(in)              :: x        !! argument
    integer,intent(inout)            :: ilo      !! an initialization parameter which must be set
                                                 !! to 1 the first time the array `xt` is
                                                 !! processed by dintrv. `ilo` contains information for
                                                 !! efficient processing after the initial call and `ilo`
                                                 !! must not be changed by the user.  each dimension
                                                 !! requires a distinct `ilo` parameter.
    integer,intent(out)              :: ileft    !! left index
    integer,intent(out)              :: iright   !! right index
    integer,intent(out)              :: mflag    !! signals when `x` lies out of bounds
    integer,intent(out),optional     :: inearest !! nearest index

    integer :: ihi, istep, imid, n

    n = size(xt)

    if (n==1) then
        ! this is only allowed for nearest interpolation
        if (present(inearest)) then
            inearest = 1
            return
        end if
    end if

    ihi = ilo + 1
    if ( ihi>=n ) then
        if ( x>=xt(n) ) then
            mflag = 1
            ileft = n-1
            iright= n
            if (present(inearest)) inearest = n
            return
        end if
        if ( n<=1 ) then
            mflag = -1
            ileft = 1
            iright= 2
            if (present(inearest)) inearest = 1
            return
        end if
        ilo = n - 1
        ihi = n
    endif

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=n ) then
                if ( x>=xt(n) ) then
                    mflag = 1
                    ileft = n-1
                    iright= n
                    if (present(inearest)) inearest = n
                    return
                end if
                ihi = n
            elseif ( x>=xt(ihi) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0
            ileft = ilo
            iright= ilo+1
            if (present(inearest)) then
                if ( abs(x-xt(ileft)) <= abs(x-xt(iright)) ) then
                    inearest = ileft
                else
                    inearest = iright
                end if
            end if
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1 ) then
                ilo = 1
                if ( x<xt(1) ) then
                    mflag = -1
                    ileft = 1
                    iright= 2
                    if (present(inearest)) inearest = 1
                    return
                end if
            elseif ( x<xt(ilo) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    endif

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        imid = (ilo+ihi)/2
        if ( imid==ilo ) then
            mflag = 0
            ileft = ilo
            iright= ilo+1
            if (present(inearest)) then
                if ( abs(x-xt(ileft)) <= abs(x-xt(iright)) ) then
                    inearest = ileft
                else
                    inearest = iright
                end if
            end if
            return
        end if
        ! note. it is assumed that imid = ilo in case ihi = ilo+1
        if ( x<xt(imid) ) then
            ihi = imid
        else
            ilo = imid
        endif
    end do

    end subroutine dintrv
!*****************************************************************************************

!*****************************************************************************************
!>
!  2D nearest neighbor interpolation routine.

    pure subroutine nearest_2d(me,x,y,f,istat)

    implicit none

    class(nearest_interp_2d),intent(inout) :: me
    real(wp),intent(in)                    :: x
    real(wp),intent(in)                    :: y
    real(wp),intent(out)                   :: f     !! Nearest \( f(x,y) \)
    integer,intent(out),optional           :: istat !! `0`  : no problems,
                                                    !! `-1` : class has not been initialized

    integer :: mflag
    integer,dimension(2) :: ix, iy
    integer :: i, j

    if (me%initialized) then

        call dintrv(me%x,x,me%ilox,ix(1),ix(2),mflag,i)
        call dintrv(me%y,y,me%iloy,iy(1),iy(2),mflag,j)

        f = me%f(i,j)
        if (present(istat)) istat = 0

    else

        if (present(istat)) istat = -1
        f = zero

    end if

    end subroutine nearest_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the validity of the inputs to the initialize routines.
!  Prints warning message if there is an error,
!  and also sets `ierr` (/=0 if there were any errors).
!
!  Supports up to 6D: x,y,z,q,r,s
!
!# History
!  * Jacob Williams, 2/24/2015 : Created this routine.
!  * Jacob Williams, 2/23/2016 : modified for linear interp module.

    pure subroutine check_inputs(me,x,y,z,q,r,s,ierr)

    implicit none

    class(linear_interp_class),intent(in)      :: me
    real(wp),dimension(:),intent(in),optional  :: x     !! `x` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: y     !! `y` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: z     !! `z` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: q     !! `q` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: r     !! `r` abscissa vector
    real(wp),dimension(:),intent(in),optional  :: s     !! `s` abscissa vector
    integer,intent(out)                        :: ierr  !! `0`   : no problems,
                                                        !! `1`   : `x` is not strictly increasing,
                                                        !! `2`   : `y` is not strictly increasing,
                                                        !! `3`   : `z` is not strictly increasing,
                                                        !! `4`   : `q` is not strictly increasing,
                                                        !! `5`   : `r` is not strictly increasing,
                                                        !! `6`   : `s` is not strictly increasing,
                                                        !! `100` : cannot use linear interpolation for only one point.

    ierr = 0  ! initialize

    if (present(x)) call check(x,1,ierr); if (ierr/=0) return
    if (present(y)) call check(y,2,ierr); if (ierr/=0) return
    if (present(z)) call check(z,3,ierr); if (ierr/=0) return
    if (present(q)) call check(q,4,ierr); if (ierr/=0) return
    if (present(r)) call check(r,5,ierr); if (ierr/=0) return
    if (present(s)) call check(s,6,ierr); if (ierr/=0) return

    if (ierr == 0) then
        select type (me)
        class is (nearest_interp_2d)
        class default
            ! need at least two points for linear interpolation:
            if (size(x)==1) ierr = 100
        end select
    end if

    contains
!*****************************************************************************************

        pure subroutine check(v,error_code,ierr)

        implicit none

        real(wp),dimension(:),intent(in) :: v          !! abcissae vector
        integer,intent(in)               :: error_code !! error code for check
        integer,intent(inout)            :: ierr       !! will be set to `error_code` if there is a problem

        integer :: i  !! counter
        integer :: n  !! size of the input `v` array

        n = size(v)
        do i=2,n
            if (v(i) <= v(i-1)) then
                ierr = error_code
                exit
            end if
        end do

        end subroutine check

    end subroutine check_inputs
!*****************************************************************************************

!*****************************************************************************************
    end module linear_interpolation_module
!*****************************************************************************************
