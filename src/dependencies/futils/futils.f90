
module futils
  implicit none
  
  integer, parameter, private :: real_kind = kind(1.d0)
  
  type Timer
    integer :: cr, cm, c1, c2
    real(real_kind) :: time
  contains
    procedure :: start => Timer_start
    procedure :: finish => Timer_finish
  end type
  
contains
  
  subroutine Timer_start(self)
    class(Timer), intent(inout) :: self
    call system_clock(count = self%c1, count_rate = self%cr, count_max = self%cm)
  end subroutine
  
  subroutine Timer_finish(self, msg)
    class(Timer), intent(inout) :: self
    character(len=*), optional, intent(in) :: msg
    call system_clock(count = self%c2)
    self%time = (self%c2-self%c1)/real(self%cr)
    if (present(msg)) then
      print"(A,1x,es10.4)",trim(msg),self%time
    endif
  end subroutine
  
  ! binning.f90 is same as binning.f 
  ! adding to module means explicit interface to function
  ! which allows for more compiler errors to be caught
  include "binning.f90"

  ! 1D linear interpolation with constant extrapolation.
  subroutine interp(ng, n, xg, x, y, yg, err)
    implicit none
    integer, intent(in) :: ng ! length of new grid (we interpolate to this new grid)
    integer, intent(in) :: n ! length of old grid
    real(real_kind), intent(in) :: xg(ng) ! new grid
    real(real_kind), intent(in) :: x(n), y(n) ! old data
    
    real(real_kind), intent(out) :: yg(ng) ! new data 
    character(len=100), intent(out) :: err 
    
    real(real_kind) :: slope
    integer :: i, j, nn
    
    do i = 1,n-1
      if (x(i+1) <= x(i)) then
        err = 'x must be sorted.'
        return
      endif
    enddo
    do i = 1,ng-1
      if (xg(i+1) <= xg(i)) then
        err = 'xg must be sorted.'
        return
      endif
    enddo
    
    nn = 1
    do i = 1,ng
      if (xg(i) < x(1)) then
        yg(i) = y(1)
      elseif ((xg(i) >= x(1)) .and. (xg(i) <= x(n))) then
        do j = nn,n
          if ((xg(i) >= x(j)) .and. (xg(i) <= x(j+1))) then
            slope = (y(j+1)-y(j))/(x(j+1)-x(j))
            yg(i) = y(j) + slope*(xg(i)-x(j))
            nn = j
            exit
          endif
        enddo
      elseif (xg(i) > x(n)) then
        yg(i) = y(n)
      endif
    enddo
    
  end subroutine
  
  pure recursive function replaceStr(string,search,substitute) result(modifiedString)
    character(len=*), intent(in)  :: string, search, substitute
    character(len=:), allocatable :: modifiedString
    integer                       :: i, stringLen, searchLen
    stringLen = len(string)
    searchLen = len(search)
    if (stringLen==0 .or. searchLen==0) then
      modifiedString = ""
      return
    elseif (stringLen<searchLen) then
      modifiedString = string
      return
    end if
    i = 1
    do
      if (string(i:i+searchLen-1)==search) then
        modifiedString = string(1:i-1) // substitute // replaceStr(string(i+searchLen:stringLen),search,substitute)
        exit
      end if
      if (i+searchLen>stringLen) then
        modifiedString = string
        exit
      end if
      i = i + 1
      cycle
    end do
  end function replaceStr
  
end module