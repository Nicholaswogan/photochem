! 1D linear interpolation
subroutine interp(ng, n, xg, x, y, yg, err)
  implicit none
  integer, intent(in) :: ng, n
  real(8), intent(in) :: xg(ng)
  real(8), intent(in) :: x(n), y(n)
  
  real(8), intent(out) :: yg(ng)
  character(len=100), intent(out) :: err 
  
  real(8) :: slope
  integer :: i, j, nn
  
  do i = 1,n-1
    if (x(i+1) < x(i)) then
      err = 'x must be sorted.'
      return
    endif
  enddo
  do i = 1,ng-1
    if (xg(i+1) < xg(i)) then
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