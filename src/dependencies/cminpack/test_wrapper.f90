module test_wrapper
  use, intrinsic :: iso_c_binding, only : c_ptr, c_loc, c_f_pointer
  implicit none
  
  integer, parameter :: real_kind = kind(1.0d0)
  
  type f_data
    real(real_kind) :: a
    real(real_kind) :: b
  end type
  
contains
  
  integer function test_func1(p, n, x, fvec, iflag) result(res)
    type(c_ptr) :: p
    integer, value :: n
    real(real_kind), intent(in) :: x(n)
    real(real_kind), intent(out) :: fvec(n)
    integer, value :: iflag
    type(f_data), pointer :: params
    
    call c_f_pointer(p, params) ! dereference pointer
    fvec(1) = x(1)**2.d0 - (params%a + params%b)
    
    res = 0
  end function
  
  integer function test_func2(p, n, x, fvec, iflag) result(res)
    type(c_ptr) :: p
    integer, value :: n
    real(real_kind), intent(in) :: x(n)
    real(real_kind), intent(out) :: fvec(n)
    integer, value :: iflag
    
    real(real_kind), pointer :: params
    
    call c_f_pointer(p, params) ! dereference pointer
    fvec(1) = x(1)**2.d0 - params
    
    res = 0
  end function
  
end module

program main
  use cminpack2fort, only: hybrd1
  use test_wrapper, only: test_func1, test_func2, real_kind, f_data
  use, intrinsic :: iso_c_binding, only : c_loc, c_ptr
  implicit none

  type(f_data), target :: params1
  real(real_kind), target :: params2
  type(c_ptr), target :: p1, p2
  integer :: n = 1
  real(real_kind) :: x(1)
  real(real_kind) :: fvec(1)
  real(real_kind) :: tol = 1.d-8
  integer :: lwa = (1*(3*1+13))/2+2
  integer :: info
  real(real_kind) :: wa((1*(3*1+13))/2+2)

  x = [10.0d0]
  params1%a = 5.d0
  params1%b = 5.d0
  p1 = c_loc(params1) ! get pointer
  ! pass pointer to hybrd1
  call hybrd1(test_func1, p1, n, x, fvec, tol, info, wa, lwa)
  
  print*,x
  
  x = [10.0d0]
  params2 = 17.d0
  p2 = c_loc(params2) ! get pointer
  call hybrd1(test_func2, p2, n, x, fvec, tol, info, wa, lwa)
  
  print*,x

end program
