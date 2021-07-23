
module cminpack2fort
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr, c_funptr, c_funloc
  implicit none
  
  interface
    subroutine hybrd1f_wrapper(funcptr, p, n, x, fvec, tol, info, wa, lwa) bind(C, name='hybrd1c_wrapper')
      use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr, c_funptr
      type(c_funptr), value :: funcptr
      type(c_ptr) :: p
      integer(c_int) :: n
      real(c_double) :: x(n)
      real(c_double) :: fvec(n)
      real(c_double) :: tol
      integer(c_int) :: info 
      real(c_double) :: wa(lwa)
      integer(c_int) :: lwa  
    end subroutine
  end interface
  
  abstract interface
    function cminpack_callback(p, n, x, fvec, iflag)
      use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr
      integer(c_int) :: cminpack_callback
      type(c_ptr) :: p
      integer(c_int), value :: n
      real(c_double), intent(in) :: x(n)
      real(c_double), intent(out) :: fvec(n)
      integer(c_int), value :: iflag
    end function    
  end interface

contains

  subroutine hybrd1(func, p, n, x, fvec, tol, info, wa, lwa)
    procedure(cminpack_callback) :: func
    type(c_ptr), intent(in) :: p
    integer(c_int), intent(in) :: n
    real(c_double), intent(inout) :: x(n)
    real(c_double), intent(inout) :: fvec(n)
    real(c_double), intent(in) :: tol
    integer(c_int), intent(inout) :: info
    real(c_double), intent(inout) :: wa(lwa)
    integer(c_int), intent(in) :: lwa
    type (c_funptr) :: funcptr
    funcptr = c_funloc(func)    
    call hybrd1f_wrapper(funcptr, p, n, x, fvec, tol, info, wa, lwa)
  end subroutine

end module





