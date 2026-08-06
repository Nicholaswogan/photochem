module block_thomas_solver
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_loc, c_ptr, c_f_pointer
  use fsundials_nvector_mod, only: N_Vector
  use fsundials_matrix_mod, only: SUNMatrix
  use fsundials_linearsolver_mod, only: SUNLinearSolver
  implicit none
  private

  public :: FSUNLinSol_BlockThomas

  interface
    function SUNLinSol_BlockThomas_c(y, A, block_size) &
        bind(c, name='SUNLinSol_BlockThomas') result(ptr)
      import :: c_int64_t, c_ptr
      type(c_ptr), value :: y, A
      integer(c_int64_t), value :: block_size
      type(c_ptr) :: ptr
    end function
  end interface

contains

  function FSUNLinSol_BlockThomas(y, A, block_size) result(solver)
    type(N_Vector), target, intent(inout) :: y
    type(SUNMatrix), target, intent(inout) :: A
    integer(c_int64_t), intent(in) :: block_size
    type(SUNLinearSolver), pointer :: solver

    type(c_ptr) :: ptr

    ptr = SUNLinSol_BlockThomas_c(c_loc(y), c_loc(A), block_size)
    call c_f_pointer(ptr, solver)

  end function

end module
