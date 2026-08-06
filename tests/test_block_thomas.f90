program test_block_thomas
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
  use block_thomas_solver, only: FSUNLinSol_BlockThomas
  use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
  use fnvector_serial_mod, only: FN_VMake_Serial
  use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
  use fsunmatrix_band_mod, only: FSUNBandMatrix, FSUNBandMatrix_Data, &
                                 FSUNBandMatrix_LDim, &
                                 FSUNBandMatrix_StoredUpperBandwidth
  use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolSetup, &
                                       FSUNLinSolSolve, FSUNLinSolFree, &
                                       SUNLS_ILL_INPUT
  use fsunlinsol_band_mod, only: FSUNLinSol_Band
  implicit none

  call test_solver(1.0_c_double, 2.0e-12_c_double)
  ! This matrix has a condition-number scale of roughly 1e8, so forward
  ! solution error near epsilon*condition_number is expected.
  call test_solver(1.0e8_c_double, 5.0e-8_c_double)
  call test_invalid_coupling()
  print *, 'test_block_thomas passed'

contains

  subroutine test_solver(scale_range, tolerance)
    integer, parameter :: block_size = 4, nblocks = 6
    integer, parameter :: n = block_size*nblocks
    real(c_double), intent(in) :: scale_range, tolerance

    real(c_double) :: diagonal(block_size,block_size,nblocks)
    real(c_double) :: upper(block_size,nblocks-1)
    real(c_double) :: lower(block_size,nblocks-1)
    real(c_double), target :: vector_storage(n), rhs_storage(n)
    real(c_double), target :: x_band(n), x_block(n)
    real(c_double) :: x_expected(n), rhs(n), residual(n)
    real(c_double) :: scales(block_size)
    real(c_double) :: relative_error, band_error, solver_difference
    real(c_double) :: relative_residual
    type(N_Vector), pointer :: vector, rhs_vector, x_band_vector, x_block_vector
    type(SUNMatrix), pointer :: matrix_band, matrix_block
    type(SUNLinearSolver), pointer :: solver_band, solver_block
    integer(c_int) :: ierr
    integer :: i, j, k

    do i = 1,block_size
      if (block_size == 1) then
        scales(i) = 1.0_c_double
      else
        scales(i) = scale_range**(real(i-1,c_double)/real(block_size-1,c_double))
      endif
    enddo

    diagonal = 0.0_c_double
    do k = 1,nblocks
      do j = 1,block_size
        do i = 1,block_size
          if (i == j) then
            diagonal(i,j,k) = scales(i)*(4.0_c_double + 0.1_c_double*k + 0.05_c_double*i)
          else
            diagonal(i,j,k) = 0.01_c_double*sqrt(scales(i)*scales(j))/real(i+j+k,c_double)
          endif
        enddo
      enddo
    enddo
    do k = 1,nblocks-1
      do i = 1,block_size
        upper(i,k) = -0.07_c_double*scales(i)/(1.0_c_double + 0.1_c_double*k)
        lower(i,k) = -0.05_c_double*scales(i)/(1.0_c_double + 0.1_c_double*k)
      enddo
    enddo

    ! Exercise partial pivoting in a local diagonal block. The first pivot is
    ! exactly zero, but the block and full matrix remain nonsingular.
    if (scale_range == 1.0_c_double) then
      diagonal(1,1,1) = 0.0_c_double
      diagonal(2,1,1) = 1.2_c_double
      diagonal(1,2,1) = 1.0_c_double
      diagonal(2,2,1) = 3.0_c_double
    endif

    do i = 1,n
      x_expected(i) = sin(0.17_c_double*i) + 0.25_c_double*cos(0.11_c_double*i)
    enddo
    call apply_blocks(diagonal, upper, lower, x_expected, rhs)

    vector_storage = 0.0_c_double
    rhs_storage = rhs
    x_band = 0.0_c_double
    x_block = 0.0_c_double
    vector => FN_VMake_Serial(int(n,c_int64_t), vector_storage)
    rhs_vector => FN_VMake_Serial(int(n,c_int64_t), rhs_storage)
    x_band_vector => FN_VMake_Serial(int(n,c_int64_t), x_band)
    x_block_vector => FN_VMake_Serial(int(n,c_int64_t), x_block)
    if (.not. associated(vector) .or. .not. associated(rhs_vector) .or. &
        .not. associated(x_band_vector) .or. .not. associated(x_block_vector)) then
      error stop 'Unable to allocate SUNDIALS vectors.'
    endif

    matrix_band => FSUNBandMatrix(int(n,c_int64_t), int(block_size,c_int64_t), &
                                 int(block_size,c_int64_t))
    matrix_block => FSUNBandMatrix(int(n,c_int64_t), int(block_size,c_int64_t), &
                                  int(block_size,c_int64_t))
    if (.not. associated(matrix_band) .or. .not. associated(matrix_block)) then
      error stop 'Unable to allocate SUNDIALS matrices.'
    endif
    call fill_band_matrix(matrix_band, diagonal, upper, lower)
    call fill_band_matrix(matrix_block, diagonal, upper, lower)

    solver_band => FSUNLinSol_Band(vector, matrix_band)
    solver_block => FSUNLinSol_BlockThomas(vector, matrix_block, int(block_size,c_int64_t))
    if (.not. associated(solver_band) .or. .not. associated(solver_block)) then
      error stop 'Unable to allocate SUNDIALS linear solvers.'
    endif

    ierr = FSUNLinSolSetup(solver_band, matrix_band)
    if (ierr /= 0) error stop 'Band solver setup failed.'
    ierr = FSUNLinSolSolve(solver_band, matrix_band, x_band_vector, rhs_vector, 0.0_c_double)
    if (ierr /= 0) error stop 'Band solve failed.'

    ierr = FSUNLinSolSetup(solver_block, matrix_block)
    if (ierr /= 0) error stop 'Block-Thomas solver setup failed.'
    ierr = FSUNLinSolSolve(solver_block, matrix_block, x_block_vector, rhs_vector, 0.0_c_double)
    if (ierr /= 0) error stop 'Block-Thomas solve failed.'

    relative_error = maxval(abs(x_block-x_expected))/max(1.0_c_double,maxval(abs(x_expected)))
    band_error = maxval(abs(x_band-x_expected))/max(1.0_c_double,maxval(abs(x_expected)))
    solver_difference = maxval(abs(x_block-x_band))/max(1.0_c_double,maxval(abs(x_band)))
    if (relative_error > tolerance) then
      print *, 'Block-Thomas relative solution error:', relative_error
      print *, 'Band relative solution error:', band_error
      error stop 'Block-Thomas solution is inaccurate.'
    endif
    if (solver_difference > tolerance) then
      print *, 'Block-Thomas difference from band solve:', solver_difference
      error stop 'Block-Thomas and band solutions disagree.'
    endif

    call apply_blocks(diagonal, upper, lower, x_block, residual)
    residual = residual-rhs
    relative_residual = maxval(abs(residual))/max(1.0_c_double,maxval(abs(rhs)))
    if (relative_residual > 2.0e-13_c_double) then
      print *, 'Block-Thomas relative residual:', relative_residual
      error stop 'Block-Thomas residual is too large.'
    endif

    ierr = FSUNLinSolFree(solver_band)
    ierr = FSUNLinSolFree(solver_block)
    call FSUNMatDestroy(matrix_band)
    call FSUNMatDestroy(matrix_block)
    call FN_VDestroy(vector)
    call FN_VDestroy(rhs_vector)
    call FN_VDestroy(x_band_vector)
    call FN_VDestroy(x_block_vector)

  end subroutine

  subroutine test_invalid_coupling()
    integer, parameter :: block_size = 3, nblocks = 2
    integer, parameter :: n = block_size*nblocks
    real(c_double), target :: vector_storage(n)
    real(c_double), pointer :: data_raw(:), data(:)
    type(N_Vector), pointer :: vector
    type(SUNMatrix), pointer :: matrix
    type(SUNLinearSolver), pointer :: solver
    integer(c_int) :: ierr
    integer :: i, ldim, stored_upper, index

    vector_storage = 0.0_c_double
    vector => FN_VMake_Serial(int(n,c_int64_t), vector_storage)
    matrix => FSUNBandMatrix(int(n,c_int64_t), int(block_size,c_int64_t), &
                            int(block_size,c_int64_t))
    if (.not. associated(vector) .or. .not. associated(matrix)) then
      error stop 'Unable to allocate invalid-structure test data.'
    endif

    ldim = int(FSUNBandMatrix_LDim(matrix))
    stored_upper = int(FSUNBandMatrix_StoredUpperBandwidth(matrix))
    data_raw => FSUNBandMatrix_Data(matrix)
    data(1:ldim*n) => data_raw
    data = 0.0_c_double
    do i = 0,n-1
      index = 1 + i*ldim + stored_upper
      data(index) = 2.0_c_double
    enddo
    ! Coupling from species 2 in layer 1 to species 1 in layer 2. This lies
    ! inside the scalar band but violates the diagonal transport-block form.
    index = 1 + block_size*ldim + stored_upper + 1-block_size
    data(index) = 0.25_c_double

    solver => FSUNLinSol_BlockThomas(vector, matrix, int(block_size,c_int64_t))
    if (.not. associated(solver)) error stop 'Unable to allocate block solver.'
    ierr = FSUNLinSolSetup(solver, matrix)
    if (ierr /= SUNLS_ILL_INPUT) then
      error stop 'Block solver accepted a non-diagonal inter-layer coupling.'
    endif

    ierr = FSUNLinSolFree(solver)
    call FSUNMatDestroy(matrix)
    call FN_VDestroy(vector)

  end subroutine

  subroutine fill_band_matrix(matrix, diagonal, upper, lower)
    type(SUNMatrix), target, intent(inout) :: matrix
    real(c_double), intent(in) :: diagonal(:,:,:), upper(:,:), lower(:,:)
    real(c_double), pointer :: data_raw(:), data(:)
    integer :: block_size, nblocks, n, ldim, stored_upper
    integer :: block, row, column, matrix_row, matrix_column, index

    block_size = size(diagonal,1)
    nblocks = size(diagonal,3)
    n = block_size*nblocks
    ldim = int(FSUNBandMatrix_LDim(matrix))
    stored_upper = int(FSUNBandMatrix_StoredUpperBandwidth(matrix))
    data_raw => FSUNBandMatrix_Data(matrix)
    data(1:ldim*n) => data_raw
    data = 0.0_c_double

    do block = 1,nblocks
      do column = 1,block_size
        matrix_column = (block-1)*block_size + column-1
        do row = 1,block_size
          matrix_row = (block-1)*block_size + row-1
          index = 1 + matrix_column*ldim + stored_upper + matrix_row-matrix_column
          data(index) = diagonal(row,column,block)
        enddo
      enddo
    enddo
    do block = 1,nblocks-1
      do row = 1,block_size
        matrix_row = (block-1)*block_size + row-1
        matrix_column = block*block_size + row-1
        index = 1 + matrix_column*ldim + stored_upper + matrix_row-matrix_column
        data(index) = upper(row,block)
        matrix_row = block*block_size + row-1
        matrix_column = (block-1)*block_size + row-1
        index = 1 + matrix_column*ldim + stored_upper + matrix_row-matrix_column
        data(index) = lower(row,block)
      enddo
    enddo

  end subroutine

  subroutine apply_blocks(diagonal, upper, lower, x, result)
    real(c_double), intent(in) :: diagonal(:,:,:), upper(:,:), lower(:,:), x(:)
    real(c_double), intent(out) :: result(:)
    integer :: block_size, nblocks, block, row, column, offset

    block_size = size(diagonal,1)
    nblocks = size(diagonal,3)
    result = 0.0_c_double
    do block = 1,nblocks
      offset = (block-1)*block_size
      do row = 1,block_size
        do column = 1,block_size
          result(offset+row) = result(offset+row) + &
              diagonal(row,column,block)*x(offset+column)
        enddo
        if (block > 1) then
          result(offset+row) = result(offset+row) + &
              lower(row,block-1)*x(offset-block_size+row)
        endif
        if (block < nblocks) then
          result(offset+row) = result(offset+row) + &
              upper(row,block)*x(offset+block_size+row)
        endif
      enddo
    enddo

  end subroutine

end program
