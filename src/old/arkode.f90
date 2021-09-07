subroutine evolve_background_atm_ark(tstart, nq, nz, usol_start, num_t_eval, t_eval, rtol, atol, &
                                 mxsteps, solution, success, err)
  use photochem_data, only: water_sat_trop, LH2O, nsp, nrT, kj, nw
  use photochem_vars, only: no_water_profile, neqs, lowerboundcond, lower_fix_mr, trop_ind, &
                         initial_dt, use_fast_jacobian, max_err_test_failures, max_order, trop_ind
  use photochem_wrk, only: cvode_mem
  
  use, intrinsic :: iso_c_binding
  use farkode_mod 
  use farkode_arkstep_mod
  
  use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
  use fnvector_serial_mod, only: FN_VMake_Serial   
  use fsunmatrix_band_mod, only: FSUNBandMatrix
  use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
  use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
  use fsunlinsol_band_mod, only: FSUNLinSol_Band
  
  ! in/out
  real(c_double), intent(in) :: tstart
  integer, intent(in) :: nq, nz
  real(real_kind), intent(in) :: usol_start(nq,nz)
  integer, intent(in) :: num_t_eval
  real(c_double), intent(in) :: t_eval(num_t_eval)
  real(c_double), intent(in) :: rtol, atol
  integer, intent(in) :: mxsteps
  real(real_kind), intent(out) :: solution(nq,nz,num_t_eval)
  logical, intent(out) :: success
  character(len=err_len), intent(out) :: err
  
  ! local variables
  real(c_double) :: tcur(1)    ! current time
  integer(c_int) :: ierr       ! error flag from C functions
  ! type(c_ptr)    :: cvode_mem  ! CVODE memory
  type(N_Vector), pointer :: sunvec_y ! sundials vector
  
  integer(c_int) :: imethod, idefault, pq    ! time step adaptivity parameters
  real(c_double) :: adapt_params(3)          ! time step adaptivity parameters
  real(c_double) :: yvec(neqs)
  integer(c_long) :: neqs_long
  integer(c_long) :: mu, ml
  integer(c_long) :: mxsteps_
  type(SUNMatrix), pointer :: sunmat
  type(SUNLinearSolver), pointer :: sunlin
  
  real(real_kind) :: fH2O(trop_ind)
  integer :: i, j, k, ii
  
  type(c_ptr)    :: user_data
  type(WrkBackgroundAtm), target :: wrk

  err = ''
  
  ! settings
  mxsteps_ = mxsteps
  neqs_long = neqs
  tcur   = tstart
  mu = nq
  ml = nq
  
  call wrk%init(nsp, nq, nz, nrT, kj, nw, trop_ind)
  user_data = c_loc(wrk)
  
  ! initialize solution vector
  do j=1,nz
    do i=1,nq
      k = i + (j-1)*nq
      yvec(k) = usol_start(i,j)
    enddo
  enddo
  do i = 1,nq
    if (lowerboundcond(i) == 1) then
      yvec(i) = lower_fix_mr(i)
    endif
  enddo
  if (water_sat_trop) then
    call water_mixing_ratio(nq, nz, trop_ind, usol_start, fH2O, err)
    if (len_trim(err) /= 0) return 
    do j = 1,trop_ind
      k = LH2O + (j-1)*nq
      yvec(k) = fH2O(j)
    enddo
    ! if there is no water profile, then we should extrapolate H2O
    ! above the tropopause
    if (no_water_profile) then
      do j = trop_ind+1,nz
        k = LH2O + (j-1)*nq
        yvec(k) = fH2O(trop_ind)
      enddo
    endif
  endif
  
  ! create SUNDIALS N_Vector
  sunvec_y => FN_VMake_Serial(neqs_long, yvec)
  if (.not. associated(sunvec_y)) then
    err = "CVODE setup error."
    return
  end if
  
  cvode_mem = FARKStepCreate(c_null_funptr, c_funloc(RhsFn), tstart, sunvec_y)
  if (.not. c_associated(cvode_mem)) print *,'ERROR: arkode_mem = NULL'
  
  ! set user data
  ierr = FARKStepSetUserData(cvode_mem, user_data)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FARKStepSStolerances(cvode_mem, rtol, atol)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  sunmat => FSUNBandMatrix(neqs_long, mu, ml)
  sunlin => FSUNLinSol_Band(sunvec_y,sunmat)
  
  ierr = FARKStepSetLinearSolver(cvode_mem, sunlin, sunmat)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  if (use_fast_jacobian) then
    ierr = FARKStepSetJacFn(cvode_mem, c_funloc(JacFn))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
  endif
  
  ierr = FARKStepSetMaxNumSteps(cvode_mem, mxsteps_)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FARKStepSetInitStep(cvode_mem, initial_dt)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FARKStepSetMaxErrTestFails(cvode_mem, max_err_test_failures)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ierr = FARKStepSetNonlinConvCoef(cvode_mem, 1.d-7)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ! ierr = FARKStepSetPredictorMethod(cvode_mem, 3)
  ! if (ierr /= 0) then
  !   err = "CVODE setup error."
  !   return
  ! end if
  
  ierr = FARKStepSetSafetyFactor(cvode_mem, 0.9d0)
  if (ierr /= 0) then
    err = "CVODE setup error."
    return
  end if
  
  ! ierr = FARKStepSetOrder(cvode_mem, 2)
  ! if (ierr /= 0) then
  !   err = "CVODE setup error."
  !   return
  ! end if
  
  imethod = 4
  idefault = 1
  pq = 0
  ierr = FARKStepSetAdaptivityMethod(cvode_mem, imethod, idefault, pq, adapt_params)
  if (ierr /= 0) then
     write(*,*) 'Error in FARKStepSetAdaptivityMethod, ierr = ', ierr, '; halting'
     stop 1
  end if
  
  do ii = 1, num_t_eval
    ierr = FARKStepEvolve(cvode_mem, t_eval(ii), sunvec_y, tcur, ARK_NORMAL)
    if (ierr /= 0) then
      success = .false.
    else
      success = .true.
      do j=1,nz
        do i=1,nq  
          k = i + (j-1)*nq
          solution(i,j,ii) = yvec(k)
        enddo
      enddo
      
      ! this will alter solution(:,:,ii) with proper fixed mixing ratios and H2O
      wrk%usol = solution(:,:,ii)
      call prep_all_background_gas(wrk, err)
      solution(:,:,ii) = wrk%usol
                                   
    endif
  enddo
      
  ! free memory
  call FN_VDestroy(sunvec_y)
  call FARKStepFree(cvode_mem)
  ierr = FSUNLinSolFree(sunlin)
  if (ierr /= 0) then
    err = "CVODE deallocation error"
    return
  end if
  call FSUNMatDestroy(sunmat)

end subroutine


subroutine jac_background_gas_wrapper(ldaa, neqs, usol_flat, ddjac, err)
  use iso_c_binding, only: c_ptr, c_loc
  use photochem_data, only: nq, nsp, nrT, kj, nw, lda
  use photochem_vars, only: nz, trop_ind
  
  integer, intent(in) :: ldaa, neqs
  real(real_kind), intent(in) :: usol_flat(neqs)
  real(real_kind), intent(out) :: ddjac(ldaa,neqs)
  character(len=err_len), intent(out) :: err
  
  type(c_ptr) :: user_data
  type(WrkBackgroundAtm), target :: wrk
  
  real(real_kind), allocatable, target  :: jac(:)
  real(real_kind), pointer :: djac(:,:)
  
  integer :: i,j
  
  call wrk%init(nsp, nq, nz, nrT, kj, nw, trop_ind)
  user_data = c_loc(wrk)
  
  allocate(jac(lda*neqs))
  djac(1:lda,1:neqs) => jac
  
  call jac_background_gas(lda*neqs, neqs, user_data, usol_flat, jac, err)
  
  do i =1,ldaa
    do j=1,neqs
      ddjac(i,j) = djac(i+nq,j)
    enddo
  enddo
  
end subroutine


subroutine jac_background_gas_wrapper(ldaa, neqs, usol_flat, ddjac, err)
  use iso_c_binding, only: c_ptr, c_loc
  use photochem_data, only: nq, nsp, nrT, kj, nw, lda
  use photochem_vars, only: nz, trop_ind
  
  integer, intent(in) :: ldaa, neqs
  real(real_kind), intent(in) :: usol_flat(neqs)
  real(real_kind), intent(out) :: ddjac(ldaa,neqs)
  character(len=err_len), intent(out) :: err
  
  type(c_ptr) :: user_data
  type(WrkBackgroundAtm), target :: wrk
  
  real(real_kind), allocatable, target  :: jac(:)
  real(real_kind), pointer :: djac(:,:)
  
  integer :: i,j
  
  call wrk%init(nsp, nq, nz, nrT, kj, nw, trop_ind)
  user_data = c_loc(wrk)
  
  allocate(jac(lda*neqs))
  djac(1:lda,1:neqs) => jac
  
  call jac_background_gas(lda*neqs, neqs, user_data, usol_flat, jac, err)
  
  do i =1,ldaa
    do j=1,neqs
      ddjac(i,j) = djac(i+nq,j)
    enddo
  enddo
  
end subroutine


subroutine rhs_background_gas_wrapper(neqs, usol_flat, rhs, err)
  use iso_c_binding, only: c_ptr, c_loc
  use photochem_data, only: nq, nsp, nrT, kj, nw
  use photochem_vars, only: nz, trop_ind
  
  integer, intent(in) :: neqs
  real(real_kind), intent(in) :: usol_flat(neqs)
  real(real_kind), intent(out) :: rhs(neqs)
  character(len=err_len), intent(out) :: err
  
  type(c_ptr) :: user_data
  type(WrkBackgroundAtm), target :: wrk
  
  call wrk%init(nsp, nq, nz, nrT, kj, nw, trop_ind)
  user_data = c_loc(wrk)
  call rhs_background_gas(neqs, user_data, usol_flat, rhs, err)
  
end subroutine

