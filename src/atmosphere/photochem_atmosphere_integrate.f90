submodule(photochem_atmosphere) photochem_atmosphere_integrate
  implicit none
  
  ! Contains routines for integrating the photochemical equations
  ! forward in time. Here, we use the CVODE BDF integrator.
  
contains
  
  integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
                                result(ierr) bind(c, name='RhsFn')
    use, intrinsic :: iso_c_binding
    use fcvode_mod
    use fsundials_nvector_mod
    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined dat
    
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)
    character(len=err_len) :: err
    integer(c_long) :: nsteps(1)
    integer(c_int) :: loc_ierr
    real(c_double) :: hcur(1)
    real(real_kind) :: tmp, mx
    integer :: k, i, j, ii
    
    type(Atmosphere), pointer :: self
  
    ierr = 0
    
    call c_f_pointer(user_data, self)
    
    ! get data arrays from SUNDIALS vectors
    yvec(1:self%var%neqs) => FN_VGetArrayPointer(sunvec_y)
    fvec(1:self%var%neqs) => FN_VGetArrayPointer(sunvec_f)
    
    ! fill RHS vector
    call self%right_hand_side(self%var%neqs, yvec, fvec, err)
    loc_ierr = FCVodeGetNumSteps(self%wrk%cvode_mem, nsteps)
    
    if (nsteps(1) /= self%wrk%nsteps_previous .and. self%var%verbose > 0) then
      loc_ierr = FCVodeGetCurrentStep(self%wrk%cvode_mem, hcur)
      
      if (self%var%verbose == 1) then
        print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,'max(dy/dt) = ',es11.5)", &
             nsteps, tn, hcur(1),maxval(abs(fvec))
             
      elseif (self%var%verbose == 2) then
        ! Find the fastest changing variable
        tmp = 0.d0
        mx = tmp
        k = 1
        do ii = 1,self%var%neqs
          tmp = abs(fvec(ii)/yvec(ii))
          if (tmp > mx .and. abs(yvec(ii)) > self%var%atol) then
            mx = tmp
            k = ii
          endif
        enddo
        j = k/self%dat%nq
        i = k-j*self%dat%nq
        
        print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,"// &
             "'dy/dt =',es12.5,3x,' y =',es12.5,3x,a8,3x,' z =',f6.2,' km')", &
             nsteps, tn, hcur(1),fvec(k),yvec(k),trim(self%dat%species_names(i)),self%var%z(j+1)/1.d5
      endif
      
      self%wrk%nsteps_previous = nsteps(1)
    endif
    
    if (len_trim(err) /= 0) then
      print*,trim(err)//". CVODE will attempt to correct the error."
      ierr = 1
    endif
    return
  end function
  
  integer(c_int) function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                                tmp1, tmp2, tmp3) &
                                result(ierr) bind(C,name='JacFn')
    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod
    use fsunmatrix_band_mod
    use fsundials_matrix_mod
    
    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f
    type(SUNMatrix)        :: sunmat_J  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined data
    type(N_Vector)        :: tmp1, tmp2, tmp3
  
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: Jmat(:)
    character(len=err_len) :: err
    
    type(Atmosphere), pointer :: self
  
    ierr = 0
    
    call c_f_pointer(user_data, self)
    yvec(1:self%var%neqs) => FN_VGetArrayPointer(sunvec_y)
    Jmat(1:self%var%neqs*self%dat%lda) => FSUNBandMatrix_Data(sunmat_J)
    call self%jacobian(self%dat%lda*self%var%neqs, self%var%neqs, yvec, Jmat, err)
    if (len_trim(err) /= 0) then
      print*,trim(err)//". CVODE will attempt to correct the error."
      ierr = 1
    endif
    return
  
  end function
  
  subroutine evolve_background_atm(self, tstart, usol_start, t_eval, &
                                   solution, success, err)
                                   
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, FCVodeInit, FCVodeSStolerances, &
                          FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                          FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                          FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd, &
                          FCVodeSetUserData
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
    use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    
    ! in/out
    class(Atmosphere), target, intent(inout) :: self
    real(c_double), intent(in) :: tstart
    real(real_kind), intent(in) :: usol_start(:,:)
    ! real(c_double), intent(in) :: t_eval(num_t_eval)
    real(c_double), intent(in) :: t_eval(:)
    ! real(real_kind), intent(out) :: solution(nq,nz,num_t_eval)
    real(real_kind), intent(out) :: solution(:,:,:)
    logical, intent(out) :: success
    character(len=err_len), intent(out) :: err
    
    ! local variables
    real(c_double) :: tcur(1)    ! current time
    integer(c_int) :: ierr       ! error flag from C functions
    type(N_Vector), pointer :: sunvec_y ! sundials vector
    
    ! solution vector, neq is set in the ode_mod module
    real(c_double) :: yvec(self%var%neqs)
    integer(c_long) :: neqs_long
    integer(c_long) :: mu, ml
    integer(c_long) :: mxsteps_
    type(SUNMatrix), pointer :: sunmat
    type(SUNLinearSolver), pointer :: sunlin
    
    ! real(real_kind) :: fH2O(self%var%trop_ind)
    integer :: i, j, k, ii
    
    type(c_ptr)    :: user_data
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    type(Atmosphere), pointer :: self_ptr
  
    err = ''
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    ! check dimensions
    if (size(usol_start,1) /= dat%nq .or. size(usol_start,2) /= var%nz) then
      err = "Problem!"
      return
    endif
    if (size(solution,1) /= dat%nq .or. size(solution,2) /= var%nz) then
      err = "Problem!"
      return
    endif
    if (size(t_eval,1) /= size(solution, 3)) then
      err = "Problem!"
      return
    endif
    
    ! settings
    mxsteps_ = var%mxsteps
    neqs_long = var%neqs
    tcur   = tstart
    mu = dat%nq
    ml = dat%nq
    
    self_ptr => self
    user_data = c_loc(self_ptr)
    
    ! initialize solution vector
    do j=1,var%nz
      do i=1,dat%nq
        k = i + (j-1)*dat%nq
        yvec(k) = usol_start(i,j)
      enddo
    enddo
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == 1) then
        yvec(i) = var%lower_fix_mr(i)
      endif
    enddo
    if (dat%fix_water_in_trop) then
      call self%prep_atmosphere(usol_start, err)
      if (len_trim(err) /= 0) return 
      do j = 1,var%trop_ind
        k = dat%LH2O + (j-1)*dat%nq
        yvec(k) = self%wrk%fH2O(j)
      enddo
      ! if there is no water profile, then we should extrapolate H2O
      ! above the tropopause
      if (var%no_water_profile) then
        do j = var%trop_ind+1,var%nz
          k = dat%LH2O + (j-1)*dat%nq
          yvec(k) = self%wrk%fH2O(var%trop_ind)
        enddo
      endif
    endif

    ! create SUNDIALS N_Vector
    sunvec_y => FN_VMake_Serial(neqs_long, yvec)
    if (.not. associated(sunvec_y)) then
      err = "CVODE setup error."
      return
    end if
    
    ! create CVode memory
    wrk%cvode_mem = FCVodeCreate(CV_BDF)
    if (.not. c_associated(wrk%cvode_mem)) then
      err = "CVODE setup error."
      return
    end if
    
    ! set user data
    ierr = FCVodeSetUserData(wrk%cvode_mem, user_data)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeInit(wrk%cvode_mem, c_funloc(RhsFn), tstart, sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSStolerances(wrk%cvode_mem, var%rtol, var%atol)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    sunmat => FSUNBandMatrix(neqs_long, mu, ml)
    sunlin => FSUNLinSol_Band(sunvec_y,sunmat)
    
    ierr = FCVodeSetLinearSolver(wrk%cvode_mem, sunlin, sunmat)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    if (var%use_fast_jacobian) then
      ierr = FCVodeSetJacFn(wrk%cvode_mem, c_funloc(JacFn))
      if (ierr /= 0) then
        err = "CVODE setup error."
        return
      end if
    endif
    
    ierr = FCVodeSetMaxNumSteps(wrk%cvode_mem, mxsteps_)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetInitStep(wrk%cvode_mem, var%initial_dt)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxErrTestFails(wrk%cvode_mem, var%max_err_test_failures)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxOrd(wrk%cvode_mem, var%max_order)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    do ii = 1, size(t_eval)
      ierr = FCVode(wrk%cvode_mem, t_eval(ii), sunvec_y, tcur, CV_NORMAL)
      if (ierr /= 0) then
        success = .false.
      else
        success = .true.
        do j=1,var%nz
          do i=1,dat%nq  
            k = i + (j-1)*dat%nq
            solution(i,j,ii) = yvec(k)
          enddo
        enddo
        
        ! this will alter solution(:,:,ii) with proper fixed mixing ratios and H2O
        call self%prep_atmosphere(solution(:,:,ii), err)
        if (len_trim(err) /= 0) return
        solution(:,:,ii) = self%wrk%usol
                                     
      endif
    enddo
    
    ! free memory
    call FN_VDestroy(sunvec_y)
    call FCVodeFree(wrk%cvode_mem)
    ierr = FSUNLinSolFree(sunlin)
    if (ierr /= 0) then
      err = "CVODE deallocation error"
      return
    end if
    call FSUNMatDestroy(sunmat)

  end subroutine
  
  subroutine photochemical_equilibrium(self, success, err)
    class(Atmosphere), target, intent(inout) :: self
    logical, intent(out) :: success
    character(len=err_len), intent(out) :: err 
    
    real(real_kind), pointer :: solution(:,:,:)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    err = ""
    
    dat => self%dat
    var => self%var
    wrk => self%wrk

    solution(1:dat%nq,1:var%nz,1:1) => var%usol_out
    call evolve_background_atm(self, 0.d0, var%usol_init, [var%equilibrium_time],  &
                               solution, success, err)
    if (len_trim(err) /= 0) return
    var%at_photo_equilibrium = success 
    
  end subroutine

  
end submodule