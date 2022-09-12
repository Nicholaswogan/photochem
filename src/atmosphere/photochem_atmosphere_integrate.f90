submodule(photochem_atmosphere) photochem_atmosphere_integrate
  implicit none
  
  ! Contains routines for integrating the photochemical equations
  ! forward in time. Here, we use the CVODE BDF integrator.
  
contains
  
  module function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
                        result(ierr) bind(c, name='RhsFn')
    use, intrinsic :: iso_c_binding
    use fcvode_mod
    use fsundials_nvector_mod
    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined dat
    integer(c_int)        :: ierr
    
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)
    character(:), allocatable :: err
    integer(c_long) :: nsteps(1)
    integer(c_int) :: loc_ierr
    real(c_double) :: hcur(1)
    real(dp) :: tmp, mx
    integer :: k, i, j, ii
    
    type(Atmosphere), pointer :: self
  
    ierr = 0
    
    call c_f_pointer(user_data, self)
    
    ! get data arrays from SUNDIALS vectors
    yvec(1:self%var%neqs) => FN_VGetArrayPointer(sunvec_y)
    fvec(1:self%var%neqs) => FN_VGetArrayPointer(sunvec_f)
    
    ! fill RHS vector
    call self%right_hand_side(self%var%neqs, yvec, fvec, err)
    loc_ierr = FCVodeGetNumSteps(self%wrk%sun%cvode_mem, nsteps)
    
    if (nsteps(1) /= self%wrk%nsteps_previous .and. self%var%verbose > 0) then
      loc_ierr = FCVodeGetCurrentStep(self%wrk%sun%cvode_mem, hcur)
      
      if (self%var%verbose == 1) then
        print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,'max(dy/dt) = ',es11.5)", &
             nsteps, tn, hcur(1),maxval(abs(fvec))
             
      elseif (self%var%verbose == 2) then
        ! Find the fastest changing variable
        tmp = 0.0_dp
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
             nsteps, tn, hcur(1),fvec(k),yvec(k),trim(self%dat%species_names(i)),self%var%z(j+1)/1.e5_dp
      endif
      
      self%wrk%nsteps_previous = nsteps(1)
    endif
    
    if (allocated(err)) then
      if (self%var%verbose > 0) then
        print*,trim(err)//". CVODE will attempt to correct the error."
      endif
      ierr = 1
    endif
    return
  end function
  
  module function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
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
    integer(c_int)        :: ierr
  
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: Jmat(:)
    character(:), allocatable :: err
    
    type(Atmosphere), pointer :: self
  
    ierr = 0
    
    call c_f_pointer(user_data, self)
    yvec(1:self%var%neqs) => FN_VGetArrayPointer(sunvec_y)
    Jmat(1:self%var%neqs*self%dat%lda) => FSUNBandMatrix_Data(sunmat_J)
    call self%jacobian(self%dat%lda*self%var%neqs, self%var%neqs, yvec, Jmat, err)
    if (allocated(err)) then
      if (self%var%verbose > 0) then
        print*,trim(err)//". CVODE will attempt to correct the error."
      endif
      ierr = 1
    endif
    return
  
  end function

  subroutine ErrHandlerFn(error_code, module_, func, msg, eh_data) bind(c, name='ErrHandlerFn')
    use iso_c_binding
    integer(c_int), value :: error_code
    character(kind=c_char) :: module_(*)
    character(kind=c_char) :: func(*)
    character(kind=c_char) :: msg(*)
    type(c_ptr), value, intent(in) :: eh_data
  end subroutine
  
  module function evolve(self, filename, tstart, usol_start, t_eval, overwrite, err) result(success)
                                   
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, CV_ONE_STEP, FCVodeInit, FCVodeSStolerances, &
                          FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                          FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                          FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd, &
                          FCVodeSetUserData, FCVodeSetErrHandlerFn
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
    use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    
    use photochem_enum, only: MixingRatioBC
    
    ! in/out
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    real(c_double), intent(in) :: tstart
    real(dp), intent(in) :: usol_start(:,:)
    ! real(c_double), intent(in) :: t_eval(num_t_eval)
    real(c_double), intent(in) :: t_eval(:)
    logical, optional, intent(in) :: overwrite
    logical :: success
    character(:), allocatable, intent(out) :: err
    
    ! local variables
    real(c_double) :: tcur(1)    ! current time
    integer(c_int) :: ierr       ! error flag from C functions
    ! type(N_Vector), pointer :: sunvec_y ! sundials vector
    
    ! solution vector, neq is set in the ode_mod module
    ! real(c_double) :: yvec(self%var%neqs)
    integer(c_int64_t) :: neqs_long
    integer(c_int64_t) :: mu, ml
    integer(c_long) :: mxsteps_
    ! type(SUNMatrix), pointer :: sunmat
    ! type(SUNLinearSolver), pointer :: sunlin
    
    integer :: i, j, k, ii, io
    logical :: overwrite_
    
    type(c_ptr)    :: user_data
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    type(Atmosphere), pointer :: self_ptr
  
    dat => self%dat
    var => self%var
    wrk => self%wrk

    if (present(overwrite)) then
      overwrite_ = overwrite
    else 
      overwrite_ = .false.
    endif
    
    call wrk%sun%finalize(err)
    if (allocated(err)) return
    
    ! check dimensions
    if (size(usol_start,1) /= dat%nq .or. size(usol_start,2) /= var%nz) then
      err = "'usol_start' has the wrong dimensions."
      return
    endif
    
    ! file prep
    if (overwrite_) then
      open(1, file = filename, status='replace', form="unformatted")
    else
      open(1, file = filename, status='new', form="unformatted",iostat=io)
      if (io /= 0) then
        err = "Unable to create file "//trim(filename)//" because it already exists"
        return
      endif
    endif
    write(1) dat%nq
    write(1) var%nz
    write(1) var%z
    write(1) dat%species_names(1:dat%nq)
    write(1) size(t_eval)
    close(1)
    
    ! settings
    mxsteps_ = var%mxsteps
    neqs_long = var%neqs
    tcur   = tstart
    mu = dat%nq
    ml = dat%nq
    
    self_ptr => self
    user_data = c_loc(self_ptr)
    
    ! initialize solution vector
    allocate(wrk%sun%yvec(var%neqs))
    do j=1,var%nz
      do i=1,dat%nq
        k = i + (j-1)*dat%nq
        wrk%sun%yvec(k) = usol_start(i,j)
      enddo
    enddo
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == MixingRatioBC) then
        wrk%sun%yvec(i) = var%lower_fix_mr(i)
      endif
    enddo

    ! create SUNDIALS N_Vector
    wrk%sun%sunvec_y => FN_VMake_Serial(neqs_long, wrk%sun%yvec)
    if (.not. associated(wrk%sun%sunvec_y)) then
      err = "CVODE setup error."
      return
    end if
    
    ! create CVode memory
    wrk%sun%cvode_mem = FCVodeCreate(CV_BDF)
    if (.not. c_associated(wrk%sun%cvode_mem)) then
      err = "CVODE setup error."
      return
    end if
    
    ! set user data
    ierr = FCVodeSetUserData(wrk%sun%cvode_mem, user_data)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeInit(wrk%sun%cvode_mem, c_funloc(RhsFn), tstart, wrk%sun%sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSStolerances(wrk%sun%cvode_mem, var%rtol, var%atol)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    wrk%sun%sunmat => FSUNBandMatrix(neqs_long, mu, ml)
    wrk%sun%sunlin => FSUNLinSol_Band(wrk%sun%sunvec_y, wrk%sun%sunmat)
    
    ierr = FCVodeSetLinearSolver(wrk%sun%cvode_mem, wrk%sun%sunlin, wrk%sun%sunmat)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetJacFn(wrk%sun%cvode_mem, c_funloc(JacFn))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxNumSteps(wrk%sun%cvode_mem, mxsteps_)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetInitStep(wrk%sun%cvode_mem, var%initial_dt)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxErrTestFails(wrk%sun%cvode_mem, var%max_err_test_failures)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxOrd(wrk%sun%cvode_mem, var%max_order)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if

    if (var%verbose == 0) then
      ierr = FCVodeSetErrHandlerFn(wrk%sun%cvode_mem, c_funloc(ErrHandlerFn), c_null_ptr)
      if (ierr /= 0) then
        err = "CVODE setup error."
        return
      end if
    endif
    
    do ii = 1, size(t_eval)
      ierr = FCVode(wrk%sun%cvode_mem, t_eval(ii), wrk%sun%sunvec_y, tcur, CV_NORMAL)
      if (ierr /= 0) then
        success = .false.
        exit
      else
        success = .true.
        do j=1,var%nz
          do i=1,dat%nq  
            k = i + (j-1)*dat%nq
            wrk%usol(i,j) = wrk%sun%yvec(k)
          enddo
        enddo
        
        call self%prep_atmosphere(wrk%usol, err)
        if (allocated(err)) return
        
        open(1, file = filename, status='old', form="unformatted",position="append")
        write(1) tcur(1)
        write(1) wrk%usol
        close(1)
                               
      endif
    enddo
    
    ! free memory
    call wrk%sun%finalize(err)
    if (allocated(err)) return

  end function
  
  module subroutine photochemical_equilibrium(self, success, err)
    class(Atmosphere), target, intent(inout) :: self
    logical, intent(out) :: success
    character(:), allocatable, intent(out) :: err 
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    real(dp) :: tn
    integer :: nsteps
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    call self%initialize_stepper(var%usol_init, err)
    if (allocated(err)) return
    
    success = .true.
    tn = 0
    nsteps = 0
    do while(tn < var%equilibrium_time)
      tn = self%step(err)
      if (allocated(err)) then
        ! If the step fails then exit
        success = .false.
        exit
      endif
      nsteps = nsteps + 1
      if (nsteps > var%mxsteps) then
        ! If we did more steps then max steps, then exit
        success = .false.
        exit
      endif
    enddo
    
    var%usol_out = wrk%usol
    call self%destroy_stepper(err)
    if (allocated(err)) return
    
    var%at_photo_equilibrium = success 
    
  end subroutine

  module subroutine initialize_stepper(self, usol_start, err)
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, CV_ONE_STEP, FCVodeInit, FCVodeSStolerances, &
                          FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                          FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                          FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd, &
                          FCVodeSetUserData, FCVodeSetErrHandlerFn
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
    use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    
    use photochem_enum, only: MixingRatioBC
    
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_start(:,:)
    character(:), allocatable, intent(out) :: err
    
    real(c_double) :: tstart
    integer(c_int) :: ierr       ! error flag from C functions
    integer(c_int64_t) :: neqs_long
    integer(c_int64_t) :: mu, ml
    integer(c_long) :: mxsteps_
    type(c_ptr)    :: user_data
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    type(Atmosphere), pointer :: self_ptr
    
    integer :: i, j, k
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(usol_start,1) /= dat%nq .or. size(usol_start,2) /= var%nz) then
      err = "Input 'usol_start' to 'initialize_stepper' is the wrong dimension"
      return
    endif
    
    mxsteps_ = var%mxsteps
    neqs_long = var%neqs
    tstart = 0
    mu = dat%nq
    ml = dat%nq
    self_ptr => self
    user_data = c_loc(self_ptr)
    
    call self%destroy_stepper(err)
    if (allocated(err)) return
    
    ! initialize solution vector
    allocate(wrk%sun%yvec(var%neqs))
    do j=1,var%nz
      do i=1,dat%nq
        k = i + (j-1)*dat%nq
        wrk%sun%yvec(k) = usol_start(i,j)
      enddo
    enddo
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == MixingRatioBC) then
        wrk%sun%yvec(i) = var%lower_fix_mr(i)
      endif
    enddo
    
    ! create SUNDIALS N_Vector
    wrk%sun%sunvec_y => FN_VMake_Serial(neqs_long, wrk%sun%yvec)
    if (.not. associated(wrk%sun%sunvec_y)) then
      err = "CVODE setup error."
      return
    end if
    
    ! create CVode memory
    wrk%sun%cvode_mem = FCVodeCreate(CV_BDF)
    if (.not. c_associated(wrk%sun%cvode_mem)) then
      err = "CVODE setup error."
      return
    end if
    
    ! set user data
    ierr = FCVodeSetUserData(wrk%sun%cvode_mem, user_data)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeInit(wrk%sun%cvode_mem, c_funloc(RhsFn), tstart, wrk%sun%sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSStolerances(wrk%sun%cvode_mem, var%rtol, var%atol)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    wrk%sun%sunmat => FSUNBandMatrix(neqs_long, mu, ml)
    wrk%sun%sunlin => FSUNLinSol_Band(wrk%sun%sunvec_y, wrk%sun%sunmat)
    
    ierr = FCVodeSetLinearSolver(wrk%sun%cvode_mem, wrk%sun%sunlin, wrk%sun%sunmat)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetJacFn(wrk%sun%cvode_mem, c_funloc(JacFn))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxNumSteps(wrk%sun%cvode_mem, mxsteps_)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetInitStep(wrk%sun%cvode_mem, var%initial_dt)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxErrTestFails(wrk%sun%cvode_mem, var%max_err_test_failures)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxOrd(wrk%sun%cvode_mem, var%max_order)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if

    if (var%verbose == 0) then
      ierr = FCVodeSetErrHandlerFn(wrk%sun%cvode_mem, c_funloc(ErrHandlerFn), c_null_ptr)
      if (ierr /= 0) then
        err = "CVODE setup error."
        return
      end if
    endif
    
  end subroutine
  
  module function step(self, err) result(tn)
    use iso_c_binding, only: c_null_ptr, c_int, c_double, c_associated
    use fcvode_mod, only: CV_ONE_STEP, FCVode
    class(Atmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    real(dp) :: tn
    
    integer(c_int) :: ierr
    real(c_double), parameter :: dum = 0.0_dp
    real(c_double) :: tcur(1)
    
    if (.not.c_associated(self%wrk%sun%cvode_mem)) then
      err = "You must first initialize the stepper with 'initialize_stepper'"
      return 
    endif
    
    ierr = FCVode(self%wrk%sun%cvode_mem, dum, self%wrk%sun%sunvec_y, tcur, CV_ONE_STEP)
    if (ierr /= 0) then
      err = "CVODE step failed"
      return
    endif
    tn = tcur(1)
  end function
  
  module subroutine destroy_stepper(self, err)
    use iso_c_binding, only: c_int, c_associated, c_null_ptr
    use fcvode_mod, only: FCVodeFree
    use fsundials_nvector_mod, only: FN_VDestroy
    use fsundials_matrix_mod, only: FSUNMatDestroy
    use fsundials_linearsolver_mod, only: FSUNLinSolFree
    
    class(Atmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    
    call self%wrk%sun%finalize(err)
    if (allocated(err)) return
    
  end subroutine
  
end submodule