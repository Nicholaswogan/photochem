submodule(photochem_evoatmosphere) photochem_evoatmosphere_integrate
  implicit none
  
  ! Contains routines for integrating the photochemical equations
  ! forward in time. Here, we use the CVODE BDF integrator.
  
contains
  
  module function RhsFn_evo(tn, sunvec_y, sunvec_f, user_data) &
                        result(ierr) bind(c, name='RhsFn_evo')
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
    
    type(EvoAtmosphere), pointer :: self
  
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
      print*,trim(err)//". CVODE will attempt to correct the error."
      ierr = 1
    endif
    return
  end function
  
  module function JacFn_evo(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                        tmp1, tmp2, tmp3) &
                        result(ierr) bind(C,name='JacFn_evo')
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
    
    type(EvoAtmosphere), pointer :: self
  
    ierr = 0
    
    call c_f_pointer(user_data, self)
    yvec(1:self%var%neqs) => FN_VGetArrayPointer(sunvec_y)
    Jmat(1:self%var%neqs*self%dat%lda) => FSUNBandMatrix_Data(sunmat_J)
    call self%jacobian(self%dat%lda*self%var%neqs, self%var%neqs, yvec, Jmat, err)
    if (allocated(err)) then
      print*,trim(err)//". CVODE will attempt to correct the error."
      ierr = 1
    endif
    return
  
  end function
  
  module function evolve(self, filename, tstart, usol_start, t_eval, overwrite, err) result(success)
                                   
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, CV_ONE_STEP, FCVodeInit, FCVodeSVtolerances, &
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
    
    use photochem_enum, only: DensityBC
    
    ! in/out
    class(EvoAtmosphere), target, intent(inout) :: self
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
    type(N_Vector), pointer :: sunvec_y ! sundials vector
    
    ! solution vector, neq is set in the ode_mod module
    real(c_double) :: yvec(self%var%neqs)
    real(c_double) :: abstol(self%var%neqs), density
    type(N_Vector), pointer :: abstol_nvec
    integer(c_int64_t) :: neqs_long
    integer(c_int64_t) :: mu, ml
    integer(c_long) :: mxsteps_
    type(SUNMatrix), pointer :: sunmat
    type(SUNLinearSolver), pointer :: sunlin
    
    ! real(dp) :: fH2O(self%var%trop_ind)
    integer :: i, j, k, ii, io
    
    type(c_ptr)    :: user_data
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    type(EvoAtmosphere), pointer :: self_ptr
  
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (c_associated(self%wrk%cvode_mem)) then
      err = "You have a time stepper initalize. To do other integrations"// &
            " you need to delete it with 'destroy_stepper' subroutine."
      return 
    endif
    
    ! check dimensions
    if (size(usol_start,1) /= dat%nq .or. size(usol_start,2) /= var%nz) then
      err = "Problem!"
      return
    endif
    
    ! file prep
    if (overwrite) then
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
    tcur = tstart
    mu = dat%nq
    ml = dat%nq
    
    self_ptr => self
    user_data = c_loc(self_ptr)
    
    ! initialize solution vector
    do j=1,var%nz
      density = sum(usol_start(:,j))
      do i=1,dat%nq
        k = i + (j-1)*dat%nq
        yvec(k) = usol_start(i,j)
        ! set abstol.
        abstol(k) = density*var%atol
      enddo
    enddo
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == DensityBC) then
        yvec(i) = var%lower_fix_den(i)
      endif
    enddo

    abstol_nvec => FN_VMake_Serial(neqs_long, abstol)
    if (.not. associated(abstol_nvec)) then
      err = "CVODE setup error."
      return
    end if

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
    
    ierr = FCVodeInit(wrk%cvode_mem, c_funloc(RhsFn_evo), tstart, sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSVtolerances(wrk%cvode_mem, var%rtol, abstol_nvec)
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
    
    ierr = FCVodeSetJacFn(wrk%cvode_mem, c_funloc(JacFn_evo))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
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
        exit
      else
        success = .true.
        do j=1,var%nz
          do i=1,dat%nq  
            k = i + (j-1)*dat%nq
            wrk%usol(i,j) = yvec(k)
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
    call FN_VDestroy(sunvec_y)
    call FCVodeFree(wrk%cvode_mem)
    ierr = FSUNLinSolFree(sunlin)
    if (ierr /= 0) then
      err = "CVODE deallocation error"
      return
    end if
    call FSUNMatDestroy(sunmat)

  end function
  
  module subroutine photochemical_equilibrium(self, success, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    logical, intent(out) :: success
    character(:), allocatable, intent(out) :: err 
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    
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
                          FCVodeSetUserData
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
    use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    
    use photochem_enum, only: DensityBC
    
    class(EvoAtmosphere), target, intent(inout) :: self
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
    type(PhotochemWrkEvo), pointer :: wrk
    type(EvoAtmosphere), pointer :: self_ptr
    
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
    wrk%tcur = tstart
    mu = dat%nq
    ml = dat%nq
    self_ptr => self
    user_data = c_loc(self_ptr)
    
    call self%destroy_stepper(err)
    if (allocated(err)) return
    
    ! initialize solution vector
    allocate(wrk%yvec(var%neqs))
    do j=1,var%nz
      do i=1,dat%nq
        k = i + (j-1)*dat%nq
        wrk%yvec(k) = usol_start(i,j)
      enddo
    enddo
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == DensityBC) then
        wrk%yvec(i) = var%lower_fix_den(i)
      endif
    enddo
    
    ! create SUNDIALS N_Vector
    wrk%sunvec_y => FN_VMake_Serial(neqs_long, wrk%yvec)
    if (.not. associated(wrk%sunvec_y)) then
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
    
    ierr = FCVodeInit(wrk%cvode_mem, c_funloc(RhsFn_evo), tstart, wrk%sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSStolerances(wrk%cvode_mem, var%rtol, var%atol)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    wrk%sunmat => FSUNBandMatrix(neqs_long, mu, ml)
    wrk%sunlin => FSUNLinSol_Band(wrk%sunvec_y,wrk%sunmat)
    
    ierr = FCVodeSetLinearSolver(wrk%cvode_mem, wrk%sunlin, wrk%sunmat)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetJacFn(wrk%cvode_mem, c_funloc(JacFn_evo))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
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
    
  end subroutine
  
  module function step(self, err) result(tn)
    use iso_c_binding, only: c_null_ptr, c_int, c_double, c_associated
    use fcvode_mod, only: CV_ONE_STEP, FCVode
    class(EvoAtmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    real(dp) :: tn
    
    integer(c_int) :: ierr
    real(c_double), parameter :: dum = 0.0_dp
    
    if (.not.c_associated(self%wrk%cvode_mem)) then
      err = "You must first initialize the stepper with 'initialize_stepper'"
      return 
    endif
    
    ierr = FCVode(self%wrk%cvode_mem, dum, self%wrk%sunvec_y, self%wrk%tcur, CV_ONE_STEP)
    if (ierr /= 0) then
      err = "CVODE step failed"
      return
    endif
    tn = self%wrk%tcur(1)
  end function
  
  module subroutine destroy_stepper(self, err)
    use iso_c_binding, only: c_int, c_associated, c_null_ptr
    use fcvode_mod, only: FCVodeFree
    use fsundials_nvector_mod, only: FN_VDestroy
    use fsundials_matrix_mod, only: FSUNMatDestroy
    use fsundials_linearsolver_mod, only: FSUNLinSolFree
    
    class(EvoAtmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    
    integer(c_int) :: ierr
    
    
    if (c_associated(self%wrk%cvode_mem)) then
      call FN_VDestroy(self%wrk%sunvec_y)
      call FCVodeFree(self%wrk%cvode_mem)
      self%wrk%cvode_mem = c_null_ptr
      ierr = FSUNLinSolFree(self%wrk%sunlin)
      if (ierr /= 0) then
        err = "CVODE deallocation error"
        return
      end if
      call FSUNMatDestroy(self%wrk%sunmat)
      deallocate(self%wrk%yvec)
    endif
    
  end subroutine
  
end submodule