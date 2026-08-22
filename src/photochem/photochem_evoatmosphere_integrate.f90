submodule(photochem_evoatmosphere) photochem_evoatmosphere_integrate
  implicit none
  
  ! Contains routines for integrating the photochemical equations
  ! forward in time. Here, we use the CVODE BDF integrator.
  
contains
  
  module function right_hand_side_callback(tn, sunvec_y, sunvec_f, user_data) &
                        result(ierr) bind(c, name='right_hand_side_callback')
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
    real(c_double), pointer :: yvec_raw(:)
    real(c_double), pointer :: fvec_raw(:)
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
    ! The SUNDIALS Fortran wrapper exposes the C array pointer as a
    ! one-element assumed-size pointer.  Do not rank-remap that pointer
    ! directly: bounds checking correctly rejects a remap larger than its
    ! declared target, and the resulting undefined behavior is compiler
    ! dependent.  Reconstruct the full view from the C address instead.
    yvec_raw => FN_VGetArrayPointer(sunvec_y)
    fvec_raw => FN_VGetArrayPointer(sunvec_f)
    call c_f_pointer(c_loc(yvec_raw(1)), yvec, [self%var%neqs])
    call c_f_pointer(c_loc(fvec_raw(1)), fvec, [self%var%neqs])
    
    ! fill RHS vector
    call self%right_hand_side(self%var%neqs, tn, yvec, fvec, err)
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
          if (abs(yvec(ii)) > self%var%atol) then
            tmp = abs(fvec(ii)/yvec(ii))
            if (tmp > mx) then
              mx = tmp
              k = ii
            endif
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
  
  module function jacobian_callback(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                            tmp1, tmp2, tmp3) &
                        result(ierr) bind(C,name='jacobian_callback')
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
    real(c_double), pointer :: yvec_raw(:)
    real(c_double), pointer :: Jmat(:)
    real(c_double), pointer :: Jmat_raw(:)
    character(:), allocatable :: err
    
    type(EvoAtmosphere), pointer :: self
  
    ierr = 0
    
    call c_f_pointer(user_data, self)
    yvec_raw => FN_VGetArrayPointer(sunvec_y)
    call c_f_pointer(c_loc(yvec_raw(1)), yvec, [self%var%neqs])
    Jmat_raw => FSUNBandMatrix_Data(sunmat_J)
    call c_f_pointer(c_loc(Jmat_raw(1)), Jmat, [self%var%neqs*self%dat%lda])
    call self%jacobian(self%dat%lda*self%var%neqs, self%var%neqs, yvec, Jmat, err)
    if (allocated(err)) then
      if (self%var%verbose > 0) then
        print*,trim(err)//". CVODE will attempt to correct the error."
      endif
      ierr = 1
    endif
    return
  
  end function

  subroutine error_handler_callback(error_code, module_, func, msg, eh_data) &
                                    bind(c, name='error_handler_callback')
    use iso_c_binding
    integer(c_int), value :: error_code
    character(kind=c_char) :: module_(*)
    character(kind=c_char) :: func(*)
    character(kind=c_char) :: msg(*)
    type(c_ptr), value, intent(in) :: eh_data
  end subroutine
  
  module function evolve(self, filename, tstart, usol_start, t_eval, overwrite, restart_from_file, err) result(success)
                                   
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use fcvode_mod, only: CV_NORMAL, FCVode, FCVodeSVtolerances, &
                          FCVodeSetInitStep, FCVodeReInit, FCVodeSetMaxStep
    use photochem_wrk, only: SundialsDataFinalizer
    
    ! in/out
    class(EvoAtmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    real(c_double), intent(inout) :: tstart
    real(dp), intent(inout) :: usol_start(:,:)
    real(c_double), intent(in) :: t_eval(:)
    logical, optional, intent(in) :: overwrite
    logical, optional, intent(in) :: restart_from_file
    logical :: success
    character(:), allocatable, intent(out) :: err
    
    ! local variables
    real(c_double) :: tcur(1)    ! current time
    integer(c_int) :: ierr       ! error flag from C functions
    
    real(c_double), pointer :: yvec_usol(:,:)
    real(dp) :: new_atol
    integer :: error_reinit_attempts
    logical :: reinitialize
    
    integer :: i, j, k, ii, io
    integer :: istart
    logical :: overwrite_, restart_from_file_
    
    type(SundialsDataFinalizer) :: sunfin
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    success = .false.
    call self%require_atmosphere_initialized('evolve', err)
    if (allocated(err)) return

    dat => self%dat
    var => self%var
    wrk => self%wrk
    ! The below association will make sure that all
    ! sundials data is destroyed after `sunfin`
    ! goes out of scope (when we leave this function)
    sunfin%sun => self%wrk%sun

    ! deal with optional arguments
    if (present(overwrite)) then
      overwrite_ = overwrite
    else
      overwrite_ = .false.
    endif
    if (present(restart_from_file)) then
      restart_from_file_ = restart_from_file
    else
      restart_from_file_ = .false.
    endif

    ! An evolution run owns its CVODE session and replaces any ordinary or
    ! robust stepper that was previously active.
    call self%destroy_stepper(err)
    if (allocated(err)) return
    
    ! check dimensions
    if (size(usol_start,1) /= dat%nq .or. size(usol_start,2) /= var%nz) then
      err = "'usol_start' has the wrong dimensions"
      return
    endif

    if (restart_from_file_) then; block
      real(dp) :: top_atmos
      integer :: restart_index

      ! read the file
      call read_end_of_evo_file(self, filename, t_eval, restart_index, tstart, top_atmos, usol_start, err)
      if (allocated(err)) return

      istart = restart_index
    endblock; else
      ! file prep
      if (overwrite_) then
        open(1, file = filename, status='replace', form="unformatted",iostat=io)
        if (io /= 0) then
          err = "Unable to replace "//trim(filename)
          return
        endif
      else
        open(1, file = filename, status='new', form="unformatted",iostat=io)
        if (io /= 0) then
          err = "Unable to create file "//trim(filename)//" because it already exists"
          return
        endif
      endif
      write(1) dat%nq
      write(1) var%nz
      write(1) dat%species_names(1:dat%nq)
      write(1) size(t_eval)
      close(1)
      istart = 1
    endif
    
    ! Construct the atmospheric integration state and common CVODE objects
    ! through the same path used by the ordinary and robust steppers.
    tcur = tstart
    new_atol = var%atol
    call self%initialize_stepper_at_time(usol_start, tstart, err)
    if (allocated(err)) return
    yvec_usol(1:dat%nq,1:var%nz) => wrk%sun%yvec

    error_reinit_attempts = 0
    do ii = istart, size(t_eval)
      success = .false.
      do
        ierr = FCVode(wrk%sun%cvode_mem, t_eval(ii), wrk%sun%sunvec_y, tcur, CV_NORMAL)

        reinitialize = .false.
        if (any(ierr == [-1, -2, -3, -4])) then; block
          use photochem_const, only: small_real
          use futils, only: linspace
          real(dp), allocatable :: atol_arr(:)
          ! -1 == The solver took mxstep internal steps but could not reach tout.
          ! -2 == The solver could not satisfy the accuracy demanded by the user
          ! for some internal step.
          ! -3 == Error test failures occurred too many times during one internal 
          ! time step or minimum step size was reached.
          ! -4 == Convergence test failures occurred too many times during one 
          ! internal time step or minimum step size was reached.
          !
          ! We will try to reinitialize a few times for all of these errors

          if (error_reinit_attempts >= var%max_error_reinit_attempts) then
            ! we give up
            return
          endif

          ! clip the results
          yvec_usol(:,:) = max(yvec_usol(:,:), small_real)

          ! Try setting a new absolute tolerance in the 
          ! vicinity of the old one.
          allocate(atol_arr(var%max_error_reinit_attempts))
          call linspace(log10(var%atol)-1.0_dp, log10(var%atol)+1.0_dp ,atol_arr)
          atol_arr = 10.0_dp**atol_arr
          new_atol = atol_arr(error_reinit_attempts+1)

          error_reinit_attempts = error_reinit_attempts + 1
          reinitialize = .true.
        endblock; elseif (ierr <= -5) then
          ! bunch of errors that we can not recover from
          return
        elseif (ierr == 0 .or. ierr == 1 .or. ierr == 99) then
          ! Successful return. Go save the results, and continue integrating.
          exit
        else
          ! in case we missed a scenario, then we assume its a failure
          err = 'Unknown CVODE return code'
          return
        endif

        if (reinitialize) then
          ierr = FCVodeReInit(wrk%sun%cvode_mem, tcur(1), wrk%sun%sunvec_y)
          if (ierr /= 0) then
            err = "CVODE reinit error."
            return
          endif

          do j=1,var%nz
            do i=1,dat%nq
              k = i + (j-1)*dat%nq
              wrk%sun%abstol(k) = wrk%density_hydro(j)*new_atol
            enddo
          enddo
          ierr = FCVodeSVtolerances(wrk%sun%cvode_mem, var%rtol, wrk%sun%abstol_nvec)
          if (ierr /= 0) then
            err = "CVODE setup error."
            return
          end if

          ierr = FCVodeSetInitStep(wrk%sun%cvode_mem, 0.0_c_double)
          if (ierr /= 0) then
            err = "CVODE setup error."
            return
          end if

          ierr = FCVodeSetMaxStep(wrk%sun%cvode_mem, var%max_dt)
          if (ierr /= 0) then
            err = "CVODE setup error."
            return
          end if

        endif
        
      enddo

      success = .true.

      call self%prep_atmosphere(yvec_usol, err)
      if (allocated(err)) return
      
      open(1, file = filename, status='old', form="unformatted",position="append")
      write(1) tcur(1)
      write(1) var%top_atmos
      write(1) var%z
      write(1) wrk%usol
      close(1)
    enddo
    
    ! free memory
    call wrk%sun%finalize(err)
    if (allocated(err)) return

  end function

  subroutine read_end_of_evo_file(self, filename, t_eval, restart_index, tcur, top_atmos, usol, err)
    use photochem_const, only: s_str_len
    use futils, only: FileCloser
    use iso_c_binding, only: c_double
    type(EvoAtmosphere), target, intent(in) :: self
    character(*), intent(in) :: filename
    real(dp), intent(in) :: t_eval(:)
    integer, intent(out) :: restart_index
    real(dp), intent(out) :: tcur, top_atmos
    real(dp), intent(out) :: usol(:,:)
    character(:), allocatable, intent(out) :: err

    integer :: io
    integer :: nq, nz
    character(s_str_len), allocatable :: species_names(:)
    integer :: nt
    real(dp), allocatable :: z(:)
    integer :: i
    real(dp) :: grid_scale
    type(FileCloser) :: file
    
    open(1, file = filename, status='old', form="unformatted",iostat=io)
    if (io /= 0) then
      err = 'Unable to open '//filename
      return
    endif
    file%unit = 1

    read(1,iostat=io) nq
    if (io /= 0) then
      err = 'Problem reading '//filename
      return
    endif
    if (nq /= self%dat%nq) then
      err = 'nq does not match EvoAtmosphere state in '//filename
      return
    endif

    read(1,iostat=io) nz
    if (io /= 0) then
      err = 'Problem reading '//filename
      return
    endif
    if (nz /= self%var%nz) then
      err = 'nz does not match EvoAtmosphere state in '//filename
      return
    endif

    allocate(species_names(nq))
    read(1,iostat=io) species_names
    if (io /= 0) then
      err = 'Problem reading '//filename
      return
    endif
    if (any(species_names /= self%dat%species_names)) then
      err = 'species_names does not match EvoAtmosphere state in '//filename
      return
    endif

    read(1,iostat=io) nt
    if (io /= 0) then
      err = 'Problem reading '//filename
      return
    endif

    ! check that there is data
    read(1,iostat=io) tcur
    if (io == -1) then
      err = 'There is no saved data in '//filename
      return
    endif
    backspace(1)

    allocate(z(nz))
    do i = 1,nt
      read(1,iostat=io) tcur
      if (io == -1) then
        ! end of file
        exit
      elseif (io < -1) then
        err = 'Problem reading '//filename
        return
      endif
      read(1,iostat=io) top_atmos
      if (io /= 0) then
        err = 'Problem reading '//filename
        return
      endif
      read(1,iostat=io) z
      if (io /= 0) then
        err = 'Problem reading '//filename
        return
      endif
      read(1,iostat=io) usol
      if (io /= 0) then
        err = 'Problem reading '//filename
        return
      endif
    enddo

    ! Evolution output is tied to the fixed grid that produced it. A restart
    ! must therefore use an atmosphere initialized with the same grid rather
    ! than silently changing the model state while reading the file.
    if (abs(top_atmos - self%var%top_atmos) > &
        1.0e-12_dp*max(1.0_dp, abs(self%var%top_atmos))) then
      err = 'The saved top-of-atmosphere altitude in '//trim(filename)// &
            ' does not match the initialized fixed grid.'
      return
    endif
    grid_scale = max(1.0_dp, maxval(abs(self%var%z)))
    if (maxval(abs(z - self%var%z)) > 1.0e-12_dp*grid_scale) then
      err = 'The saved altitude grid in '//trim(filename)// &
            ' does not match the initialized fixed grid.'
      return
    endif

    ! tcur is the end of the file
    restart_index = -1
    do i = 1,size(t_eval)
      if (t_eval(i) > tcur) then
        restart_index = i
        exit
      endif
    enddo

    if (restart_index == -1) then
      err = 'Was unable to find a time in t_eval that is greater than tcur in '//filename
      return
    endif
    
  end subroutine

  module function check_for_convergence(self, err) result(converged)
    use, intrinsic :: iso_c_binding
    class(EvoAtmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    logical :: converged

    integer :: i,j,ind
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    converged = .false.
    call self%require_atmosphere_initialized('check_for_convergence', err)
    if (allocated(err)) return

    dat => self%dat
    var => self%var
    wrk => self%wrk

    if (.not.c_associated(wrk%sun%cvode_mem)) then
      err = "You must first initialize the stepper with 'initialize_stepper'"
      return
    endif

    ! If we reach equilibrium time, then converged
    if (wrk%tn > var%equilibrium_time) then
      converged = .true.
      return
    endif

    ! Now consider step history.

    ! Can't do analysis on step 0
    if (wrk%nsteps == 0) return 

    ! Find index in history closest to time of interest. Note that this will only
    ! consider a limited step history. We can not save all history.
    ind = minloc(abs(wrk%t_history - var%conv_hist_factor*wrk%t_history(1)),1)

    ! Can't be current time, so we will check the previous step if needed.
    if (ind == 1) ind = 2 

    ! Compute difference between current mixing ratios, and mixing ratios
    ! at our index of interest.
    do j = 1,var%nz
      do i = 1,dat%nq
        if (wrk%mix_history(i,j,1) > var%conv_min_mix) then
          wrk%dmix(i,j) = abs(wrk%mix_history(i,j,1) - wrk%mix_history(i,j,ind))
        else 
          ! Ignore small mixing ratios
          wrk%dmix(i,j) = 0.0_dp
        endif
      enddo
    enddo

    ! Maximum normalized change
    wrk%longdy = maxval(abs(wrk%dmix/wrk%mix_history(:,:,1)))
    ! Also consider that change over time
    wrk%longdydt = wrk%longdy/(wrk%t_history(1) - wrk%t_history(ind))

    ! Check for convergence
    if (wrk%longdy < var%conv_longdy .and. wrk%longdydt < var%conv_longdydt) then
      converged = .true.
      return
    endif

  end function

  module subroutine initialize_stepper(self, usol_start, err)
    use, intrinsic :: iso_c_binding, only: c_associated
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_start(:,:)
    character(:), allocatable, intent(out) :: err

    call self%initialize_stepper_at_time(usol_start, 0.0_dp, err)
    if (.not.allocated(err) .or. .not.c_associated(self%wrk%sun%cvode_mem)) then
      self%wrk%robust_stepper_initialized = .false.
    endif

  end subroutine

  module subroutine initialize_stepper_at_time(self, usol_start, tstart, err, initial_step)
    use, intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fcvode_mod, only: CV_BDF, FCVodeInit, FCVodeSetLinearSolver, &
                          FCVodeCreate, FCVodeSetJacFn, FCVodeSetUserData
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_start(:,:)
    real(dp), intent(in) :: tstart
    character(:), allocatable, intent(out) :: err
    real(dp), optional, intent(in) :: initial_step
    
    integer(c_int) :: ierr       ! error flag from C functions
    integer(c_int64_t) :: neqs_long
    integer(c_int64_t) :: mu, ml
    type(c_ptr)    :: user_data
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    type(EvoAtmosphere), pointer :: self_ptr
    real(dp) :: initial_step_
    
    call self%require_atmosphere_initialized('initialize_stepper', err)
    if (allocated(err)) return

    dat => self%dat
    var => self%var
    wrk => self%wrk

    if (size(usol_start,1) /= dat%nq .or. size(usol_start,2) /= var%nz) then
      err = "Input 'usol_start' to 'initialize_stepper' is the wrong dimension"
      return
    endif
    if (.not.ieee_is_finite(tstart) .or. tstart < 0.0_dp) then
      err = "Input 'tstart' to the internal stepper initializer must be finite and nonnegative"
      return
    endif

    ! settings
    neqs_long = var%neqs
    mu = dat%nq
    ml = dat%nq
    self_ptr => self
    user_data = c_loc(self_ptr)
    initial_step_ = var%initial_dt
    if (present(initial_step)) initial_step_ = initial_step

    call wrk%sun%finalize(err)
    if (allocated(err)) return

    ! Allocate storage, then prepare the common atmospheric integration state.
    allocate(wrk%sun%yvec(var%neqs))
    allocate(wrk%sun%abstol(var%neqs))
    call self%prepare_stepper_state(usol_start, tstart, err)
    if (allocated(err)) then
      call cleanup_after_setup_failure()
      return
    endif

    wrk%sun%abstol_nvec => FN_VMake_Serial(neqs_long, wrk%sun%abstol)
    if (.not. associated(wrk%sun%abstol_nvec)) then
      err = "CVODE setup error while creating the absolute-tolerance vector."
      call cleanup_after_setup_failure()
      return
    end if

    ! create SUNDIALS N_Vector
    wrk%sun%sunvec_y => FN_VMake_Serial(neqs_long, wrk%sun%yvec)
    if (.not. associated(wrk%sun%sunvec_y)) then
      err = "CVODE setup error while creating the solution vector."
      call cleanup_after_setup_failure()
      return
    end if

    ! create CVode memory
    wrk%sun%cvode_mem = FCVodeCreate(CV_BDF)
    if (.not. c_associated(wrk%sun%cvode_mem)) then
      err = "CVODE setup error while creating CVODE memory."
      call cleanup_after_setup_failure()
      return
    end if
    
    ! set user data
    ierr = FCVodeSetUserData(wrk%sun%cvode_mem, user_data)
    if (ierr /= 0) then
      err = "CVODE setup error while setting user data."
      call cleanup_after_setup_failure()
      return
    end if
    
    ierr = FCVodeInit(wrk%sun%cvode_mem, c_funloc(right_hand_side_callback), tstart, wrk%sun%sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error while initializing CVODE."
      call cleanup_after_setup_failure()
      return
    end if
    
    wrk%sun%sunmat => FSUNBandMatrix(neqs_long, mu, ml)
    if (.not.associated(wrk%sun%sunmat)) then
      err = "CVODE setup error while creating the band matrix."
      call cleanup_after_setup_failure()
      return
    endif
    wrk%sun%sunlin => FSUNLinSol_Band(wrk%sun%sunvec_y, wrk%sun%sunmat)
    if (.not.associated(wrk%sun%sunlin)) then
      err = "CVODE setup error while creating the band linear solver."
      call cleanup_after_setup_failure()
      return
    endif
    
    ierr = FCVodeSetLinearSolver(wrk%sun%cvode_mem, wrk%sun%sunlin, wrk%sun%sunmat)
    if (ierr /= 0) then
      err = "CVODE setup error while attaching the linear solver."
      call cleanup_after_setup_failure()
      return
    end if
    
    ierr = FCVodeSetJacFn(wrk%sun%cvode_mem, c_funloc(jacobian_callback))
    if (ierr /= 0) then
      err = "CVODE setup error while setting the Jacobian function."
      call cleanup_after_setup_failure()
      return
    end if
    
    call self%configure_stepper(initial_step_, err)
    if (allocated(err)) then
      call cleanup_after_setup_failure()
      return
    endif

  contains

    subroutine cleanup_after_setup_failure()
      character(:), allocatable :: cleanup_err
      call wrk%sun%finalize(cleanup_err)
      if (allocated(cleanup_err)) then
        err = err//" Cleanup also failed: "//cleanup_err
      endif
    end subroutine

  end subroutine

  module subroutine prepare_stepper_state(self, usol_start, tstart, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_start(:,:)
    real(dp), intent(in) :: tstart
    character(:), allocatable, intent(out) :: err

    real(dp), pointer :: yvec_usol(:,:)
    integer :: i, j, k
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    dat => self%dat
    var => self%var
    wrk => self%wrk

    if (.not.allocated(wrk%sun%yvec) .or. .not.allocated(wrk%sun%abstol)) then
      err = "Internal stepper state storage is not allocated."
      return
    endif
    if (size(wrk%sun%yvec) /= var%neqs .or. size(wrk%sun%abstol) /= var%neqs) then
      err = "Internal stepper state storage has the wrong dimensions."
      return
    endif

    yvec_usol(1:dat%nq,1:var%nz) => wrk%sun%yvec
    yvec_usol = usol_start
    call self%apply_lower_boundary_conditions(var%temperature(1), yvec_usol(:,1), err)
    if (allocated(err)) return

    call self%prepare_atmosphere_structure( &
      yvec_usol, wrk%usol, wrk%molecules_per_particle, wrk%pressure, &
      wrk%density, wrk%mix, wrk%mubar, wrk%pressure_hydro, &
      wrk%density_hydro, err=err &
    )
    if (allocated(err)) return
    do j=1,var%nz
      do i=1,dat%nq
        k = i + (j-1)*dat%nq
        wrk%sun%abstol(k) = wrk%density_hydro(j)*var%atol
      enddo
    enddo

    ! A restart begins a new CVODE and convergence-history segment without
    ! changing the integration's logical time or robust-session totals.
    wrk%nsteps = 0
    wrk%nsteps_previous = -10
    wrk%tn = tstart
    wrk%t_history = -1.0_dp
    wrk%t_history(1) = tstart
    wrk%mix_history = -1.0_dp
    wrk%mix_history(:,:,1) = wrk%mix

  end subroutine

  module subroutine configure_stepper(self, initial_step, err)
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long, &
                                           c_funloc, c_null_ptr
    use fcvode_mod, only: FCVodeSVtolerances, FCVodeSetMaxNumSteps, &
                          FCVodeSetInitStep, FCVodeSetMaxStep, &
                          FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd, &
                          FCVodeSetErrHandlerFn
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: initial_step
    character(:), allocatable, intent(out) :: err

    integer(c_int) :: ierr
    integer(c_long) :: mxsteps_

    mxsteps_ = self%var%mxsteps
    ierr = FCVodeSVtolerances(self%wrk%sun%cvode_mem, self%var%rtol, &
                             self%wrk%sun%abstol_nvec)
    if (ierr /= 0) then
      err = "CVODE setup error while setting tolerances."
      return
    endif
    ierr = FCVodeSetMaxNumSteps(self%wrk%sun%cvode_mem, mxsteps_)
    if (ierr /= 0) then
      err = "CVODE setup error while setting the maximum number of steps."
      return
    endif
    ierr = FCVodeSetInitStep(self%wrk%sun%cvode_mem, real(initial_step,c_double))
    if (ierr /= 0) then
      err = "CVODE setup error while setting the initial step."
      return
    endif
    ierr = FCVodeSetMaxStep(self%wrk%sun%cvode_mem, self%var%max_dt)
    if (ierr /= 0) then
      err = "CVODE setup error while setting the maximum step."
      return
    endif
    ierr = FCVodeSetMaxErrTestFails(self%wrk%sun%cvode_mem, self%var%max_err_test_failures)
    if (ierr /= 0) then
      err = "CVODE setup error while setting error-test failures."
      return
    endif
    ierr = FCVodeSetMaxOrd(self%wrk%sun%cvode_mem, self%var%max_order)
    if (ierr /= 0) then
      err = "CVODE setup error while setting the maximum order."
      return
    endif

    if (self%var%verbose == 0) then
      ierr = FCVodeSetErrHandlerFn(self%wrk%sun%cvode_mem, &
                                   c_funloc(error_handler_callback), c_null_ptr)
      if (ierr /= 0) then
        err = "CVODE setup error while setting the error handler."
        return
      endif
    endif

  end subroutine

  module subroutine restart_robust_stepper(self, usol_restart, tstart, err)
    use, intrinsic :: iso_c_binding, only: c_associated, c_int
    use fcvode_mod, only: FCVodeReInit
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_restart(:,:)
    real(dp), intent(in) :: tstart
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: usol_clipped(:,:)
    real(dp) :: restart_initial_step
    character(:), allocatable :: reinit_err
    integer(c_int) :: ierr
    logical :: can_reinit, attempted_reinit
    type(PhotochemWrk), pointer :: wrk

    wrk => self%wrk
    usol_clipped = max(usol_restart, self%var%reinit_min_density)
    ! Preserve the configured restart behavior unless initial_dt is too small
    ! to advance floating-point time at the current absolute time.
    restart_initial_step = max(self%var%initial_dt, 2.0_dp*spacing(tstart))
    attempted_reinit = .false.
    can_reinit = c_associated(wrk%sun%cvode_mem) .and. &
                 allocated(wrk%sun%yvec) .and. allocated(wrk%sun%abstol) .and. &
                 associated(wrk%sun%sunvec_y) .and. associated(wrk%sun%abstol_nvec) .and. &
                 associated(wrk%sun%sunmat) .and. associated(wrk%sun%sunlin)
    if (can_reinit) then
      can_reinit = size(wrk%sun%yvec) == self%var%neqs .and. &
                   size(wrk%sun%abstol) == self%var%neqs
    endif

    if (can_reinit) then
      attempted_reinit = .true.
      call self%prepare_stepper_state(usol_clipped, tstart, err)
      if (.not.allocated(err)) then
        ierr = FCVodeReInit(wrk%sun%cvode_mem, tstart, wrk%sun%sunvec_y)
        if (ierr /= 0) err = "CVodeReInit returned an error."
      endif
      if (.not.allocated(err)) call self%configure_stepper(restart_initial_step, err)
      if (.not.allocated(err)) return
      reinit_err = err
      deallocate(err)
    endif

    ! Missing infrastructure or an unsuccessful in-place restart requires a
    ! clean reconstruction. The robust-session counters remain untouched.
    call self%initialize_stepper_at_time(usol_clipped, tstart, err, &
                                         initial_step=restart_initial_step)
    if (allocated(err) .and. attempted_reinit) then
      err = "In-place CVODE restart failed ("//reinit_err// &
            "); full reconstruction also failed: "//err
    endif

  end subroutine

  module function step(self, err) result(tn)
    use iso_c_binding, only: c_null_ptr, c_int, c_double, c_associated, c_long
    use fcvode_mod, only: CV_ONE_STEP, FCVode, FCVodeGetNumSteps
    class(EvoAtmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    real(dp) :: tn
    
    integer(c_int) :: ierr
    integer(c_long) :: nsteps_(1)
    integer :: i, k
    real(c_double) :: tout
    real(c_double) :: tcur(1)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    tn = 0.0_dp
    call self%require_atmosphere_initialized('step', err)
    if (allocated(err)) return

    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (.not.c_associated(self%wrk%sun%cvode_mem)) then
      err = "You must first initialize the stepper with 'initialize_stepper'"
      return 
    endif

    ! CV_ONE_STEP still uses tout to establish integration direction on its
    ! first call after CVodeInit/CVodeReInit. Keep it safely forward of the
    ! committed time, including after a restart at nonzero time.
    if (wrk%t_history(1) <= 0.5_dp*huge(1.0_dp)) then
      tout = wrk%t_history(1) + max(1.0_dp, abs(wrk%t_history(1)))
    else
      tout = huge(1.0_dp)
    endif
    ierr = FCVode(self%wrk%sun%cvode_mem, tout, self%wrk%sun%sunvec_y, tcur, CV_ONE_STEP)
    if (ierr /= 0) then
      err = "CVODE step failed"
      return
    endif
    tn = tcur(1)

    block
      real(c_double), pointer :: usol_tmp(:,:)
      usol_tmp(1:dat%nq,1:var%nz) => wrk%sun%yvec
      ! RHS and Jacobian callbacks may have prepared var for trial states.
      ! Prepare the accepted CVODE solution before this public call returns.
      call self%prepare_atmosphere_structure(usol_tmp, wrk%usol, &
           wrk%molecules_per_particle, wrk%pressure, wrk%density, wrk%mix, wrk%mubar, &
           wrk%pressure_hydro, wrk%density_hydro, err=err)
      if (allocated(err)) return
    endblock

    ! Commit counters and convergence history only after both CVODE and the
    ! atmospheric-state preparation have succeeded.
    ierr = FCVodeGetNumSteps(wrk%sun%cvode_mem, nsteps_)
    if (ierr /= 0) then
      err = "Unable to obtain the CVODE step count"
      return
    endif
    wrk%nsteps = nsteps_(1)
    wrk%tn = tn
    k = min(wrk%nsteps+1,size(wrk%t_history))
    do i = k,2,-1
      wrk%t_history(i) = wrk%t_history(i-1)
      wrk%mix_history(:,:,i) = wrk%mix_history(:,:,i-1)
    enddo
    wrk%t_history(1) = tn
    wrk%mix_history(:,:,1) = wrk%mix

  end function

  module subroutine destroy_stepper(self, err)
    use iso_c_binding, only: c_int, c_associated, c_null_ptr
    use fcvode_mod, only: FCVodeFree
    use fsundials_nvector_mod, only: FN_VDestroy
    use fsundials_matrix_mod, only: FSUNMatDestroy
    use fsundials_linearsolver_mod, only: FSUNLinSolFree
    
    class(EvoAtmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    
    call self%wrk%sun%finalize(err)
    self%wrk%robust_stepper_initialized = .false.
    if (allocated(err)) return
    
  end subroutine

  module subroutine validate_robust_stepper_settings(self, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), intent(in) :: self
    character(:), allocatable, intent(out) :: err

    if (self%var%nerrors_before_giveup < 1) then
      err = "`nerrors_before_giveup` must be positive"
    elseif (self%var%nsteps_before_conv_check < 0) then
      err = "`nsteps_before_conv_check` must be nonnegative"
    elseif (self%var%nsteps_before_reinit < 1) then
      err = "`nsteps_before_reinit` must be positive"
    elseif (self%var%nsteps_before_giveup < 1) then
      err = "`nsteps_before_giveup` must be positive"
    elseif (self%var%nsteps_before_conv_check >= self%var%nsteps_before_reinit) then
      err = "`nsteps_before_conv_check` must be less than `nsteps_before_reinit`"
    elseif (.not.ieee_is_finite(self%var%reinit_min_density) .or. &
            self%var%reinit_min_density <= 0.0_dp) then
      err = "`reinit_min_density` must be finite and positive"
    elseif (self%var%toa_pressure_maintenance%enabled .and. &
            .not.self%var%press_temp_edd_profile%enabled) then
      err = "TOA-pressure maintenance requires an enabled persistent pressure-based temperature and eddy-diffusion profile"
    elseif (self%var%toa_pressure_maintenance%enabled .and. &
            (.not.ieee_is_finite(self%var%toa_pressure_maintenance%target_pressure) .or. &
             self%var%toa_pressure_maintenance%target_pressure <= 0.0_dp)) then
      err = "`toa_pressure_maintenance%target_pressure` must be finite and positive"
    elseif (self%var%toa_pressure_maintenance%enabled .and. &
            (.not.ieee_is_finite(self%var%toa_pressure_maintenance%pressure_factor) .or. &
             self%var%toa_pressure_maintenance%pressure_factor < 1.0_dp)) then
      err = "`toa_pressure_maintenance%pressure_factor` must be finite and at least one"
    elseif (self%var%toa_pressure_maintenance%enabled .and. &
            self%var%toa_pressure_maintenance%nsteps_between_updates < 1) then
      err = "`toa_pressure_maintenance%nsteps_between_updates` must be positive"
    elseif (self%var%toa_pressure_maintenance%enabled .and. &
            self%var%toa_pressure_maintenance%max_failures < 0) then
      err = "`toa_pressure_maintenance%max_failures` must be nonnegative"
    endif

  end subroutine

  module subroutine initialize_robust_stepper(self, usol_start, err)
    use, intrinsic :: iso_c_binding, only: c_associated
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_start(:,:)
    character(:), allocatable, intent(out) :: err
    real(dp) :: current_pressure, pressure_ratio
    real(dp), allocatable :: usol_preflight(:,:)
    logical :: initial_toa_update

    call self%require_atmosphere_initialized('initialize_robust_stepper', err)
    if (allocated(err)) return

    call self%validate_robust_stepper_settings(err)
    if (allocated(err)) return

    if (size(usol_start,1) /= self%dat%nq .or. size(usol_start,2) /= self%var%nz) then
      err = "Input 'usol_start' to 'initialize_robust_stepper' is the wrong dimension"
      return
    endif

    ! Prepare the supplied starting composition before checking the TOA. This
    ! makes the preflight pressure correspond to the state that will actually
    ! enter CVODE, including persistent-profile and boundary-condition logic.
    initial_toa_update = .false.
    if (self%var%toa_pressure_maintenance%enabled) then
      ! Keep the caller's array independent of the work allocation: a
      ! successful vertical-grid transaction moves that allocation.
      usol_preflight = usol_start
      call self%prep_atmosphere(usol_preflight, err)
      if (allocated(err)) return

      current_pressure = self%wrk%pressure(self%var%nz)
      if (.not.ieee_is_finite(current_pressure) .or. current_pressure <= 0.0_dp) then
        err = 'Initial TOA-pressure maintenance failed: the current TOA pressure '// &
              'was not finite and positive.'
        return
      endif
      pressure_ratio = current_pressure / &
                       self%var%toa_pressure_maintenance%target_pressure
      if (pressure_ratio < 1.0_dp / self%var%toa_pressure_maintenance%pressure_factor .or. &
          pressure_ratio > self%var%toa_pressure_maintenance%pressure_factor) then
        call self%update_vertical_grid( &
             TOA_pressure=self%var%toa_pressure_maintenance%target_pressure, err=err)
        if (allocated(err)) then
          err = 'Initial TOA-pressure maintenance failed: '//err
          return
        endif
        initial_toa_update = .true.
      endif
    endif

    ! A successful preflight regrid remaps the supplied starting state into
    ! the new grid and may move the original work allocation. Use the
    ! committed remapped state in that case; otherwise preserve the caller's
    ! exact initial array as the CVODE starting state.
    if (initial_toa_update) then
      call self%initialize_stepper_at_time(self%wrk%usol, 0.0_dp, err)
    else
      call self%initialize_stepper_at_time(usol_start, 0.0_dp, err)
    endif
    if (allocated(err)) then
      if (.not.c_associated(self%wrk%sun%cvode_mem)) then
        self%wrk%robust_stepper_initialized = .false.
      endif
      return
    endif

    self%wrk%nsteps_total = 0
    self%wrk%nerrors_total = 0
    self%wrk%n_toa_pressure_updates = merge(1, 0, initial_toa_update)
    self%wrk%n_toa_pressure_failures = 0
    self%wrk%nsteps_since_toa_pressure_update = 0
    self%wrk%robust_stepper_initialized = .true.

  end subroutine

  module subroutine robust_step(self, give_up, converged, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    logical, intent(out) :: give_up
    logical, intent(out) :: converged
    character(:), allocatable, intent(out) :: err

    real(dp) :: tn
    real(dp) :: t_committed
    real(dp) :: usol_committed(self%dat%nq,self%var%nz)
    character(:), allocatable :: cleanup_err
    logical :: chemistry_converged, toa_updated, toa_failed

    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    converged = .false.
    give_up = .false.

    call self%require_atmosphere_initialized('robust_step', err)
    if (allocated(err)) return

    var => self%var
    wrk => self%wrk

    if (.not.wrk%robust_stepper_initialized) then
      err = "You must first initialize a robust stepper with 'initialize_robust_stepper'"
      return
    endif
    call self%validate_robust_stepper_settings(err)
    if (allocated(err)) return

    ! RHS and atmospheric preparation use shared workspace, so retain an
    ! explicit snapshot until the whole attempted step has committed.
    usol_committed = wrk%usol
    t_committed = wrk%t_history(1)
    tn = self%step(err)
    if (.not.allocated(err)) then
      ! If step worked, then we add it to counter
      wrk%nsteps_total = wrk%nsteps_total + 1
      wrk%nsteps_since_toa_pressure_update = wrk%nsteps_since_toa_pressure_update + 1
    else
      ! There was an error
      deallocate(err)
      wrk%nerrors_total = wrk%nerrors_total + 1

      ! If there are too many errors, then give up
      if (wrk%nerrors_total > var%nerrors_before_giveup) then
        wrk%usol = usol_committed
        wrk%tn = t_committed
        call self%destroy_stepper(cleanup_err)
        if (allocated(cleanup_err)) err = cleanup_err
        give_up = .true.
        return
      endif

      ! Recover from the last committed state and time. Do not use the failed
      ! call's returned time and do not run convergence logic on this call.
      call self%restart_robust_stepper(usol_committed, t_committed, err)
      if (allocated(err)) then
        wrk%robust_stepper_initialized = .false.
        return
      endif
      return

    endif

    ! Determine the chemistry-only convergence result first. TOA maintenance
    ! below may invalidate that result by changing the vertical grid.
    chemistry_converged = .false.
    if (tn > var%equilibrium_time) then
      chemistry_converged = .true.
    elseif (self%wrk%nsteps > var%nsteps_before_conv_check) then
      chemistry_converged = self%check_for_convergence(err)
      if (allocated(err)) return
    endif

    ! An accepted step can be chemically converged while the model top has
    ! drifted outside the requested pressure range. Maintenance must happen
    ! before reporting convergence, and a successful regrid starts a fresh
    ! segment-local convergence history.
    call self%maybe_maintain_toa_pressure(chemistry_converged, toa_updated, toa_failed, err)
    if (allocated(err)) return
    if (toa_failed) return
    if (toa_updated) return
    if (chemistry_converged) then
      converged = .true.
      return
    endif

    ! The total-step ceiling counts accepted steps exactly. Check it before a
    ! scheduled restart so a terminal call does not rebuild/reset CVODE state.
    if (wrk%nsteps_total >= var%nsteps_before_giveup) then
      give_up = .true.
      return
    endif

    ! Reinitialize after exactly this many accepted steps in the current
    ! segment. Restarting resets local CVODE and convergence history only.
    if (self%wrk%nsteps >= var%nsteps_before_reinit) then
      call self%restart_robust_stepper(wrk%usol, wrk%t_history(1), err)
      if (allocated(err)) then
        wrk%robust_stepper_initialized = .false.
        return
      endif
    endif

  end subroutine

  module subroutine maybe_maintain_toa_pressure(self, chemistry_converged, updated, failed, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(inout) :: self
    logical, intent(in) :: chemistry_converged
    logical, intent(out) :: updated
    logical, intent(out) :: failed
    character(:), allocatable, intent(out) :: err

    real(dp) :: current_pressure, pressure_ratio, t_current
    character(:), allocatable :: failure_message
    integer :: nsteps_total, nerrors_total
    integer :: nupdates, nfailures

    updated = .false.
    failed = .false.
    if (.not.self%var%toa_pressure_maintenance%enabled) return

    nfailures = self%wrk%n_toa_pressure_failures

    current_pressure = self%wrk%pressure(self%var%nz)
    if (.not.ieee_is_finite(current_pressure) .or. current_pressure <= 0.0_dp) then
      self%wrk%n_toa_pressure_failures = nfailures + 1
      failed = .true.
      if (self%wrk%n_toa_pressure_failures > &
          self%var%toa_pressure_maintenance%max_failures) then
        err = 'TOA-pressure maintenance failed (failure limit exceeded): '// &
              'the current TOA pressure was not finite and positive.'
      endif
      return
    endif
    pressure_ratio = current_pressure / &
                     self%var%toa_pressure_maintenance%target_pressure
    if (pressure_ratio >= 1.0_dp / self%var%toa_pressure_maintenance%pressure_factor .and. &
        pressure_ratio <= self%var%toa_pressure_maintenance%pressure_factor) return
    if (self%wrk%nsteps_since_toa_pressure_update < &
        self%var%toa_pressure_maintenance%nsteps_between_updates .and. &
        .not.chemistry_converged) return

    ! update_vertical_grid intentionally installs a fresh work structure, so
    ! preserve robust-session totals and maintenance counters across it.
    nsteps_total = self%wrk%nsteps_total
    nerrors_total = self%wrk%nerrors_total
    nupdates = self%wrk%n_toa_pressure_updates
    t_current = self%wrk%tn

    call self%update_vertical_grid( &
         TOA_pressure=self%var%toa_pressure_maintenance%target_pressure, err=err)
    if (allocated(err)) then
      ! The vertical-grid transaction is failure atomic, so the active solver
      ! and committed state remain available for the caller.
      self%wrk%n_toa_pressure_failures = nfailures + 1
      failed = .true.
      if (self%wrk%n_toa_pressure_failures > &
          self%var%toa_pressure_maintenance%max_failures) then
        failure_message = err
        deallocate(err)
        err = 'TOA-pressure maintenance failed (failure limit exceeded): '// &
              failure_message
      else
        deallocate(err)
      endif
      return
    endif

    self%wrk%nsteps_total = nsteps_total
    self%wrk%nerrors_total = nerrors_total
    self%wrk%n_toa_pressure_updates = nupdates + 1
    self%wrk%n_toa_pressure_failures = nfailures
    self%wrk%nsteps_since_toa_pressure_update = 0

    ! The successful regrid invalidates the old CVODE infrastructure. The
    ! restart helper reconstructs compatible resources and preserves t_current.
    call self%restart_robust_stepper(self%wrk%usol, t_current, err)
    if (allocated(err)) then
      self%wrk%robust_stepper_initialized = .false.
      err = 'TOA-pressure maintenance regrid succeeded, but CVODE restart failed: '//err
      return
    endif
    self%wrk%robust_stepper_initialized = .true.
    updated = .true.

  end subroutine

  module function find_steady_state(self, err) result(converged)
    class(EvoAtmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    logical :: converged

    logical :: give_up

    converged = .false.

    call self%require_atmosphere_initialized('find_steady_state', err)
    if (allocated(err)) return

    call self%initialize_robust_stepper(self%wrk%usol, err)
    if (allocated(err)) return

    do
      call self%robust_step(give_up, converged, err)
      if (allocated(err)) return

      if (give_up) then
        converged = .false.
        return
      endif

      if (converged) return
    enddo

  end function
  
end submodule
