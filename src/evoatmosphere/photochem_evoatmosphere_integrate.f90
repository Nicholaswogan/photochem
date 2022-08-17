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
    loc_ierr = FCVodeGetNumSteps(self%wrk%sun%cvode_mem, nsteps)
    
    if (nsteps(1) /= self%wrk%nsteps_previous .and. self%var%verbose > 0) then
      loc_ierr = FCVodeGetCurrentStep(self%wrk%sun%cvode_mem, hcur)
      
      if (self%var%verbose == 1) then
        if (self%evolve_climate) then
          print"(1x,'N =',i6,2x,'Time = ',es11.5,2x,'dt = ',es11.5,"// &
               "2x,'max(dy/dt) = ',es11.5,2x,'T_surf = ',f11.5,2x,'trop_ind = ',i3,2x,'Ptop = ',es11.5)", &
              nsteps, tn, hcur(1),maxval(abs(fvec)),self%T_surf,self%var%trop_ind,self%wrk%pressure_hydro(self%var%nz)/1.0e6_dp
        else
          print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,'max(dy/dt) = ',es11.5)", &
              nsteps, tn, hcur(1),maxval(abs(fvec))
        endif
             
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

  function RootFn(tn, sunvec_y, gvec, user_data) result(ierr) bind(c,name='RootFn')
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    real(c_double), value :: tn        ! current time
    type(N_Vector) :: sunvec_y  ! solution N_Vector
    real(c_double) :: gvec(*)
    type(c_ptr), value :: user_data ! user-defined data
    integer(c_int) :: ierr

    real(c_double), pointer :: yvec(:)
    real(dp), pointer :: usol_in(:,:)
    type(EvoAtmosphere), pointer :: self
    character(:), allocatable :: err
    real(dp), parameter :: tol = 1.0e-2

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    
    ierr = 0

    call c_f_pointer(user_data, self)
    dat => self%dat
    var => self%var
    wrk => self%wrk
    yvec(1:var%neqs) => FN_VGetArrayPointer(sunvec_y)
    usol_in(1:dat%nq,1:var%nz) => yvec(1:var%neqs)
    
    call self%prep_atm_evo_gas(usol_in, wrk%usol, &
                               wrk%molecules_per_particle, wrk%pressure, wrk%density, wrk%mix, wrk%mubar, &
                               wrk%pressure_hydro, wrk%density_hydro, err)

    ! tropopause
    ! We only dynamically adjust the tropopause when climate is evolving
    if (self%evolve_climate) then
      gvec(1) = var%trop_alt - (var%z(var%trop_ind+1) + (0.5_dp+tol)*var%dz(var%trop_ind+1))
      gvec(2) = var%trop_alt - (var%z(var%trop_ind+1) - (0.5_dp+tol)*var%dz(var%trop_ind+1))
    else
      gvec(1) = 1.0_dp
      gvec(2) = 1.0_dp
    endif

    ! pressure at the top of the atmosphere
    gvec(3) = wrk%pressure_hydro(var%nz)/1.0e6_dp - self%P_top_min
    gvec(4) = wrk%pressure_hydro(var%nz)/1.0e6_dp - self%P_top_max

  end function
  
  module function evolve(self, filename, tstart, usol_start, t_eval, overwrite, err) result(success)
                                   
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, CV_ONE_STEP, FCVodeInit, FCVodeSVtolerances, &
                          FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                          FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                          FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd, &
                          FCVodeSetUserData, FCVodeRootInit, FCVodeReInit, FCVodeGetRootInfo
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
    use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    
    use photochem_enum, only: DensityBC
    use photochem_types, only: SundialsDataFinalizer
    
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
    
    ! solution vector, neq is set in the ode_mod module
    integer(c_int), parameter :: nrtfn = 4
    integer(c_int) :: rootsfound(nrtfn)
    real(c_double) :: usol_new(self%dat%nq,self%var%nz)
    real(c_double), pointer :: yvec_usol(:,:)
    integer(c_int64_t) :: neqs_long
    integer(c_int64_t) :: mu, ml
    integer(c_long) :: mxsteps_
    real(dp) :: new_atol
    integer :: error_reinit_attempts
    logical :: reinitialize
    
    integer :: i, j, k, ii, io
    
    type(SundialsDataFinalizer) :: sunfin
    type(c_ptr)    :: user_data
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    type(EvoAtmosphere), pointer :: self_ptr
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    ! The below association will make sure that all
    ! sundials data is destroyed after `sunfin`
    ! goes out of scope (when we leave this function)
    sunfin%sun => self%wrk%sun

    ! free memory if possible
    call wrk%sun%finalize(err)
    if (allocated(err)) return
    
    ! check dimensions
    if (size(usol_start,1) /= dat%nq .or. size(usol_start,2) /= var%nz) then
      err = "'usol_start' has the wrong dimensions"
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
    write(1) dat%species_names(1:dat%nq)
    write(1) size(t_eval)
    close(1)
    
    ! settings
    mxsteps_ = var%mxsteps
    neqs_long = var%neqs
    tcur = tstart
    mu = dat%nq
    ml = dat%nq
    new_atol = var%atol
    
    self_ptr => self
    user_data = c_loc(self_ptr)
    
    ! initialize solution vector
    allocate(wrk%sun%yvec(var%neqs))
    yvec_usol(1:dat%nq,1:var%nz) => wrk%sun%yvec
    do j=1,var%nz
      do i=1,dat%nq
        k = i + (j-1)*dat%nq
        wrk%sun%yvec(k) = usol_start(i,j)
      enddo
    enddo
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == DensityBC) then
        wrk%sun%yvec(i) = var%lower_fix_den(i)
      endif
    enddo
    ! set abstol
    allocate(wrk%sun%abstol(var%neqs))
    call self%set_trop_ind(yvec_usol, err)
    if (allocated(err)) return
    do j=1,var%nz
      do i=1,dat%nq
        k = i + (j-1)*dat%nq
        wrk%sun%abstol(k) = self%wrk%density_hydro(j)*var%atol
      enddo
    enddo

    wrk%sun%abstol_nvec => FN_VMake_Serial(neqs_long, wrk%sun%abstol)
    if (.not. associated(wrk%sun%abstol_nvec)) then
      err = "CVODE setup error."
      return
    end if

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
    
    ierr = FCVodeInit(wrk%sun%cvode_mem, c_funloc(RhsFn_evo), tstart, wrk%sun%sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSVtolerances(wrk%sun%cvode_mem, var%rtol, wrk%sun%abstol_nvec)
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
    
    ierr = FCVodeSetJacFn(wrk%sun%cvode_mem, c_funloc(JacFn_evo))
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

    ierr = FCVodeRootInit(wrk%sun%cvode_mem, nrtfn, c_funloc(RootFn))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if

    error_reinit_attempts = 0
    do ii = 1, size(t_eval)
      success = .false.
      do
        ierr = FCVode(wrk%sun%cvode_mem, t_eval(ii), wrk%sun%sunvec_y, tcur, CV_NORMAL)

        reinitialize = .false.
        if (any(ierr == [-1, -2, -3, -4])) then; block
          use photochem_const, only: small_real
          real(dp) :: rnd_num
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

          ! Try setting a new absolute tolerance. We select
          ! a number from a 2 order of magnitude log-normal 
          ! distribution centered at var%atol.
          call random_number(rnd_num)
          rnd_num = rnd_num*2.0_dp - 1.0_dp + log10(var%atol)
          new_atol = 10.0_dp**rnd_num

          error_reinit_attempts = error_reinit_attempts + 1
          reinitialize = .true.
        endblock; elseif (ierr <= -5) then
          ! bunch of errors that we can not recover from
          return
        elseif (ierr == 0 .or. ierr == 1 .or. ierr == 99) then
          ! Successful return. Go save the results, and continue integrating.
          exit
        elseif (ierr == 2) then; block
          real(dp) :: new_top_atmos
          ! root was found

          ierr = FCVodeGetRootInfo(wrk%sun%cvode_mem, rootsfound)
          if (ierr /= 0) then
            err = 'CVODE roots error.'
            return
          endif

          if (rootsfound(3) == -1) then
            ! pressure at the top of the atmosphere is going down
            ! we must decrease the top of the atmosphere
            new_top_atmos = (1.0_dp-self%top_atmos_adjust_frac)*var%top_atmos
            call self%rebin_update_vertical_grid(yvec_usol, new_top_atmos, usol_new, err)
            if (allocated(err)) return
            yvec_usol = usol_new

          elseif (rootsfound(4) == 1) then
            ! pressure at the top of the atmosphere is going up
            ! we must increase the top of the atmosphere
            new_top_atmos = (1.0_dp+self%top_atmos_adjust_frac)*var%top_atmos
            call self%rebin_update_vertical_grid(yvec_usol, new_top_atmos, usol_new, err)
            if (allocated(err)) return
            yvec_usol = usol_new

          elseif (rootsfound(1) /= 0 .or. rootsfound(2) /= 0) then
            ! tropopause index needs changing
          
            ! set the tropopause index
            call self%set_trop_ind(yvec_usol, err)
            if (allocated(err)) return
            
          endif

          reinitialize = .true.
        endblock; else
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
  
end submodule