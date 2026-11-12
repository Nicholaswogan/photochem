program ptc_benchmark
  use photochem, only: EvoAtmosphere, EvoAtmosphereState, ProductionLoss, dp
  use photochem_evoatmosphere_ptc, only: EvoAtmospherePTC
  use photochem_enum, only: PeriodicPressTempEdd
  use fcvode_mod, only: FCVodeGetNumRhsEvals, FCVodeGetNumJacEvals, FCVodeGetNumErrTestFails
  use iso_c_binding, only: c_long,c_int
  implicit none
  type(EvoAtmosphere), target :: pc
  type(EvoAtmosphereState) :: state
  type(EvoAtmospherePTC) :: ptc
  type(ProductionLoss) :: budget
  character(2048) :: inputs,data,mode,output,arg
  character(:), allocatable :: err
  character(512) :: reason='step_limit'
  real(dp) :: start,now,limit,dt,time,budget_error,rhs_column,normalizer,relative_error
  real(dp), allocatable :: atom_change(:)
  real(dp), allocatable :: y(:,:),f(:)
  integer :: n,steps,max_steps,unit,budget_unit,i,j,rejected,nrhs,njac,maint
  integer(c_long) :: crhs(1),cjac(1),crej(1)
  integer(c_int) :: status
  logical :: finished,give_up,is_ptc
  integer :: nsolves,clips
  real(dp) :: max_clip
  real(dp) :: atol,rtol,equilibrium_time,smooth_factor,k_cond,k_evap
  logical :: upwind_molec_diff
  integer :: jacobian_method
  namelist /benchmark_controls/ atol,rtol,equilibrium_time,smooth_factor,k_cond,k_evap, &
                                 upwind_molec_diff,jacobian_method
  call get_command_argument(1,inputs);call get_command_argument(2,data)
  call get_command_argument(3,mode);call get_command_argument(4,output)
  call get_command_argument(5,arg);read(arg,*)max_steps
  call get_command_argument(6,arg);read(arg,*)limit
  pc=EvoAtmosphere(trim(inputs)//'/mechanism.yaml',trim(inputs)//'/settings.yaml', &
                   trim(inputs)//'/flux.txt',trim(data),err);call check()
  call state%initialize_from_h5(trim(inputs)//'/initial.h5',err);call check()
  call pc%initialize_from_state(state,err);call check()
  ! Solver settings and physical controls are not part of an atmosphere state.
  ! Reapply the case-owned controls identically for every backend.
  open(newunit=unit,file=trim(inputs)//'/controls.nml',status='old')
  read(unit,nml=benchmark_controls);close(unit)
  pc%var%atol=atol;pc%var%rtol=rtol
  pc%var%conv_longdy=0;pc%var%equilibrium_time=equilibrium_time
  pc%var%verbose=0;pc%var%nsteps_before_giveup=max_steps
  pc%var%cond_params(:)%smooth_factor=smooth_factor
  pc%var%cond_params(:)%k_cond=k_cond
  pc%var%cond_params(:)%k_evap=k_evap
  pc%var%upwind_molec_diff=upwind_molec_diff
  pc%var%jacobian_method=jacobian_method
  if(index(mode,'cvode1')==1)pc%var%max_order=1
  if(trim(mode)=='cvode1-tight')then
    pc%var%rtol=1e-6_dp;pc%var%atol=1e-23_dp
  endif
  is_ptc=index(mode,'ptc')==1
  if(trim(mode)=='ptc-periodic')then
    pc%var%press_temp_edd_profile%mode=PeriodicPressTempEdd
    ! Force a real profile mismatch after initialization, not a synthetic call.
  endif
  if(trim(mode)=='ptc-toa')then
    pc%var%toa_pressure_maintenance%target_pressure=0.8_dp*pc%wrk%pressure(pc%var%nz)
    pc%var%toa_pressure_maintenance%pressure_factor=1.02_dp
    pc%var%toa_pressure_maintenance%extreme_pressure_factor=1.1_dp
  endif
  call cpu_time(start)
  if(is_ptc)then
    if(index(mode,'floor12')>0)ptc%abundance_floor=1e-12_dp
    call ptc%initialize(pc,err);call check()
    ptc%solver%residual_policy=trim(mode)=='ptc-residual'
    if(trim(mode)=='ptc-periodic') &
      pc%var%press_temp_edd_profile%temperature=1.1_dp*pc%var%press_temp_edd_profile%temperature
  else
    call pc%initialize_robust_stepper(pc%wrk%usol,err);call check()
  endif
  open(newunit=unit,file=trim(output)//'.csv',status='replace')
  write(unit,'(A)')'step,time,dt,change,residual,limiting_species,layer,rejects,maintenance'
  finished=.false.;give_up=.false.;steps=0;time=0;dt=0
  do n=1,max_steps
    if(is_ptc)then
      call ptc%step(pc,finished,err)
      if(allocated(err))exit
      steps=ptc%solver%steps;time=ptc%solver%time;dt=ptc%solver%last_dt
      i=mod(ptc%solver%limiting_index-1,pc%dat%nq)+1
      j=(ptc%solver%limiting_index-1)/pc%dat%nq+1
      write(unit,'(I0,4(A,ES23.15),A,A,A,I0,A,I0,A,I0)') &
        steps,',',time,',',dt,',',ptc%solver%change,',',ptc%solver%residual, &
        ',',trim(pc%dat%species_names(i)),',',j,',',ptc%solver%rejects,',',ptc%maintenance_updates
    else
      time=pc%wrk%tn
      call pc%robust_step(give_up,finished,err)
      if(allocated(err))exit
      steps=pc%wrk%nsteps_total;dt=pc%wrk%tn-time;time=pc%wrk%tn
      write(unit,'(I0,A,ES23.15,A,ES23.15)')steps,',',time,',',dt
    endif
    call cpu_time(now)
    if(finished)then
      reason='time_threshold';exit
    elseif(give_up)then
      reason='give_up';exit
    elseif(now-start>limit)then
      reason='cpu_limit';exit
    endif
  enddo
  close(unit)
  call cpu_time(now)
  if(allocated(err))then
    reason='failure: '//err
    print *, 'FAILURE: ',err
    deallocate(err)
  endif
  if(is_ptc)then
    rejected=ptc%solver%rejects;nrhs=ptc%solver%rhs_calls;njac=ptc%solver%jac_calls
    maint=ptc%maintenance_updates;nsolves=ptc%solver%solves
    clips=ptc%solver%clipped_components;max_clip=ptc%solver%maximum_clipped_scaled
  else
    status=FCVodeGetNumRhsEvals(pc%wrk%sun%cvode_mem,crhs)
    status=FCVodeGetNumJacEvals(pc%wrk%sun%cvode_mem,cjac)
    status=FCVodeGetNumErrTestFails(pc%wrk%sun%cvode_mem,crej)
    rejected=int(crej(1));nrhs=int(crhs(1));njac=int(cjac(1));maint=0
    nsolves=-1;clips=0;max_clip=0
    call pc%destroy_stepper(err);call check()
  endif
  ! Reconcile column-integrated RHS against process-resolved budgets locally.
  allocate(f(pc%var%neqs));y=pc%wrk%usol
  call pc%right_hand_side(size(f),0.0_dp,reshape(y,[size(f)]),f,err);call check()
  budget_error=0
  open(newunit=budget_unit,file=trim(output)//'.budgets.csv',status='replace')
  write(budget_unit,'(A)')'species,rhs_column,diagnostic_net,gross_turnover,normalized_difference'
  do i=1,pc%dat%nq
    call pc%production_and_loss(trim(pc%dat%species_names(i)),y,budget,err);call check()
    rhs_column=sum(f(i::pc%dat%nq)*pc%var%dz)
    normalizer=max(sum(abs(budget%integrated_production))+sum(abs(budget%integrated_loss)),1e-100_dp)
    relative_error=abs(rhs_column-budget%integrated_net)/normalizer
    budget_error=max(budget_error,relative_error)
    write(budget_unit,'(A,4(A,ES23.15))')trim(pc%dat%species_names(i)),',',rhs_column, &
      ',',budget%integrated_net,',',normalizer,',',relative_error
  enddo
  close(budget_unit)
  allocate(atom_change(pc%dat%natoms));atom_change=0
  do i=1,pc%dat%nq
    atom_change=atom_change+pc%dat%species_composition(:,i)*sum(f(i::pc%dat%nq)*pc%var%dz)
  enddo
  open(newunit=unit,file=trim(output)//'.summary',status='replace')
  write(unit,*)'jacobian_method ',pc%var%jacobian_method
  write(unit,*)'reason ',trim(reason)
  write(unit,*)'finished ',finished
  write(unit,*)'steps ',steps
  write(unit,*)'time ',time
  write(unit,*)'cpu_seconds ',now-start
  write(unit,*)'rejects ',rejected
  write(unit,*)'rhs_calls ',nrhs
  write(unit,*)'jacobian_calls ',njac
  write(unit,*)'maintenance ',maint
  write(unit,*)'linear_solves ',nsolves
  write(unit,*)'clipped_components ',clips
  write(unit,*)'maximum_clipped_scaled ',max_clip
  write(unit,*)'budget_reconciliation_error ',budget_error
  write(unit,*)'maximum_element_column_tendency ',maxval(abs(atom_change))
  write(unit,*)'minimum_density ',minval(pc%wrk%usol)
  write(unit,*)'maximum_mix_tendency ',maxval(abs(reshape(f,shape(y)))/spread(pc%wrk%density,1,pc%dat%nq))
  close(unit)
  state=pc%model_state(err);call check()
  call state%save_to_h5(trim(output)//'.h5',err);call check()
  print *,trim(mode),' finished=',finished,' steps=',steps,' time=',time,' CPU seconds=',now-start
contains
  subroutine check()
    if(allocated(err))then
      print *,err
      stop 1
    endif
  end subroutine
end program
