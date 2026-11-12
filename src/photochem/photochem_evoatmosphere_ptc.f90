module photochem_evoatmosphere_ptc
  ! Experimental Fortran-only adapter. No Python API or CVODE default changes.
  use photochem_const, only: dp
  use photochem_ptc, only: PTCSolver
  use photochem_evoatmosphere, only: EvoAtmosphere, maintain_ptc_atmosphere
  use iso_c_binding, only: c_associated
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: EvoAtmospherePTC
  type :: EvoAtmospherePTC
    type(PTCSolver) :: solver
    integer :: maintenance_updates=0
    real(dp) :: abundance_floor=1e-20_dp
  contains
    procedure :: initialize
    procedure :: step
  end type
  ! Context avoids storing a model pointer beyond the duration of a call.
  type :: AtmosphereContext
    type(EvoAtmosphere), pointer :: pc
    real(dp) :: floor
  end type
contains
  subroutine initialize(self, pc, err)
    class(EvoAtmospherePTC), intent(inout) :: self
    type(EvoAtmosphere), target, intent(inout) :: pc
    character(:), allocatable, intent(out) :: err
    integer :: i,j
    logical :: finished, updated
    self%solver%initialized=.false.
    if (.not.pc%atmosphere_initialized) then
      err='PTC requires an initialized atmosphere'; return
    endif
    if (c_associated(pc%wrk%sun%cvode_mem)) then
      err='Destroy the CVODE stepper before initializing PTC'; return
    endif
    if (associated(pc%var%photon_flux_fcn)) then
      err='PTC prototype does not support time-dependent forcing'; return
    endif
    do i=1,size(pc%var%rate_fcns)
      if (associated(pc%var%rate_fcns(i)%fcn)) then
        err='PTC prototype does not support rate callbacks'; return
      endif
    enddo
    if (.not.ieee_is_finite(self%abundance_floor).or.self%abundance_floor<=0) then
      err='PTC abundance floor must be positive'; return
    endif
    finished=.false.
    call maintain_ptc_atmosphere(pc,finished,updated,initializing=.true.,err=err)
    if (allocated(err)) return
    do j=1,pc%var%nz
      if (any(pc%wrk%usol(:,j)<-1e-4_dp*self%abundance_floor*pc%wrk%density(j))) then
        err='PTC initial atmosphere has significant negative densities'; return
      endif
    enddo
    call self%solver%initialize(reshape(max(pc%wrk%usol,0.0_dp),[pc%var%neqs]),pc%dat%nq,pc%dat%nq, &
                               pc%var%initial_dt,err)
    if (allocated(err)) return
    self%solver%nonnegative=.true.
    self%maintenance_updates=0
    pc%wrk%tn=0; pc%wrk%t_history(1)=0
  end subroutine

  subroutine step(self, pc, finished, err)
    class(EvoAtmospherePTC), intent(inout) :: self
    type(EvoAtmosphere), target, intent(inout) :: pc
    logical, intent(out) :: finished
    character(:), allocatable, intent(out) :: err
    type(AtmosphereContext) :: context
    real(dp) :: time, dt, last_dt
    logical :: updated
    if (.not.self%solver%initialized.or..not.pc%atmosphere_initialized) then
      finished=.false.; err='PTC atmosphere is not initialized'; return
    endif
    if (size(self%solver%y)/=pc%var%neqs.or.c_associated(pc%wrk%sun%cvode_mem)) then
      finished=.false.; err='PTC model changed or another stepper is active; reinitialize PTC'; return
    endif
    context%pc=>pc; context%floor=self%abundance_floor
    finished=.false.
    call self%solver%step(context,atmosphere_rhs,atmosphere_jac,err)
    pc%wrk%tn=self%solver%time
    pc%wrk%t_history(1)=self%solver%time
    if (allocated(err)) return
    finished=self%solver%time>pc%var%equilibrium_time
    call maintain_ptc_atmosphere(pc,finished,updated,err=err)
    if (allocated(err)) then
      self%solver%initialized=.false.
      return
    endif
    if (updated) then
      ! Regridding can change all work allocations. Reinitialize numerical
      ! state while retaining cumulative cost counters and pseudo-time.
      time=self%solver%time; dt=self%solver%dt; last_dt=self%solver%last_dt
      block
        integer :: steps,rejects,nrhs,njac,nsolves,nclipped
        real(dp) :: max_clip
        steps=self%solver%steps; rejects=self%solver%rejects
        nclipped=self%solver%clipped_components; max_clip=self%solver%maximum_clipped_scaled
        nrhs=self%solver%rhs_calls; njac=self%solver%jac_calls; nsolves=self%solver%solves
        call self%solver%initialize(reshape(pc%wrk%usol,[pc%var%neqs]),pc%dat%nq,pc%dat%nq,dt,err)
        if(allocated(err))return
        self%solver%steps=steps; self%solver%rejects=rejects
        self%solver%clipped_components=nclipped; self%solver%maximum_clipped_scaled=max_clip
        self%solver%rhs_calls=nrhs; self%solver%jac_calls=njac; self%solver%solves=nsolves
      end block
      self%solver%time=time; self%solver%last_dt=last_dt
      self%maintenance_updates=self%maintenance_updates+1
    endif
  end subroutine

  subroutine atmosphere_rhs(context, y, f, scale, err)
    class(*), intent(inout) :: context
    real(dp), target, intent(in) :: y(:)
    real(dp), intent(out) :: f(:), scale(:)
    character(:), allocatable, intent(out) :: err
    integer :: j, first,last
    select type(context)
    type is(AtmosphereContext)
      ! Autonomous prototype: time-dependent callbacks are rejected at setup.
      call context%pc%right_hand_side(size(y),0.0_dp,y,f,err)
      if(allocated(err))return
      do j=1,context%pc%var%nz
        first=(j-1)*context%pc%dat%nq+1; last=j*context%pc%dat%nq
        scale(first:last)=max(abs(y(first:last)),context%floor*context%pc%wrk%density(j))
      enddo
    class default
      err='Invalid PTC atmosphere context'
    end select
  end subroutine

  subroutine atmosphere_jac(context,y,jac,err)
    class(*), intent(inout) :: context
    real(dp), target, intent(in) :: y(:)
    real(dp), target, intent(out) :: jac(:,:)
    character(:), allocatable, intent(out) :: err
    real(dp), allocatable :: flat(:)
    select type(context)
    type is(AtmosphereContext)
      ! Photochem and LAPACK both place the diagonal at 2*nq+1.
      allocate(flat(size(jac)))
      call context%pc%jacobian(size(flat),size(y),y,flat,err)
      if(allocated(err))return
      jac=reshape(flat,shape(jac))
    class default
      err='Invalid PTC atmosphere context'
    end select
  end subroutine
end module
