program test_ptc
  use photochem_const, only: dp
  use photochem_ptc, only: PTCSolver
  use photochem, only: EvoAtmosphere
  use photochem_evoatmosphere_ptc, only: EvoAtmospherePTC
  use photochem_test_paths, only: test_file, data_file, data_dir
  implicit none
  type :: TestContext
    real(dp) :: last(3)=0
    integer :: mode=1
  end type
  type(TestContext) :: context
  type(PTCSolver) :: solver
  type(EvoAtmosphere), target :: pc
  type(EvoAtmospherePTC) :: atmosphere
  character(:), allocatable :: err
  real(dp) :: eigen(3),expected(3),lambda,time
  logical :: finished
  integer :: i
  eigen=[1.0_dp,sqrt(2.0_dp),1.0_dp];lambda=-2+sqrt(2.0_dp)
  call solver%initialize(3+eigen,1,1,0.1_dp,err);call check()
  call solver%step(context,rhs,jac,err);call check()
  expected=3+eigen/(1-0.1_dp*lambda)
  if(maxval(abs(solver%y-expected))>1e-13_dp)stop 1
  ! A large implicit step damps the stiff oscillation about a positive state.
  context%mode=2
  call solver%initialize([1001.0_dp,1000.0_dp],1,1,10.0_dp,err);call check()
  call solver%step(context,rhs,jac,err);call check()
  if(sqrt(sum((solver%y-1000)**2))>0.002_dp)stop 2
  ! Nonlinear relaxation exercises rejection, then converges without tuning.
  context%mode=3
  call solver%initialize([2.0_dp],0,0,1.0_dp,err);call check()
  do i=1,100
    call solver%step(context,rhs,jac,err);call check()
  enddo
  if(solver%rejects<1.or.abs(solver%y(1)-1)>1e-12_dp)stop 3
  ! Exhausting retries must restore the callback context and accepted iterate.
  call solver%initialize([2.0_dp],0,0,1.0_dp,err);call check()
  solver%max_retries=1
  call solver%step(context,rhs,jac,err)
  if(.not.allocated(err))stop 4
  deallocate(err)
  if(solver%steps/=0.or.solver%time/=0.or.context%last(1)/=2.or.solver%y(1)/=2)stop 5
  ! A failed trial callback must also restore the callback context.
  context%mode=5
  call solver%initialize([2.0_dp],0,0,1.0_dp,err);call check()
  solver%change_limit=0.9_dp
  call solver%step(context,rhs,jac,err)
  if(.not.allocated(err))stop 11
  deallocate(err)
  if(context%last(1)/=2.or.solver%steps/=0)stop 12
  ! Only explicitly bounded trace negatives may be projected to zero.
  context%mode=4;solver%nonnegative=.true.;solver%max_retries=100
  call solver%initialize([0.0_dp],0,0,1e-6_dp,err);call check()
  call solver%step(context,rhs,jac,err);call check()
  if(solver%y(1)/=0.or.solver%clipped_components/=1.or.solver%maximum_clipped_scaled>1e-4_dp)stop 13
  ! Bounded real-atmosphere smoke; numerical storage may be reinitialized.
  pc=EvoAtmosphere(data_file('reaction_mechanisms/zahnle_earth.yaml'), &
       test_file('settings.yaml'),test_file('sun.txt'),test_file('atmosphere.txt'),data_dir,err)
  call check()
  call atmosphere%initialize(pc,err);call check()
  call atmosphere%step(pc,finished,err);call check()
  if(atmosphere%solver%steps/=1.or.pc%wrk%tn<=0)stop 6
  time=pc%wrk%tn
  ! Force shared TOA maintenance through an actual atmospheric grid update.
  pc%var%toa_pressure_maintenance%enabled=.true.
  pc%var%toa_pressure_maintenance%target_pressure=0.8_dp*pc%wrk%pressure(pc%var%nz)
  pc%var%toa_pressure_maintenance%pressure_factor=1.02_dp
  pc%var%toa_pressure_maintenance%extreme_pressure_factor=1.1_dp
  ! TOA transactions require a persistent pressure profile; use the current one.
  call pc%set_press_temp_edd_profile(pc%wrk%pressure,pc%var%temperature,pc%var%edd,trop_p=pc%wrk%pressure_hydro(pc%var%trop_ind), &
       target_pressure=pc%var%toa_pressure_maintenance%target_pressure,err=err)
  call check()
  pc%var%toa_pressure_maintenance%pressure_factor=1.02_dp
  pc%var%toa_pressure_maintenance%extreme_pressure_factor=1.1_dp
  call atmosphere%step(pc,finished,err);call check()
  if(atmosphere%maintenance_updates/=1.or.pc%wrk%tn<=time)stop 7
  print *, 'PTC numerical and atmospheric smoke checks passed'
contains
  subroutine check()
    if(allocated(err))then
      print *,err
      stop 10
    endif
  end subroutine
  subroutine rhs(raw,y,f,scale,err)
    class(*),intent(inout)::raw
    real(dp),target,intent(in)::y(:)
    real(dp),intent(out)::f(:),scale(:)
    character(:),allocatable,intent(out)::err
    real(dp)::u(size(y))
    select type(raw)
    type is(TestContext)
      raw%last(1:size(y))=y
      select case(raw%mode)
      case(1)
        u=y-3;f=-2*u
        f(1:2)=f(1:2)+u(2:3);f(2:3)=f(2:3)+u(1:2)
      case(2)
        u=y-1000;f=[-u(1)-100*u(2),100*u(1)-u(2)]
      case(3,5)
        f=1000*(1-y*y)
        if(raw%mode==5.and.y(1)<1.99_dp)err='Synthetic recoverable trial failure'
      case(4)
        f=-1e-30_dp
      end select
    end select
    scale=max(abs(y),1e-20_dp)
  end subroutine
  subroutine jac(raw,y,a,err)
    class(*),intent(inout)::raw
    real(dp),target,intent(in)::y(:)
    real(dp),target,intent(out)::a(:,:)
    character(:),allocatable,intent(out)::err
    a=0
    select type(raw)
    type is(TestContext)
      select case(raw%mode)
      case(1)
        a(3,:)=-2;a(2,2:3)=1;a(4,1:2)=1
      case(2)
        a(3,:)=-1;a(2,2)=-100;a(4,1)=100
      case(3,5)
        a(1,:)=-2000*y
      end select
    end select
  end subroutine
end program
