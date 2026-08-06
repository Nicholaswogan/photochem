program benchmark_block_thomas
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long
  use, intrinsic :: iso_fortran_env, only: int64
  use photochem, only: EvoAtmosphere, dp
  implicit none

  real(dp), parameter :: default_final_time = 1.0e17_dp
  real(dp), parameter :: ch4_flux = 1.0e10_dp

  type :: BenchmarkResult
    character(12) :: solver
    real(dp) :: elapsed_seconds
    real(dp) :: final_time
    integer(c_long) :: nsteps
    integer(c_long) :: nrhs
    integer(c_long) :: nlinsetups
    integer(c_long) :: nerrfails
    integer(c_long) :: nniters
    integer(c_long) :: nncfails
    integer(c_long) :: njevals
    real(dp), allocatable :: state(:)
  end type

  type(BenchmarkResult) :: band, block
  character(:), allocatable :: err
  real(dp) :: final_time

  call parse_final_time(final_time, err)
  if (allocated(err)) then
    print *, trim(err)
    stop 1
  endif

  print '(a)', 'Block-Thomas CVODE benchmark'
  print '(a,es12.4)', 'Target integration time (s): ', final_time
  print '(a,es12.4)', 'Perturbed lower CH4 flux (molecules/cm2/s): ', ch4_flux

  call run_case(.false., final_time, band, err)
  if (allocated(err)) then
    print *, trim(err)
    stop 1
  endif
  call run_case(.true., final_time, block, err)
  if (allocated(err)) then
    print *, trim(err)
    stop 1
  endif

  call print_results(band, block)

contains

  subroutine parse_final_time(final_time, err)
    real(dp), intent(out) :: final_time
    character(:), allocatable, intent(out) :: err

    character(128) :: argument
    integer :: io

    final_time = default_final_time
    if (command_argument_count() > 1) then
      err = 'Usage: ./tests/benchmark_block_thomas [final_time_seconds]'
      return
    endif
    if (command_argument_count() == 1) then
      call get_command_argument(1, argument)
      read(argument, *, iostat=io) final_time
      if (io /= 0 .or. final_time <= 0.0_dp) then
        err = 'The optional final time must be a positive number.'
      endif
    endif

  end subroutine

  subroutine run_case(use_block_thomas, final_time, result, err)
    use fcvode_mod, only: CV_NORMAL, FCVode, FCVodeGetIntegratorStats, &
                          FCVodeGetNonlinSolvStats, FCVodeGetNumJacEvals
    logical, intent(in) :: use_block_thomas
    real(dp), intent(in) :: final_time
    type(BenchmarkResult), intent(out) :: result
    character(:), allocatable, intent(out) :: err

    type(EvoAtmosphere) :: pc
    integer(c_int) :: ierr
    integer(c_long) :: nsteps(1), nrhs(1), nlinsetups(1), nerrfails(1)
    integer(c_long) :: nniters(1), nncfails(1), njevals(1)
    integer(c_int) :: qlast(1), qcur(1)
    integer(int64) :: count_start, count_end, count_rate
    real(c_double) :: hinused(1), hlast(1), hcur(1), tcur(1)

    pc = EvoAtmosphere('../data/reaction_mechanisms/zahnle_earth.yaml', &
                       '../examples/ModernEarth/settings.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    if (allocated(err)) return

    pc%var%verbose = 0
    pc%var%use_block_thomas = use_block_thomas
    call pc%set_lower_bc('CH4', 'flux', flux=ch4_flux, err=err)
    if (allocated(err)) return
    call pc%initialize_stepper(pc%var%usol_init, err)
    if (allocated(err)) return

    call system_clock(count_start, count_rate=count_rate)
    ierr = FCVode(pc%wrk%sun%cvode_mem, final_time, pc%wrk%sun%sunvec_y, &
                  tcur, CV_NORMAL)
    call system_clock(count_end)
    if (ierr /= 0) then
      if (use_block_thomas) then
        err = 'CVODE integration failed with the block-Thomas solver.'
      else
        err = 'CVODE integration failed with the band solver.'
      endif
      return
    endif

    ierr = FCVodeGetIntegratorStats(pc%wrk%sun%cvode_mem, nsteps, nrhs, &
                                    nlinsetups, nerrfails, qlast, qcur, &
                                    hinused, hlast, hcur, tcur)
    if (ierr /= 0) then
      err = 'Unable to retrieve CVODE integrator statistics.'
      return
    endif
    ierr = FCVodeGetNonlinSolvStats(pc%wrk%sun%cvode_mem, nniters, nncfails)
    if (ierr /= 0) then
      err = 'Unable to retrieve CVODE nonlinear solver statistics.'
      return
    endif
    ierr = FCVodeGetNumJacEvals(pc%wrk%sun%cvode_mem, njevals)
    if (ierr /= 0) then
      err = 'Unable to retrieve the CVODE Jacobian evaluation count.'
      return
    endif

    if (use_block_thomas) then
      result%solver = 'block-thomas'
    else
      result%solver = 'band'
    endif
    result%elapsed_seconds = real(count_end-count_start,dp)/real(count_rate,dp)
    result%final_time = tcur(1)
    result%nsteps = nsteps(1)
    result%nrhs = nrhs(1)
    result%nlinsetups = nlinsetups(1)
    result%nerrfails = nerrfails(1)
    result%nniters = nniters(1)
    result%nncfails = nncfails(1)
    result%njevals = njevals(1)
    result%state = pc%wrk%sun%yvec

  end subroutine

  subroutine print_results(band, block)
    type(BenchmarkResult), intent(in) :: band, block

    real(dp) :: state_scale, significant_cutoff
    real(dp) :: normalized_state_difference, componentwise_difference
    logical, allocatable :: significant(:)

    print '(a)', ''
    print '(a)', 'Metric                         Band          Block-Thomas'
    print '(a)', '----------------------------------------------------------'
    print '(a28,2(1x,es14.6))', 'Wall time (s)', &
          band%elapsed_seconds, block%elapsed_seconds
    print '(a28,2(1x,i14))', 'Internal steps', band%nsteps, block%nsteps
    print '(a28,2(1x,i14))', 'RHS evaluations', band%nrhs, block%nrhs
    print '(a28,2(1x,i14))', 'Linear solver setups', &
          band%nlinsetups, block%nlinsetups
    print '(a28,2(1x,i14))', 'Jacobian evaluations', band%njevals, block%njevals
    print '(a28,2(1x,i14))', 'Nonlinear iterations', band%nniters, block%nniters
    print '(a28,2(1x,i14))', 'Nonlinear convergence fails', &
          band%nncfails, block%nncfails
    print '(a28,2(1x,i14))', 'Error-test failures', &
          band%nerrfails, block%nerrfails
    print '(a28,2(1x,es14.6))', 'Reached time (s)', &
          band%final_time, block%final_time

    print '(a)', ''
    print '(a,f10.4)', 'Band/block wall-time speedup: ', &
          band%elapsed_seconds/block%elapsed_seconds

    state_scale = max(1.0_dp, maxval(abs(band%state)), maxval(abs(block%state)))
    normalized_state_difference = maxval(abs(block%state-band%state))/state_scale
    significant_cutoff = 1.0e-20_dp*state_scale
    significant = max(abs(band%state),abs(block%state)) > significant_cutoff
    if (any(significant)) then
      componentwise_difference = maxval(abs(block%state-band%state) / &
          max(abs(band%state),abs(block%state),tiny(1.0_dp)), &
          mask=significant)
    else
      componentwise_difference = 0.0_dp
    endif
    print '(a,es12.4)', 'Normalized final-state max difference: ', &
          normalized_state_difference
    print '(a,es12.4)', 'Max relative difference for significant entries: ', &
          componentwise_difference

  end subroutine

end program
