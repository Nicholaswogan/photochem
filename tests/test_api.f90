program test_api
  ! This executable contains focused correctness and error-behavior checks for
  ! public EvoAtmosphere operations. See tests/README.md.
  use photochem, only: EvoAtmosphere, dp
  implicit none

  call test()
  print *, 'test_api passed'

contains

  subroutine test()
    call test_initialization_state()
    call test_initialize_atmosphere_z()
    call test_initialize_atmosphere_p()
    call test_legacy_file_grid()
    call test_robust_stepper_initialization()
    call test_robust_stepper_restarts()
    call test_robust_stepper_limits()
    call test_set_press_temp_edd_nonmonotonic()
    call test_update_vertical_grid_inputs()
    call test_update_vertical_grid_pressure()
    call test_update_vertical_grid_repeated()
    call test_update_vertical_grid_particles()
    call test_update_vertical_grid_atomicity()
    call test_methods('../data/reaction_mechanisms/zahnle_earth.yaml')
    call test_methods('../tests/no_particle_test.yaml')
  end subroutine

  subroutine test_update_vertical_grid_inputs()
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use photochem_const, only: k_boltz
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: z_before(:), usol_before(:,:), gas_mix_top(:)
    real(dp) :: top_before, invalid_alt, z_top, density_top, pressure_top, temperature_top

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    top_before = pc%var%top_atmos
    z_before = pc%var%z
    usol_before = pc%wrk%usol
    z_top = pc%var%z(pc%var%nz)
    density_top = sum(pc%wrk%usol(pc%dat%ng_1:,pc%var%nz))
    pressure_top = density_top*k_boltz*pc%var%temperature(pc%var%nz)
    temperature_top = pc%var%temperature(pc%var%nz)
    gas_mix_top = pc%wrk%usol(pc%dat%ng_1:,pc%var%nz)/density_top

    call pc%update_vertical_grid(err=err)
    if (.not.allocated(err)) then
      print *, 'A vertical-grid update without a target was accepted'
      stop 1
    endif
    deallocate(err)
    call pc%update_vertical_grid(TOA_alt=top_before, TOA_pressure=1.0_dp, err=err)
    if (.not.allocated(err)) then
      print *, 'A vertical-grid update with two targets was accepted'
      stop 1
    endif
    deallocate(err)

    invalid_alt = ieee_value(0.0_dp, ieee_quiet_nan)
    call pc%update_vertical_grid(TOA_alt=invalid_alt, err=err)
    if (.not.allocated(err)) then
      print *, 'A nonfinite TOA altitude was accepted'
      stop 1
    endif
    deallocate(err)
    call pc%update_vertical_grid(TOA_alt=pc%var%bottom_atmos, err=err)
    if (.not.allocated(err)) then
      print *, 'A TOA altitude at the model bottom was accepted'
      stop 1
    endif
    deallocate(err)
    call pc%update_vertical_grid(TOA_pressure=0.0_dp, err=err)
    if (.not.allocated(err)) then
      print *, 'A zero TOA pressure was accepted'
      stop 1
    endif
    deallocate(err)
    call pc%update_vertical_grid(TOA_pressure=invalid_alt, err=err)
    if (.not.allocated(err)) then
      print *, 'A nonfinite TOA pressure was accepted'
      stop 1
    endif
    deallocate(err)

    if (pc%var%top_atmos /= top_before .or. any(pc%var%z /= z_before) .or. &
        any(pc%wrk%usol /= usol_before)) then
      print *, 'A rejected vertical-grid target changed model state'
      stop 1
    endif

    ! Exercise both directions through the candidate-remapping path and check
    ! the hydrostatic policy used above the old model domain.
    call pc%update_vertical_grid(TOA_alt=1.20_dp*top_before, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (pc%var%top_atmos /= 1.20_dp*top_before .or. &
        .not.all(ieee_is_finite(pc%var%z)) .or. any(pc%var%dz <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(pc%wrk%usol)) .or. &
        any(sum(pc%wrk%usol(pc%dat%ng_1:,:),dim=1) <= 0.0_dp)) then
      print *, 'Raising the TOA produced an invalid candidate grid'
      stop 1
    endif
    call check_hydrostatic_extension(pc, z_top, pressure_top, temperature_top, gas_mix_top)
    call pc%update_vertical_grid(TOA_alt=0.98_dp*top_before, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (pc%var%top_atmos /= 0.98_dp*top_before .or. &
        .not.all(ieee_is_finite(pc%var%z)) .or. any(pc%var%dz <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(pc%wrk%usol))) then
      print *, 'Lowering the TOA produced an invalid candidate grid'
      stop 1
    endif

  end subroutine

  subroutine test_update_vertical_grid_repeated()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp) :: top_initial, target_pressure, target_altitude
    integer :: i

    pc = make_pressure_test_model(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    top_initial = pc%var%top_atmos
    target_pressure = pc%wrk%pressure(pc%var%nz)

    ! Alternate explicit-altitude remaps with pressure-targeted remaps. This
    ! catches cumulative ownership, allocation, and remapping failures without
    ! requiring a long integration in routine tests.
    do i = 1,3
      if (mod(i,2) == 1) then
        target_altitude = 1.08_dp*top_initial
      else
        target_altitude = 0.94_dp*top_initial
      endif
      call pc%update_vertical_grid(TOA_alt=target_altitude, err=err)
      if (allocated(err)) then
        print *, trim(err)
        stop 1
      endif
      call check_repeated_grid_state(pc)

      call pc%update_vertical_grid(TOA_pressure=target_pressure, err=err)
      if (allocated(err)) then
        print *, trim(err)
        stop 1
      endif
      if (abs(log10(pc%wrk%pressure(pc%var%nz)/target_pressure)) > 2.0e-8_dp) then
        print *, 'Repeated regrid missed its TOA pressure target'
        stop 1
      endif
      call check_repeated_grid_state(pc)
    enddo

  end subroutine

  subroutine check_repeated_grid_state(model)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(EvoAtmosphere), intent(in) :: model

    if (.not.all(ieee_is_finite(model%var%z)) .or. &
        .not.all(ieee_is_finite(model%wrk%usol)) .or. &
        .not.all(ieee_is_finite(model%wrk%pressure)) .or. &
        any(model%var%dz <= 0.0_dp) .or. &
        any(model%wrk%pressure <= 0.0_dp) .or. &
        any(sum(model%wrk%usol(model%dat%ng_1:,:),dim=1) <= 0.0_dp)) then
      print *, 'Repeated vertical regridding produced an invalid atmosphere'
      stop 1
    endif

  end subroutine

  subroutine test_update_vertical_grid_pressure()
    type(EvoAtmosphere) :: reference_high, reference_low
    type(EvoAtmosphere) :: pc_endpoint, pc_high, pc_low, pc_failure
    character(:), allocatable :: err
    real(dp), allocatable :: z_before(:), usol_before(:,:)
    real(dp) :: top_initial, top_high, top_low
    real(dp) :: pressure_initial, pressure_high, pressure_low

    reference_high = make_pressure_test_model(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    top_initial = reference_high%var%top_atmos
    top_high = 1.40_dp*top_initial
    call reference_high%update_vertical_grid(TOA_alt=top_high, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    pressure_high = reference_high%wrk%pressure(reference_high%var%nz)

    reference_low = make_pressure_test_model(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    top_low = 0.80_dp*top_initial
    call reference_low%update_vertical_grid(TOA_alt=top_low, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    pressure_low = reference_low%wrk%pressure(reference_low%var%nz)
    if (pressure_high >= pressure_low) then
      print *, 'Reference TOA pressure was not monotonic with model-top altitude'
      stop 1
    endif

    pc_endpoint = make_pressure_test_model(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    pressure_initial = pc_endpoint%wrk%pressure(pc_endpoint%var%nz)
    call pc_endpoint%update_vertical_grid(TOA_pressure=pressure_initial, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call check_pressure_target(pc_endpoint, top_initial, pressure_initial, 'endpoint')

    pc_high = make_pressure_test_model(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call pc_high%update_vertical_grid(TOA_pressure=pressure_high, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call check_pressure_target(pc_high, top_high, pressure_high, 'raised top')

    pc_low = make_pressure_test_model(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call pc_low%update_vertical_grid(TOA_pressure=pressure_low, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call check_pressure_target(pc_low, top_low, pressure_low, 'lowered top')

    pc_failure = make_pressure_test_model(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    z_before = pc_failure%var%z
    usol_before = pc_failure%wrk%usol
    call pc_failure%update_vertical_grid(TOA_pressure=1.0e100_dp, err=err)
    if (.not.allocated(err)) then
      print *, 'An unreachable TOA pressure was accepted'
      stop 1
    endif
    if (index(err,'maximum reachable pressure') == 0 .or. &
        any(pc_failure%var%z /= z_before) .or. &
        any(pc_failure%wrk%usol /= usol_before)) then
      print *, 'Failed TOA-pressure bracketing was unclear or changed model state'
      stop 1
    endif
    deallocate(err)

  end subroutine

  function make_pressure_test_model(err) result(pc)
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphere) :: pc

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
  end function

  subroutine check_pressure_target(pc, expected_top, target_pressure, label)
    type(EvoAtmosphere), intent(in) :: pc
    real(dp), intent(in) :: expected_top, target_pressure
    character(*), intent(in) :: label

    real(dp) :: pressure_result

    pressure_result = pc%wrk%pressure(pc%var%nz)
    if (abs(pc%var%top_atmos/expected_top-1.0_dp) > 2.0e-8_dp .or. &
        abs(log10(pressure_result/target_pressure)) > 2.0e-8_dp) then
      print *, 'Bracketed TOA-pressure solve missed the ',trim(label),' target'
      stop 1
    endif
  end subroutine

  subroutine test_update_vertical_grid_particles()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: particle_mix_top(:), particle_radius_top(:)
    real(dp) :: top_before, z_top, density_top, density_layer
    integer :: i, first_extended

    pc = EvoAtmosphere('../data/reaction_mechanisms/zahnle_earth.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    top_before = pc%var%top_atmos
    z_top = pc%var%z(pc%var%nz)
    density_top = sum(pc%wrk%usol(pc%dat%ng_1:,pc%var%nz))
    particle_mix_top = pc%wrk%usol(1:pc%dat%npq,pc%var%nz)/density_top
    particle_radius_top = pc%var%particle_radius(:,pc%var%nz)

    call pc%update_vertical_grid(TOA_alt=1.20_dp*top_before, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    first_extended = findloc(pc%var%z > z_top, .true., dim=1)
    if (first_extended == 0) then
      print *, 'Particle regrid test did not extend beyond the old top center'
      stop 1
    endif
    do i = first_extended,pc%var%nz
      density_layer = sum(pc%wrk%usol(pc%dat%ng_1:,i))
      if (maxval(abs(pc%wrk%usol(1:pc%dat%npq,i)/density_layer- &
                     particle_mix_top)) > 1.0e-12_dp .or. &
          maxval(abs(pc%var%particle_radius(:,i)/particle_radius_top-1.0_dp)) > &
          1.0e-12_dp) then
        print *, 'Particle abundance or radius was not constant above the old model top'
        stop 1
      endif
    enddo

  end subroutine

  subroutine test_update_vertical_grid_atomicity()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: z_before(:), usol_before(:,:), temperature_before(:)
    real(dp), allocatable :: particle_radius_before(:,:)
    real(dp) :: top_before, surface_pressure_before, tn
    integer :: i, optical_particle

    pc = EvoAtmosphere('../data/reaction_mechanisms/zahnle_earth.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    call pc%initialize_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    optical_particle = 0
    do i = 1,pc%dat%np
      if (pc%dat%part_xs_file(i)%ThereIsData) then
        optical_particle = i
        exit
      endif
    enddo
    if (optical_particle == 0) then
      print *, 'Atomic regrid test requires a particle with optical data'
      stop 1
    endif

    ! The interpolation contract rejects radii exactly on its tabulated lower
    ! boundary. Introduce that condition only after CVODE is initialized so
    ! candidate preparation fails after the original solver has been parked.
    pc%var%particle_radius(optical_particle,pc%var%nz) = &
        pc%dat%radii_file(1,optical_particle)
    top_before = pc%var%top_atmos
    surface_pressure_before = pc%var%surface_pressure
    z_before = pc%var%z
    usol_before = pc%wrk%usol
    temperature_before = pc%var%temperature
    particle_radius_before = pc%var%particle_radius

    call pc%update_vertical_grid(TOA_alt=1.20_dp*top_before, err=err)
    if (.not.allocated(err)) then
      print *, 'Invalid candidate optical properties were accepted during regridding'
      stop 1
    endif
    deallocate(err)
    if (pc%var%top_atmos /= top_before .or. &
        pc%var%surface_pressure /= surface_pressure_before .or. &
        any(pc%var%z /= z_before) .or. any(pc%wrk%usol /= usol_before) .or. &
        any(pc%var%temperature /= temperature_before) .or. &
        any(pc%var%particle_radius /= particle_radius_before)) then
      print *, 'Failed candidate preparation changed committed model state'
      stop 1
    endif

    tn = pc%step(err)
    if (allocated(err) .or. tn <= 0.0_dp) then
      if (allocated(err)) print *, trim(err)
      print *, 'Failed regridding did not preserve the active stepper'
      stop 1
    endif
    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

  end subroutine

  subroutine test_set_press_temp_edd_nonmonotonic()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp) :: P(2), T(2), edd(2)
    real(dp), allocatable :: temperature_before(:), edd_before(:)

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    ! Construct a deterministic actual-pressure inversion between the bottom
    ! two layers. Constant input temperature makes pressure proportional to
    ! gas density in actual-pressure mode.
    pc%wrk%usol(:,2) = 2.0_dp*pc%wrk%usol(:,1)
    P = [1.0e7_dp, 1.0_dp]
    T = [200.0_dp, 200.0_dp]
    edd = [1.0e7_dp, 1.0e7_dp]
    allocate(temperature_before(pc%var%nz), edd_before(pc%var%nz))
    temperature_before = pc%var%temperature
    edd_before = pc%var%edd

    call pc%set_press_temp_edd(P, T, edd, 1.0e5_dp, &
                               hydro_pressure=.false., err=err)
    if (.not.allocated(err)) then
      print *, 'A tropopause was mapped through nonmonotonic actual pressure'
      stop 1
    endif
    if (index(err, 'actual pressure is not strictly decreasing between layers 1 and 2') == 0) then
      print *, 'Nonmonotonic actual pressure returned an unclear error'
      stop 1
    endif
    if (any(pc%var%temperature /= temperature_before) .or. &
        any(pc%var%edd /= edd_before)) then
      print *, 'Failed pressure-profile mapping changed model state'
      stop 1
    endif
    deallocate(err)

  end subroutine

  subroutine test_robust_stepper_limits()
    use iso_c_binding, only: c_associated
    use fcvode_mod, only: FCVodeFree
    type(EvoAtmosphere) :: pc_errors, pc_steps
    character(:), allocatable :: err
    logical :: give_up, converged

    pc_errors = EvoAtmosphere('../tests/no_particle_test.yaml', &
                              '../tests/test_settings_minimal.yaml', &
                              '../examples/ModernEarth/Sun_now.txt', &
                              '../examples/ModernEarth/atmosphere.txt', &
                              '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    ! A limit of one permits one recovery restart and gives up on the next
    ! failed step. Forced missing CVODE memory makes both failures immediate.
    pc_errors%var%nerrors_before_giveup = 1
    call pc_errors%initialize_robust_stepper(pc_errors%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call FCVodeFree(pc_errors%wrk%sun%cvode_mem)
    call pc_errors%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (give_up .or. converged .or. pc_errors%wrk%nerrors_total /= 1 .or. &
        .not.pc_errors%wrk%robust_stepper_initialized .or. &
        .not.c_associated(pc_errors%wrk%sun%cvode_mem)) then
      print *, 'robust error limit did not permit exactly one recovery'
      stop 1
    endif
    call FCVodeFree(pc_errors%wrk%sun%cvode_mem)
    call pc_errors%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (.not.give_up .or. converged .or. &
        pc_errors%wrk%nerrors_total /= 2 .or. &
        pc_errors%wrk%robust_stepper_initialized .or. &
        c_associated(pc_errors%wrk%sun%cvode_mem)) then
      print *, 'robust error limit did not terminate on the next failure'
      stop 1
    endif

    pc_steps = EvoAtmosphere('../tests/no_particle_test.yaml', &
                             '../tests/test_settings_minimal.yaml', &
                             '../examples/ModernEarth/Sun_now.txt', &
                             '../examples/ModernEarth/atmosphere.txt', &
                             '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    ! The accepted-step ceiling is exact and is checked before an otherwise
    ! coincident scheduled restart.
    pc_steps%var%nsteps_before_conv_check = 0
    pc_steps%var%nsteps_before_reinit = 1
    pc_steps%var%nsteps_before_giveup = 1
    pc_steps%var%equilibrium_time = huge(1.0_dp)
    call pc_steps%initialize_robust_stepper(pc_steps%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call pc_steps%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (.not.give_up .or. converged .or. &
        pc_steps%wrk%nsteps_total /= 1 .or. pc_steps%wrk%nsteps /= 1) then
      print *, 'robust accepted-step limit was not exact'
      stop 1
    endif

    call pc_steps%destroy_stepper(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

  end subroutine

  subroutine test_robust_stepper_initialization()
    use iso_c_binding, only: c_associated
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: usol_wrong(:,:)
    logical :: give_up, converged
    integer :: max_order_original
    real(dp) :: reinit_min_density_original

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    ! Predictable validation failures must not disturb an existing ordinary
    ! stepper or claim that a robust session was initialized.
    call pc%initialize_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    pc%var%nerrors_before_giveup = 0
    call pc%initialize_robust_stepper(pc%wrk%usol, err)
    if (.not.allocated(err)) then
      print *, 'invalid robust error limit was accepted'
      stop 1
    endif
    if (.not.c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'robust-setting validation destroyed the existing stepper'
      stop 1
    endif
    if (pc%wrk%robust_stepper_initialized) then
      print *, 'failed robust validation initialized a robust session'
      stop 1
    endif
    deallocate(err)
    pc%var%nerrors_before_giveup = 10

    pc%var%nsteps_before_conv_check = pc%var%nsteps_before_reinit
    call pc%initialize_robust_stepper(pc%wrk%usol, err)
    if (.not.allocated(err)) then
      print *, 'invalid robust convergence interval was accepted'
      stop 1
    endif
    if (.not.c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'robust-interval validation destroyed the existing stepper'
      stop 1
    endif
    deallocate(err)
    pc%var%nsteps_before_conv_check = 300

    reinit_min_density_original = pc%var%reinit_min_density
    pc%var%reinit_min_density = ieee_value(0.0_dp, ieee_quiet_nan)
    call pc%initialize_robust_stepper(pc%wrk%usol, err)
    if (.not.allocated(err)) then
      print *, 'nonfinite robust clipping density was accepted'
      stop 1
    endif
    if (.not.c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'robust-density validation destroyed the existing stepper'
      stop 1
    endif
    deallocate(err)
    pc%var%reinit_min_density = reinit_min_density_original

    allocate(usol_wrong(pc%dat%nq-1,pc%var%nz))
    call pc%initialize_robust_stepper(usol_wrong, err)
    if (.not.allocated(err)) then
      print *, 'wrong-sized robust initial state was accepted'
      stop 1
    endif
    if (.not.c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'robust dimension validation destroyed the existing stepper'
      stop 1
    endif
    deallocate(err)

    ! A successful robust initialization commits state and counters together.
    call pc%initialize_robust_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (.not.pc%wrk%robust_stepper_initialized .or. &
        pc%wrk%nsteps_total /= 0 .or. pc%wrk%nerrors_total /= 0) then
      print *, 'robust initialization did not commit its state and counters'
      stop 1
    endif
    if (pc%wrk%nsteps /= 0 .or. pc%wrk%tn /= 0.0_dp .or. &
        pc%wrk%t_history(1) /= 0.0_dp) then
      print *, 'public robust initialization did not begin at time zero'
      stop 1
    endif

    ! A rejected ordinary initialization retains the robust session, while a
    ! successful ordinary initialization explicitly replaces it.
    call pc%initialize_stepper(usol_wrong, err)
    if (.not.allocated(err)) then
      print *, 'wrong-sized ordinary initial state was accepted'
      stop 1
    endif
    if (.not.pc%wrk%robust_stepper_initialized .or. &
        .not.c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'rejected ordinary initialization changed the robust session'
      stop 1
    endif
    deallocate(err)
    call pc%initialize_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (pc%wrk%robust_stepper_initialized) then
      print *, 'ordinary initialization retained stale robust-session state'
      stop 1
    endif
    call pc%robust_step(give_up, converged, err)
    if (.not.allocated(err)) then
      print *, 'robust_step accepted an ordinary stepper'
      stop 1
    endif
    deallocate(err)

    ! An unexpected CVODE setup failure must clean up all partial resources
    ! and must not reset robust counters as if initialization had succeeded.
    pc%wrk%nsteps_total = 7
    pc%wrk%nerrors_total = 8
    max_order_original = pc%var%max_order
    pc%var%max_order = 0
    call pc%initialize_robust_stepper(pc%wrk%usol, err)
    if (.not.allocated(err)) then
      print *, 'invalid CVODE maximum order did not fail setup'
      stop 1
    endif
    if (c_associated(pc%wrk%sun%cvode_mem) .or. &
        allocated(pc%wrk%sun%yvec) .or. allocated(pc%wrk%sun%abstol) .or. &
        associated(pc%wrk%sun%sunvec_y) .or. &
        associated(pc%wrk%sun%abstol_nvec) .or. &
        associated(pc%wrk%sun%sunmat) .or. associated(pc%wrk%sun%sunlin)) then
      print *, 'failed CVODE setup retained partial resources'
      stop 1
    endif
    if (pc%wrk%robust_stepper_initialized .or. &
        pc%wrk%nsteps_total /= 7 .or. pc%wrk%nerrors_total /= 8) then
      print *, 'failed CVODE setup committed robust-session state'
      stop 1
    endif
    deallocate(err)
    pc%var%max_order = max_order_original

    call pc%initialize_robust_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (pc%wrk%robust_stepper_initialized) then
      print *, 'destroy_stepper retained robust-session state'
      stop 1
    endif

  end subroutine

  subroutine test_robust_stepper_restarts()
    use iso_c_binding, only: c_associated, c_ptr
    use fcvode_mod, only: FCVodeFree
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    logical :: give_up, converged
    real(dp) :: restart_time
    type(c_ptr) :: cvode_mem_before

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    ! Force a scheduled restart after one accepted step. The CVODE allocation
    ! should be reused, while local CVODE/history counters reset at the same
    ! nonzero logical integration time.
    pc%var%nsteps_before_conv_check = 0
    pc%var%nsteps_before_reinit = 1
    call pc%initialize_robust_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    cvode_mem_before = pc%wrk%sun%cvode_mem
    call pc%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (give_up .or. converged) then
      print *, 'scheduled robust restart returned a terminal result'
      stop 1
    endif
    if (.not.c_associated(cvode_mem_before, pc%wrk%sun%cvode_mem)) then
      print *, 'scheduled robust restart rebuilt compatible CVODE memory'
      stop 1
    endif
    if (pc%wrk%nsteps /= 0 .or. pc%wrk%nsteps_total /= 1 .or. &
        pc%wrk%nerrors_total /= 0 .or. pc%wrk%t_history(1) <= 0.0_dp .or. &
        pc%wrk%tn /= pc%wrk%t_history(1)) then
      print *, 'scheduled robust restart did not preserve time and total counters'
      stop 1
    endif

    ! The first step after a nonzero-time ReInit must retain the forward
    ! direction. Keeping the one-step threshold also exercises a second ReInit.
    restart_time = pc%wrk%t_history(1)
    call pc%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (give_up .or. converged .or. pc%wrk%nerrors_total /= 0 .or. &
        pc%wrk%nsteps_total /= 2 .or. pc%wrk%t_history(1) <= restart_time) then
      print *, 'first step after scheduled robust restart did not advance'
      stop 1
    endif

    ! A missing CVODE object forces the same recovery path to reconstruct all
    ! infrastructure from the last committed state. Recovery itself must not
    ! run convergence checks, even when the committed time exceeds the limit.
    restart_time = pc%wrk%t_history(1)
    pc%var%nsteps_before_reinit = 1000
    pc%var%equilibrium_time = 0.0_dp
    call FCVodeFree(pc%wrk%sun%cvode_mem)
    call pc%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (give_up .or. converged) then
      print *, 'robust recovery used failed-step convergence state'
      stop 1
    endif
    if (.not.pc%wrk%robust_stepper_initialized .or. &
        .not.c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'robust recovery did not reconstruct missing CVODE memory'
      stop 1
    endif
    if (pc%wrk%nsteps /= 0 .or. pc%wrk%nsteps_total /= 2 .or. &
        pc%wrk%nerrors_total /= 1 .or. pc%wrk%tn /= restart_time .or. &
        pc%wrk%t_history(1) /= restart_time) then
      print *, 'robust recovery changed committed time or counters'
      stop 1
    endif
    ! Full reconstruction at nonzero time must permit the next forward step.
    pc%var%equilibrium_time = huge(1.0_dp)
    call pc%robust_step(give_up, converged, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (give_up .or. converged .or. pc%wrk%nerrors_total /= 1 .or. &
        pc%wrk%nsteps_total /= 3 .or. pc%wrk%t_history(1) <= restart_time) then
      print *, 'first step after robust fallback reconstruction did not advance'
      stop 1
    endif

    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

  end subroutine

  subroutine test_initialize_atmosphere_z()
    use iso_c_binding, only: c_associated
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    integer, parameter :: nprofile = 3
    real(dp), parameter :: surface_pressure = 1.0e6_dp
    real(dp), parameter :: fixed_H2_density = 2.5e15_dp
    real(dp) :: z(nprofile), temperature(nprofile), edd(nprofile)
    real(dp), allocatable :: mix(:,:), particle_radius(:,:)
    real(dp), allocatable :: temperature_before(:)
    real(dp) :: profile_pressure(2), profile_temperature(2), profile_edd(2)
    real(dp) :: fraction, expected_temperature, expected_log10edd
    integer :: i, ind_H2

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    allocate(mix(pc%dat%nq,nprofile))
    allocate(particle_radius(pc%dat%npq,nprofile))
    do i = 1,nprofile
      mix(:,i) = pc%wrk%mix(1:pc%dat%nq,1)
    enddo

    z = [0.0_dp, 5.0e6_dp, 1.0e7_dp]
    temperature = [300.0_dp, 240.0_dp, 180.0_dp]
    edd = [1.0e5_dp, 1.0e6_dp, 1.0e7_dp]

    ! A failed replacement must retain both atmospheric and CVODE state.
    profile_pressure = [2.0_dp*maxval(pc%wrk%pressure_hydro), &
                        0.5_dp*minval(pc%wrk%pressure_hydro)]
    profile_temperature = [290.0_dp, 180.0_dp]
    profile_edd = [1.0e5_dp, 1.0e7_dp]
    call pc%set_press_temp_edd_profile(profile_pressure, profile_temperature, &
                                       profile_edd, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call pc%initialize_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    temperature_before = pc%var%temperature
    z(1) = 1.0_dp
    call pc%initialize_atmosphere_z(z, temperature, edd, surface_pressure, &
                                    mix, particle_radius, err)
    if (.not. allocated(err)) then
      print *, 'invalid altitude initialization did not return an error'
      stop 1
    endif
    if (any(pc%var%temperature /= temperature_before)) then
      print *, 'failed altitude initialization changed atmospheric state'
      stop 1
    endif
    if (.not. c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'failed altitude initialization destroyed the active stepper'
      stop 1
    endif
    if (.not. pc%var%press_temp_edd_profile%enabled) then
      print *, 'failed altitude initialization changed the persistent profile'
      stop 1
    endif

    ! Fixed lower boundary conditions override the supplied bottom mixing
    ! ratio and are reflected in both prepared and canonical initial state.
    deallocate(err)
    call pc%set_lower_bc('H2', 'den', den=fixed_H2_density, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    ind_H2 = findloc(pc%dat%species_names(1:pc%dat%nq), 'H2', 1)
    z(1) = 0.0_dp
    call pc%initialize_atmosphere_z(z, temperature, edd, surface_pressure, &
                                    mix, particle_radius, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    if (.not. pc%atmosphere_initialized) then
      print *, 'altitude initialization did not set lifecycle state'
      stop 1
    endif
    if (c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'successful altitude initialization retained stale CVODE state'
      stop 1
    endif
    if (pc%var%top_atmos /= z(nprofile)) then
      print *, 'altitude initialization did not set the model top'
      stop 1
    endif
    if (pc%var%press_temp_edd_profile%enabled) then
      print *, 'altitude initialization retained a pressure-based profile'
      stop 1
    endif
    if (pc%wrk%usol(ind_H2,1) /= fixed_H2_density .or. &
        pc%var%usol_init(ind_H2,1) /= fixed_H2_density) then
      print *, 'altitude initialization did not apply the fixed-density boundary condition'
      stop 1
    endif
    if (any(pc%var%usol_init /= pc%wrk%usol)) then
      print *, 'canonical and prepared altitude-initialization states differ'
      stop 1
    endif
    if (any(pc%wrk%pressure_hydro <= 0.0_dp) .or. any(pc%wrk%density_hydro <= 0.0_dp)) then
      print *, 'altitude initialization produced invalid hydrostatic state'
      stop 1
    endif

    do i = 1,pc%var%nz
      fraction = pc%var%z(i)/z(nprofile)
      expected_temperature = temperature(1) + fraction*(temperature(nprofile)-temperature(1))
      expected_log10edd = log10(edd(1)) + fraction*(log10(edd(nprofile))-log10(edd(1)))
      if (abs(pc%var%temperature(i)-expected_temperature) > 1.0e-10_dp) then
        print *, 'altitude initialization mapped temperature incorrectly'
        stop 1
      endif
      if (abs(log10(pc%var%edd(i))-expected_log10edd) > 1.0e-10_dp) then
        print *, 'altitude initialization mapped eddy diffusion incorrectly'
        stop 1
      endif
    enddo

    call test_initialize_atmosphere_z_particles(z, temperature, edd, surface_pressure)
  end subroutine

  subroutine test_initialize_atmosphere_p()
    use iso_c_binding, only: c_associated
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    integer, parameter :: nprofile = 4
    real(dp) :: pressure(nprofile), pressure_bad(nprofile)
    real(dp) :: temperature(nprofile), edd(nprofile)
    real(dp), allocatable :: mix(:,:), particle_radius(:,:)
    real(dp), allocatable :: temperature_before(:)
    real(dp) :: top_before
    integer :: i

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    pressure = [1.0e6_dp, 1.0e5_dp, 1.0e4_dp, 1.0e2_dp]
    temperature = [300.0_dp, 260.0_dp, 220.0_dp, 180.0_dp]
    edd = [1.0e5_dp, 3.0e5_dp, 1.0e6_dp, 1.0e7_dp]
    allocate(mix(pc%dat%nq,nprofile))
    allocate(particle_radius(pc%dat%npq,nprofile))
    do i = 1,nprofile
      mix(:,i) = pc%wrk%mix(1:pc%dat%nq,1)
    enddo

    call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                    particle_radius, persistent=.true., err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (.not. pc%atmosphere_initialized .or. pc%var%top_atmos <= 0.0_dp .or. &
        pc%var%bottom_atmos /= 0.0_dp) then
      print *, 'pressure initialization produced an invalid model domain'
      stop 1
    endif
    if (any(pc%wrk%pressure_hydro(2:) >= &
            pc%wrk%pressure_hydro(:pc%var%nz-1))) then
      print *, 'pressure initialization did not produce decreasing hydrostatic pressure'
      stop 1
    endif
    if (pc%wrk%pressure_hydro(1) >= pressure(1) .or. &
        pc%wrk%pressure_hydro(pc%var%nz) <= pressure(nprofile)) then
      print *, 'pressure initialization did not treat input endpoints as domain boundaries'
      stop 1
    endif
    if (.not. pc%var%press_temp_edd_profile%enabled .or. &
        any(pc%var%press_temp_edd_profile%pressure /= pressure) .or. &
        any(pc%var%press_temp_edd_profile%temperature /= temperature) .or. &
        any(pc%var%press_temp_edd_profile%edd /= edd)) then
      print *, 'pressure initialization did not retain the persistent profile'
      stop 1
    endif
    if (any(pc%var%usol_init /= pc%wrk%usol)) then
      print *, 'canonical and prepared pressure-initialization states differ'
      stop 1
    endif

    ! A failed pressure reinitialization must retain atmospheric, profile, and
    ! CVODE state. Equal adjacent pressure points violate strict ordering.
    call pc%initialize_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    top_before = pc%var%top_atmos
    temperature_before = pc%var%temperature
    pressure_bad = pressure
    pressure_bad(2) = pressure_bad(1)
    call pc%initialize_atmosphere_p(pressure_bad, temperature, edd, mix, &
                                    particle_radius, err=err)
    if (.not. allocated(err)) then
      print *, 'invalid pressure initialization did not return an error'
      stop 1
    endif
    if (pc%var%top_atmos /= top_before .or. &
        any(pc%var%temperature /= temperature_before)) then
      print *, 'failed pressure initialization changed atmospheric state'
      stop 1
    endif
    if (.not. c_associated(pc%wrk%sun%cvode_mem) .or. &
        .not. pc%var%press_temp_edd_profile%enabled) then
      print *, 'failed pressure initialization changed integrator or profile state'
      stop 1
    endif

    ! Persistence is opt-in. A successful default initialization replaces the
    ! profile and invalidates CVODE state associated with the old atmosphere.
    deallocate(err)
    call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                    particle_radius, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (pc%var%press_temp_edd_profile%enabled .or. &
        c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'default pressure initialization retained profile or integrator state'
      stop 1
    endif

    call test_initialize_atmosphere_p_particles(pressure, temperature, edd)

  end subroutine

  subroutine test_initialize_atmosphere_p_particles(pressure, temperature, edd)
    real(dp), intent(in) :: pressure(:), temperature(:), edd(:)
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: mix(:,:), particle_radius(:,:)
    integer :: i, nprofile

    pc = EvoAtmosphere('../data/reaction_mechanisms/zahnle_earth.yaml', &
                       '../examples/ModernEarth/settings.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    nprofile = size(pressure)
    allocate(mix(pc%dat%nq,nprofile))
    allocate(particle_radius(pc%dat%npq,nprofile))
    do i = 1,nprofile
      mix(:,i) = pc%wrk%mix(1:pc%dat%nq,1)
      particle_radius(:,i) = pc%var%particle_radius(:,1)
    enddo

    call pc%initialize_atmosphere_p(pressure, temperature, edd, mix, &
                                    particle_radius, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (pc%var%press_temp_edd_profile%enabled .or. &
        any(pc%var%particle_radius <= 0.0_dp) .or. &
        any(pc%var%usol_init(1:pc%dat%npq,:) <= 0.0_dp)) then
      print *, 'particle-bearing pressure initialization produced invalid state'
      stop 1
    endif

  end subroutine

  subroutine test_initialize_atmosphere_z_particles(z, temperature, edd, surface_pressure)
    real(dp), intent(in) :: z(:), temperature(:), edd(:), surface_pressure
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: mix(:,:), particle_radius(:,:)
    integer :: i, nprofile

    pc = EvoAtmosphere('../data/reaction_mechanisms/zahnle_earth.yaml', &
                       '../examples/ModernEarth/settings.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    nprofile = size(z)
    allocate(mix(pc%dat%nq,nprofile))
    allocate(particle_radius(pc%dat%npq,nprofile))
    do i = 1,nprofile
      mix(:,i) = pc%wrk%mix(1:pc%dat%nq,1)
      particle_radius(:,i) = pc%var%particle_radius(:,1)
    enddo

    call pc%initialize_atmosphere_z(z, temperature, edd, surface_pressure, &
                                    mix, particle_radius, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (any(pc%var%particle_radius <= 0.0_dp) .or. &
        any(pc%var%usol_init(1:pc%dat%npq,:) <= 0.0_dp)) then
      print *, 'particle-bearing altitude initialization produced invalid particle state'
      stop 1
    endif
  end subroutine

  subroutine test_initialization_state()
    use iso_c_binding, only: c_associated
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: temperature_before(:)
    real(dp) :: profile_pressure(2), profile_temperature(2), profile_edd(2)

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (.not. pc%atmosphere_initialized) then
      print *, 'successful initialization did not set atmosphere_initialized'
      stop 1
    endif

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../tests/does_not_exist.txt', &
                       '../data', &
                       err)
    if (.not. allocated(err)) then
      print *, 'missing atmosphere file did not produce an initialization error'
      stop 1
    endif
    if (pc%atmosphere_initialized) then
      print *, 'failed initialization left atmosphere_initialized set'
      stop 1
    endif

    ! The static-only constructor succeeds without reading an atmosphere and
    ! leaves mechanism-sized state available for inspection and configuration.
    deallocate(err)
    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (pc%atmosphere_initialized) then
      print *, 'static-only construction initialized an atmosphere'
      stop 1
    endif
    if (.not. allocated(pc%dat) .or. .not. allocated(pc%var) .or. &
        .not. allocated(pc%wrk) .or. pc%dat%nq <= 0 .or. pc%var%nz <= 0) then
      print *, 'static-only construction did not complete static setup'
      stop 1
    endif

    ! Static configuration methods remain available before atmosphere
    ! initialization.
    call pc%set_lower_bc('H2', 'flux', flux=1.0e8_dp, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    ! Atmosphere-dependent methods consistently reject the same object before
    ! touching uninitialized atmospheric arrays.
    call pc%prep_atmosphere(pc%var%usol_init, err)
    call check_not_initialized_error(err, 'prep_atmosphere')
    deallocate(err)

    call pc%set_temperature(pc%var%temperature, err=err)
    call check_not_initialized_error(err, 'set_temperature')
    deallocate(err)

    call pc%initialize_stepper(pc%var%usol_init, err)
    call check_not_initialized_error(err, 'initialize_stepper')

    ! A statically configured object can be initialized later.
    deallocate(err)
    call pc%initialize_from_atmosphere_file('../examples/ModernEarth/atmosphere.txt', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (.not. pc%atmosphere_initialized) then
      print *, 'explicit atmosphere initialization did not set lifecycle state'
      stop 1
    endif
    if (any(pc%var%usol_init /= pc%wrk%usol)) then
      print *, 'legacy-file canonical and prepared initial states differ'
      stop 1
    endif

    ! Failed reinitialization is atomic: the initialized atmosphere and active
    ! CVODE state are retained.
    profile_pressure = [2.0_dp*maxval(pc%wrk%pressure_hydro), &
                        0.5_dp*minval(pc%wrk%pressure_hydro)]
    profile_temperature = [290.0_dp, 180.0_dp]
    profile_edd = [1.0e5_dp, 1.0e7_dp]
    call pc%set_press_temp_edd_profile(profile_pressure, profile_temperature, &
                                       profile_edd, err=err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    call pc%initialize_stepper(pc%var%usol_init, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    temperature_before = pc%var%temperature

    call pc%initialize_from_atmosphere_file('../tests/does_not_exist.txt', err)
    if (.not. allocated(err)) then
      print *, 'failed reinitialization did not return an error'
      stop 1
    endif
    if (.not. pc%atmosphere_initialized) then
      print *, 'failed reinitialization cleared atmosphere_initialized'
      stop 1
    endif
    if (any(pc%var%temperature /= temperature_before)) then
      print *, 'failed reinitialization changed atmospheric state'
      stop 1
    endif
    if (.not. c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'failed reinitialization destroyed the active stepper'
      stop 1
    endif
    if (.not. pc%var%press_temp_edd_profile%enabled) then
      print *, 'failed reinitialization changed the persistent profile'
      stop 1
    endif

    ! Successful reinitialization commits the replacement atmosphere and
    ! invalidates CVODE state associated with the old atmosphere.
    deallocate(err)
    call pc%initialize_from_atmosphere_file('../examples/ModernEarth/atmosphere.txt', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    if (.not. pc%atmosphere_initialized) then
      print *, 'successful reinitialization did not set atmosphere_initialized'
      stop 1
    endif
    if (c_associated(pc%wrk%sun%cvode_mem)) then
      print *, 'successful reinitialization retained stale CVODE state'
      stop 1
    endif
    if (pc%var%press_temp_edd_profile%enabled) then
      print *, 'successful legacy-file initialization retained a pressure-based profile'
      stop 1
    endif
  end subroutine

  subroutine check_not_initialized_error(err, operation)
    character(:), allocatable, intent(in) :: err
    character(len=*), intent(in) :: operation

    if (.not. allocated(err)) then
      print *, trim(operation)//' did not reject an uninitialized atmosphere'
      stop 1
    endif
    if (index(err, 'atmosphere is not initialized') == 0 .or. &
        index(err, '"'//trim(operation)//'"') == 0) then
      print *, trim(operation)//' returned an unexpected lifecycle error: '//err
      stop 1
    endif
  end subroutine

  subroutine test_legacy_file_grid()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), parameter :: expected_top = 1.0e7_dp

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    if (.not. pc%atmosphere_initialized) then
      print *, 'legacy file initialization did not set lifecycle state'
      stop 1
    endif

    if (abs(pc%var%top_atmos - expected_top) > 1.0e-6_dp) then
      print *, 'legacy file initialization did not use the atmosphere-file edge'
      stop 1
    endif
  end subroutine

  subroutine test_methods(filename)
    character(*), intent(in) :: filename
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err

    pc = EvoAtmosphere(filename, &
                       "../tests/test_settings_minimal.yaml", &
                       "../examples/ModernEarth/Sun_now.txt", &
                       "../examples/ModernEarth/atmosphere.txt", &
                       "../data", &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    ! Exercise repeated legacy-file loading for both particle-bearing and
    ! particle-free mechanisms.
    call pc%initialize_from_atmosphere_file('../examples/ModernEarth/atmosphere.txt', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    call test_set_press_temp_edd(pc)
    call test_press_temp_edd_profile(pc)
  end subroutine

  subroutine test_set_press_temp_edd(pc)
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    real(dp) :: P(2), T(2), edd(2)
    real(dp) :: T_expected(pc%var%nz), log10edd_expected(pc%var%nz)
    real(dp) :: fraction
    integer :: i

    ! Two input points are sufficient because interpolation is linear in log10(P).
    ! First test the default hydrostatic-pressure mode with profiles that differ
    ! substantially from the current atmosphere.
    P = [2.0_dp*pc%var%surface_pressure*1.0e6_dp, &
         0.5_dp*pc%wrk%pressure_hydro(pc%var%nz)]
    T = [310.0_dp, 160.0_dp]
    edd = [1.0e7_dp, 2.0e5_dp]

    call pc%set_press_temp_edd(P, T, edd, pc%wrk%pressure_hydro(10), &
                               hydro_pressure=.true., err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    do i = 1,pc%var%nz
      fraction = (log10(pc%wrk%pressure_hydro(i))-log10(P(2)))/ &
                 (log10(P(1))-log10(P(2)))
      fraction = min(max(fraction,0.0_dp),1.0_dp)
      T_expected(i) = T(2) + fraction*(T(1)-T(2))
      log10edd_expected(i) = log10(edd(2)) + fraction*(log10(edd(1))-log10(edd(2)))
    enddo
    if (maxval(abs(pc%var%temperature-T_expected)) > 1.0e-6_dp) then
      print*,'set_press_temp_edd did not reproduce the hydrostatic P-T profile'
      stop 1
    endif
    if (maxval(abs(log10(pc%var%edd)-log10edd_expected)) > 1.0e-8_dp) then
      print*,'set_press_temp_edd did not reproduce the hydrostatic P-Kzz profile'
      stop 1
    endif

    ! Then test actual-pressure mode. This pressure is density*k*T and each
    ! layer can therefore be solved independently.
    P = [2.0_dp*maxval(pc%wrk%pressure), 0.5_dp*minval(pc%wrk%pressure)]
    T = [280.0_dp, 140.0_dp]
    edd = [2.0e8_dp, 8.0e4_dp]

    call pc%set_press_temp_edd(P, T, edd, pc%wrk%pressure(10), &
                               hydro_pressure=.false., err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    do i = 1,pc%var%nz
      fraction = (log10(pc%wrk%pressure(i))-log10(P(2)))/ &
                 (log10(P(1))-log10(P(2)))
      fraction = min(max(fraction,0.0_dp),1.0_dp)
      T_expected(i) = T(2) + fraction*(T(1)-T(2))
      log10edd_expected(i) = log10(edd(2)) + fraction*(log10(edd(1))-log10(edd(2)))
    enddo
    if (maxval(abs(pc%var%temperature-T_expected)) > 1.0e-6_dp) then
      print*,'set_press_temp_edd did not reproduce the actual P-T profile'
      stop 1
    endif
    if (maxval(abs(log10(pc%var%edd)-log10edd_expected)) > 1.0e-8_dp) then
      print*,'set_press_temp_edd did not reproduce the actual P-Kzz profile'
      stop 1
    endif

  end subroutine

  subroutine test_press_temp_edd_profile(pc)
    use photochem_const, only: k_boltz
    type(EvoAtmosphere), intent(inout) :: pc
    character(:), allocatable :: err
    real(dp) :: P(2), T(2), edd(2)
    real(dp) :: P_bad(2)
    real(dp) :: usol_original(pc%dat%nq,pc%var%nz)
    real(dp) :: usol_trial(pc%dat%nq,pc%var%nz)
    real(dp), target :: usol_flat(pc%var%neqs)
    real(dp) :: rhs(pc%var%neqs)
    real(dp) :: temperature_before(pc%var%nz)
    real(dp) :: temperature_disabled(pc%var%nz)
    real(dp), allocatable :: gas_mix_top(:)
    real(dp) :: top_atmos_original, z_top, density_top, pressure_top, temperature_top
    real(dp) :: tn

    usol_original = pc%wrk%usol
    P = [2.0_dp*pc%var%surface_pressure*1.0e6_dp, &
         0.5_dp*pc%wrk%pressure_hydro(pc%var%nz)]
    T = [300.0_dp, 180.0_dp]
    edd = [3.0e7_dp, 4.0e5_dp]

    ! Invalid input must fail without enabling or otherwise changing the mode.
    P_bad = [P(2), P(1)]
    call pc%set_press_temp_edd_profile(P_bad, T, edd, &
         pc%wrk%pressure_hydro(pc%var%nz/2), hydro_pressure=.true., err=err)
    if (.not.allocated(err)) then
      print*,'An increasing persistent pressure profile was accepted'
      stop 1
    endif
    deallocate(err)
    if (pc%var%press_temp_edd_profile%enabled) then
      print*,'An invalid persistent profile left the mode enabled'
      stop 1
    endif

    call pc%set_press_temp_edd_profile(P, T, edd, &
         pc%wrk%pressure_hydro(pc%var%nz/2), hydro_pressure=.true., err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    call check_press_temp_edd_profile(pc, P, T, edd)

    ! Altitude-based setters must not be silently overridden by persistent
    ! mode. The caller must first make the mode change explicit by clearing it.
    call pc%set_temperature(pc%var%temperature, err=err)
    if (.not.allocated(err)) then
      print*,'set_temperature was accepted while a persistent profile was enabled'
      stop 1
    endif
    deallocate(err)
    call pc%set_press_temp_edd(P, T, edd, &
         pc%wrk%pressure_hydro(pc%var%nz/2), hydro_pressure=.true., err=err)
    if (.not.allocated(err)) then
      print*,'set_press_temp_edd was accepted while a persistent profile was enabled'
      stop 1
    endif
    deallocate(err)

    ! Exercise repeated RHS and Jacobian preparation through CVODE while the
    ! persistent mode is enabled.
    call pc%initialize_stepper(usol_original, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

    ! Changing persistent mode would invalidate CVODE's current history and
    ! must therefore require an explicit stepper destruction/reinitialization.
    call pc%set_press_temp_edd_profile(P, T, edd, &
         pc%wrk%pressure_hydro(pc%var%nz/2), hydro_pressure=.true., err=err)
    if (.not.allocated(err)) then
      print*,'Persistent profile was replaced while a stepper was initialized'
      stop 1
    endif
    deallocate(err)
    call pc%clear_press_temp_edd_profile(err)
    if (.not.allocated(err)) then
      print*,'Persistent profile was cleared while a stepper was initialized'
      stop 1
    endif
    deallocate(err)
    if (.not.pc%var%press_temp_edd_profile%enabled) then
      print*,'Rejected profile clear changed persistent mode'
      stop 1
    endif

    tn = pc%step(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    call check_press_temp_edd_profile(pc, P, T, edd)

    ! Public vertical regridding must preserve the pressure-based inputs and
    ! remap them onto both the changed grid and the restored grid.
    top_atmos_original = pc%var%top_atmos
    call pc%update_vertical_grid(TOA_alt=0.98_dp*top_atmos_original, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    if (.not.pc%var%press_temp_edd_profile%enabled) then
      print*,'Vertical regridding disabled the persistent profile'
      stop 1
    endif
    call check_press_temp_edd_profile(pc, P, T, edd)
    call pc%update_vertical_grid(TOA_alt=top_atmos_original, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    call check_press_temp_edd_profile(pc, P, T, edd)

    ! Raising the top while persistence is active must use the mapped
    ! temperature in the hydrostatic continuation, not a stale altitude-based
    ! extrapolation.
    z_top = pc%var%z(pc%var%nz)
    density_top = sum(pc%wrk%usol(pc%dat%ng_1:,pc%var%nz))
    pressure_top = density_top*k_boltz*pc%var%temperature(pc%var%nz)
    temperature_top = pc%var%temperature(pc%var%nz)
    gas_mix_top = pc%wrk%usol(pc%dat%ng_1:,pc%var%nz)/density_top
    call pc%initialize_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    call pc%update_vertical_grid(TOA_alt=1.20_dp*top_atmos_original, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    call check_press_temp_edd_profile(pc, P, T, edd)
    call check_hydrostatic_extension(pc, z_top, pressure_top, temperature_top, gas_mix_top)
    tn = pc%step(err)
    if (.not.allocated(err)) then
      print*,'A successful vertical regrid retained the stale CVODE stepper'
      stop 1
    endif
    deallocate(err)
    ! Match the gas-giant maintenance sequence after the successful regrid has
    ! invalidated the old solver resources.
    call pc%initialize_stepper(pc%wrk%usol, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    call pc%destroy_stepper(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    call pc%update_vertical_grid(TOA_pressure=pressure_top, err=err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    if (abs(log10(pc%wrk%pressure(pc%var%nz)/pressure_top)) > 2.0e-8_dp) then
      print*,'Persistent-profile TOA-pressure solve missed its pressure target'
      stop 1
    endif
    call check_press_temp_edd_profile(pc, P, T, edd)

    ! A different trial composition must produce a different altitude-based
    ! temperature while still reproducing the prescribed pressure profiles.
    temperature_before = pc%var%temperature
    usol_trial = usol_original
    usol_trial(pc%dat%ng_1:,:) = 1.5_dp*usol_trial(pc%dat%ng_1:,:)
    usol_flat = reshape(usol_trial, [pc%var%neqs])
    call pc%right_hand_side(pc%var%neqs, 0.0_dp, usol_flat, rhs, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    if (maxval(abs(pc%var%temperature-temperature_before)) < 1.0e-3_dp) then
      print*,'Persistent P-T-Kzz profile did not respond to a changed trial composition'
      stop 1
    endif
    call check_press_temp_edd_profile(pc, P, T, edd)

    ! Clearing the profile leaves the latest mapped altitude profiles in
    ! place but prevents subsequent RHS calls from remapping them.
    call pc%clear_press_temp_edd_profile(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    temperature_disabled = pc%var%temperature
    usol_trial(pc%dat%ng_1:,:) = 0.5_dp*usol_trial(pc%dat%ng_1:,:)
    usol_flat = reshape(usol_trial, [pc%var%neqs])
    call pc%right_hand_side(pc%var%neqs, 0.0_dp, usol_flat, rhs, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    if (any(pc%var%temperature /= temperature_disabled)) then
      print*,'Cleared P-T-Kzz profile continued to alter temperature'
      stop 1
    endif

    call pc%prep_atmosphere(usol_original, err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif

  end subroutine

  subroutine check_hydrostatic_extension(pc, z_anchor, pressure_anchor, &
                                         temperature_anchor, gas_mix_anchor)
    use photochem_const, only: k_boltz, N_avo
    type(EvoAtmosphere), intent(in) :: pc
    real(dp), intent(in) :: z_anchor, pressure_anchor, temperature_anchor
    real(dp), intent(in) :: gas_mix_anchor(:)
    real(dp) :: density_layer, pressure_layer, pressure_expected
    real(dp) :: pressure_previous, temperature_previous, delta_z, mubar
    real(dp) :: gas_mix(size(gas_mix_anchor))
    integer :: i, first_extended

    first_extended = findloc(pc%var%z > z_anchor, .true., dim=1)
    if (first_extended == 0) then
      print *, 'Hydrostatic regrid test did not extend beyond the old top center'
      stop 1
    endif

    pressure_previous = pressure_anchor
    temperature_previous = temperature_anchor
    delta_z = pc%var%z(first_extended)-z_anchor
    do i = first_extended,pc%var%nz
      if (i > first_extended) delta_z = pc%var%z(i)-pc%var%z(i-1)
      density_layer = sum(pc%wrk%usol(pc%dat%ng_1:,i))
      gas_mix = pc%wrk%usol(pc%dat%ng_1:,i)/density_layer
      if (abs(sum(gas_mix)-1.0_dp) > 1.0e-14_dp .or. &
          maxval(abs(gas_mix-gas_mix_anchor)) > 1.0e-12_dp) then
        print *, 'Gas mixing ratios were not normalized and constant above the old model top'
        stop 1
      endif
      mubar = sum(pc%dat%species_mass(pc%dat%ng_1:pc%dat%nq)*gas_mix)
      pressure_expected = pressure_previous*exp( &
          -(mubar*pc%var%grav(i)*delta_z)/ &
           (N_avo*k_boltz*0.5_dp*(temperature_previous+pc%var%temperature(i))))
      pressure_layer = density_layer*k_boltz*pc%var%temperature(i)
      if (pressure_layer >= pressure_previous .or. &
          abs(pressure_layer/pressure_expected-1.0_dp) > 2.0e-9_dp) then
        print *, 'Density above the old model top is not a hydrostatic continuation'
        stop 1
      endif
      pressure_previous = pressure_layer
      temperature_previous = pc%var%temperature(i)
    enddo

  end subroutine

  subroutine check_press_temp_edd_profile(pc, P, T, edd)
    type(EvoAtmosphere), intent(in) :: pc
    real(dp), intent(in) :: P(:), T(:), edd(:)
    real(dp) :: T_expected(pc%var%nz), log10edd_expected(pc%var%nz)
    real(dp) :: fraction
    integer :: i

    do i = 1,pc%var%nz
      fraction = (log10(pc%wrk%pressure_hydro(i))-log10(P(2)))/ &
                 (log10(P(1))-log10(P(2)))
      fraction = min(max(fraction,0.0_dp),1.0_dp)
      T_expected(i) = T(2) + fraction*(T(1)-T(2))
      log10edd_expected(i) = log10(edd(2)) + &
          fraction*(log10(edd(1))-log10(edd(2)))
    enddo
    if (maxval(abs(pc%var%temperature-T_expected)) > 1.0e-6_dp) then
      print*,'Persistent profile did not reproduce the prescribed P-T profile'
      stop 1
    endif
    if (maxval(abs(log10(pc%var%edd)-log10edd_expected)) > 1.0e-8_dp) then
      print*,'Persistent profile did not reproduce the prescribed P-Kzz profile'
      stop 1
    endif

  end subroutine

end program
