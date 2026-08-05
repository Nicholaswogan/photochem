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
    call test_top_from_atmosphere_file()
    call test_methods('../data/reaction_mechanisms/zahnle_earth.yaml')
    call test_methods('../tests/no_particle_test.yaml')
  end subroutine

  subroutine test_initialization_state()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_top_atmospherefile.yaml', &
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
                       '../tests/test_settings_top_atmospherefile.yaml', &
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

    ! Static configuration methods remain available before atmosphere
    ! initialization.
    deallocate(err)
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

  subroutine test_top_from_atmosphere_file()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), parameter :: expected_top = 1.0e7_dp

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_top_atmospherefile.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', &
                       err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    if (.not. pc%atmosphere_initialized) then
      print *, 'top: atmospherefile initialization did not set lifecycle state'
      stop 1
    endif

    if (.not. pc%var%top_atmos_from_file) then
      print *, 'top: atmospherefile was not retained during static setup'
      stop 1
    endif
    if (abs(pc%var%top_atmos - expected_top) > 1.0e-6_dp) then
      print *, 'top: atmospherefile did not resolve to the atmosphere-file edge'
      stop 1
    endif
  end subroutine

  subroutine test_methods(filename)
    character(*), intent(in) :: filename
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err

    pc = EvoAtmosphere(filename, &
                       "../examples/ModernEarth/settings.yaml", &
                       "../examples/ModernEarth/Sun_now.txt", &
                       "../examples/ModernEarth/atmosphere.txt", &
                       "../data", &
                       err)
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
    real(dp) :: top_atmos_original
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
