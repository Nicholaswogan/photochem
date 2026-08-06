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
    call test_methods('../data/reaction_mechanisms/zahnle_earth.yaml')
    call test_methods('../tests/no_particle_test.yaml')
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
                       "../examples/ModernEarth/settings.yaml", &
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
