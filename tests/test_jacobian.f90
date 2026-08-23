program test_jacobian
  ! Focused numerical validation of Photochem Jacobian implementations.
  use photochem, only: AnalyticalJacobian, AutodiffJacobian, &
                       EvoAtmosphere, dp
  implicit none

  call test_analytical_jacobian_core()
  call test_analytical_jacobian_short_lived()
  call test_analytical_jacobian_rainout()
  print *, 'test_jacobian passed'

contains

  subroutine test_analytical_jacobian_core()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: usol_flat(:)
    integer :: h_index, j

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    if (pc%dat%nsl /= 0 .or. pc%dat%there_are_particles .or. &
        pc%dat%gas_rainout) then
      print *, 'Core analytical Jacobian test unexpectedly requires fallback'
      stop 1
    endif

    usol_flat = reshape(pc%wrk%usol, [pc%var%neqs])
    call compare_analytical_jacobian(pc, usol_flat, &
                                     'ordinary reaction state')

    ! Exercise a repeated-reactant derivative at zero reactant density.
    h_index = findloc(pc%dat%species_names, 'H', dim=1)
    if (h_index < pc%dat%ng_1 .or. h_index > pc%dat%nq) then
      print *, 'Could not locate evolved H for analytical Jacobian test'
      stop 1
    endif
    do j = 1,pc%var%nz
      usol_flat(h_index+pc%dat%nq*(j-1)) = 0.0_dp
    enddo
    call compare_analytical_jacobian( &
        pc, usol_flat, 'zero-density repeated-reactant state')

  end subroutine

  subroutine test_analytical_jacobian_short_lived()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: usol_flat(:)

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_short_lived.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    if (pc%dat%nsl /= 2 .or. pc%dat%there_are_particles .or. &
        pc%dat%gas_rainout) then
      print *, 'Short-lived Jacobian test has unexpected model features'
      stop 1
    endif

    usol_flat = reshape(pc%wrk%usol, [pc%var%neqs])
    call compare_analytical_jacobian(pc, usol_flat, &
                                     'short-lived species state')

  end subroutine

  subroutine test_analytical_jacobian_rainout()
    type(EvoAtmosphere) :: pc
    character(:), allocatable :: err
    real(dp), allocatable :: usol_flat(:)

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_rainout.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif

    if (pc%dat%nsl /= 0 .or. pc%dat%there_are_particles .or. &
        .not.pc%dat%gas_rainout) then
      print *, 'Rainout Jacobian test has unexpected model features'
      stop 1
    endif
    if (pc%var%trop_ind < 1 .or. &
        maxval(pc%wrk%rainout_rates(:,1:pc%var%trop_ind)) <= 0.0_dp) then
      print *, 'Rainout Jacobian test did not produce positive rainout rates'
      stop 1
    endif

    usol_flat = reshape(pc%wrk%usol, [pc%var%neqs])
    call compare_analytical_jacobian(pc, usol_flat, 'gas-rainout state')
    call check_analytical_rainout_contribution(pc, usol_flat)

  end subroutine

  subroutine compare_analytical_jacobian(pc, usol_flat, label)
    type(EvoAtmosphere), intent(inout) :: pc
    real(dp), intent(in) :: usol_flat(:)
    character(*), intent(in) :: label

    character(:), allocatable :: err
    real(dp), allocatable :: jac(:), jac_autodiff(:,:), jac_analytical(:,:)
    real(dp) :: jacobian_scale, relative_error

    allocate(jac(pc%dat%lda*pc%var%neqs))

    pc%var%jacobian_method = AutodiffJacobian
    call pc%jacobian(size(jac), pc%var%neqs, usol_flat, jac, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    jac_autodiff = pc%wrk%djac_chem

    pc%var%jacobian_method = AnalyticalJacobian
    call pc%jacobian(size(jac), pc%var%neqs, usol_flat, jac, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    jac_analytical = pc%wrk%djac_chem

    jacobian_scale = max(1.0_dp, maxval(abs(jac_autodiff)), &
                         maxval(abs(jac_analytical)))
    relative_error = &
        maxval(abs(jac_analytical-jac_autodiff))/jacobian_scale
    if (relative_error > 1.0e-12_dp) then
      print *, 'Core analytical Jacobian differs from autodiff for ', &
               trim(label), relative_error
      stop 1
    endif

  end subroutine

  subroutine check_analytical_rainout_contribution(pc, usol_flat)
    type(EvoAtmosphere), intent(inout) :: pc
    real(dp), intent(in) :: usol_flat(:)

    character(:), allocatable :: err
    real(dp), allocatable :: jac(:), chemistry_with(:,:)
    real(dp), allocatable :: chemistry_without(:,:), expected(:,:)
    real(dp) :: scale, relative_error
    integer :: i, j, k

    allocate(jac(pc%dat%lda*pc%var%neqs))
    allocate(expected(pc%dat%nq,pc%var%neqs))
    expected = 0.0_dp
    do j = 1,pc%var%trop_ind
      do i = 1,pc%dat%nq
        k = i + pc%dat%nq*(j-1)
        expected(i,k) = -pc%wrk%rainout_rates(i,j)
      enddo
    enddo

    pc%var%jacobian_method = AnalyticalJacobian
    call pc%jacobian(size(jac), pc%var%neqs, usol_flat, jac, err)
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    chemistry_with = pc%wrk%djac_chem

    pc%dat%gas_rainout = .false.
    call pc%jacobian(size(jac), pc%var%neqs, usol_flat, jac, err)
    pc%dat%gas_rainout = .true.
    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
    chemistry_without = pc%wrk%djac_chem

    scale = maxval(abs(expected))
    relative_error = maxval(abs(chemistry_with-chemistry_without-expected))/ &
                     scale
    if (relative_error > 1.0e-5_dp) then
      print *, 'Analytical rainout contribution is incorrect', relative_error
      stop 1
    endif

  end subroutine

end program
