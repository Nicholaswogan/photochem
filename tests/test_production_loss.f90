program test_production_loss
  ! Characterization tests for the current production-and-loss diagnostic.
  ! These checks intentionally cover only reactions and rainout. Later 8.8
  ! passes will extend reconciliation to every term in the full RHS.
  use photochem, only: EvoAtmosphere, ProductionLoss, dp
  implicit none

  call test_reaction_accounting()
  call test_rainout_accounting()
  print *, 'test_production_loss passed'

contains

  subroutine test_reaction_accounting()
    type(EvoAtmosphere) :: pc
    type(ProductionLoss) :: pl
    character(:), allocatable :: err

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    call check_error(err)

    call pc%production_and_loss('HCN', pc%wrk%usol, pl, err)
    call check_error(err)
    call check_result_structure(pc, pl)
    call check_chemistry_reconciliation(pc, 'HCN', pl)
  end subroutine

  subroutine test_rainout_accounting()
    type(EvoAtmosphere) :: pc
    type(ProductionLoss) :: pl
    character(:), allocatable :: err
    real(dp), allocatable :: rainout_column(:), expected(:)
    integer :: i, j, rainout_ind, species_ind

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_rainout.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    call check_error(err)

    allocate(rainout_column(pc%dat%nq), source=0.0_dp)
    do i = 1,pc%dat%nq
      rainout_column(i) = sum(pc%wrk%rainout_rates(i,1:pc%var%trop_ind)* &
                              pc%wrk%usol(i,1:pc%var%trop_ind)* &
                              pc%var%dz(1:pc%var%trop_ind))
    enddo
    species_ind = maxloc(rainout_column, dim=1)
    if (rainout_column(species_ind) <= 0.0_dp) then
      print *, 'Production-loss rainout case has no positive rainout loss'
      stop 1
    endif

    call pc%production_and_loss(trim(pc%dat%species_names(species_ind)), &
                                pc%wrk%usol, pl, err)
    call check_error(err)
    call check_result_structure(pc, pl)

    rainout_ind = 0
    do i = 1,size(pl%loss_rx)
      if (trim(pl%loss_rx(i)) == 'rainout') then
        rainout_ind = i
        exit
      endif
    enddo
    if (rainout_ind == 0) then
      print *, 'Production-loss result omitted configured rainout'
      stop 1
    endif

    allocate(expected(pc%var%nz), source=0.0_dp)
    do j = 1,pc%var%trop_ind
      expected(j) = pc%wrk%rainout_rates(species_ind,j)* &
                    pc%wrk%usol(species_ind,j)
    enddo
    call check_close(pl%loss(:,rainout_ind), expected, &
                     'Reported rainout profile does not match prepared rates')
    call check_chemistry_reconciliation( &
        pc, trim(pc%dat%species_names(species_ind)), pl)
  end subroutine

  subroutine check_result_structure(pc, pl)
    type(EvoAtmosphere), intent(in) :: pc
    type(ProductionLoss), intent(in) :: pl
    real(dp) :: expected
    integer :: i

    if (size(pl%production,1) /= pc%var%nz .or. &
        size(pl%loss,1) /= pc%var%nz .or. &
        size(pl%production,2) /= size(pl%production_rx) .or. &
        size(pl%loss,2) /= size(pl%loss_rx) .or. &
        size(pl%integrated_production) /= size(pl%production_rx) .or. &
        size(pl%integrated_loss) /= size(pl%loss_rx)) then
      print *, 'Production-loss result has inconsistent dimensions'
      stop 1
    endif
    if (any(pl%production < 0.0_dp) .or. any(pl%loss < 0.0_dp)) then
      print *, 'Production-loss result contains a negative contribution'
      stop 1
    endif

    do i = 1,size(pl%integrated_production)
      expected = sum(pl%production(:,i)*pc%var%dz)
      call check_scalar_close(pl%integrated_production(i), expected, &
                              'Integrated production is inconsistent')
      if (i > 1 .and. pl%integrated_production(i) > &
                     pl%integrated_production(i-1)) then
        print *, 'Production contributions are not sorted'
        stop 1
      endif
    enddo
    do i = 1,size(pl%integrated_loss)
      expected = sum(pl%loss(:,i)*pc%var%dz)
      call check_scalar_close(pl%integrated_loss(i), expected, &
                              'Integrated loss is inconsistent')
      if (i > 1 .and. pl%integrated_loss(i) > pl%integrated_loss(i-1)) then
        print *, 'Loss contributions are not sorted'
        stop 1
      endif
    enddo
  end subroutine

  subroutine check_chemistry_reconciliation(pc, species, pl)
    type(EvoAtmosphere), intent(inout) :: pc
    character(len=*), intent(in) :: species
    type(ProductionLoss), intent(in) :: pl
    character(:), allocatable :: err
    real(dp), allocatable :: rhs(:), reported(:), expected(:)
    integer :: j, species_ind

    species_ind = findloc(pc%dat%species_names(1:pc%dat%nq), species, 1)
    if (species_ind == 0) then
      print *, 'Reconciliation species is not evolved: ', trim(species)
      stop 1
    endif

    allocate(rhs(pc%var%neqs), reported(pc%var%nz), expected(pc%var%nz))
    call pc%chemistry_right_hand_side(pc%wrk%usol, rhs, err)
    call check_error(err)

    reported = sum(pl%production, dim=2) - sum(pl%loss, dim=2)
    do j = 1,pc%var%nz
      expected(j) = rhs(species_ind + (j - 1)*pc%dat%nq)
    enddo
    call check_close(reported, expected, &
        'Reported reactions and rainout do not reconstruct chemistry RHS')
  end subroutine

  subroutine check_close(actual, expected, message)
    real(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: message
    real(dp) :: error, scale

    scale = max(maxval(abs(actual)), maxval(abs(expected)), 1.0e-100_dp)
    error = maxval(abs(actual - expected))/scale
    if (error > 5.0e-13_dp) then
      print *, trim(message), error
      stop 1
    endif
  end subroutine

  subroutine check_scalar_close(actual, expected, message)
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: message
    real(dp) :: scale

    scale = max(abs(actual), abs(expected), 1.0e-100_dp)
    if (abs(actual - expected)/scale > 5.0e-14_dp) then
      print *, trim(message), actual, expected
      stop 1
    endif
  end subroutine

  subroutine check_error(err)
    character(:), allocatable, intent(in) :: err

    if (allocated(err)) then
      print *, trim(err)
      stop 1
    endif
  end subroutine

end program
