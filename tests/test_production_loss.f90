program test_production_loss
  ! Characterization tests for the current production-and-loss diagnostic.
  ! These checks cover local chemistry and internal vertical transport. Later
  ! passes will extend reconciliation to every term in the full RHS.
  use photochem, only: EvoAtmosphere, ProductionLoss, dp
  implicit none

  call test_reaction_accounting()
  call test_rainout_accounting()
  call test_condensation_and_evaporation_accounting()
  call test_custom_rate_accounting()
  call test_zahnle_escape_accounting()
  call test_vertical_transport_accounting()
  call test_short_lived_accounting()
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
    call check_reaction_consolidation(pc, 'HCN', pl)
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

  subroutine test_condensation_and_evaporation_accounting()
    use photochem_enum, only: CondensingParticle
    type(EvoAtmosphere) :: pc
    type(ProductionLoss) :: pl
    character(:), allocatable :: err
    real(dp), allocatable :: usol(:,:)
    integer :: gas_ind, i, j

    pc = EvoAtmosphere('../data/reaction_mechanisms/zahnle_earth.yaml', &
                       '../tests/test_settings_condensation_jacobian.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    call check_error(err)
    i = findloc(pc%dat%particle_formation_method, CondensingParticle, 1)
    if (i == 0) then
      print *, 'Production-loss condensation case has no condensing particle'
      stop 1
    endif
    gas_ind = pc%dat%particle_gas_phase_ind(i)
    usol = pc%wrk%usol
    do j = 1,pc%var%nz
      if (j <= pc%var%nz/2) then
        usol(gas_ind,j) = 1.25_dp*pc%var%cond_params(i)%RHc* &
                          pc%wrk%gas_sat_den(i,j)
      else
        usol(gas_ind,j) = 0.75_dp*pc%var%cond_params(i)%RHc* &
                          pc%wrk%gas_sat_den(i,j)
      endif
      usol(i,j) = max(usol(i,j), 1.0e5_dp)
    enddo

    call pc%production_and_loss(trim(pc%dat%species_names(gas_ind)), &
                                usol, pl, err)
    call check_error(err)
    call check_result_structure(pc, pl)
    call require_label(pl%loss_rx, 'condensation')
    call require_label(pl%production_rx, 'evaporation')
    call check_chemistry_reconciliation( &
        pc, trim(pc%dat%species_names(gas_ind)), pl)

    call pc%production_and_loss(trim(pc%dat%species_names(i)), usol, pl, err)
    call check_error(err)
    call check_result_structure(pc, pl)
    call require_label(pl%production_rx, 'condensation')
    call require_label(pl%loss_rx, 'evaporation')
    call check_chemistry_reconciliation( &
        pc, trim(pc%dat%species_names(i)), pl)
  end subroutine

  subroutine test_custom_rate_accounting()
    use photochem_vars, only: time_dependent_rate_fcn
    type(EvoAtmosphere) :: pc
    type(ProductionLoss) :: pl
    procedure(time_dependent_rate_fcn), pointer :: rate_fcn
    character(:), allocatable :: err

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    call check_error(err)
    rate_fcn => signed_custom_rate
    call pc%set_rate_fcn('HCN', rate_fcn, err)
    call check_error(err)
    pc%wrk%tn = 7.0_dp

    call pc%production_and_loss('HCN', pc%wrk%usol, pl, err)
    call check_error(err)
    call check_result_structure(pc, pl)
    call require_label(pl%production_rx, 'custom rate')
    call require_label(pl%loss_rx, 'custom rate')
    call check_chemistry_reconciliation(pc, 'HCN', pl)
  end subroutine

  subroutine test_zahnle_escape_accounting()
    type(EvoAtmosphere) :: pc
    type(ProductionLoss) :: pl
    character(:), allocatable :: err

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_zahnle_escape.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    call check_error(err)

    call pc%production_and_loss('H2', pc%wrk%usol, pl, err)
    call check_error(err)
    call check_result_structure(pc, pl)
    call require_label(pl%loss_rx, 'hydrogen escape')
    call check_chemistry_reconciliation(pc, 'H2', pl)
  end subroutine

  subroutine test_short_lived_accounting()
    type(EvoAtmosphere) :: pc
    type(ProductionLoss) :: pl
    character(:), allocatable :: err
    real(dp) :: balance_scale

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_short_lived.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    call check_error(err)

    call pc%production_and_loss('O1D', pc%wrk%usol, pl, err)
    call check_error(err)
    call check_result_structure(pc, pl)
    balance_scale = max(maxval(pl%total_production), &
                        maxval(pl%total_loss), 1.0e-100_dp)
    if (maxval(abs(pl%net))/balance_scale > 5.0e-13_dp) then
      print *, 'Short-lived production and loss are not algebraically balanced'
      stop 1
    endif
  end subroutine

  subroutine test_vertical_transport_accounting()
    type(EvoAtmosphere) :: pc
    type(ProductionLoss) :: pl
    character(:), allocatable :: err
    real(dp) :: integrated_transport, transport_scale
    integer :: i

    pc = EvoAtmosphere('../tests/no_particle_test.yaml', &
                       '../tests/test_settings_minimal.yaml', &
                       '../examples/ModernEarth/Sun_now.txt', &
                       '../examples/ModernEarth/atmosphere.txt', &
                       '../data', err)
    call check_error(err)

    call pc%production_and_loss('HCN', pc%wrk%usol, pl, err)
    call check_error(err)
    call check_result_structure(pc, pl)

    integrated_transport = 0.0_dp
    transport_scale = 0.0_dp
    do i = 1,size(pl%production_rx)
      if (trim(pl%production_rx(i)) == 'vertical transport') then
        integrated_transport = integrated_transport + &
                               pl%integrated_production(i)
        transport_scale = transport_scale + pl%integrated_production(i)
      endif
    enddo
    do i = 1,size(pl%loss_rx)
      if (trim(pl%loss_rx(i)) == 'vertical transport') then
        integrated_transport = integrated_transport - pl%integrated_loss(i)
        transport_scale = transport_scale + pl%integrated_loss(i)
      endif
    enddo
    if (transport_scale == 0.0_dp) then
      print *, 'Production-loss result omitted vertical transport'
      stop 1
    endif
    if (abs(integrated_transport)/transport_scale > 5.0e-13_dp) then
      print *, 'Internal vertical transport is not column conserving', &
               integrated_transport
      stop 1
    endif
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
    if (size(pl%total_production) /= pc%var%nz .or. &
        size(pl%total_loss) /= pc%var%nz .or. &
        size(pl%net) /= pc%var%nz) then
      print *, 'Production-loss totals have inconsistent dimensions'
      stop 1
    endif
    call check_close(pl%total_production, sum(pl%production, dim=2), &
                     'Total production is inconsistent')
    call check_close(pl%total_loss, sum(pl%loss, dim=2), &
                     'Total loss is inconsistent')
    call check_close(pl%net, pl%total_production-pl%total_loss, &
                     'Net production-loss tendency is inconsistent')

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
    call check_scalar_close(pl%integrated_total_production, &
                            sum(pl%integrated_production), &
                            'Integrated total production is inconsistent')
    call check_scalar_close(pl%integrated_total_loss, &
                            sum(pl%integrated_loss), &
                            'Integrated total loss is inconsistent')
    call check_scalar_close(pl%integrated_net, &
                            pl%integrated_total_production- &
                            pl%integrated_total_loss, &
                            'Integrated net tendency is inconsistent')
  end subroutine

  subroutine check_reaction_consolidation(pc, species, pl)
    type(EvoAtmosphere), intent(in) :: pc
    character(len=*), intent(in) :: species
    type(ProductionLoss), intent(in) :: pl
    integer :: expected, actual, i, j, species_ind

    species_ind = findloc(pc%dat%species_names(1:pc%dat%nsp), species, 1)
    expected = 0
    do i = 1,pc%dat%pl(species_ind)%nump
      if (.not.any(pc%dat%pl(species_ind)%iprod(1:i-1) == &
                  pc%dat%pl(species_ind)%iprod(i))) expected = expected + 1
    enddo
    actual = 0
    do i = 1,size(pl%production_rx)
      if (any(pl%production_rx(i) == &
              pc%dat%reaction_equations( &
                  pc%dat%pl(species_ind)%iprod))) actual = actual + 1
    enddo
    if (actual /= expected) then
      print *, 'Repeated production stoichiometry was not consolidated'
      stop 1
    endif
    do i = 1,size(pl%production_rx)
      do j = i+1,size(pl%production_rx)
        if (pl%production_rx(i) == pl%production_rx(j)) then
          print *, 'Production result contains a duplicate reaction label'
          stop 1
        endif
      enddo
    enddo
  end subroutine

  subroutine require_label(labels, label)
    character(len=*), intent(in) :: labels(:)
    character(len=*), intent(in) :: label
    integer :: i

    do i = 1,size(labels)
      if (trim(labels(i)) == label) return
    enddo
    print *, 'Production-loss result omitted process: ', trim(label)
    stop 1
  end subroutine

  subroutine check_chemistry_reconciliation(pc, species, pl)
    type(EvoAtmosphere), intent(inout) :: pc
    character(len=*), intent(in) :: species
    type(ProductionLoss), intent(in) :: pl
    character(:), allocatable :: err
    real(dp), allocatable :: rhs(:), reported(:), expected(:), transport(:)
    integer :: i, j, species_ind

    species_ind = findloc(pc%dat%species_names(1:pc%dat%nq), species, 1)
    if (species_ind == 0) then
      print *, 'Reconciliation species is not evolved: ', trim(species)
      stop 1
    endif

    allocate(rhs(pc%var%neqs), reported(pc%var%nz), expected(pc%var%nz), &
             transport(pc%var%nz), source=0.0_dp)
    call pc%chemistry_right_hand_side(pc%wrk%usol, rhs, err)
    call check_error(err)

    reported = sum(pl%production, dim=2) - sum(pl%loss, dim=2)
    do i = 1,size(pl%production_rx)
      if (trim(pl%production_rx(i)) == 'vertical transport') then
        transport = transport + pl%production(:,i)
      endif
    enddo
    do i = 1,size(pl%loss_rx)
      if (trim(pl%loss_rx(i)) == 'vertical transport') then
        transport = transport - pl%loss(:,i)
      endif
    enddo
    do j = 1,pc%var%nz
      expected(j) = rhs(species_ind + (j - 1)*pc%dat%nq)
    enddo
    expected = expected + transport
    call check_close(reported, expected, &
        'Reported local and transport terms do not reconstruct their RHS')
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

  subroutine signed_custom_rate(tn, nz, rate)
    use iso_c_binding, only: c_double, c_int
    real(c_double), value, intent(in) :: tn
    integer(c_int), value, intent(in) :: nz
    real(c_double), intent(out) :: rate(nz)

    rate(1:nz/2) = tn*2.0e4_dp
    rate(nz/2+1:nz) = -tn*3.0e4_dp
  end subroutine

end program
