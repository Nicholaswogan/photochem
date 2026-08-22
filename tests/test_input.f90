program test_input
  ! Focused input parsing and validation tests. See tests/README.md.
  implicit none
  call test_parsing()
  call test_removed_evolve_climate_setting()
  call test_removed_conserving_initialization()
  call test_removed_water_settings()
  print *, 'test_input passed'
contains

  subroutine test_parsing()
    call expect_valid('H2 + O => H2O', .false., [character(len=2) :: 'H2', 'O'], ['H2O'])
    call expect_valid('H + OH <=> H2O', .true., [character(len=2) :: 'H', 'OH'], ['H2O'])
    call expect_valid('H+ + e- => H', .false., ['H+', 'e-'], ['H'])
    call expect_valid('H2 + O (+M) => H2O (+M)', .false., &
                      [character(len=2) :: 'H2', 'O', 'M'], [character(len=3) :: 'H2O', 'M'])
    call expect_valid('H2 + O ( + M ) => H2O ( + M )', .false., &
                      [character(len=2) :: 'H2', 'O', 'M'], [character(len=3) :: 'H2O', 'M'])

    call expect_invalid('', 'empty')
    call expect_invalid('H2 + O', 'no valid arrow')
    call expect_invalid('H2 + O =>', 'no products')
    call expect_invalid('=> H2O', 'no reactants')
    call expect_invalid('H2 + + O => H2O', 'without a species')
    call expect_invalid('H2 + O + => H2O', 'ends with a plus sign')
    call expect_invalid('H2 O => H2O', 'adjacent species')
    call expect_invalid('H2 => H2O => O', 'more than one reaction arrow')
    call expect_invalid('H2 + O) => H2O', 'unmatched closing parenthesis')
    call expect_invalid('H2 + O => (H2O', 'unmatched opening parenthesis')
    call expect_invalid('ABCDEFGHIJKLMNOPQRSTU => H2O', 'longer than the supported')

  end subroutine

  subroutine test_removed_evolve_climate_setting()
    use photochem_settings, only: PhotoSettings
    type(PhotoSettings) :: settings
    character(:), allocatable :: err

    settings = PhotoSettings('../tests/test_settings1.yaml', err)
    if (.not. allocated(err)) then
      call fail('The obsolete evolve-climate setting was accepted')
    endif
    if (index(err, 'evolve-climate') == 0) then
      call fail('Unexpected error for obsolete evolve-climate setting: '//trim(err))
    endif
  end subroutine

  subroutine test_removed_conserving_initialization()
    use photochem_settings, only: PhotoSettings
    type(PhotoSettings) :: settings
    character(:), allocatable :: err

    settings = PhotoSettings('../tests/test_settings_conserving_init.yaml', err)
    if (.not. allocated(err)) then
      call fail('The obsolete conserving-initialization setting was accepted')
    endif
    if (index(err, 'conserving-initialization') == 0) then
      call fail('Unexpected error for obsolete conserving-initialization setting: '//trim(err))
    endif
  end subroutine

  subroutine test_removed_water_settings()
    use photochem_settings, only: PhotoSettings
    type(PhotoSettings) :: settings
    character(:), allocatable :: err

    settings = PhotoSettings('../tests/test_settings_water_compat.yaml', err)
    if (allocated(err)) then
      call fail('Disabled legacy water settings were rejected: '//trim(err))
    endif

    settings = PhotoSettings('../tests/test_settings_fix_water.yaml', err)
    if (.not. allocated(err)) then
      call fail('The obsolete fix-water-in-troposphere setting was accepted')
    endif
    if (index(err, 'fix-water-in-troposphere: true') == 0) then
      call fail('Unexpected error for fix-water-in-troposphere: '//trim(err))
    endif

    settings = PhotoSettings('../tests/test_settings_water_condensation.yaml', err)
    if (.not. allocated(err)) then
      call fail('The obsolete water-condensation setting was accepted')
    endif
    if (index(err, 'water-condensation: true') == 0) then
      call fail('Unexpected error for water-condensation: '//trim(err))
    endif
  end subroutine

  subroutine expect_valid(equation, expected_reverse, expected_reactants, expected_products)
    use photochem_input, only: parse_reaction
    character(len=*), intent(in) :: equation
    logical, intent(in) :: expected_reverse
    character(len=*), intent(in) :: expected_reactants(:), expected_products(:)
    character(len=20), allocatable :: reactants(:), products(:)
    character(:), allocatable :: err
    logical :: reverse
    integer :: i

    call parse_reaction(equation, reverse, reactants, products, err)
    if (allocated(err)) call fail('Valid reaction rejected: '//trim(equation)//' ('//trim(err)//')')
    if (reverse .neqv. expected_reverse) call fail('Wrong arrow direction for: '//trim(equation))
    if (size(reactants) /= size(expected_reactants)) then
      call fail('Wrong number of reactants for: '//trim(equation))
    endif
    if (size(products) /= size(expected_products)) then
      call fail('Wrong number of products for: '//trim(equation))
    endif
    do i = 1,size(reactants)
      if (trim(reactants(i)) /= trim(expected_reactants(i))) then
        call fail('Wrong reactant token for: '//trim(equation))
      endif
    enddo
    do i = 1,size(products)
      if (trim(products(i)) /= trim(expected_products(i))) then
        call fail('Wrong product token for: '//trim(equation))
      endif
    enddo
  end subroutine

  subroutine expect_invalid(equation, expected_message)
    use photochem_input, only: parse_reaction
    character(len=*), intent(in) :: equation, expected_message
    character(len=20), allocatable :: reactants(:), products(:)
    character(:), allocatable :: err
    logical :: reverse

    call parse_reaction(equation, reverse, reactants, products, err)
    if (.not. allocated(err)) then
      call fail('Malformed reaction was accepted: '//trim(equation))
    endif
    if (index(err, expected_message) == 0) then
      call fail('Unexpected error for '//trim(equation)//': '//trim(err))
    endif
  end subroutine

  subroutine fail(message)
    character(len=*), intent(in) :: message
    error stop trim(message)
  end subroutine

end program
