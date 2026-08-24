program test_input
  ! Focused input parsing and validation tests. See tests/README.md.
  use photochem_test_paths, only: test_file, data_file, data_dir
  implicit none
  call test_parsing()
  call test_data_construction()
  call test_vars_construction()
  call test_wrk_construction()
  call test_removed_evolve_climate_setting()
  call test_removed_conserving_initialization()
  call test_removed_water_settings()
  print *, 'test_input passed'
contains

  subroutine test_data_construction()
    use photochem_data, only: PhotochemData
    use photochem_settings, only: PhotoSettings
    type(PhotoSettings) :: settings
    type(PhotochemData) :: dat
    character(:), allocatable :: err

    settings = PhotoSettings(test_file('settings.yaml'), err)
    if (allocated(err)) call fail('Could not construct test settings: '//trim(err))

    dat = PhotochemData(data_file('reaction_mechanisms/zahnle_earth.yaml'), &
                        settings, data_dir, err)
    if (allocated(err)) call fail('Could not construct PhotochemData: '//trim(err))

    if (dat%nsl /= settings%nsl) call fail('PhotochemData has the wrong short-lived species count')
    if (.not. allocated(dat%SL_names)) call fail('PhotochemData did not store short-lived species names')
    if (any(dat%SL_names /= settings%SL_names)) call fail('PhotochemData has the wrong short-lived species names')
    if (dat%planet_mass /= settings%planet_mass) call fail('PhotochemData has the wrong planet mass')
    if (dat%planet_radius /= settings%planet_radius) call fail('PhotochemData has the wrong planet radius')
    if (dat%gas_rainout .neqv. settings%gas_rainout) call fail('PhotochemData has the wrong rainout setting')

    if (dat%nsp <= 0 .or. dat%nq <= 0 .or. dat%nrT <= 0) then
      call fail('PhotochemData mechanism dimensions were not initialized')
    endif
    if (.not. allocated(dat%species_names)) call fail('PhotochemData species were not initialized')
    if (.not. allocated(dat%reaction_equations)) call fail('PhotochemData reactions were not initialized')
    if (.not. allocated(dat%wavl)) call fail('PhotochemData wavelength grid was not initialized')
    if (.not. allocated(dat%photolysis_xs)) call fail('PhotochemData cross sections were not initialized')
    if (size(dat%wavl) /= dat%nw + 1) call fail('PhotochemData wavelength dimensions are inconsistent')
    if (size(dat%photolysis_xs, 1) /= dat%kj .or. &
        size(dat%photolysis_xs, 2) /= dat%nw) then
      call fail('PhotochemData photolysis dimensions are inconsistent')
    endif

    if (dat%kd /= 2*dat%nq + 1) call fail('PhotochemData has the wrong Jacobian diagonal width')
    if (dat%kl /= dat%kd + dat%nq) call fail('PhotochemData has the wrong lower bandwidth')
    if (dat%ku /= dat%kd - dat%nq) call fail('PhotochemData has the wrong upper bandwidth')
    if (dat%lda /= 3*dat%nq + 1) call fail('PhotochemData has the wrong Jacobian leading dimension')
  end subroutine

  subroutine test_vars_construction()
    use photochem_data, only: PhotochemData
    use photochem_settings, only: PhotoSettings
    use photochem_vars, only: PhotochemVars
    type(PhotoSettings) :: settings
    type(PhotochemData) :: dat
    type(PhotochemVars) :: var
    character(:), allocatable :: err

    settings = PhotoSettings(test_file('settings.yaml'), err)
    if (allocated(err)) call fail('Could not construct test settings: '//trim(err))

    dat = PhotochemData(data_file('reaction_mechanisms/zahnle_earth.yaml'), &
                        settings, data_dir, err)
    if (allocated(err)) call fail('Could not construct test data: '//trim(err))

    var = PhotochemVars(dat, settings, test_file('sun.txt'), err)
    if (allocated(err)) call fail('Could not construct PhotochemVars: '//trim(err))

    if (var%nz /= settings%nz) call fail('PhotochemVars has the wrong grid size')
    if (var%neqs /= dat%nq*var%nz) call fail('PhotochemVars has the wrong equation count')
    if (.not. allocated(var%lowerboundcond)) call fail('PhotochemVars boundary conditions were not initialized')
    if (size(var%lowerboundcond) /= dat%nq) call fail('PhotochemVars boundary dimensions are inconsistent')
    if (.not. allocated(var%photon_flux)) call fail('PhotochemVars photon flux was not initialized')
    if (size(var%photon_flux) /= dat%nw) call fail('PhotochemVars photon-flux dimensions are inconsistent')
    if (.not. allocated(var%temperature)) call fail('PhotochemVars temperature storage was not allocated')
    if (.not. allocated(var%z)) call fail('PhotochemVars altitude storage was not allocated')
    if (.not. allocated(var%xs_x_qy)) call fail('PhotochemVars cross-section storage was not allocated')
    if (size(var%temperature) /= var%nz .or. size(var%z) /= var%nz) then
      call fail('PhotochemVars atmospheric storage dimensions are inconsistent')
    endif
    if (size(var%xs_x_qy,1) /= var%nz .or. &
        size(var%xs_x_qy,2) /= dat%kj .or. &
        size(var%xs_x_qy,3) /= dat%nw) then
      call fail('PhotochemVars cross-section dimensions are inconsistent')
    endif
  end subroutine

  subroutine test_wrk_construction()
    use photochem_const, only: nsteps_save
    use photochem_wrk, only: PhotochemWrk
    type(PhotochemWrk) :: wrk
    integer, parameter :: nsp = 12, np = 2, nq = 8, nz = 20
    integer, parameter :: nrT = 15, kj = 3, nw = 25

    wrk = PhotochemWrk(nsp, np, nq, nz, nrT, kj, nw)

    if (size(wrk%t_history) /= nsteps_save) call fail('PhotochemWrk has the wrong history size')
    if (any(shape(wrk%mix_history) /= [nq,nz,nsteps_save])) then
      call fail('PhotochemWrk has inconsistent mixing-ratio history dimensions')
    endif
    if (any(shape(wrk%usol) /= [nq,nz])) call fail('PhotochemWrk has inconsistent solution dimensions')
    if (any(shape(wrk%densities) /= [nsp+1,nz])) call fail('PhotochemWrk has inconsistent density dimensions')
    if (any(shape(wrk%rx_rates) /= [nz,nrT])) call fail('PhotochemWrk has inconsistent reaction-rate dimensions')
    if (any(shape(wrk%prates) /= [nz,kj])) call fail('PhotochemWrk has inconsistent photolysis-rate dimensions')
    if (any(shape(wrk%amean_grd) /= [nz,nw])) call fail('PhotochemWrk has inconsistent radiation dimensions')
    if (any(shape(wrk%mix) /= [nq,nz])) call fail('PhotochemWrk has inconsistent mixing-ratio dimensions')
    if (size(wrk%pressure_hydro) /= nz .or. size(wrk%density_hydro) /= nz) then
      call fail('PhotochemWrk has inconsistent hydrostatic dimensions')
    endif
    if (wrk%robust_stepper_initialized) call fail('PhotochemWrk robust state was initialized prematurely')
    if (wrk%nsteps_total /= -1 .or. wrk%nerrors_total /= -1) then
      call fail('PhotochemWrk robust counters have the wrong initial values')
    endif
  end subroutine

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

    settings = PhotoSettings(test_file('test_settings1.yaml'), err)
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

    settings = PhotoSettings(test_file('test_settings_conserving_init.yaml'), err)
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

    settings = PhotoSettings(test_file('test_settings_water_compat.yaml'), err)
    if (allocated(err)) then
      call fail('Disabled legacy water settings were rejected: '//trim(err))
    endif

    settings = PhotoSettings(test_file('test_settings_fix_water.yaml'), err)
    if (.not. allocated(err)) then
      call fail('The obsolete fix-water-in-troposphere setting was accepted')
    endif
    if (index(err, 'fix-water-in-troposphere: true') == 0) then
      call fail('Unexpected error for fix-water-in-troposphere: '//trim(err))
    endif

    settings = PhotoSettings(test_file('test_settings_water_condensation.yaml'), err)
    if (.not. allocated(err)) then
      call fail('The obsolete water-condensation setting was accepted')
    endif
    if (index(err, 'water-condensation: true') == 0) then
      call fail('Unexpected error for water-condensation: '//trim(err))
    endif
  end subroutine

  subroutine expect_valid(equation, expected_reverse, expected_reactants, expected_products)
    use photochem_data, only: parse_reaction
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
    use photochem_data, only: parse_reaction
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
