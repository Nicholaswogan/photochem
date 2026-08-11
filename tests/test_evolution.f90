program test_evolution
  ! This longer-running end-to-end evolution exercise is intended for manual
  ! use and is compiled, but not run, in CI. See tests/README.md.
  implicit none
  call main_()
contains

  subroutine main_()
    use photochem, only: EvoAtmosphere, version, dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(:), allocatable :: err
    type(EvoAtmosphere) :: pc
    logical :: converged

    ! Print version
    print*,'photochem version == ',trim(version)

    ! Initialize code
    pc = EvoAtmosphere(&
                      "../data/reaction_mechanisms/zahnle_earth.yaml", &
                      "../examples/ModernEarth/settings.yaml", &
                      "../examples/ModernEarth/Sun_now.txt", &
                      "../examples/ModernEarth/atmosphere.txt", &
                      "../data", &
                      err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    
    converged = pc%find_steady_state(err)
    if (allocated(err)) then
      print*,trim(err)
      stop 1
    endif
    if (.not.converged) then
      print*,'Modern Earth integration did not converge'
      stop 1
    endif
    if (pc%wrk%nsteps_total <= 0) then
      print*,'Modern Earth integration converged without an accepted step'
      stop 1
    endif
    if (.not.ieee_is_finite(pc%wrk%tn) .or. pc%wrk%tn <= 0.0_dp .or. &
        .not.ieee_is_finite(pc%wrk%longdy) .or. &
        .not.ieee_is_finite(pc%wrk%longdydt) .or. &
        .not.all(ieee_is_finite(pc%wrk%usol))) then
      print*,'Modern Earth integration produced an invalid final state'
      stop 1
    endif

    print*,'Modern Earth converged in accepted steps: ',pc%wrk%nsteps_total
    print*,'Final integration time: ',pc%wrk%tn
    print*,'Final longdy and longdydt: ',pc%wrk%longdy,pc%wrk%longdydt

  end subroutine

end program
