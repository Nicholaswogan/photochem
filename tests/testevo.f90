program main
  implicit none
  call main_()
contains

  subroutine main_()
    use photochem, only: EvoAtmosphere, version, dp
    character(:), allocatable :: err
    type(EvoAtmosphere) :: pc
    integer :: ind
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

  end subroutine

end program