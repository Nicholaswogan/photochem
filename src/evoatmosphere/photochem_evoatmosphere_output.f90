submodule(photochem_evoatmosphere) photochem_evoatmosphere_output
  implicit none

contains

  module subroutine out2atmosphere_txt(self, filename, number_of_decimals, overwrite, clip, err)
    use photochem_common, only: out2atmosphere_txt_base
    class(EvoAtmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: number_of_decimals
    logical, intent(in) :: overwrite
    logical, intent(in) :: clip
    character(:), allocatable, intent(out) :: err

    real(dp) :: rhs(self%var%neqs)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    call self%require_atmosphere_initialized('out2atmosphere_txt', err)
    if (allocated(err)) return

    dat => self%dat
    var => self%var
    wrk => self%wrk

    ! Update work variables needed for output.
    call self%right_hand_side_chem(wrk%usol, rhs, err)
    if (allocated(err)) return

    call out2atmosphere_txt_base(dat, var, &
                                 wrk%pressure, wrk%density, wrk%densities, &
                                 wrk%molecules_per_particle, filename, &
                                 number_of_decimals, overwrite, clip, err)
    if (allocated(err)) return

  end subroutine

end submodule
