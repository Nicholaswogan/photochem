module photochem_evoatmosphere
  use photochem_const, only: dp
  use photochem_types, only : PhotochemData, PhotochemVars, PhotochemWrkEvo
  implicit none

  private
  public :: EvoAtmosphere

  type :: EvoAtmosphere
    type(PhotochemData), allocatable :: dat
    type(PhotochemVars), allocatable :: var
    type(PhotochemWrkEvo), allocatable :: wrk

  contains
  !!! photochem_evoatmosphere_init.f90 !!!
    procedure :: init => EvoAtmosphere_init

    !!! photochem_evoatmosphere_rhs.f90 !!!
    procedure :: prep_atmosphere => prep_all_evo_gas

  end type

  interface
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_evoatmosphere_init.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module subroutine EvoAtmosphere_init(self, data_dir, mechanism_file, &
                                         settings_file, flux_file, atmosphere_txt, err)
      class(EvoAtmosphere), intent(inout) :: self
      character(len=*), intent(in) :: data_dir
      character(len=*), intent(in) :: mechanism_file
      character(len=*), intent(in) :: settings_file
      character(len=*), intent(in) :: flux_file
      character(len=*), intent(in) :: atmosphere_txt
      character(:), allocatable, intent(out) :: err
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_atmosphere_rhs.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module subroutine prep_all_evo_gas(self, dsol_in, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: dsol_in(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine

  end interface

end module