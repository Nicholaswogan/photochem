module photochem_object
  use photochem_const, only: err_len, real_kind
  use photochem_types, only : PhotochemData, PhotochemVars, PhotochemWrk
  implicit none
  
  private
  
  public :: Photochem
  
  type :: Photochem
    type(PhotochemData) :: dat
    type(PhotochemVars) :: var
    type(PhotochemWrk) :: wrk
  contains
    procedure :: init => Photochem_init
    procedure :: rhs => rhs_background_gas
  end type
  
  
  interface 
    module subroutine Photochem_init(self, data_dir, mechanism_file, &
                                     settings_file, flux_file, atmosphere_txt, err)
      class(Photochem), intent(inout) :: self
      character(len=*), intent(in) :: data_dir
      character(len=*), intent(in) :: mechanism_file
      character(len=*), intent(in) :: settings_file
      character(len=*), intent(in) :: flux_file
      character(len=*), intent(in) :: atmosphere_txt
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine rhs_background_gas(self, neqs, usol_flat, rhs, err)
      class(Photochem), target, intent(inout) :: self
      integer, intent(in) :: neqs
      real(real_kind), target, intent(in) :: usol_flat(neqs)
      real(real_kind), intent(out) :: rhs(neqs)
      character(len=err_len), intent(out) :: err
    end subroutine
    
  end interface
  
end module