module photochem_atmosphere
  use photochem_const, only: err_len, real_kind
  use photochem_types, only : PhotochemData, PhotochemVars, PhotochemWrk
  implicit none
  
  ! Contains the Atmosphere type. This type computes everything
    
  private
  public :: Atmosphere
  
  type :: Atmosphere
    type(PhotochemData), allocatable :: dat
    type(PhotochemVars), allocatable :: var
    type(PhotochemWrk), allocatable :: wrk
  contains
    procedure :: init => Atmosphere_init
    procedure, private :: prep_atmosphere => prep_all_background_gas
    procedure, private :: dochem_implicit => dochem_implicit
    procedure :: right_hand_side => rhs_background_gas
    procedure :: jacobian => jac_background_gas
    procedure :: photochemical_equilibrium
    procedure :: out2atmosphere_txt
    procedure :: out2in
    procedure :: surface_fluxes
  end type
  
  
  interface
    ! photochem_atmosphere_init.f90
    module subroutine Atmosphere_init(self, data_dir, mechanism_file, &
                                     settings_file, flux_file, atmosphere_txt, err)
      class(Atmosphere), intent(inout) :: self
      character(len=*), intent(in) :: data_dir
      character(len=*), intent(in) :: mechanism_file
      character(len=*), intent(in) :: settings_file
      character(len=*), intent(in) :: flux_file
      character(len=*), intent(in) :: atmosphere_txt
      character(len=err_len), intent(out) :: err
    end subroutine
    
    ! photochem_atmosphere_rhs.f90
    module subroutine prep_all_background_gas(self, usol_in, err)
      class(Atmosphere), target, intent(inout) :: self
      real(real_kind), intent(in) :: usol_in(:,:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    ! photochem_atmosphere_rhs.f90
    module subroutine dochem_implicit(self, usol, rhs, err)
      class(Atmosphere), target, intent(inout) :: self
      real(real_kind), intent(in) :: usol(:,:)
      real(real_kind), intent(out) :: rhs(:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    ! photochem_atmosphere_rhs.f90
    module subroutine rhs_background_gas(self, neqs, usol_flat, rhs, err)
      class(Atmosphere), target, intent(inout) :: self
      integer, intent(in) :: neqs
      real(real_kind), target, intent(in) :: usol_flat(neqs)
      real(real_kind), intent(out) :: rhs(neqs)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    ! photochem_atmosphere_rhs.f90
    module subroutine jac_background_gas(self, lda_neqs, neqs, usol_flat, jac, err)
      class(Atmosphere), target, intent(inout) :: self
      integer, intent(in) :: lda_neqs, neqs
      real(real_kind), target, intent(in) :: usol_flat(neqs)
      real(real_kind), intent(out), target :: jac(lda_neqs)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    ! photochem_atmosphere_integrate.f90
    module subroutine photochemical_equilibrium(self, success, err)
      class(Atmosphere), target, intent(inout) :: self
      logical, intent(out) :: success
      character(len=err_len), intent(out) :: err 
    end subroutine
    
    ! utilities
    module subroutine out2atmosphere_txt(self, filename, overwrite, clip, err)
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename
      logical, intent(in) :: overwrite, clip
      character(len=1024), intent(out) :: err
    end subroutine
    
    module subroutine out2in(self, err)
      class(Atmosphere), intent(inout) :: self
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine surface_fluxes(self, fluxes, err)
      class(Atmosphere), target, intent(inout) :: self
      real(real_kind), intent(out) :: fluxes(:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
  end interface
  
end module