module photochem_atmosphere
  use photochem_const, only: err_len, real_kind
  use photochem_types, only : PhotochemData, PhotochemVars, PhotochemWrk
  implicit none
  
  ! The Atmosphere type. It contains all data describing the atmosphere
  ! such as temperature profile, photon flux, chemical species, 
  ! particles, and chemical reactions. It also contains several
  ! procedures which compute atmospheric photochemistry.
    
  private
  public :: Atmosphere
  
  type :: Atmosphere
    type(PhotochemData), allocatable :: dat
    type(PhotochemVars), allocatable :: var
    type(PhotochemWrk), allocatable :: wrk
  contains
    ! public
    procedure :: init => Atmosphere_init
    procedure :: right_hand_side_chem
    procedure :: production_and_loss
    procedure :: right_hand_side => rhs_background_gas
    procedure :: jacobian => jac_background_gas
    procedure :: photochemical_equilibrium
    procedure :: initialize_stepper
    procedure :: step
    procedure :: destroy_stepper
    procedure :: out2atmosphere_txt
    procedure :: out2in
    procedure :: gas_fluxes
    procedure :: redox_conservation
    procedure :: set_lower_bc
    procedure :: set_upper_bc
    
    ! private
    procedure, private :: prep_atmosphere => prep_all_background_gas
  end type
  
  
  interface
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_atmosphere_init.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_atmosphere_rhs.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module subroutine prep_all_background_gas(self, usol_in, err)
      class(Atmosphere), target, intent(inout) :: self
      real(real_kind), intent(in) :: usol_in(:,:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine right_hand_side_chem(self, usol, rhs, err)
      class(Atmosphere), target, intent(inout) :: self
      real(real_kind), intent(in) :: usol(:,:)
      real(real_kind), intent(out) :: rhs(:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine production_and_loss(self, species, usol, pl, err)     
      use photochem_types, only: ProductionLoss
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: species
      real(real_kind), intent(in) :: usol(:,:)
      class(ProductionLoss), intent(out) :: pl
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine rhs_background_gas(self, neqs, usol_flat, rhs, err)
      class(Atmosphere), target, intent(inout) :: self
      integer, intent(in) :: neqs
      real(real_kind), target, intent(in) :: usol_flat(neqs)
      real(real_kind), intent(out) :: rhs(neqs)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine jac_background_gas(self, lda_neqs, neqs, usol_flat, jac, err)
      class(Atmosphere), target, intent(inout) :: self
      integer, intent(in) :: lda_neqs, neqs
      real(real_kind), target, intent(in) :: usol_flat(neqs)
      real(real_kind), intent(out), target :: jac(lda_neqs)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module function fcn_fH2O(ptr, n, x, fvec, iflag) result(res) bind(c)
      use, intrinsic :: iso_c_binding, only : c_ptr  
      type(c_ptr) :: ptr
      integer, value :: n, iflag
      real(real_kind), intent(in) :: x(n)
      real(real_kind), intent(out) :: fvec(n)
      integer :: res
    end function  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_atmosphere_integrate.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module subroutine photochemical_equilibrium(self, success, err)
      class(Atmosphere), target, intent(inout) :: self
      logical, intent(out) :: success
      character(len=err_len), intent(out) :: err 
    end subroutine
    
    module subroutine initialize_stepper(self, usol_start, err)      
      class(Atmosphere), target, intent(inout) :: self
      real(real_kind), intent(in) :: usol_start(:,:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module function step(self, err) result(tn)
      class(Atmosphere), target, intent(inout) :: self
      character(len=err_len), intent(out) :: err
      real(real_kind) :: tn
    end function
    
    module subroutine destroy_stepper(self, err)
      class(Atmosphere), target, intent(inout) :: self
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
                          result(ierr) bind(c, name='RhsFn')
      use, intrinsic :: iso_c_binding
      use fcvode_mod
      use fsundials_nvector_mod
      real(c_double), value :: tn
      type(N_Vector)        :: sunvec_y
      type(N_Vector)        :: sunvec_f
      type(c_ptr), value    :: user_data
      integer(c_int)        :: ierr
    end function
    
    module function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                          tmp1, tmp2, tmp3) &
                          result(ierr) bind(C,name='JacFn')
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      use fnvector_serial_mod
      use fsunmatrix_band_mod
      use fsundials_matrix_mod
      real(c_double), value :: tn
      type(N_Vector)        :: sunvec_y 
      type(N_Vector)        :: sunvec_f
      type(SUNMatrix)       :: sunmat_J 
      type(c_ptr), value    :: user_data 
      type(N_Vector)        :: tmp1, tmp2, tmp3
      integer(c_int)        :: ierr
    end function
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_atmosphere_utils.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    
    module subroutine gas_fluxes(self, surf_fluxes, top_fluxes, err)
      class(Atmosphere), target, intent(inout) :: self
      real(real_kind), intent(out) :: surf_fluxes(:)
      real(real_kind), intent(out) :: top_fluxes(:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module function redox_conservation(self, err) result(redox_factor)
      class(Atmosphere), target, intent(inout) :: self
      character(len=err_len), intent(out) :: err
      real(real_kind) :: redox_factor
    end function
    
    module subroutine set_lower_bc(self, species, bc_type, vdep, mix, flux, height, err)
      class(Atmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species
      character(len=*), intent(in) :: bc_type
      real(real_kind), optional, intent(in) :: vdep
      real(real_kind), optional, intent(in) :: mix
      real(real_kind), optional, intent(in) :: flux
      real(real_kind), optional, intent(in) :: height
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
      class(Atmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species
      character(len=*), intent(in) :: bc_type
      real(real_kind), optional, intent(in) :: veff
      real(real_kind), optional, intent(in) :: flux
      character(len=err_len), intent(out) :: err
    end subroutine
    
  end interface
  
end module