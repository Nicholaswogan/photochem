module photochem_atmosphere
  use photochem_const, only: err_len, dp
  use photochem_types, only : PhotochemData, PhotochemVars, PhotochemWrk
  implicit none
  
  ! The Atmosphere type. It contains all data describing the atmosphere
  ! such as temperature profile, photon flux, chemical species, 
  ! particles, and chemical reactions. It also contains
  ! procedures which compute atmospheric photochemistry.
    
  private
  public :: Atmosphere
  
  type :: Atmosphere
    type(PhotochemData), allocatable :: dat
    type(PhotochemVars), allocatable :: var
    type(PhotochemWrk), allocatable :: wrk
  contains
    !!! photochem_atmosphere_init.f90 !!!
    procedure :: init => Atmosphere_init
    
    !!! photochem_atmosphere_rhs.f90 !!!
    procedure :: prep_atmosphere => prep_all_background_gas
    procedure :: right_hand_side_chem
    procedure :: production_and_loss
    procedure :: right_hand_side => rhs_background_gas
    procedure :: jacobian => jac_background_gas
    
    !!! photochem_atmosphere_integrate.f90 !!!
    procedure :: evolve
    procedure :: photochemical_equilibrium
    procedure :: initialize_stepper
    procedure :: step
    procedure :: destroy_stepper
    
    !!! photochem_atmosphere_utils.f90 !!!
    procedure :: out2atmosphere_txt
    procedure :: out2in
    procedure :: gas_fluxes
    procedure :: atom_conservation
    procedure :: redox_conservation
    procedure :: set_lower_bc
    procedure :: set_upper_bc
    procedure :: set_temperature
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
      ! initializes the Atmosphere object. Reads all files, 
      ! allocates memry, and prepares for photochemical calculations.
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_atmosphere_rhs.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module subroutine prep_all_background_gas(self, usol_in, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:)
      character(len=err_len), intent(out) :: err
      ! Given usol_in, the mixing ratios of each species in the atmosphere,
      ! this subroutine calculates reaction rates, photolysis rates, etc.
      ! and puts this information into self%wrk. self%wrk contains all the
      ! information needed for `dochem` to compute chemistry
    end subroutine
    
    module subroutine right_hand_side_chem(self, usol, rhs, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol(:,:)
      real(dp), intent(out) :: rhs(:)
      character(len=err_len), intent(out) :: err
      ! The right-hand-side of the ODEs describing atmospheric chemistry
      ! but does not include transport
    end subroutine
    
    module subroutine production_and_loss(self, species, usol, pl, err)     
      use photochem_types, only: ProductionLoss
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: species
      real(dp), intent(in) :: usol(:,:)
      type(ProductionLoss), intent(out) :: pl
      character(len=err_len), intent(out) :: err
      ! Computes the production and loss of input `species`.
      ! See ProductionLoss object in photochem_types.
    end subroutine
    
    module subroutine rhs_background_gas(self, neqs, usol_flat, rhs, err)
      class(Atmosphere), target, intent(inout) :: self
      integer, intent(in) :: neqs
      real(dp), target, intent(in) :: usol_flat(neqs)
      real(dp), intent(out) :: rhs(neqs)
      character(len=err_len), intent(out) :: err
      ! Computes the right-hand-side of the ODEs describing atmospheric chemistry
      ! and transport.
    end subroutine
    
    module subroutine jac_background_gas(self, lda_neqs, neqs, usol_flat, jac, err)
      class(Atmosphere), target, intent(inout) :: self
      integer, intent(in) :: lda_neqs, neqs
      real(dp), target, intent(in) :: usol_flat(neqs)
      real(dp), intent(out), target :: jac(lda_neqs)
      character(len=err_len), intent(out) :: err
      ! The jacobian of the rhs_background_gas.
    end subroutine
    
    module function fcn_fH2O(ptr, n, x, fvec, iflag) result(res) bind(c)
      use, intrinsic :: iso_c_binding, only : c_ptr  
      type(c_ptr) :: ptr
      integer, value :: n, iflag
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: fvec(n)
      integer :: res
    end function  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_atmosphere_integrate.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module subroutine evolve(self, filename, tstart, usol_start, t_eval, success, err)
      use, intrinsic :: iso_c_binding
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename
      real(c_double), intent(in) :: tstart
      real(dp), intent(in) :: usol_start(:,:)
      real(c_double), intent(in) :: t_eval(:)
      logical, intent(out) :: success
      character(len=err_len), intent(out) :: err
      ! Evolve atmosphere through time, and saves output in a binary Fortran file.
    end subroutine
    
    module subroutine photochemical_equilibrium(self, success, err)
      class(Atmosphere), target, intent(inout) :: self
      logical, intent(out) :: success
      character(len=err_len), intent(out) :: err 
      ! Integrates to photochemical equilibrium starting from self%var%usol_init
    end subroutine
    
    module subroutine initialize_stepper(self, usol_start, err)      
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_start(:,:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module function step(self, err) result(tn)
      class(Atmosphere), target, intent(inout) :: self
      character(len=err_len), intent(out) :: err
      real(dp) :: tn
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
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_atmosphere_utils.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      real(dp), intent(out) :: surf_fluxes(:)
      real(dp), intent(out) :: top_fluxes(:)
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module function atom_conservation(self, atom, err) result(con)
      use photochem_types, only: AtomConservation
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: atom
      character(len=err_len), intent(out) :: err
      type(AtomConservation) :: con
    end function
    
    module function redox_conservation(self, err) result(redox_factor)
      class(Atmosphere), target, intent(inout) :: self
      character(len=err_len), intent(out) :: err
      real(dp) :: redox_factor
    end function
    
    module subroutine set_lower_bc(self, species, bc_type, vdep, mix, flux, height, err)
      class(Atmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species
      character(len=*), intent(in) :: bc_type
      real(dp), optional, intent(in) :: vdep
      real(dp), optional, intent(in) :: mix
      real(dp), optional, intent(in) :: flux
      real(dp), optional, intent(in) :: height
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
      class(Atmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species
      character(len=*), intent(in) :: bc_type
      real(dp), optional, intent(in) :: veff
      real(dp), optional, intent(in) :: flux
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine set_temperature(self, temperature, trop_alt, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: temperature(:)
      real(dp), optional, intent(in) :: trop_alt
      character(len=err_len), intent(out) :: err
    end subroutine
    
  end interface
  
end module