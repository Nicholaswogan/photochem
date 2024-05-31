module photochem_atmosphere
  use photochem_const, only: dp
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
    
    ! photochem_atmosphere_rhs.f90
    procedure, private :: prep_atm_background_gas
    procedure :: prep_atmosphere => prep_all_background_gas
    procedure :: right_hand_side_chem
    procedure :: production_and_loss
    procedure :: right_hand_side => rhs_background_gas
    procedure :: jacobian => jac_background_gas
    
    ! photochem_atmosphere_integrate.f90
    procedure :: evolve
    procedure :: check_for_convergence
    procedure :: photochemical_equilibrium
    procedure :: initialize_stepper
    procedure :: step
    procedure :: destroy_stepper
    
    ! photochem_atmosphere_utils.f90
    procedure :: out2atmosphere_txt
    procedure :: out2in
    procedure :: gas_fluxes
    procedure :: atom_conservation
    procedure :: redox_conservation
    procedure :: set_lower_bc
    procedure :: set_upper_bc
    procedure :: set_temperature
    procedure :: set_press_temp_edd
    procedure :: set_rate_fcn
    procedure :: update_vertical_grid
  end type
  interface Atmosphere
    module procedure :: create_Atmosphere
  end interface
  
  
  interface

    !~~ photochem_atmosphere_init.f90 ~~!

    !> Initializes the Atmosphere object by reading input files.
    module function create_Atmosphere(mechanism_file, settings_file, flux_file, atmosphere_txt, data_dir, err) result(self)
      character(len=*), intent(in) :: mechanism_file !! Path to reaction mechanism file.
      character(len=*), intent(in) :: settings_file !! Path to setting file
      character(len=*), intent(in) :: flux_file !! Path to file that specifies stellar flux
      !> Path to file that specifies temperature profile and initial conditions
      character(len=*), intent(in) :: atmosphere_txt
      character(len=*), intent(in) :: data_dir !! Directory data is contained in.
      character(:), allocatable, intent(out) :: err
      type(Atmosphere) :: self
    end function

    !~~ photochem_atmosphere_rhs.f90 ~~!

    module subroutine prep_atm_background_gas(self, usol_in, usol, molecules_per_particle)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:)
      real(dp), intent(out) :: usol(:,:)
      real(dp), intent(out) :: molecules_per_particle(:,:)
    end subroutine

    !> Given usol_in, the mixing ratios of each species in the atmosphere,
    !> this subroutine calculates reaction rates, photolysis rates, etc.
    !> and puts this information into self%wrk. self%wrk contains all the
    !> information needed for `dochem` to compute chemistry.
    module subroutine prep_all_background_gas(self, usol_in, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:) !! Input mixing ratios (nq,nz)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> The right-hand-side of the ODEs describing atmospheric chemistry
    !> but does not include vertical transport terms.
    module subroutine right_hand_side_chem(self, usol, rhs, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol(:,:) !! Input mixing ratios (nq,nz)
      real(dp), intent(out) :: rhs(:) !! rate of change of each species in 
                                      !! mixing ratio per second (neqs)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Computes the production and loss of input `species`.
    !> See ProductionLoss object in photochem_types.f90.
    module subroutine production_and_loss(self, species, usol, pl, err)     
      use photochem_types, only: ProductionLoss
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: species !! name of species
      real(dp), intent(in) :: usol(:,:) !! Mixing ratios (nq,nz)
      type(ProductionLoss), intent(out) :: pl !! Type describing production and loss of species
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Computes the right-hand-side of the ODEs describing atmospheric chemistry
    !> and transport.
    module subroutine rhs_background_gas(self, neqs, tn, usol_flat, rhs, err)
      class(Atmosphere), target, intent(inout) :: self
      integer, intent(in) :: neqs !! number of equations
      real(dp), intent(in) :: tn !! time in seconds
      real(dp), target, intent(in) :: usol_flat(neqs) !! mixing ratios, flattened to 1D array.
      real(dp), intent(out) :: rhs(neqs) !! rate of change of each species in 
                                         !! mixing ratio per second (neqs)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Computes the banded jacobian approximation of the system of ODEs.
    module subroutine jac_background_gas(self, lda_neqs, neqs, tn, usol_flat, jac, err)
      class(Atmosphere), target, intent(inout) :: self
      integer, intent(in) :: lda_neqs !! total length of jacobian
      integer, intent(in) :: neqs !! number of equations
      real(dp), intent(in) :: tn !! time in seconds
      real(dp), target, intent(in) :: usol_flat(neqs) !! mixing ratios, flattened to 1D array.
      !> Jacobian of ODEs, with length lda*neqs. We assume the Jacobian is banded, with bandwith
      !> 3*nq + 1 = lda.
      real(dp), intent(out), target :: jac(lda_neqs)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !~~ photochem_atmosphere_integrate.f90 ~~!

    !> Evolve atmosphere through time, and saves output in a binary Fortran file.
    module function evolve(self, filename, tstart, usol_start, t_eval, overwrite, err) result(success)
      use, intrinsic :: iso_c_binding
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename !! Filename to save results.
      real(c_double), intent(in) :: tstart !! start time in seconds
      real(dp), intent(in) :: usol_start(:,:) !! Initial mixing ratios
      real(c_double), intent(in) :: t_eval(:) !! times to evaluate the solution
      logical, optional, intent(in) :: overwrite !! If true, then overwrites pre-existing files with `filename`
      logical :: success !! If True, then integration was successful.
      character(:), allocatable, intent(out) :: err
    end function

    !> Determines if integration has converged to photochemical steady-state.
    module function check_for_convergence(self, err) result(converged)
      class(Atmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      logical :: converged
    end function
    
    !> Integrates to photochemical equilibrium starting from self%var%usol_init
    module subroutine photochemical_equilibrium(self, success, err)
      class(Atmosphere), target, intent(inout) :: self
      logical, intent(out) :: success
      character(:), allocatable, intent(out) :: err 
    end subroutine
    
    !> Initializes an integration starting at `usol_start` mixing ratios
    module subroutine initialize_stepper(self, usol_start, err)      
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_start(:,:) !! Initial mixing ratios
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Takes one internal integration step. Function `initialize_stepper`
    !> must have been called befe this
    module function step(self, err) result(tn)
      class(Atmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      real(dp) :: tn
    end function
    
    !> Deallocates memory created during `initialize_stepper`
    module subroutine destroy_stepper(self, err)
      class(Atmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
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
    
    !~~ photochem_atmosphere_utils.f90 ~~!

    !> Saves state of the atmosphere, using the mixing ratios in self%wrk%usol.
    module subroutine out2atmosphere_txt(self, filename, number_of_decimals, overwrite, clip, err)
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename !! Output filename
      integer, intent(in) :: number_of_decimals
      logical, intent(in) :: overwrite !! If true, then output file can be overwritten
      logical, intent(in) :: clip !! If true, then mixing ratios are 
                                  !! clipped at a very small positive number
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Copies self%var%usol_out to self%var%usol_init
    module subroutine out2in(self, err)
      class(Atmosphere), intent(inout) :: self
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Computes gas fluxes at model boundaries in order to maintain
    !> current atmospheric concentrations. Uses the mixing ratios stored in
    !> self%wrk%usol.
    module subroutine gas_fluxes(self, surf_fluxes, top_fluxes, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(out) :: surf_fluxes(:) !! surface fluxes (molecules/cm^2/s)
      real(dp), intent(out) :: top_fluxes(:) !! top-of-atmosphere fluxes (molecules/cm^2/s)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Computes atom conservation. This is useful for determining if a model
    !> is in steady-state or not. Uses the mixing ratios in self%wrk%usol.
    module function atom_conservation(self, atom, err) result(con)
      use photochem_types, only: AtomConservation
      class(Atmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: atom !! Name of atom
      character(:), allocatable, intent(out) :: err
      type(AtomConservation) :: con !! Type containing conservation information
    end function
    
    !> Computes redox conservation. This is useful for determining if a model
    !> is in steady-state or not. Uses the mixing ratios in self%wrk%usol.
    module function redox_conservation(self, err) result(redox_factor)
      class(Atmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      real(dp) :: redox_factor !! should be small (< 1e-5) at equilibrium.
    end function
    
    !> Sets a lower boundary condition.
    module subroutine set_lower_bc(self, species, bc_type, vdep, mix, flux, height, err)
      class(Atmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species !! Species to set boundary condition
      character(len=*), intent(in) :: bc_type !! Boundary condition type
      real(dp), optional, intent(in) :: vdep !! Deposition velocity (cm/s)
      real(dp), optional, intent(in) :: mix !! Mixing ratio
      real(dp), optional, intent(in) :: flux !! Flux (molecules/cm^2/s)
      real(dp), optional, intent(in) :: height !! Height in atmosphere (km)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Sets upper boundary condition
    module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
      class(Atmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species !! Species to set boundary condition
      character(len=*), intent(in) :: bc_type !! Boundary condition type
      real(dp), optional, intent(in) :: veff !! effusion velocity (cm/s)
      real(dp), optional, intent(in) :: flux !! Flux (molecules/cm^2/s)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Changes the temperature profile.
    module subroutine set_temperature(self, temperature, trop_alt, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: temperature(:) !! new temperature at each atomspheric layer
      real(dp), optional, intent(in) :: trop_alt !! Tropopause altitude (cm). Only necessary if
                                                 !! rainout == True, or fix_water_in_trop == True.
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Given an input P, T, and edd, the code will find the temperature and eddy diffusion profile
    !> on the current altitude-grid that matches the inputs.
    module subroutine set_press_temp_edd(self, P, T, edd, trop_p, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: P(:) !! Pressure (dynes/cm^2)
      real(dp), intent(in) :: T(:) !! Temperature (K)
      real(dp), intent(in) :: edd(:) !! Eddy diffusion (cm^2/s)
      real(dp), optional, intent(in) :: trop_p !! Tropopause pressure (dynes/cm^2)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Sets a function describing a custom rate for a species.
    !> This could be useful for modeling external processes not in the
    !> model.
    module subroutine set_rate_fcn(self, species, fcn, err)
      use photochem_types, only: time_dependent_rate_fcn
      class(Atmosphere), target, intent(inout) :: self
      character(*), intent(in) :: species !! Species name
      procedure(time_dependent_rate_fcn), pointer :: fcn
      character(:), allocatable, intent(inout) :: err
    end subroutine

    !> Re-does the vertical grid so that the pressure at the top of the
    !> atmosphere is a `TOA_alt` or `TOA_pressure`. If the TOA needs to be raised above the current
    !> TOA, then the function constantly extrapolates mixing ratios, temperature,
    !> eddy diffusion, and particle radii.
    module subroutine update_vertical_grid(self, TOA_alt, TOA_pressure, err)
      class(Atmosphere), target, intent(inout) :: self
      real(dp), optional, intent(in) :: TOA_alt !! New top of atmosphere altitude (cm)
      real(dp), optional, intent(in) :: TOA_pressure !! New top of atmosphere pressure (dynes/cm^2)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
  end interface
  
end module