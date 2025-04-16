module photochem_evoatmosphere
  use photochem_const, only: dp
  use photochem_types, only : PhotochemData, PhotochemVars, PhotochemWrkEvo
  use clima_radtran, only: Radtran
  implicit none

  private
  public :: EvoAtmosphere
  public :: temp_dependent_albedo_fcn

  abstract interface
    !> A temperature dependent surface albedo
    function temp_dependent_albedo_fcn(T_surf) result(albedo)
      use iso_c_binding, only: c_double
      real(c_double), value, intent(in) :: T_surf !! K
      real(c_double) :: albedo
    end function
  end interface

  type :: EvoAtmosphere
    type(PhotochemData), allocatable :: dat
    type(PhotochemVars), allocatable :: var
    type(PhotochemWrkEvo), allocatable :: wrk

    logical :: evolve_climate
    ! Below are only relevant for evolve_climate = .true.
    !> The surface temperature (K). Only relevent when doing time-dependent
    !> photochemical-climate simulation.
    real(dp) :: T_surf
    !> Assumed tropopause temperature for climate calculations (K).
    real(dp) :: T_trop = 200.0_dp
    !> Callback function to set a temperature dependent albedo (e.g. ice-albedo feedback).
    procedure(temp_dependent_albedo_fcn), nopass, pointer :: albedo_fcn => null()
    type(Radtran), allocatable :: rad
    ! Above are only relevant for evolve_climate = .true.

    !> When running the `evolve` routine, this determines
    !> the minimum pressure of the top of the model domain (bars). 
    !> If the pressure gets smaller than this value, then the integration will stop
    !> and re-grid the model domain before continuing integration, so that the 
    !> top of the atmosphere has a bigger pressure than `P_top_min`.
    real(dp) :: P_top_min = 1.0e-50_dp
    !> When running the `evolve` routine, this determines
    !> the maximum pressure of the top of the model domain (bars). 
    !> If the pressure gets larger than this value, then the integration will stop
    !> and re-grid the model domain before continuing integration, so that the 
    !> top of the atmosphere has a smaller pressure than `P_top_max`.
    real(dp) :: P_top_max = 1.0e50_dp
    !> Sets the fractional amount that the top of the model domain changes
    !> when integration is haulted by `P_top_min` or `P_top_max`
    real(dp) :: top_atmos_adjust_frac = 0.02

  contains

    !~~ photochem_evoatmosphere_rhs.f90 ~~!
    procedure :: set_trop_ind
    procedure :: prep_atm_evo_gas
    procedure :: prep_atmosphere => prep_all_evo_gas
    procedure :: right_hand_side_chem
    procedure :: production_and_loss
    procedure :: right_hand_side => rhs_evo_gas
    procedure :: jacobian => jac_evo_gas

    !~~ photochem_evoatmosphere_integrate.f90 ~~!
    procedure :: evolve
    procedure :: check_for_convergence
    procedure :: initialize_stepper
    procedure :: step
    procedure :: destroy_stepper
    procedure :: initialize_robust_stepper
    procedure :: robust_step
    procedure :: find_steady_state

    !~~ photochem_evoatmosphere_utils.f90 ~~!
    procedure :: out2atmosphere_txt
    procedure :: gas_fluxes
    procedure :: set_lower_bc
    procedure :: set_upper_bc
    procedure :: set_rate_fcn
    procedure :: set_temperature
    procedure :: set_press_temp_edd
    procedure :: update_vertical_grid
    procedure :: rebin_update_vertical_grid
    procedure :: regrid_prep_atmosphere

  end type
  interface EvoAtmosphere
    module procedure :: create_EvoAtmosphere
  end interface

  interface

    !~~ photochem_evoatmosphere_init.f90 ~~!

    !> Initializes the photochemical model.
    module function create_EvoAtmosphere(mechanism_file, &
                                         settings_file, flux_file, atmosphere_txt, data_dir, err) result(self)
      character(len=*), intent(in) :: mechanism_file !! Path to the reaction mechanism file (yaml format).
      character(len=*), intent(in) :: settings_file !! Path to the settings file (yaml format).
      character(len=*), intent(in) :: flux_file !! Path to the file describing the stellar flux.
      !> Path to the file containing altitude, total number density, temperature, 
      !> eddy diffusion, initial concentrations of each gas (mixing ratios), 
      !> and particle radii.
      character(len=*), intent(in) :: atmosphere_txt
      !> Path to the data directory containing photolysis cross sections and other data
      !> needed to run the model
      character(len=*), intent(in) :: data_dir
      character(:), allocatable, intent(out) :: err
      type(EvoAtmosphere) :: self
    end function

    !~~ photochem_atmosphere_rhs.f90 ~~!

    module subroutine set_trop_ind(self, usol_in, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine prep_atm_evo_gas(self, usol_in, usol, &
                                      molecules_per_particle, pressure, density, mix, mubar, &
                                      pressure_hydro, density_hydro, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:)
      real(dp), intent(out) :: usol(:,:)
      real(dp), intent(out) :: molecules_per_particle(:,:)
      real(dp), intent(out) :: pressure(:), density(:), mix(:,:), mubar(:)
      real(dp), intent(out) :: pressure_hydro(:), density_hydro(:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Given `usol`, the densities of each species in the atmosphere,
    !> this subroutine calculates reaction rates, photolysis rates, etc.
    !> and puts this information into self.wrk. self.wrk contains all the
    !> information needed for `dochem` to compute chemistry.
    module subroutine prep_all_evo_gas(self, usol_in, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:) !! Number densities (molecules/cm^3)
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine right_hand_side_chem(self, usol, rhs, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol(:,:)
      real(dp), intent(out) :: rhs(:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Computes the right-hand-side of the ODEs describing atmospheric chemistry
    !> and transport.
    module subroutine rhs_evo_gas(self, neqs, tn, usol_flat, rhs, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      integer, intent(in) :: neqs
      real(dp), intent(in) :: tn
      real(dp), target, intent(in) :: usol_flat(neqs)
      real(dp), intent(out) :: rhs(neqs)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> The jacobian of the rhs_background_gas.
    module subroutine jac_evo_gas(self, lda_neqs, neqs, usol_flat, jac, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      integer, intent(in) :: lda_neqs, neqs
      real(dp), target, intent(in) :: usol_flat(neqs)
      real(dp), intent(out), target :: jac(lda_neqs)
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine production_and_loss(self, species, usol, pl, err)  
      use photochem_types, only: ProductionLoss   
      class(EvoAtmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: species
      real(dp), intent(in) :: usol(:,:)
      type(ProductionLoss), intent(out) :: pl
      character(:), allocatable, intent(out) :: err
    end subroutine

    !~~ photochem_evoatmosphere_integrate.f90 ~~!

    !> Evolve atmosphere through time, and saves output in a binary Fortran file.
    module function evolve(self, filename, tstart, usol_start, t_eval, overwrite, restart_from_file, err) result(success)
      use, intrinsic :: iso_c_binding
      class(EvoAtmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename !! Filename to save results.
      real(c_double), intent(inout) :: tstart !! start time in seconds
      real(dp), intent(inout) :: usol_start(:,:) !! Initial number densities (molecules/cm^3)
      real(c_double), intent(in) :: t_eval(:) !! times to evaluate the solution
      logical, optional, intent(in) :: overwrite !! If true, then overwrites pre-existing files with `filename`
      logical, optional, intent(in) :: restart_from_file !! If true, then the integration restarts from the input file.
      logical :: success !! If True, then integration was successful.
      character(:), allocatable, intent(out) :: err
    end function

    !> Determines if integration has converged to photochemical steady-state.
    module function check_for_convergence(self, err) result(converged)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      logical :: converged
    end function

    !> Initializes an integration starting at `usol_start`
    module subroutine initialize_stepper(self, usol_start, err)      
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_start(:,:) !! Initial number densities (molecules/cm^3)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Takes one internal integration step. Function `initialize_stepper`
    !> must have been called befe this
    module function step(self, err) result(tn)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      real(dp) :: tn !! Current time in the integration.
    end function
    
    !> Deallocates memory created during `initialize_stepper`
    module subroutine destroy_stepper(self, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Initializes a robust integration starting at `usol_start`
    module subroutine initialize_robust_stepper(self, usol_start, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_start(:,:) !! Initial number densities (molecules/cm^3)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Takes one robust integration step. Function `initialize_robust_stepper`
    !> must have been called before this.
    module subroutine robust_step(self, give_up, converged, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      logical, intent(out) :: give_up !! If .true., then the algorithm thinks it is time to give up.
      logical, intent(out) :: converged !! If .true., then the integration has converged to a steady state.
      character(:), allocatable, intent(out) :: err  
    end subroutine
    
    !> Integrates using a robust stepper until a steady state has been achieved.
    module function find_steady_state(self, err) result(converged)
      class(EvoAtmosphere), target, intent(inout) :: self
      character(:), allocatable, intent(out) :: err
      logical :: converged !! If .true., then the integration has converged to a steady state.
    end function
    
    module function RhsFn_evo(tn, sunvec_y, sunvec_f, user_data) &
                          result(ierr) bind(c, name='RhsFn_evo')
      use, intrinsic :: iso_c_binding
      use fcvode_mod
      use fsundials_nvector_mod
      real(c_double), value :: tn
      type(N_Vector)        :: sunvec_y
      type(N_Vector)        :: sunvec_f
      type(c_ptr), value    :: user_data
      integer(c_int)        :: ierr
    end function
    
    module function JacFn_evo(tn, sunvec_y, sunvec_f, sunmat_J, user_data, &
                          tmp1, tmp2, tmp3) &
                          result(ierr) bind(C,name='JacFn_evo')
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

    !~~ photochem_evoatmosphere_utils.f90 ~~!

    !> Saves state of the atmosphere using the concentrations in self%wrk%usol.
    module subroutine out2atmosphere_txt(self, filename, number_of_decimals, overwrite, clip, err)
      use photochem_common, only: out2atmosphere_txt_base
      class(EvoAtmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename !! Output filename
      integer, intent(in) :: number_of_decimals !! Number of decimals
      logical, intent(in) :: overwrite !! If true, then output file can be overwritten, by default False
      !> If true, then mixing ratios are clipped at a very small 
      !> positive number, by default False
      logical, intent(in) :: clip
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Computes gas fluxes at model boundaries in order to maintain
    !> current atmospheric concentrations. Uses the densities stored in
    !> self%wrk%usol.
    module subroutine gas_fluxes(self, surf_fluxes, top_fluxes, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(out) :: surf_fluxes(:) !! surface fluxes (molecules/cm^2/s)
      real(dp), intent(out) :: top_fluxes(:) !! top-of-atmosphere fluxes (molecules/cm^2/s)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Sets a lower boundary condition.
    module subroutine set_lower_bc(self, species, bc_type, vdep, den, press, flux, height, err)
      class(EvoAtmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species !! Species to set boundary condition
      character(len=*), intent(in) :: bc_type !! Boundary condition type
      real(dp), optional, intent(in) :: vdep !! Deposition velocity (cm/s)
      real(dp), optional, intent(in) :: den !! density (molecules/cm^3)
      real(dp), optional, intent(in) :: press !! pressure (dynes/cm^2)
      real(dp), optional, intent(in) :: flux !! Flux (molecules/cm^2/s)
      real(dp), optional, intent(in) :: height !! Height in atmosphere (km)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Sets upper boundary condition
    module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
      class(EvoAtmosphere), intent(inout) :: self
      character(len=*), intent(in) :: species !! Species to set boundary condition
      character(len=*), intent(in) :: bc_type !! Boundary condition type
      real(dp), optional, intent(in) :: veff !! effusion velocity (cm/s)
      real(dp), optional, intent(in) :: flux !! Flux (molecules/cm^2/s)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Sets a function describing a custom rate for a species.
    !> This could be useful for modeling external processes not in the
    !> model.
    module subroutine set_rate_fcn(self, species, fcn, err)
      use photochem_types, only: time_dependent_rate_fcn
      class(EvoAtmosphere), target, intent(inout) :: self
      character(*), intent(in) :: species !! Species name
      procedure(time_dependent_rate_fcn), pointer :: fcn
      character(:), allocatable, intent(inout) :: err
    end subroutine

    !> Changes the temperature profile.
    module subroutine set_temperature(self, temperature, trop_alt, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: temperature(:) !! new temperature at each atomspheric layer
      real(dp), optional, intent(in) :: trop_alt !! Tropopause altitude (cm). Only necessary if
                                                 !! rainout == True, or fix_water_in_trop == True.
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Given an input P, T, and edd, the code will find the temperature and eddy diffusion profile
    !> on the current altitude-grid that matches the inputs.
    module subroutine set_press_temp_edd(self, P, T, edd, trop_p, hydro_pressure, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: P(:) !! Pressure (dynes/cm^2)
      real(dp), intent(in) :: T(:) !! Temperature (K)
      real(dp), intent(in) :: edd(:) !! Eddy diffusion (cm^2/s)
      real(dp), optional, intent(in) :: trop_p !! Tropopause pressure (dynes/cm^2)
      !> If .true., then use hydrostatic pressure. If .false. then use the
      !> actual pressure in the atmosphere. Default is .true..
      logical, optional, intent(in) :: hydro_pressure
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Re-does the vertical grid so that the pressure at the top of the
    !> atmosphere is a `TOA_alt` or `TOA_pressure`. If the TOA needs to be raised above the current
    !> TOA, then the function constantly extrapolates mixing ratios, temperature,
    !> eddy diffusion, and particle radii.
    module subroutine update_vertical_grid(self, TOA_alt, TOA_pressure, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), optional, intent(in) :: TOA_alt !! New top of atmosphere altitude (cm)
      real(dp), optional, intent(in) :: TOA_pressure !! New top of atmosphere pressure (dynes/cm^2)
      character(:), allocatable, intent(out) :: err
    end subroutine

    ! Both below are needed only for `evolve` routine and probably should not be used most of the time.

    module subroutine rebin_update_vertical_grid(self, usol_old, top_atmos, usol_new, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_old(:,:)
      real(dp), intent(in) :: top_atmos
      real(dp), intent(out) :: usol_new(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> This subroutine calculates re-grids the model so that the top of the model domain 
    !> is at `top_atmos` the computes reaction rates, photolysis rates, etc.
    !> and puts this information into self.wrk. self.wrk contains all the
    !> information needed for `dochem` to compute chemistry.
    module subroutine regrid_prep_atmosphere(self, usol_new, top_atmos, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_new(:,:) !! The number densities (molecules/cm^3)
      real(dp), intent(in) :: top_atmos !! The top of the model domain (cm)
      character(:), allocatable, intent(out) :: err
    end subroutine

  end interface

end module