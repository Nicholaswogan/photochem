module photochem_evoatmosphere
  use photochem_const, only: dp
  use photochem_types, only : PhotochemData, PhotochemVars, PhotochemWrkEvo
  use clima_radtran, only: Radtran
  implicit none

  private
  public :: EvoAtmosphere

  type :: EvoAtmosphere
    type(PhotochemData), allocatable :: dat
    type(PhotochemVars), allocatable :: var
    type(PhotochemWrkEvo), allocatable :: wrk

    logical :: evolve_climate
    real(dp) :: T_surf
    real(dp) :: T_trop = 200.0_dp
    type(Radtran), allocatable :: rad

    ! Modern Earth has a pressure of 4e-7 at 100 km
    ! so it makes sense to try to keep pressure between these values
    real(dp) :: P_top_min = 1.0e-7_dp ! bars
    real(dp) :: P_top_max = 1.0e-6_dp ! bars
    real(dp) :: top_atmos_adjust_frac = 0.02 ! fraction

  contains
  !!! photochem_evoatmosphere_init.f90 !!!
    procedure :: init => EvoAtmosphere_init

    !!! photochem_evoatmosphere_rhs.f90 !!!
    procedure :: set_trop_ind
    procedure :: prep_atm_evo_gas
    procedure :: prep_atmosphere => prep_all_evo_gas
    procedure :: right_hand_side => rhs_evo_gas
    procedure :: jacobian => jac_evo_gas
    procedure :: production_and_loss

    !!! photochem_evoatmosphere_integrate.f90 !!!
    procedure :: evolve

    !!! photochem_evoatmosphere_utils.f90 !!!
    procedure :: rebin_update_vertical_grid
    procedure :: regrid_prep_atmosphere

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

    module subroutine prep_all_evo_gas(self, usol_in, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_in(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine rhs_evo_gas(self, neqs, usol_flat, rhs, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      integer, intent(in) :: neqs
      real(dp), target, intent(in) :: usol_flat(neqs)
      real(dp), intent(out) :: rhs(neqs)
      character(:), allocatable, intent(out) :: err
      ! Computes the right-hand-side of the ODEs describing atmospheric chemistry
      ! and transport.
    end subroutine
    
    module subroutine jac_evo_gas(self, lda_neqs, neqs, usol_flat, jac, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      integer, intent(in) :: lda_neqs, neqs
      real(dp), target, intent(in) :: usol_flat(neqs)
      real(dp), intent(out), target :: jac(lda_neqs)
      character(:), allocatable, intent(out) :: err
      ! The jacobian of the rhs_background_gas.
    end subroutine

    module subroutine production_and_loss(self, species, usol, top_atmos, pl, err)  
      use photochem_types, only: ProductionLoss   
      class(EvoAtmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: species
      real(dp), intent(in) :: usol(:,:)
      real(dp), intent(in) :: top_atmos
      type(ProductionLoss), intent(out) :: pl
      character(:), allocatable, intent(out) :: err
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_evoatmosphere_integrate.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module function evolve(self, filename, tstart, usol_start, t_eval, overwrite, restart_from_file, err) result(success)
      use, intrinsic :: iso_c_binding
      class(EvoAtmosphere), target, intent(inout) :: self
      character(len=*), intent(in) :: filename
      real(c_double), intent(inout) :: tstart
      real(dp), intent(inout) :: usol_start(:,:)
      real(c_double), intent(in) :: t_eval(:)
      logical, optional, intent(in) :: overwrite
      logical, optional, intent(in) :: restart_from_file
      logical :: success
      character(:), allocatable, intent(out) :: err
      ! Evolve atmosphere through time, and saves output in a binary Fortran file.
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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! photochem_evoatmosphere_utils.f90 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module subroutine rebin_update_vertical_grid(self, usol_old, top_atmos, usol_new, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_old(:,:)
      real(dp), intent(in) :: top_atmos
      real(dp), intent(out) :: usol_new(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine regrid_prep_atmosphere(self, usol_new, top_atmos, err)
      class(EvoAtmosphere), target, intent(inout) :: self
      real(dp), intent(in) :: usol_new(:,:)
      real(dp), intent(in) :: top_atmos
      character(:), allocatable, intent(out) :: err
    end subroutine

  end interface

end module