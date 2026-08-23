module photochem_wrk
  use, intrinsic :: iso_c_binding, only: c_double, c_long, c_ptr, c_null_ptr
  use photochem_const, only: dp
  use fsundials_nvector_mod, only: N_Vector
  use fsundials_matrix_mod, only: SUNMatrix
  use fsundials_linearsolver_mod, only: SUNLinearSolver
  implicit none
  private

  public :: PhotochemWrk
  public :: SundialsDataFinalizer

  type :: SundialsData
    !> cvode memory
    type(c_ptr) :: cvode_mem = c_null_ptr
    ! solution vector
    real(c_double), allocatable :: yvec(:)
    type(N_Vector), pointer :: sunvec_y => NULL()
    ! absolute tolerance (Used in EvoAtmosphere only)
    real(c_double), allocatable :: abstol(:)
    type(N_Vector), pointer :: abstol_nvec => NULL()
    ! matrix and linear solver
    type(SUNMatrix), pointer :: sunmat => NULL()
    type(SUNLinearSolver), pointer :: sunlin => NULL()
  contains
    procedure :: finalize => SundialsData_finalize
    final :: SundialsData_final
  end type

  !> Scoped finalizer for SUNDIALS resources owned by a work object.
  type :: SundialsDataFinalizer
    type(SundialsData), pointer :: sun => NULL()
  contains
    final :: SundialsDataFinalizer_final
  end type

  !> Runtime work arrays and integration state for an evolving atmosphere.
  type :: PhotochemWrk

    ! Total step counter for robust_step method
    !> True while the CVODE stepper belongs to an initialized robust
    !> integration session.
    logical :: robust_stepper_initialized = .false.
    !> Total number of steps in a robust integration.
    integer :: nsteps_total = -1
    !> Total number of errors experienced in the robust integration.
    integer :: nerrors_total = -1

    ! used in cvode
    integer(c_long) :: nsteps_previous = -10 !! For printing
    type(SundialsData) :: sun !! CVODE data

    ! All for determining convergence
    integer :: nsteps = 0 !! Number of integration steps excuted. Updated
                          !! after every successful step.
    !> History of times at previous integration steps. Index 1 is current,
    !> while index 2, 3, 4 are previous steps. Updated after every successful step.
    real(dp), allocatable :: t_history(:)
    !> History of mixing ratios at previous integration steps. Index 1 is
    !> current, while index 2, 3, 4 are previous steps. Updated after
    !> every successful step.
    real(dp), allocatable :: mix_history(:,:,:)
    !> Change in mixing ratio over some number of integrations steps. Updated
    !> during convergence checking.
    real(dp), allocatable :: dmix(:,:)
    !> Normalized change in mixing ratios over some number of integrations steps
    real(dp) :: longdy = 0.0_dp
    !> Normalized change in mixing ratios divided by change in time
    !> over some number of integrations steps.
    real(dp) :: longdydt = 0.0_dp
    ! end stuff for determining convergence

    !> The current time (seconds). The is updated with each call to the
    !> right hand side, and jacobian. It is only important if
    !> var%photon_flux_fcn is set.
    real(dp) :: tn = 0.0_dp

    ! Used in prep_all_background_gas
    ! work arrays
    real(dp), allocatable :: usol(:,:) !! (nq,nz)
    real(dp), allocatable :: densities(:,:) !! (nsp+1,nz)
    real(dp), allocatable :: density(:) !! (nz)
    real(dp), allocatable :: rx_rates(:,:) !! (nz,nrT)
    real(dp), allocatable :: mubar(:) !! (nz)
    real(dp), allocatable :: pressure(:) !! (nz)
    real(dp), allocatable :: prates(:,:) !! (nz,kj)
    real(dp), allocatable :: surf_radiance(:) !! (nw)
    real(dp), allocatable :: amean_grd(:,:) !! (nz,nw)
    real(dp), allocatable :: optical_depth(:,:) !! (nz,nw)
    real(dp), allocatable :: upper_veff_copy(:) !! (nq)
    real(dp), allocatable :: lower_vdep_copy(:) !! (nq)
    real(dp), allocatable :: xp(:) !! (nz)
    real(dp), allocatable :: xl(:) !! (nz)
    ! diffusion and H escape
    real(dp), allocatable :: DU(:,:) !! (nq,nz)
    real(dp), allocatable :: DD(:,:) !! (nq,nz)
    real(dp), allocatable :: DL(:,:) !! (nq,nz)
    real(dp), allocatable :: ADU(:,:) !! (nq,nz)
    real(dp), allocatable :: ADL(:,:) !! (nq,nz)
    real(dp), allocatable :: ADD(:,:) !! (nq,nz)
    real(dp) :: VH2_esc
    real(dp) :: VH_esc
    ! other
    real(dp), allocatable :: scale_height(:)
    real(dp), allocatable :: wfall(:,:)
    real(dp), allocatable :: gas_sat_den(:,:)
    real(dp), allocatable :: molecules_per_particle(:,:)
    real(dp), allocatable :: rainout_rates(:,:)
    ! end used in prep_all_background_gas

    ! Work space for block-diagonal chemistry Jacobians
    !> A sparse representation of the block diagonal Jacobian.
    !> chemistry terms only.
    real(dp), allocatable :: djac_chem(:,:)
    ! end chemistry Jacobian work space

    !> Surface pressure derived from the current atmospheric column (bars).
    real(dp) :: surface_pressure = 0.0_dp
    real(dp), allocatable :: mix(:,:) !! (nq,nz) mixing ratio.
    real(dp), allocatable :: pressure_hydro(:) !! (nz)
    real(dp), allocatable :: density_hydro(:) !! (nz)

    ! Runtime bookkeeping for optional robust-stepper TOA maintenance.
    integer :: n_toa_pressure_updates = 0
    integer :: n_toa_pressure_failures = 0
    integer :: nsteps_since_toa_pressure_update = 0
  end type

  interface PhotochemWrk
    module procedure :: create_PhotochemWrk
  end interface

contains

  !> Allocate a work object for the supplied model dimensions.
  function create_PhotochemWrk(nsp, np, nq, nz, nrT, kj, nw) result(wrk)
    use photochem_const, only: nsteps_save
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw
    type(PhotochemWrk) :: wrk

    allocate(wrk%t_history(nsteps_save))
    allocate(wrk%mix_history(nq,nz,nsteps_save))
    allocate(wrk%dmix(nq,nz))
    allocate(wrk%usol(nq,nz))
    allocate(wrk%mubar(nz))
    allocate(wrk%pressure(nz))
    allocate(wrk%density(nz))
    allocate(wrk%densities(nsp+1,nz))
    allocate(wrk%rx_rates(nz,nrT))
    allocate(wrk%prates(nz,kj))
    allocate(wrk%surf_radiance(nw))
    allocate(wrk%amean_grd(nz,nw))
    allocate(wrk%optical_depth(nz,nw))
    allocate(wrk%xp(nz))
    allocate(wrk%xl(nz))
    allocate(wrk%DU(nq,nz))
    allocate(wrk%DD(nq,nz))
    allocate(wrk%DL(nq,nz))
    allocate(wrk%ADU(nq,nz))
    allocate(wrk%ADL(nq,nz))
    allocate(wrk%ADD(nq,nz))
    allocate(wrk%upper_veff_copy(nq))
    allocate(wrk%lower_vdep_copy(nq))
    allocate(wrk%scale_height(nz))
    allocate(wrk%wfall(np,nz))
    allocate(wrk%gas_sat_den(np,nz))
    allocate(wrk%molecules_per_particle(np,nz))
    allocate(wrk%rainout_rates(nq,nz))
    allocate(wrk%djac_chem(nq,nz*nq))
    allocate(wrk%mix(nq,nz))
    allocate(wrk%pressure_hydro(nz))
    allocate(wrk%density_hydro(nz))
  end function

  subroutine SundialsData_finalize(self, err)
    use iso_c_binding, only: c_associated, c_null_ptr, c_int
    use fcvode_mod, only: FCVodeFree
    use fsundials_nvector_mod, only: FN_VDestroy
    use fsundials_matrix_mod, only: FSUNMatDestroy
    use fsundials_linearsolver_mod, only: FSUNLinSolFree
    class(SundialsData), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    integer(c_int) :: ierr

    if (associated(self%sunvec_y)) then
      call FN_VDestroy(self%sunvec_y)
      nullify(self%sunvec_y)
    endif
    if (allocated(self%yvec)) then
      deallocate(self%yvec)
    endif

    if (associated(self%abstol_nvec)) then
      call FN_VDestroy(self%abstol_nvec)
      nullify(self%abstol_nvec)
    endif
    if (allocated(self%abstol)) then
      deallocate(self%abstol)
    endif

    if (c_associated(self%cvode_mem)) then
      call FCVodeFree(self%cvode_mem)
      self%cvode_mem = c_null_ptr
    endif

    if (associated(self%sunlin)) then
      ierr = FSUNLinSolFree(self%sunlin)
      if (ierr /= 0) then
        err = "Sundials failed to deallocated linear solver"
      end if
      nullify(self%sunlin)
    endif

    if (associated(self%sunmat)) then
      call FSUNMatDestroy(self%sunmat)
      nullify(self%sunmat)
    endif
  end subroutine

  subroutine SundialsData_final(self)
    type(SundialsData), intent(inout) :: self
    character(:), allocatable :: err
    call SundialsData_finalize(self, err)
  end subroutine

  subroutine SundialsDataFinalizer_final(self)
    type(SundialsDataFinalizer), intent(inout) :: self
    character(:), allocatable :: err
    if (associated(self%sun)) then
      call self%sun%finalize(err)
    endif
  end subroutine

end module
