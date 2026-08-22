
module photochem_types ! make a giant IO object
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long, c_ptr, c_null_ptr
  use photochem_const, only: dp, m_str_len
  use photochem_data, only: PhotochemData, XsectionData, ParticleXsections, &
                            ThermodynamicData, Reaction, Efficiencies, BaseRate, &
                            ReverseRate, PhotolysisRate, ElementaryRate, &
                            ThreeBodyRate, FalloffRate, PressDependentRate, &
                            MultiArrheniusRate, ProdLoss
  use photochem_settings, only: PhotoSettings, SettingsBC, CondensationParameters
  use photochem_vars, only: PhotochemVars, PressureTempEddProfile, &
                            TOAPressureMaintenance, time_dependent_flux_fcn, &
                            time_dependent_rate_fcn, binary_diffusion_fcn
  use fsundials_nvector_mod, only: N_Vector
  use fsundials_matrix_mod, only: SUNMatrix
  use fsundials_linearsolver_mod, only: SUNLinearSolver
  implicit none
  private
  
  public :: PhotoSettings, SettingsBC
  public :: XsectionData, ParticleXsections
  public :: PhotochemData, PhotochemVars, PhotochemWrk, PhotochemWrkEvo
  public :: AtmosphereState
  public :: ProductionLoss, ThermodynamicData, CondensationParameters
  public :: PressureTempEddProfile
  public :: Reaction, Efficiencies, BaseRate, PhotolysisRate, PressDependentRate, MultiArrheniusRate
  public :: TOAPressureMaintenance
  public :: ElementaryRate, ThreeBodyRate, FalloffRate, ProdLoss, ReverseRate
  public :: SundialsDataFinalizer
  public :: time_dependent_flux_fcn, time_dependent_rate_fcn, binary_diffusion_fcn
  
  type :: ProductionLoss
    real(dp), allocatable :: production(:,:)
    real(dp), allocatable :: loss(:,:)
    real(dp), allocatable :: integrated_production(:)
    real(dp), allocatable :: integrated_loss(:)
    character(len=m_str_len), allocatable :: production_rx(:)
    character(len=m_str_len), allocatable :: loss_rx(:)
  end type
  
  !> Mutable atmospheric state assembled by an initialization or mapping path.
  !! This is deliberately separate from PhotochemVars so a candidate can be
  !! validated before its fields are committed to the live model state.
  type :: AtmosphereState

    ! The state below MUST be assigned by routines that build an
    ! AtmosphereState type.
    real(dp) :: bottom_atmos
    real(dp) :: top_atmos
    real(dp) :: trop_alt
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: dz(:)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: edd(:)
    real(dp), allocatable :: particle_radius(:,:)
    real(dp), allocatable :: usol(:,:)
    real(dp), allocatable :: grav(:) ! derived
    integer :: trop_ind ! derived
    real(dp), allocatable :: xs_x_qy(:,:,:) ! derived
    type(ParticleXsections), allocatable :: particle_xs(:) ! derived
    real(dp), allocatable :: gibbs_energy(:,:) ! derived, and only allocated if `dat%reverse`

    ! Persistent atmospheric-profile configuration associated with this state.
    type(PressureTempEddProfile) :: press_temp_edd_profile
    type(TOAPressureMaintenance) :: toa_pressure_maintenance

  contains
    procedure :: allocate => AtmosphereState_allocate
  end type

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

  type :: SundialsDataFinalizer
    type(SundialsData), pointer :: sun => NULL()
  contains
    final :: SundialsDataFinalizer_final
  end type
  
  type :: PhotochemWrk
    ! PhotochemWrk are work variables that change
    ! through the course of an integration

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

    ! Work space for autodiff jacobian
    !> A sparse representation of the block diagonal Jacobian.
    !> chemistry terms only.
    real(dp), allocatable :: djac_chem(:,:)
    ! end work space for autodiff jacobian
    
  contains
    procedure :: init => PhotochemWrk_init
  end type

  type, extends(PhotochemWrk) :: PhotochemWrkEvo
    !> Surface pressure derived from the current atmospheric column (bars).
    real(dp) :: surface_pressure = 0.0_dp
    real(dp), allocatable :: mix(:,:) !! (nq,nz) mixing ratio.
    real(dp), allocatable :: pressure_hydro(:) !! (nz)
    real(dp), allocatable :: density_hydro(:) !! (nz)

    ! Runtime bookkeeping for optional robust-stepper TOA maintenance.
    integer :: n_toa_pressure_updates = 0
    integer :: n_toa_pressure_failures = 0
    integer :: nsteps_since_toa_pressure_update = 0

  contains
    procedure :: init => PhotochemWrkEvo_init

  end type
  
contains

  subroutine AtmosphereState_allocate(self, dat, nz)
    class(AtmosphereState), intent(inout) :: self
    type(PhotochemData), intent(in) :: dat
    integer, intent(in) :: nz
    integer :: i

    allocate(self%z(nz), self%dz(nz), self%temperature(nz))
    allocate(self%edd(nz), self%particle_radius(dat%npq, nz), self%usol(dat%nq, nz))

    allocate(self%grav(nz))
    allocate(self%xs_x_qy(nz, dat%kj, dat%nw))
    allocate(self%particle_xs(dat%np))

    if (dat%reverse) allocate(self%gibbs_energy(nz, dat%ng))

    do i = 1,dat%np
      ! only allocate space if there is data
      if (dat%part_xs_file(i)%ThereIsData) then
        self%particle_xs(i)%ThereIsData = .true.
        allocate(self%particle_xs(i)%w0(nz,dat%nw))
        allocate(self%particle_xs(i)%qext(nz,dat%nw))
        allocate(self%particle_xs(i)%gt(nz,dat%nw))
      else
        self%particle_xs(i)%ThereIsData = .false.
      endif
    enddo

  end subroutine

  subroutine PhotochemWrkEvo_init(self, nsp, np, nq, nz, nrT, kj, nw)
    class(PhotochemWrkEvo), intent(inout) :: self
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw

    call PhotochemWrk_init(self, nsp, np, nq, nz, nrT, kj, nw)

    self%surface_pressure = 0.0_dp
    self%n_toa_pressure_updates = 0
    self%n_toa_pressure_failures = 0
    self%nsteps_since_toa_pressure_update = 0

    if (allocated(self%mix)) then
      deallocate(self%mix)
      deallocate(self%pressure_hydro)
      deallocate(self%density_hydro)
    endif

    allocate(self%mix(nq,nz))
    allocate(self%pressure_hydro(nz))
    allocate(self%density_hydro(nz))

  end subroutine
 
  subroutine PhotochemWrk_init(self, nsp, np, nq, nz, nrT, kj, nw)
    use photochem_const, only: nsteps_save
    class(PhotochemWrk), intent(inout) :: self
    integer, intent(in) :: nsp, np, nq, nz, nrT, kj, nw

    self%robust_stepper_initialized = .false.
    self%nsteps_total = -1
    self%nerrors_total = -1
    
    if (allocated(self%usol)) then
      deallocate(self%t_history)
      deallocate(self%mix_history)
      deallocate(self%dmix)
      deallocate(self%usol)
      deallocate(self%mubar)
      deallocate(self%pressure)
      deallocate(self%density)
      deallocate(self%densities)
      deallocate(self%rx_rates)
      deallocate(self%prates)
      deallocate(self%surf_radiance)
      deallocate(self%amean_grd)
      deallocate(self%optical_depth)
      deallocate(self%xp)
      deallocate(self%xl)
      deallocate(self%DU)
      deallocate(self%DD)
      deallocate(self%DL)
      deallocate(self%ADU)
      deallocate(self%ADL)
      deallocate(self%ADD)
      deallocate(self%upper_veff_copy)
      deallocate(self%lower_vdep_copy)
      deallocate(self%scale_height)
      deallocate(self%wfall)
      deallocate(self%gas_sat_den)
      deallocate(self%molecules_per_particle)
      deallocate(self%rainout_rates)
      deallocate(self%djac_chem)
    endif
    
    allocate(self%t_history(nsteps_save))
    allocate(self%mix_history(nq,nz,nsteps_save))
    allocate(self%dmix(nq,nz))
    allocate(self%usol(nq,nz))
    allocate(self%mubar(nz))
    allocate(self%pressure(nz))
    allocate(self%density(nz))
    allocate(self%densities(nsp+1,nz))
    allocate(self%rx_rates(nz,nrT))
    allocate(self%prates(nz,kj))
    allocate(self%surf_radiance(nw))
    allocate(self%amean_grd(nz,nw))
    allocate(self%optical_depth(nz,nw))
    allocate(self%xp(nz))
    allocate(self%xl(nz))
    allocate(self%DU(nq,nz))
    allocate(self%DD(nq,nz))
    allocate(self%DL(nq,nz))
    allocate(self%ADU(nq,nz))
    allocate(self%ADL(nq,nz))
    allocate(self%ADD(nq,nz))
    allocate(self%upper_veff_copy(nq))
    allocate(self%lower_vdep_copy(nq))
    allocate(self%scale_height(nz))
    allocate(self%wfall(np,nz))
    allocate(self%gas_sat_den(np,nz))
    allocate(self%molecules_per_particle(np,nz))
    allocate(self%rainout_rates(nq,nz))
    allocate(self%djac_chem(nq,nz*nq))
  end subroutine

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
