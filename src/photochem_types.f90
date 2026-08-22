
module photochem_types ! make a giant IO object
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
  use photochem_wrk, only: PhotochemWrk, SundialsDataFinalizer
  implicit none
  private
  
  public :: PhotoSettings, SettingsBC
  public :: XsectionData, ParticleXsections
  public :: PhotochemData, PhotochemVars, PhotochemWrk
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

end module
