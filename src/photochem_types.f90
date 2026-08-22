
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
  
end module
