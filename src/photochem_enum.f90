module photochem_enum
  implicit none
  public
  
  enum, bind(c)
  
  ! dtype
  enumerator :: &
    ShomatePolynomial = 1, &
    Nasa9Polynomial = 2, &
    Nasa7Polynomial = 3

  ! rxtype
  enumerator :: &
    ReverseRateType = -1, &
    PhotolysisRateType = 0, &
    ElementaryRateType = 1, &
    ThreeBodyRateType = 2, &
    FalloffRateType = 3, &
    PressDependentRateType = 4

  ! falloff_type
  enumerator :: &
    NoFalloff = 0, &
    TroeWithoutT2Falloff = 1, &
    TroeWithT2Falloff = 2, &
    JPLFalloff = 3

  ! particle_formation_method
  enumerator :: &
    CondensingParticle = 1, &
    ReactionParticle = 2

  ! particle_optical_type
  enumerator ::  &
    MieParticle = 0, &
    FractalParticle = 1
    
  ! lowerboundcond
  ! upperboundcond
  enumerator :: &
    MosesBC = -1, &
    VelocityBC = 0, &
    FluxBC = 2, &
    VelocityDistributedFluxBC = 3, &
    DensityBC = 4, &
    PressureBC = 5

  enumerator :: &
    DiffusionLimHydrogenEscape, &
    ZahnleHydrogenEscape, &
    NoHydrogenEscape

  end enum
end module