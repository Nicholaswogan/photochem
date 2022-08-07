module photochem_enum
  implicit none
  public
  
  enum, bind(c)
  
  ! dtype
  enumerator :: &
    ShomatePolynomial = 1, &
    Nasa9Polynomial = 2
  
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
    
  ! particle_sat_type
  enumerator :: &
    ArrheniusSaturation = 1, &
    H2SO4Saturation = 2
  
  ! particle_optical_type
  enumerator ::  &
    MieParticle = 0, &
    FractalParticle = 1
    
  ! lowerboundcond
  ! upperboundcond
  enumerator :: &
    MosesBC = -1, &
    VelocityBC = 0, &
    MixingRatioBC = 1, &
    FluxBC = 2, &
    VelocityDistributedFluxBC = 3, &
    DensityBC = 4, &
    MixDependentFluxBC = 5

  end enum
end module