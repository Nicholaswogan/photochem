module photochem
  use photochem_const, only: err_len
  use photochem_const, only: str_len
  use photochem_const, only: real_kind
  use photochem_atmosphere, only: Atmosphere
  use photochem_types, only: ProductionLoss
  use photochem_types, only: AtomConservation
  implicit none
  
  private
  
  public :: Atmosphere, ProductionLoss, AtomConservation, err_len, str_len, real_kind

end module