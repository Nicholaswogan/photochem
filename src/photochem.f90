module photochem
  use photochem_const, only: err_len_tmp => err_len
  use photochem_const, only: str_len_tmp => str_len
  use photochem_const, only: real_kind_tmp => real_kind
  use photochem_atmosphere, only: Atmosphere_tmp => Atmosphere
  use photochem_types, only: ProductionLoss_tmp => ProductionLoss
  implicit none
  
  private
  
  public :: Atmosphere, ProductionLoss, err_len, str_len, real_kind
  
  integer, parameter :: err_len = err_len_tmp
  integer, parameter :: str_len = str_len_tmp
  integer, parameter :: real_kind = real_kind_tmp
  
  type, extends(Atmosphere_tmp) :: Atmosphere
  end type
  
  type, extends(ProductionLoss_tmp) :: ProductionLoss
  end type
  
end module