module photochem
  use photochem_const, only: err_len_tmp => err_len
  use photochem_const, only: str_len_tmp => str_len
  use photochem_atmosphere, only: Atmosphere_tmp => Atmosphere
  implicit none
  
  private
  
  public :: Atmosphere, err_len, str_len
  
  integer, parameter :: err_len = err_len_tmp
  integer, parameter :: str_len = str_len_tmp
  
  type, extends(Atmosphere_tmp) :: Atmosphere
  end type
  
end module