module Atmos
  use photochem_const, only: err_len_tmp => err_len
  use photochem_const, only: str_len_tmp => str_len
  use photochem_object, only: Photochem_tmp => Photochem
  implicit none
  
  private
  
  public :: Photochem, err_len, str_len
  
  integer, parameter :: err_len = err_len_tmp
  integer, parameter :: str_len = str_len_tmp
  
  type, extends(Photochem_tmp) :: Photochem
  end type
  
end module