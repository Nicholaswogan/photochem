module photochem_c_api
  use iso_c_binding
  use clima_saturationdata, only: SaturationData
  use photochem, only: EvoAtmosphere
  use photochem_types, only: PhotochemData
  use photochem_types, only: PhotochemVars
  use photochem_types, only: PhotochemWrk, PhotochemWrkEvo
  use photochem_types, only: ProductionLoss
  use photochem_types, only: CondensationParameters
  use photochem, only: err_len
  use photochem_const, only: s_str_len, m_str_len
  implicit none

contains

  include "EvoAtmosphere_wrapper.f90"

  include "PhotochemData_wrapper.f90"
  include "PhotochemVars_wrapper.f90"
  include "PhotochemWrk_wrapper.f90"

  include "AtomConservation_wrapper.f90"
  include "ProductionLoss_wrapper.f90"

  subroutine photochem_version_get(version_c) bind(c)
    use photochem, only: version 
    character(kind=c_char), intent(out) :: version_c(100+1)
    call copy_string_ftoc(version, version_c)
  end subroutine

  !!!!!!!!!!!!!!!!!!
  !!! Utilities  !!!
  !!!!!!!!!!!!!!!!!!
  
  function len_cstring(stringc) result (length)
    ! DOES NOT include the null character terminating c string
    character(kind=c_char), intent(in) :: stringc(*)
    integer(c_int) :: length
    integer, parameter :: max_len = 10000
    integer :: j  
    j = 1
    do
      if (stringc(j)==c_null_char) then
        length = j - 1
        exit
      endif
      if (j == max_len) then
        print*,"'len_cstring' tried to determine the length of an invalid C string"
        stop 1
      endif
      j = j + 1
    end do
  end function
  
  subroutine copy_string_ctof(stringc,stringf)
    ! utility function to convert c string to fortran string
    character(len=*), intent(out) :: stringf
    character(c_char), intent(in) :: stringc(*)
    integer j
    stringf = ''
    char_loop: do j=1,len(stringf)
       if (stringc(j)==c_null_char) exit char_loop
       stringf(j:j) = stringc(j)
    end do char_loop
  end subroutine copy_string_ctof

  subroutine copy_string_ftoc(stringf,stringc)
    ! utility function to convert c string to fortran string
    character(len=*), intent(in) :: stringf
    character(c_char), intent(out) :: stringc(:)
    integer j, n, n1, n2
    n1 = len_trim(stringf)  
    n2 = size(stringc) - 1
    n = min(n1, n2)
    do j=1,n    
      stringc(j) = stringf(j:j)   
    end do
    stringc(n+1) = c_null_char
  end subroutine copy_string_ftoc

end module