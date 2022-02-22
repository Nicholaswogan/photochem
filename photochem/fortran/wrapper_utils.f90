module wrapper_utils
  use iso_c_binding
  implicit none
  
contains
  
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