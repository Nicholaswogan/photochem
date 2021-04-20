module mod
contains
  subroutine test1(mapping)
    use yaml_types
    class (type_dictionary),intent(in)  :: mapping  
    type (type_error), pointer :: config_error
    
    print*,mapping%get_real('gravity',error=config_error)
  end subroutine
  
  subroutine list2realarray(mapping,listlength,arr)
    use yaml_types
    class (type_list), intent(in) :: mapping
    type (type_list_item), pointer :: item
    real(8) :: arr(listlength)
    integer :: i
    logical :: first
    i = 1
    first = .true.
    item => mapping%first
    do while (associated(item))
      if (first) then
        select type (temp => item%node)
        class is (type_scalar)
          arr(i) = temp%to_real(0.d0)
        class default
        end select
        first = .false.
        item => item%next
        i = i + 1
      else
        select type (temp => item%node)
        class is (type_scalar)
          arr(i) = temp%to_real(0.d0)
        class default
        end select
        item => item%next
        i = i + 1
      end if
    end do
  end subroutine
  
  subroutine findlistlength(mapping, length)
    use yaml_types
    class (type_list), intent(in) :: mapping
    type (type_list_item), pointer :: item
    integer, intent(out) :: length
    logical :: first
    
    length = 0
    first = .true.
    item => mapping%first
    do while (associated(item))
      if (first) then
        first = .false.
        item => item%next
        length = length + 1
      else
        item => item%next
        length = length + 1
      end if
    enddo
  end subroutine


  subroutine test(infile)
    use yaml_types
    use yaml
    implicit none

    character(error_length) :: error
    class (type_node), pointer :: root
    type (type_error), pointer :: config_error
    class (type_dictionary), pointer :: lightning
    class (type_list), pointer :: species
    type (type_list_item), pointer :: item
    logical :: first
    integer :: listlength
    real(8), allocatable :: arr(:)
    
    character(len=*), intent(in) :: infile
    
    real(8) :: gravity, pressure

    root => parse(infile,unit=100,error=error)
    if (error/='') then
      write (*,*) 'PARSE ERROR: '//trim(error)
    end if
     
    select type (root)
    class is (type_dictionary)
      ! call test1(root)
      ! gravity = root%get_real('gravity',error=config_error)
      ! pressure = root%get_real('surface-pressure',error=config_error)
      lightning => root%get_dictionary('lightning',.true.,error=config_error)
      if (associated(config_error)) call handle_yaml_error(infile,config_error%message)  
      
      species => root%get_list('species',.true.,error=config_error)
      if (associated(config_error)) call handle_yaml_error(infile,config_error%message) 
      
      ! species%append
      ! 
      call findlistlength(species, listlength)
      print*,listlength
      allocate(arr(listlength))
      call list2realarray(species,listlength,arr)
      print*,arr
      
      print*,lightning%get_logical("on-off",.true.,error=config_error)
      
    class is (type_list)
      print*, infile//" must contain dictionaries at root level"
      stop
      ! item1 => root%first%node
      ! 
      ! select type (item1)
      ! class is (type_scalar)
      !   print*, item1%to_real(1)
      ! class default
      ! end select
    class default
      print*, infile//" must contain dictionaries at root level"
      stop
    end select

  end subroutine

  subroutine handle_yaml_error(infile, message)
    implicit none
    character(len=*), intent(in) :: message, infile
    print*,infile,trim(message)
    stop
  end subroutine


end module


  
  
program main
  use mod
  call test("../templates/zahnle/photosettings.yaml")
end program





