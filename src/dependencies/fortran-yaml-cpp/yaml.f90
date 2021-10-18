module yaml
  use iso_c_binding
  use yaml_types
  implicit none
  
  private
  
  integer,parameter :: string_length = 1024
  integer,parameter :: error_length = 1024
  
  public :: parse, error_length
  
  type type_node_c
    
    ! Node
    integer(c_int) :: T
    ! character(len=string_length, kind=c_char) :: path 
    
    ! Scalar
    type(c_ptr) :: string = c_null_ptr
    
    ! Dictionary
    type(c_ptr) :: first_keyvaluepair = c_null_ptr
    type(c_ptr) :: key = c_null_ptr
    type(c_ptr) :: value = c_null_ptr
    type(c_ptr) :: next_keyvaluepair = c_null_ptr
    
    ! List
    type(c_ptr) :: first_listitem = c_null_ptr
    type(c_ptr) :: node = c_null_ptr
    type(c_ptr) :: next_listitem = c_null_ptr
  end type
  
  interface
    subroutine LoadFile_c(filename, ptr, error) bind(C, name="LoadFile_c")
      use, intrinsic :: iso_c_binding
      character(len=1, kind = c_char), intent(in) :: filename
      type(c_ptr), intent(out) :: ptr
      character(len=1, kind = c_char), intent(in) :: error
    end subroutine
    
    subroutine DestroyNode(root) bind(C, name="DestroyNode")
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(inout) :: root
    end subroutine
  end interface
  
contains
  
  function LoadFile(filename, error) result(root)
    character(len=*), intent(in) :: filename
    character(len=string_length, kind=c_char), intent(out) :: error
    type(c_ptr) :: root
    character(len=:, kind=c_char), allocatable :: filename_copy
    filename_copy = filename//char(0)
    call LoadFile_c(filename_copy, root, error)
    deallocate(filename_copy)
  end function
  
  recursive subroutine read_value(node_c_ptr, node)
    type(c_ptr), intent(in) :: node_c_ptr
    class(type_node), intent(inout), pointer :: node
    
    type(c_ptr) :: pair_c, item_c
    type(type_node_c), pointer :: node_c, pair, item
    character(len=string_length, kind=c_char), pointer :: key
    character(len=string_length, kind=c_char), pointer :: string
    
    class (type_node), pointer :: list_item
    class (type_node), pointer :: value
    
    call c_f_pointer(node_c_ptr, node_c)
    
    if (node_c%T == 1) then
      ! is map
      allocate(type_dictionary::node)      
      pair_c = node_c%first_keyvaluepair
      do while(c_associated(pair_c))
        call c_f_pointer(pair_c, pair)
        call c_f_pointer(pair%key, key)
        call read_value(pair%value, value)
        select type (node)
          class is (type_dictionary)
            call node%set(key, value)
        end select
        pair_c = pair%next_keyvaluepair
      enddo
    elseif (node_c%T == 2) then
      ! is sequence
      allocate(type_list::node)  
      item_c = node_c%first_listitem
      do while(c_associated(item_c))
        call c_f_pointer(item_c, item)
        call read_value(item%node, list_item)
        select type (node)
          class is (type_list)
            call node%append(list_item)
        end select
        item_c = item%next_listitem
      enddo
    elseif (node_c%T == 3) then
      ! is scalar
      allocate(type_scalar::node) 
      call c_f_pointer(node_c%string, string)
      select type (node)
        class is (type_scalar)
          node%string = string
      end select
    elseif (node_c%T == 4) then
      ! is null
      allocate(type_null::node) 
    else
      print*,'Problem!!!'
    endif
    
  end subroutine  
  
  function parse(path, error) result(root)
    character(len=*), intent(in) :: path
    character(len=string_length), intent(out) :: error
    class (type_node), pointer :: root
    type(c_ptr) :: root_c
    
    nullify(root)
    
    if (len_trim(path) >= string_length) then
      error = "The path can not be longer than 1024 characters."
      return
    endif
    
    root_c = LoadFile(trim(path), error)
    if (c_associated(root_c)) then
      call read_value(root_c, root)
      call root%set_path("")
      call DestroyNode(root_c)
    endif
    
  end function
  
end module