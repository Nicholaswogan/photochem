
program test
  use, intrinsic :: iso_fortran_env, only:  output_unit
  use yaml, only: parse, error_length
  use yaml_types, only: type_node
  
  class(type_node), pointer :: root
  character(len=error_length) :: error
  
  root => parse("../test.yaml", error = error)
  if (error/='') then
    print*,trim(error)
    stop 1
  endif
  
  call root%dump(unit=output_unit,indent=0)
  
  call root%finalize()
  deallocate(root)
end program


