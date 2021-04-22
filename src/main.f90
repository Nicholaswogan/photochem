
program main
  use photochem_io, only: get_photodata, print_reaction
  use photochem_types, only: PhotoMechanism
  implicit none
  type(PhotoMechanism) :: photomech
  integer i, j, k

  call get_photodata("../zahnle.yaml", photomech)

  do i=1,photomech%nrT
    call print_reaction(photomech,i)
  enddo
  
end program