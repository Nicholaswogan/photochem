
program main
  use photochem_io, only: get_photomech, print_reaction
  use photochem_types, only: PhotoMechanism
  implicit none
  type(PhotoMechanism) :: photomech
  integer i

  call get_photomech("../zahnle.yaml", photomech)

  do i=1,photomech%nrT
    call print_reaction(photomech,i)
  enddo
  
  do i=1,photomech%nsp
    print*,photomech%species_mass(i), photomech%species_names(i)
  enddo
  
end program