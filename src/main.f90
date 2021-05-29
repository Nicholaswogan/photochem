
program main
  use photochem_input, only: read_all_files, reaction_string
  use photochem_types, only: PhotoMechanism, PhotoRadTran, PhotoSettings, PhotoInitAtm
  use photochem_data, only: read_input_files
  implicit none
  type(PhotoSettings) :: photoset
  type(PhotoMechanism) :: photomech
  type(PhotoRadTran) :: photorad
  type(PhotoInitAtm) :: photoinit
  character(len=1024) :: err
  character(len=:), allocatable :: rxstring
  integer i
  
  ! read(*,*) i
  
  ! call read_all_files("../zahnle.yaml", "../settings.yaml", "../Sun_4.0Ga.txt", "../atmosphere.txt", &
  !                     photomech, photoset, photorad, photoinit, err)
  ! if (len(trim(err)) > 0) then
  !   print*,trim(err)
  !   print*,'error worked rad'
  !   stop
  ! endif
  
  call read_input_files("../zahnle.yaml", "../settings.yaml", "../Sun_4.0Ga.txt", "../atmosphere.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    print*,'error worked rad'
    stop
  endif
  
  ! read(*,*) i
  
  ! print*,photomech%reactants_sp_inds
 
  ! do i=1,photomech%nrT
  !   call reaction_string(photomech,i,rxstring)
  !   print*,i,rxstring
  ! enddo
  ! deallocate(rxstring)
  
  ! do i=1,photomech%nsp
  !   print*,photomech%species_mass(i), photomech%species_names(i)
  ! enddo
  
end program






subroutine test(filename,rxn,rxstring1,err)
  use photochem_input, only: get_photomech, reaction_string
  use photochem_types, only: PhotoMechanism
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: rxn
  character(len=1024), intent(out) :: err
  character(len=1024), intent(out) :: rxstring1
  
  
  type(PhotoMechanism) :: photomech
  character(len=:), allocatable :: rxstring
  ! integer i
  
  ! call get_photomech(filename, photomech, err)
  ! if (len_trim(err) > 0) return 
  
  ! get
  ! if ((rxn>photomech%nrT) .or. (rxn<1)) then
  !   err = 'Index is not in range'
  !   return
  ! endif
  ! do i=1,photomech%nrT
  ! call reaction_string(photomech,rxn,rxstring)
  ! rxstring1 = rxstring
    ! print*,rxstring
  ! enddo
  
end subroutine