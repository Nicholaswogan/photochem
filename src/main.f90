
program main
  use photochem_io, only: get_photomech, get_photorad, get_photoset, reaction_string
  use photochem_types, only: PhotoMechanism, PhotoRadTran, PhotoSettings
  implicit none
  type(PhotoSettings) :: photoset
  type(PhotoMechanism) :: photomech
  type(PhotoRadTran) :: photorad
  character(len=1000) :: err
  character(len=:), allocatable :: rxstring
  integer i
  
  call get_photomech("../zahnle_rx.yaml", photomech, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    print*,'error worked'
    stop
  endif
  
  call get_photoset("../settings.yaml", photomech, photoset, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    print*,'error worked!'
    stop
  endif
  
  call get_photorad(photomech, photoset, photorad, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    print*,'error worked rad'
    stop
  endif
  
  ! print*,photomech%reactants_sp_inds
 
  ! do i=1,photomech%nrT
  !   call reaction_string(photomech,i,rxstring)
  !   print*,rxstring
  !   if (i <= photomech%nrF) then
  !     print*,photomech%reactants_sp_inds(:,i)
  !   endif
  ! enddo
  ! deallocate(rxstring)
  
  ! do i=1,photomech%nsp
  !   print*,photomech%species_mass(i), photomech%species_names(i)
  ! enddo
  
end program






subroutine test(filename,rxn,rxstring1,err)
  use photochem_io, only: get_photomech, reaction_string
  use photochem_types, only: PhotoMechanism
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: rxn
  character(len=1000), intent(out) :: err
  character(len=1000), intent(out) :: rxstring1
  
  
  type(PhotoMechanism) :: photomech
  character(len=:), allocatable :: rxstring
  ! integer i
  
  call get_photomech(filename, photomech, err)
  if (len_trim(err) > 0) return 
  
  ! get
  if ((rxn>photomech%nrT) .or. (rxn<1)) then
    err = 'Index is not in range'
    return
  endif
  ! do i=1,photomech%nrT
  call reaction_string(photomech,rxn,rxstring)
  rxstring1 = rxstring
    ! print*,rxstring
  ! enddo
  
end subroutine