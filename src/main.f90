
program main
  use Photochem, only: Atmosphere, err_len, real_kind, ProductionLoss, AtomConservation
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  type(ProductionLoss) :: pl
  logical :: success
  character(len=:), allocatable :: reaction
  integer :: i
  real(real_kind) :: tn, redox
  type(AtomConservation) :: con
  
  err = ""
  call pc%init("../Photochem/data", &
               "../Photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif


  
  
  
  ! stop 1
  
  call pc%initialize_stepper(pc%wrk%usol,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  tn = 0.d0
  do while(tn < 1d17)
    do i = 1,10
      tn = pc%step(err)
      if (len(trim(err)) > 0) then
        print*,trim(err)
        stop 1
      endif
      
    enddo
    redox = pc%redox_conservation(err)
  
  enddo
  stop 1
  call pc%out2atmosphere_txt("../atmosphere_tiny.txt",.true.,.true.,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  

end program
