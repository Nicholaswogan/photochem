
program main
  use Photochem, only: Atmosphere, err_len, real_kind, ProductionLoss
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  type(ProductionLoss) :: pl
  logical :: success
  character(len=:), allocatable :: reaction
  integer :: i
  real(real_kind) :: tn
  
  err = ""
  call pc%init("../Photochem/data", &
               "../Photochem/data/reaction_mechanisms/apl_new.yaml", &
               "../templates/Venus/settings_Venus_apl.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/Venus/atmosphere.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  call pc%photochemical_equilibrium(success,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  ! call pc%out2atmosphere_txt("../templates/ModernEarth/atmosphere_ModernEarth.txt",.true.,.true.,err)
  ! if (len(trim(err)) > 0) then
  !   print*,trim(err)
  !   stop 1
  ! endif

end program
