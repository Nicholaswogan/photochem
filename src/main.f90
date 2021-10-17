
program main
  use photochem_types
  use photochem_input, only: read_all_files
  implicit none
  character(len=1024) :: err
  type(PhotochemData) :: photodata
  type(PhotochemVars) :: photovars
  
  photovars%data_dir = "../data"
  
  call read_all_files("../data/reaction_mechanisms/zahnle_earth.yaml", &
                      "../templates/ModernEarth/settings_ModernEarth.yaml", &
                      "../templates/ModernEarth/Sun_now.txt", &
                      "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
                      photodata, photovars, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif


end program
