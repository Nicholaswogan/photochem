
program main
  use photochem_setup, only: setup, out2atmosphere_txt
  use photochem_data, only: nq, species_names
  use photochem_vars, only: data_dir, max_order, initial_dt, equilibrium_time, &
                            verbose, usol_init, nz, usol_out, xs_folder_name
  use photochem, only: photo_equilibrium, compute_surface_fluxes
  implicit none
  character(len=1024) :: err
  real(8) :: rtol, atol
  logical :: success
  integer :: i
  real(8), allocatable :: surface_flux(:)

  data_dir = "../data"
  xs_folder_name = "xsections"

  call setup("../data/reaction_mechanisms/zahnle_earth.yaml", &
             "../templates/ModernEarth/settings_ModernEarth.yaml", &
             "../templates/ModernEarth/Sun_now.txt", &
             "../templates/ModernEarth/atmosphere_ModernEarth.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif

  equilibrium_time = 1.d17
  max_order = 5
  rtol = 1.d-3
  atol = 1.d-25
  initial_dt = 1.d-9
  verbose = 1
  call photo_equilibrium(100000, rtol, atol, success, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  allocate(surface_flux(nq))
  call compute_surface_fluxes(nq, nz, usol_out, surface_flux, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  ! do i = 1,nq
    ! print"(A10,' = ',es10.2)",species_names(i),surface_flux(i)
  ! enddo
  
  
  call out2atmosphere_txt("../atmosphere_ModernEarth.txt",.true.,.true.,err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  

end program
