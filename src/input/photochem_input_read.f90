submodule (photochem_input) photochem_input_read
  implicit none

contains

  module subroutine read_static_files(mechanism_file, s, flux_file, data_dir, &
                                      dat, var, err)
    character(len=*), intent(in) :: mechanism_file
    type(PhotoSettings), intent(in) :: s
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: data_dir
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    call read_photochem_mechanism(mechanism_file, s, dat, err)
    if (allocated(err)) return

    call apply_data_settings(s, dat, err)
    if (allocated(err)) return

    call apply_vars_settings(dat, s, var, err)
    if (allocated(err)) return

    call read_photochem_supporting_data(dat, data_dir, s, err)
    if (allocated(err)) return

    allocate(var%photon_flux(dat%nw))
    call read_stellar_flux(flux_file, dat%nw, dat%wavl, var%photon_flux, err)
  end subroutine

  subroutine apply_data_settings(s, dat, err)
    use photochem_enum, only: ZahnleHydrogenEscape, NoHydrogenEscape
    use photochem_eqns, only: zahnle_Hescape_coeff
    type(PhotoSettings), intent(in) :: s
    type(PhotochemData), intent(inout) :: dat
    character(:), allocatable, intent(out) :: err

    integer :: ind(1)

    dat%planet_mass = s%planet_mass
    dat%planet_radius = s%planet_radius
    dat%H_escape_type = s%H_escape_type
    if (dat%H_escape_type /= NoHydrogenEscape) then
      ind = findloc(dat%species_names, 'H2')
      dat%LH2 = ind(1)
      if (ind(1) == 0) then
        err = 'IOError: H2 must be a species if hydrogen-escape is on.'
        return
      endif
      ind = findloc(dat%species_names, 'H')
      dat%LH = ind(1)
      if (ind(1) == 0) then
        err = 'IOError: H must be a species if hydrogen-escape is on.'
        return
      endif

      ind = findloc(dat%species_names(dat%nq+1:dat%nq+dat%nsl), 'H2')
      if (ind(1) /= 0) then
        err = 'H2 must not be short lived if hydrogen escape is on'
        return
      endif
      ind = findloc(dat%species_names(dat%nq+1:dat%nq+dat%nsl), 'H')
      if (ind(1) /= 0) then
        err = 'H must not be short lived if hydrogen escape is on'
        return
      endif

      if (dat%H_escape_type == ZahnleHydrogenEscape) then
        dat%H_escape_coeff = zahnle_Hescape_coeff(s%H_escape_S1)
      endif
    endif

    dat%gas_rainout = s%gas_rainout
    ind = findloc(dat%species_names, 'H2O')
    dat%LH2O = ind(1)
    if (ind(1) == 0 .and. dat%gas_rainout) then
      err = 'IOError: H2O must be a species if gas-rainout = True.'
      return
    endif
  end subroutine

  module subroutine read_atmosphere_file(atmosphere_txt, dat, profile, err)
    use futils, only: FileCloser
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(in) :: dat
    type(AtmosphereFileProfile), intent(out) :: profile
    character(:), allocatable, intent(out) :: err

    character(len=10000) :: line
    character(len=s_str_len) :: arr1(1000)
    character(len=s_str_len) :: arr11(1000)
    character(len=s_str_len),allocatable, dimension(:) :: labels
    integer :: ind(1)
    real(dp), allocatable :: temp(:,:)
    integer :: io, i, n, nn, ii
    type(FileCloser) :: file

    open(4, file=trim(atmosphere_txt),status='old',iostat=io)
    file%unit = 4
    if (io /= 0) then
      err = 'Can not open file '//trim(atmosphere_txt)
      return
    endif
    read(4,'(A)') line

    profile%nlayer = -1
    io = 0
    do while (io == 0)
      read(4,*,iostat=io)
      profile%nlayer = profile%nlayer + 1
    enddo

    allocate(profile%z(profile%nlayer))
    allocate(profile%temperature(profile%nlayer))
    allocate(profile%edd(profile%nlayer))
    allocate(profile%density(profile%nlayer))
    allocate(profile%mix(dat%nq, profile%nlayer))
    profile%z = 0.0_dp
    profile%temperature = 0.0_dp
    profile%edd = 0.0_dp
    profile%mix = 1.0e-40_dp
    if (dat%there_are_particles) then
      allocate(profile%particle_radius(dat%npq, profile%nlayer))
    endif

    rewind(4)
    read(4,'(A)') line
    n = 0
    nn = 0
    do i=1,1000
      read(line,*,iostat=io) arr1(1:i)
      if (io==-1) exit
      n = n+1
    enddo
    read(4,'(A)') line
    do i=1,1000
      read(line,*,iostat=io) arr11(1:i)
      if (io==-1) exit
      nn = nn+1
    enddo
    if (n /= nn) then
      err = 'There is a missing column label in the file '//trim(atmosphere_txt)
      return
    endif

    ! allocate memory
    allocate(labels(n))
    allocate(temp(n,profile%nlayer))
    rewind(4)
    read(4,'(A)') line
    read(line,*) (labels(i),i=1,n)

    ! First read in all the data into big array
    do i = 1,profile%nlayer
      read(4,*,iostat=io) (temp(ii,i),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in initial atmosphere in '//trim(atmosphere_txt)
        return
      endif
    enddo

    ! reads in mixing ratios
    do i=1,dat%nq
      ind = findloc(labels,dat%species_names(i))
      if (ind(1) /= 0) then
        profile%mix(i,:) = temp(ind(1),:)
      endif
    enddo

    if (dat%there_are_particles) then
      do i=1,dat%npq
        ind = findloc(labels,trim(dat%species_names(i))//"_r")
        if (ind(1) /= 0) then
          profile%particle_radius(i,:) = temp(ind(1),:)
        else
          ! did not find the data
          ! will set to 0.1 micron
          profile%particle_radius(i,:) = 1.0e-5_dp
        endif
      enddo
    endif

    ! reads in temperature
    ind = findloc(labels,'temp')
    if (ind(1) /= 0) then
      profile%temperature(:) = temp(ind(1),:)
    else
      err = '"temp" was not found in input file '//trim(atmosphere_txt)
      return
    endif

    ! reads in alt
    ind = findloc(labels,'alt')
    if (ind(1) /= 0) then
      profile%z(:) = temp(ind(1),:)*1.e5_dp ! convert to cm
    else
      err = '"alt" was not found in input file '//trim(atmosphere_txt)
      return
    endif

    ! reads in eddy diffusion
    ind = findloc(labels,'eddy')
    if (ind(1) /= 0) then
      profile%edd(:) = temp(ind(1),:)
    else
      err = '"eddy" was not found in input file '//trim(atmosphere_txt)
      return
    endif

    ! reads in density.
    ind = findloc(labels,'den')
    if (ind(1) /= 0) then
      profile%density(:) = temp(ind(1),:)
    else
      err = '"den" was not found in input file '//trim(atmosphere_txt)
      return
    endif

  end subroutine

end submodule
