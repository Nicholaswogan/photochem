module photochem_mie
  implicit none
  integer, private, parameter :: real_kind = 8
  integer, private, parameter :: err_len = 1024
  
contains
  

  ! Reads mie binary data file, then interpolates optical data to wavelength grid.
  ! also returns the radii of particles for that file.
  subroutine read_mie_data_file(filename, nw, wavl, &
                                 nrad_file, radii_file, w0_file, qext_file, g_file, err)
    use interp_tools, only: addpnt, inter2
    
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nw
    real(real_kind), intent(in) :: wavl(nw+1)
    
    integer, intent(out) :: nrad_file
    real(real_kind), allocatable, intent(out) :: radii_file(:)
    real(real_kind), allocatable, intent(out) :: w0_file(:,:), qext_file(:,:), g_file(:,:)
    character(len=err_len), intent(out) :: err
    
    real(real_kind), allocatable :: wavl_tmp(:)
    real(real_kind), allocatable :: w0_tmp(:,:), qext_tmp(:,:), g_tmp(:,:)
    real(real_kind), allocatable :: temp_data(:), temp_wavelength(:)
    
    integer :: nw_tmp
    real(real_kind) :: dum
    integer :: i, j, io, ierr
    
    err = ''
    
    open(2,file=filename,form="unformatted",status='old',iostat=io)
    if (io /= 0) then
      err = "Was unable to open mie data file "//trim(filename)
      return
    endif
    
    read(2, iostat=io) nw_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    allocate(wavl_tmp(nw_tmp))
    read(2, iostat=io) wavl_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    
    read(2, iostat=io) nrad_file
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    allocate(radii_file(nrad_file))
    read(2, iostat=io) radii_file
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    
    allocate(w0_tmp(nw_tmp,nrad_file))
    allocate(qext_tmp(nw_tmp,nrad_file))
    allocate(g_tmp(nw_tmp,nrad_file))
    
    read(2, iostat=io) w0_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    read(2, iostat=io) qext_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    read(2, iostat=io) g_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    
    ! this next read should be usuccessful
    read(2, iostat=io) dum
    if (io == 0) then
      err = "Problem reading mie data file "//trim(filename)// &
            ". Should have reached the end of the file, but did not."
      return
    endif
    close(2)
    
    ! now lets interpolate a bunch
    allocate(w0_file(nrad_file,nw))
    allocate(qext_file(nrad_file,nw))
    allocate(g_file(nrad_file,nw))
    
    allocate(temp_data(nw_tmp+2)) ! for interpolation
    allocate(temp_wavelength(nw_tmp+2))
    
    ! We do a constant extrapolation
    ! beyond the aerosol data, if it is necessary.
    ! If extrapolation happens, we print a warning
    ! if (wavl(1) < wavl_tmp(1)) then
    !   print*,"Warning: Constantly extrapolating aerosol data file"//trim(filename)// &
    !           " to lower wavelengths than there is data."
    ! endif
    ! if (wavl(nw+1) > wavl_tmp(nw_tmp)) then
    !   print*,"Warning: Constantly extrapolating aerosol data file"//trim(filename)// &
    !           " to larger wavelengths than there is data."
    ! endif
    do i = 1, nrad_file
      
      ! w0 (single scattering albedo)
      temp_data(1:nw_tmp) = w0_tmp(1:nw_tmp,i)
      temp_wavelength(1:nw_tmp) =  wavl_tmp(1:nw_tmp)
      j = nw_tmp
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, 0.d0, w0_tmp(1,i), ierr)
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, huge(1.d0), w0_tmp(nw_tmp,i), ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      call inter2(nw+1, wavl, w0_file(i,:), &
                  nw_tmp+2, temp_wavelength, temp_data, ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      
      ! qext (The extinction efficiency)
      temp_data(1:nw_tmp) = qext_tmp(1:nw_tmp,i)
      temp_wavelength(1:nw_tmp) =  wavl_tmp(1:nw_tmp)
      j = nw_tmp
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, 0.d0, qext_tmp(1,i), ierr)
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, huge(1.d0), qext_tmp(nw_tmp,i), ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      call inter2(nw+1, wavl, qext_file(i,:), &
                  nw_tmp+2, temp_wavelength, temp_data, ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      
      ! g (The scattering anisotropy or asymmetry factor)
      temp_data(1:nw_tmp) = g_tmp(1:nw_tmp,i)
      temp_wavelength(1:nw_tmp) =  wavl_tmp(1:nw_tmp)
      j = nw_tmp
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, 0.d0, g_tmp(1,i), ierr)
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, huge(1.d0), g_tmp(nw_tmp,i), ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      call inter2(nw+1, wavl, g_file(i,:), &
                  nw_tmp+2, temp_wavelength, temp_data, ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      
    enddo
  
  end subroutine
  
  
end module