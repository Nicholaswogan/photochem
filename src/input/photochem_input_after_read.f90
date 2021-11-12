submodule (photochem_input) photochem_input_after_read
  implicit none
  
contains
  
  module subroutine after_read_setup(photodata, photovars, err)
    use photochem_eqns, only: vertical_grid, gravity
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(inout) :: photovars
    character(len=err_len), intent(out) :: err
    err = ""
    
    photodata%kd = 2*photodata%nq + 1
    photodata%kl = photodata%kd + photodata%nq
    photodata%ku = photodata%kd - photodata%nq
    photodata%lda = 3*photodata%nq + 1
    
    call allocate_nz_vars(photodata, photovars)
    ! set up the atmosphere grid
    call vertical_grid(photovars%bottom_atmos, photovars%top_atmos, &
                       photovars%nz, photovars%z, photovars%dz)
    if (photodata%fix_water_in_trop) then
      photovars%trop_ind = minloc(photovars%z,1, &
                          photovars%z .ge. photovars%trop_alt) - 1
    endif
    call gravity(photodata%planet_radius, photodata%planet_mass, &
                 photovars%nz, photovars%z, photovars%grav)
    call interp2atmosfile(photodata, photovars, err)
    if (len_trim(err) /= 0) return
    
    ! lets do xsections
    call interp2xsdata(photodata, photovars, err)
    if (len_trim(err) /= 0) return
    
  end subroutine
  
  subroutine interp2atmosfile(dat, var, err)
    use futils, only: interp
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(len=err_len), intent(out) :: err
    
    integer :: i
    
    err = ''
    
    call interp(var%nz, dat%nzf, var%z, dat%z_file, dat%T_file, var%Temperature, err)
    if (len_trim(err) /= 0) return
    
    call interp(var%nz, dat%nzf, var%z, dat%z_file, dlog10(dabs(dat%edd_file)), var%edd, err)
    if (len_trim(err) /= 0) return
    var%edd = 10.d0**var%edd
    
    do i = 1,dat%nq
      call interp(var%nz, dat%nzf, var%z, dat%z_file,&
                  dlog10(dabs(dat%usol_file(i,:))), var%usol_init(i,:), err)
      if (len_trim(err) /= 0) return
    enddo
    var%usol_init = 10.d0**var%usol_init
    
    if (dat%there_are_particles) then
      do i = 1,dat%npq
        call interp(var%nz, dat%nzf, var%z, dat%z_file, &
                    log10(abs(dat%particle_radius_file(i,:))), var%particle_radius(i,:), err)
        if (len_trim(err) /= 0) return
      enddo
      var%particle_radius = 10.d0**var%particle_radius
    endif
    
    if (var%z(1) < dat%z_file(1)) then
      print*,'Warning: vertical grid is being extrapolated below where there is input data.'
    endif
    
    if (var%z(var%nz) > dat%z_file(dat%nzf)) then
      print*,'Warning: vertical grid is being extrapolated above where there is input data.'
    endif

  end subroutine
  
  subroutine interp2xsdata(dat, var, err)
    use futils, only: interp
    use photochem_const, only: smaller_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var

    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k, jj
    real(real_kind) :: val(1), T_temp(1)
    real(real_kind) :: dr, slope, intercept
    
    err = ''

    do k = 1, dat%nw
      do i = 1,dat%kj
        do j = 1,var%nz
          T_temp(1) = var%temperature(j)
    
          call interp(1, dat%xs_data(i)%n_temps, T_temp, dat%xs_data(i)%xs_temps, &
                      log10(abs(dat%xs_data(i)%xs(:,k))+smaller_real), val, err)
                      
          var%xs_x_qy(j,i,k) = 10.d0**val(1)
        enddo
      enddo
    enddo
    
    ! particles
    if (dat%there_are_particles) then
      do j = 1,var%nz
        do k = 1,dat%np
          if (var%particle_radius(k,j) <= dat%radii_file(1,k)) then
            err = "There is not any optical data for the "// &
                  "particle radii specified in the atmosphere."
            return
          endif
          if (var%particle_radius(k,j) >= dat%radii_file(dat%nrad_file,k)) then
            err = "There is not any optical data for the "// &
                  "particle radii specified in the atmosphere."
            return
          endif
        enddo
      enddo
      do i = 1,dat%nw
        do j = 1,var%nz
          do k = 1,dat%np
            do jj = 1,dat%nrad_file-1
              if (var%particle_radius(k,j) >= dat%radii_file(jj,k) .and. &
                  var%particle_radius(k,j) < dat%radii_file(jj+1,k)) then
    
                dr = dat%radii_file(jj+1,k) - dat%radii_file(jj,k)
    
                slope = (dat%w0_file(jj+1,k,i) - dat%w0_file(jj,k,i))/dr
                intercept = dat%w0_file(jj,k,i) - dat%radii_file(jj,k)*slope
                var%w0_particles(k,j,i) = slope*var%particle_radius(k,j) + intercept
    
                slope = (dat%qext_file(jj+1,k,i) - dat%qext_file(jj,k,i))/dr
                intercept = dat%qext_file(jj,k,i) - dat%radii_file(jj,k)*slope
                var%qext_particles(k,j,i) = slope*var%particle_radius(k,j) + intercept
    
                slope = (dat%g_file(jj+1,k,i) - dat%g_file(jj,k,i))/dr
                intercept = dat%g_file(jj,k,i) - dat%radii_file(jj,k)*slope
                var%gt_particles(k,j,i) = slope*var%particle_radius(k,j) + intercept
              endif
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine
  
  subroutine allocate_nz_vars(photodata, vars)
    type(PhotochemData), intent(in) :: photodata
    type(PhotochemVars), intent(inout) :: vars
    
    vars%neqs = photodata%nq*vars%nz

    allocate(vars%temperature(vars%nz))
    allocate(vars%z(vars%nz))
    allocate(vars%dz(vars%nz))
    allocate(vars%edd(vars%nz))
    allocate(vars%grav(vars%nz))
    allocate(vars%usol_init(photodata%nq,vars%nz))
    allocate(vars%particle_radius(photodata%npq,vars%nz))
    allocate(vars%xs_x_qy(vars%nz,photodata%kj,photodata%nw))
    allocate(vars%usol_out(photodata%nq,vars%nz))
    allocate(vars%w0_particles(photodata%np,vars%nz,photodata%nw))
    allocate(vars%qext_particles(photodata%np,vars%nz,photodata%nw))
    allocate(vars%gt_particles(photodata%np,vars%nz,photodata%nw))

  end subroutine
  
end submodule