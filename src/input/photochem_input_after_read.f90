submodule (photochem_input) photochem_input_after_read
  implicit none
  
contains
  
  module subroutine after_read_setup(photodata, photovars, err)
    use photochem_eqns, only: vertical_grid, gravity
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(inout) :: photovars
    character(:), allocatable, intent(out) :: err
    
    photodata%kd = 2*photodata%nq + 1
    photodata%kl = photodata%kd + photodata%nq
    photodata%ku = photodata%kd - photodata%nq
    photodata%lda = 3*photodata%nq + 1
    
    call allocate_nz_vars(photodata, photovars)
    ! set up the atmosphere grid
    call vertical_grid(photovars%bottom_atmos, photovars%top_atmos, &
                       photovars%nz, photovars%z, photovars%dz)
    call gravity(photodata%planet_radius, photodata%planet_mass, &
                 photovars%nz, photovars%z, photovars%grav)
    call interp2atmosfile(photodata, photovars, err)
    if (allocated(err)) return
    
    ! all below depends on Temperature
    call interp2xsdata(photodata, photovars, err)
    if (allocated(err)) return
    
    if (photodata%reverse) then
      call compute_gibbs_energy(photodata, photovars, err)
      if (allocated(err)) return
    endif
    
    if (photodata%fix_water_in_trop .or. photodata%gas_rainout) then
      ! we have a tropopause
      photovars%trop_ind = minloc(photovars%z,1, &
                          photovars%z .ge. photovars%trop_alt) - 1
    else
      photovars%trop_ind = 0
    endif
    
  end subroutine
  
  subroutine interp2atmosfile(dat, var, err)
    use futils, only: interp
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    integer :: i
    
    
    call interp(var%nz, dat%nzf, var%z, dat%z_file, dat%T_file, var%Temperature, err)
    if (allocated(err)) return
    
    call interp(var%nz, dat%nzf, var%z, dat%z_file, dlog10(dabs(dat%edd_file)), var%edd, err)
    if (allocated(err)) return
    var%edd = 10.0_dp**var%edd
    
    do i = 1,dat%nq
      call interp(var%nz, dat%nzf, var%z, dat%z_file,&
                  dlog10(dabs(dat%usol_file(i,:))), var%usol_init(i,:), err)
      if (allocated(err)) return
    enddo
    var%usol_init = 10.0_dp**var%usol_init
    
    if (dat%there_are_particles) then
      do i = 1,dat%npq
        call interp(var%nz, dat%nzf, var%z, dat%z_file, &
                    log10(abs(dat%particle_radius_file(i,:))), var%particle_radius(i,:), err)
        if (allocated(err)) return
      enddo
      var%particle_radius = 10.0_dp**var%particle_radius
    endif
    
    if (var%z(1) < dat%z_file(1)) then
      print*,'Warning: vertical grid is being extrapolated below where there is input data.'
    endif
    
    if (var%z(var%nz) > dat%z_file(dat%nzf)) then
      print*,'Warning: vertical grid is being extrapolated above where there is input data.'
    endif

  end subroutine
  
  module subroutine compute_gibbs_energy(dat, var, err)
    use photochem_eqns, only: gibbs_energy_eval
    
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    integer :: i, j
    logical :: found

    
    do i = 1,dat%ng
      do j = 1,var%nz
        call gibbs_energy_eval(dat%thermo_data(i), var%temperature(j), &
                               found, var%gibbs_energy(j,i))
        if (.not. found) then
          err = 'The temperature is not within the ranges '// &
                'given for the thermodynamic data for '//trim(dat%species_names(i+dat%npq))
          return
        endif
      enddo
    enddo

  end subroutine
  
  module subroutine interp2xsdata(dat, var, err)
    use futils, only: interp
    use photochem_const, only: smaller_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var

    character(:), allocatable, intent(out) :: err
    
    integer :: i, j, k, jj
    real(dp) :: val(1), T_temp(1)
    real(dp) :: dr, slope, intercept
    

    do k = 1, dat%nw
      do i = 1,dat%kj
        do j = 1,var%nz
          T_temp(1) = var%temperature(j)
    
          call interp(1, dat%xs_data(i)%n_temps, T_temp, dat%xs_data(i)%xs_temps, &
                      log10(abs(dat%xs_data(i)%xs(:,k))+smaller_real), val, err)
                      
          var%xs_x_qy(j,i,k) = 10.0_dp**val(1)
        enddo
      enddo
    enddo
    
    ! particles
    if (dat%there_are_particles) then
      do j = 1,var%nz
        do k = 1,dat%np
          ! if there is optical data, then check that the
          ! data covers the particle radii in the atmosphere.
          if (dat%part_xs_file(k)%ThereIsData) then
          
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
            
          endif
          
        enddo
      enddo
      do i = 1,dat%nw
        do j = 1,var%nz
          do k = 1,dat%np
            
          ! if there is particle optical data, then linearly interpolate
          ! it to to the particle radii in the atmosphere.
          if (dat%part_xs_file(k)%ThereIsData) then
            
            do jj = 1,dat%nrad_file-1
              if (var%particle_radius(k,j) >= dat%radii_file(jj,k) .and. &
                  var%particle_radius(k,j) < dat%radii_file(jj+1,k)) then
    
                dr = dat%radii_file(jj+1,k) - dat%radii_file(jj,k)
    
                slope = (dat%part_xs_file(k)%w0(jj+1,i) - dat%part_xs_file(k)%w0(jj,i))/dr
                intercept = dat%part_xs_file(k)%w0(jj,i) - dat%radii_file(jj,k)*slope
                var%particle_xs(k)%w0(j,i) = slope*var%particle_radius(k,j) + intercept
                
                slope = (dat%part_xs_file(k)%qext(jj+1,i) - dat%part_xs_file(k)%qext(jj,i))/dr
                intercept = dat%part_xs_file(k)%qext(jj,i) - dat%radii_file(jj,k)*slope
                var%particle_xs(k)%qext(j,i) = slope*var%particle_radius(k,j) + intercept
    
                slope = (dat%part_xs_file(k)%gt(jj+1,i) - dat%part_xs_file(k)%gt(jj,i))/dr
                intercept = dat%part_xs_file(k)%gt(jj,i) - dat%radii_file(jj,k)*slope
                var%particle_xs(k)%gt(j,i) = slope*var%particle_radius(k,j) + intercept
              endif
            enddo
          
          endif 
            
          enddo
        enddo
      enddo
    endif
    
    

  end subroutine
  
  subroutine allocate_nz_vars(photodata, vars)
    type(PhotochemData), intent(in) :: photodata
    type(PhotochemVars), intent(inout) :: vars
    
    integer :: i
    
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
    
    allocate(vars%particle_xs(photodata%np))
    do i = 1,photodata%np
      ! only allocate space if there is data
      if (photodata%part_xs_file(i)%ThereIsData) then
        vars%particle_xs(i)%ThereIsData = .true.
        allocate(vars%particle_xs(i)%w0(vars%nz,photodata%nw))
        allocate(vars%particle_xs(i)%qext(vars%nz,photodata%nw))
        allocate(vars%particle_xs(i)%gt(vars%nz,photodata%nw))
      else
        vars%particle_xs(i)%ThereIsData = .false.
      endif
    enddo
    
    if (photodata%reverse) then
      allocate(vars%gibbs_energy(vars%nz,photodata%ng))
    endif

  end subroutine
  
end submodule