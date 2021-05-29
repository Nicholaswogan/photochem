module photochem_setup
  implicit none
  private
  integer, private, parameter :: real_kind = kind(1.0d0)
  integer, private, parameter :: err_len = 1024
  
  public :: setup
  
contains
  
  subroutine setup(mechanism_file, settings_file, flux_file, atmosphere_txt, err)
    use photochem_data, only: setup_files, &
                              planet_radius, planet_mass, nq, kj, nw
    use photochem_vars, only: bottom_atmos, top_atmos, nz, &
                              z, dz, grav, temperature, edd, usol_init, &
                              xs_x_qy
    
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(len=err_len), intent(out) :: err
    
    integer :: i
    
    call setup_files(mechanism_file, settings_file, flux_file, atmosphere_txt, err)
    if (len_trim(err) /= 0) return
    
    ! allocate nz vars
    call allocate_nz_vars()
    ! set up the atmosphere grid
    call vertical_grid(bottom_atmos, top_atmos, nz, z, dz)
    call gravity(planet_radius, planet_mass, nz, z, grav)
    call interp2atmosfile(nz, nq, z, temperature, edd, usol_init, err)
    if (len_trim(err) /= 0) return
    
    ! lets do xsections
    call interp2xsdata(nz, kj, nw, temperature, xs_x_qy, err)
    if (len_trim(err) /= 0) return

  end subroutine
  
  subroutine interp2xsdata(nz,kj,nw,temperature,xs_x_qy, err)
    use photochem_const, only: small_real
    use photochem_data, only: num_temp_cols, sum_temp_cols, &
                              xs_data, xs_data_temps
    
    integer, intent(in) :: nz, kj, nw
    real(real_kind), intent(in) :: temperature(nz)
    real(real_kind), intent(out) :: xs_x_qy(nz,kj,nw)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k, l, m, mm, ncol
    real(real_kind) :: val(1), T_temp(1)
    real(real_kind) ,allocatable :: tmp(:)
    
    allocate(tmp(size(xs_data_temps,1)))
    err = ''
    
    do j = 1,nz
      T_temp(1) = temperature(j)
      do i = 1,kj
        ncol = num_temp_cols(i)
    
        do k = 1, nw
          do l = 1, ncol
            m = ((l-1)*nw + 1) + (sum_temp_cols(i)*nw)
            tmp(l) = xs_data(m+k-1+(l-1)*nw)
          enddo
          ! print*,dlog10(tmp(1:ncol))
          call interp(1, ncol, T_temp, xs_data_temps(1:ncol,i), dlog10(dabs(tmp(1:ncol)+small_real)), val, err)
          if (len_trim(err) /= 0) return
          xs_x_qy(j,i,k) = 10.d0**val(1)
        enddo
      enddo
    enddo

  end subroutine
  
  subroutine interp2atmosfile(nz, nq, z, T, edd, usol, err)
    use photochem_data, only: nzf, z_file, &
                              T_file, edd_file, usol_file
    integer, intent(in) :: nz, nq
    real(real_kind), intent(in) :: z(nz)
    real(real_kind), intent(out) :: T(nz), edd(nz), usol(nq,nz)
    character(len=err_len), intent(out) :: err
    
    integer :: i
    
    err = ''
    
    call interp(nz, nzf, z, z_file, T_file, T, err)
    if (len_trim(err) /= 0) return
    
    call interp(nz, nzf, z, z_file, dlog10(dabs(edd_file)), edd, err)
    if (len_trim(err) /= 0) return
    edd = 10.d0**edd
    
    do i = 1,nq
      call interp(nz, nzf, z, z_file, dlog10(dabs(usol_file(i,:))), usol(i,:), err)
      if (len_trim(err) /= 0) return
    enddo
    usol = 10.d0**usol
    
    if (z(1) < z_file(1)) then
      print*,'Warning: vertical grid is being extrapolated below where there is input data.'
    endif
    
    if (z(nz) > z_file(nzf)) then
      print*,'Warning: vertical grid is being extrapolated above where there is input data.'
    endif
    
    
  end subroutine
  
  
  
  subroutine allocate_nz_vars()
    use photochem_data, only: nq, kj, nw
    use photochem_vars
    use photochem_wrk
    neqs = nq*nz
    
    if (allocated(temperature)) then
      deallocate(temperature)
      deallocate(z)
      deallocate(dz)
      deallocate(edd)
      deallocate(grav)
      deallocate(usol_init)
      deallocate(xs_x_qy)
    endif
    
    allocate(temperature(nz))
    allocate(z(nz))
    allocate(dz(nz))
    allocate(edd(nz))
    allocate(grav(nz))
    allocate(usol_init(nq,nz))
    allocate(xs_x_qy(nz,kj,nw))

  end subroutine
  

  subroutine gravity(radius, mass, nz, z, grav)
    use photochem_const, only: G_grav
    real(real_kind), intent(in) :: radius, mass ! radius in cm, mass in grams
    integer, intent(in) :: nz
    real(real_kind), intent(in) :: z(nz) ! cm
    real(real_kind), intent(out) :: grav(nz) ! cm/s2

    integer :: i
    
    do i = 1, nz              
      grav(i) = G_grav * (mass/1.d3) / ((radius + z(i))/1.d2)**2.d0
      grav(i) = grav(i)*1.d2 ! convert to cgs
    enddo 
    
  end subroutine
  
  subroutine vertical_grid(bottom, top, nz, z, dz)
    real(real_kind), intent(in) :: bottom, top
    integer, intent(in) :: nz
    real(real_kind), intent(out) :: z(nz), dz(nz)
  
    integer :: i
  
    dz = (top - bottom)/nz
    z(1) = dz(1)/2.d0
    do i = 2,nz
      z(i) = z(i-1) + dz(i)
    enddo
  end subroutine
  
  
end module