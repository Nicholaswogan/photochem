module photochem_setup
  implicit none
  private
  
contains
  
  subroutine setup_initial_atmosphere()
    
    
    ! first allocate memory

    ! set up the atmosphere grid
    allocate(z(nz))
    allocate(photoset%dz(photoset%nz))
    call vertical_grid(photoset%bottom_atmos,photoset%top_atmos,photoset%nz, &
                       photoset%z, photoset%dz)
    
    ! compute the gravity
    allocate(photoset%grav(photoset%nz))
    do i =1,photoset%nz
      call compute_gravity((photoset%planet_radius + photoset%z(i))/1.d2, &
                            photoset%planet_mass/1.d3, grav)
      photoset%grav(i) = grav*1.d2 ! convert to cgs
    enddo    
    
  end subroutine
  
  subroutine compute_gravity(radius, mass, grav)
    use photochem_const, only: G_grav
    real(real_kind), intent(in) :: radius, mass
    real(real_kind), intent(out) :: grav
    grav = G_grav * mass / radius**2.d0
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