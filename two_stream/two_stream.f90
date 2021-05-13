
module radtran
  implicit none
  private
  public two_stream
  integer,parameter :: real_kind = kind(1.0d0)

contains
  
  subroutine two_stream(nz, tau, w0, u0, Rsfc, amean, surface_radiance, method)
    real(real_kind), parameter :: pi = 3.14159
    integer, intent(in) :: nz
    real(8), intent(in) :: tau(nz)
    real(8), intent(in) :: w0(nz)
    real(8), intent(in) :: u0, Rsfc
    integer, intent(in), optional :: method
    real(8), intent(out) :: amean(nz+1)
    real(8), intent(out) :: surface_radiance
    
    ! real(real_kind) :: intensity(nz)
    
    ! local
    real(real_kind) :: gt(nz)
    real(real_kind) :: gam1(nz), gam2(nz), gam3(nz), gam4(nz), u1
    real(real_kind) :: lambda(nz), cap_gam(nz)
    real(real_kind) :: e1(nz), e2(nz), e3(nz), e4(nz)
    real(real_kind) :: tauc(nz+1), direct(nz+1)
    real(real_kind) :: cp0(nz), cpb(nz), cm0(nz), cmb(nz)
    real(real_kind) :: A(nz*2), B(nz*2), D(nz*2), E(nz*2)
    real(real_kind) :: y1(nz), y2(nz)! fup(nz+1), fdn(nz+1), net_flux(nz+1)
    
    integer :: i, l
    real(real_kind) :: wrk_real, facp, facm, et0, etb, denom, Ssfc, fs_pi
    integer ::  info

    gt = 0.0d0 ! asymetry factor. Zero for now
    
    i = 0
    if (present(method)) i = method
    if (i == 1) then
      ! Eddington Two-Stream coefficients (Table 1)
      gam1 = (7.d0 - w0*(4.d0 + 3.d0*gt))/4.d0
      gam2 = -(1.d0- w0*(4.d0 - 3.d0*gt))/4.d0
      gam3 = (2.d0 - 3.d0*gt*u0)/4.d0
      gam4 = 1.d0 - gam3
      u1 = 0.5d0
    else ! default
      ! Quadrature Two-Stream coefficients (Table 1)
      gam1 = dsqrt(3.d0)*(2.d0-w0*(1+gt))/2.d0
      gam2 = dsqrt(3.d0)*w0*(1.d0-gt)/2.d0
      gam3 = (1.d0-dsqrt(3.d0)*gt*u0)/2.d0
      gam4 = 1.d0 - gam3
      u1 = 1.d0/dsqrt(3.d0)
    endif

    ! lambda, and capital Gamma (Equations 21, 22)
    lambda = (gam1**2.d0 - gam2**2.d0)**(0.5d0)
    cap_gam = gam2 / (gam1 + lambda)
    
    ! e's (Equation 44)
    do i = 1,nz
      wrk_real = dexp(-lambda(i)*tau(i))
      e1(i) = 1.d0 + cap_gam(i)*wrk_real
      e2(i) = 1.d0 - cap_gam(i)*wrk_real
      e3(i) = cap_gam(i) + wrk_real
      e4(i) = cap_gam(i) - wrk_real
    enddo
    
    ! tauc - cumulative optical depth at the top of each layer
    tauc(1) = 0.d0
    do i = 2,nz+1
      tauc(i) = tauc(i-1) + tau(i-1)
    enddo
    
    ! C+ and C- (Equation 23, 24)
    Fs_pi = 1.d0 ! We take the solar flux at the top of the atmosphere to be = 1
    ! Note: Toon et al. 1989 defines Fs * pi = solar flux. So Fs = solar flux / pi
    direct(1) = u0*Fs_pi
    do i=1,nz
      facp = w0(i)*Fs_pi*((gam1(i)-1.d0/u0)*gam3(i)+gam4(i)*gam2(i))
      facm = w0(i)*Fs_pi*((gam1(i)+1.d0/u0)*gam4(i)+gam2(i)*gam3(i)) 
      et0 = dexp(-tauc(i)/u0)
      etb = et0*dexp(-tau(i)/u0)!*pi
      denom = lambda(i)**2.d0 - 1.d0/u0**2.d0
      
      direct(i+1) = u0*Fs_pi*etb
      cp0(i) = et0*facp/denom
      cpb(i) = etb*facp/denom
      cm0(i) = et0*facm/denom
      cmb(i) = etb*facm/denom
    enddo
    Ssfc = Rsfc*direct(nz+1)
      
    ! Coefficients of tridiagonal linear system (Equations 39 - 43)
    ! Odd coeficients (Equation 41)
    A(1) = 0.d0
    B(1) = e1(1)
    D(1) = -e2(1)
    E(1) = 0.d0 - cm0(1) ! assumes no downward diffuse flux at top-of-atmosphere
    do i = 1, nz-1
      l = 2*i + 1 ! is this right?
      A(l) = e2(i)*e3(i) - e4(i)*e1(i)
      B(l) = e1(i)*e1(i+1) - e3(i)*e3(i+1)
      D(l) = e3(i)*e4(i+1) - e1(i)*e2(i+1)
      E(l) = e3(i)*(cp0(i+1) - cpb(i)) + e1(i)*(cmb(i) - cm0(i+1))
    enddo
    
    ! Even coefficients (Equation 42)
    do i = 1, nz-1
      l = 2*i ! is this right?
      A(l) = e2(i+1)*e1(i) - e3(i)*e4(i+1)
      B(l) = e2(i)*e2(i+1) - e4(i)*e4(i+1)
      D(l) = e1(i+1)*e4(i+1) - e2(i+1)*e3(i+1)
      E(l) = e2(i+1)*(cp0(i+1) - cpb(i)) + e4(i+1)*(cm0(i+1) - cmb(i))
    enddo
    l = 2*nz
    A(l) = e1(nz) - Rsfc*e3(nz)
    B(l) = e2(nz) - Rsfc*e4(nz)
    D(l) = 0.d0
    E(l) = Ssfc - cpb(nz) + Rsfc*cmb(nz)
    
    ! Solve tridiagonal system. e is solution
    call dgtsv(nz*2, 1, a(2:nz*2), b, d(1:nz*2-1), e, nz*2, info)
    if (info /= 0) then
      print*,info
      ! ierr = 1
    endif
    
    ! unpack solution
    do i = 1, nz
      l = 2*i
      y1(i) = e(l-1)
      y2(i) = e(l)
    enddo

    ! amean = integral(J_n d_Omega) = J_n*4*pi
    ! Above is an integration of J_n over a complete sphere, 
    ! and is thus has units of flux density, or irradiance.
    ! amean(i) is for at top of layer i. amean(nz + 1) is the ground.
    
    ! very top edge of atmosphere (Not in paper. Derive from Equation 17 and 31. bit confusing)
    amean(1) = (1.d0/u1)*(y1(1)*e3(1)-y2(1)*e4(1)+cp0(1)) + direct(1)/u0
    do i = 1, nz
      ! J_n*4*pi = amean (Equation 49)
      amean(i+1) = (1.d0/u1)*(y1(i)*(e1(i)+e3(i))+y2(i)*(e2(i)+e4(i)) &
                   + cpb(i)+cmb(i)) + direct(i+1)/u0
    enddo
    
    ! net_flux(1) = y1(1)*e3(1) - y2(1)*e4(1)+cp0(1) - direct(1)
    ! do i = 1, nz
    !   net_flux(i+1) = y1(i)*(e1(i) - e3(i)) + y2(i)*(e2(i) - e4(i)) &
    !                   + cpb(i) - cmb(i) - direct(i+1)
    ! enddo
    
    ! Equations 31 and 32. But direct flux added
    ! fup(1) = ((y1(1)*e3(1)-y2(1)*e4(1))+cp0(1))
    ! fdn(1) = direct(1)
    ! do i = 1 , NZ
    !    fup(i+1) = (y1(i)*e1(i)+y2(i)*e2(i)+cpb(i))
    !    fdn(i+1) = (y1(i)*e3(i)+y2(i)*e4(i)+cmb(i)) + direct(i+1)
    ! enddo
    
    ! surface radiance (Ranjan and Sasselov 2017). photons actually hitting the ground (i think)
    i = nz
    surface_radiance = (y1(i)*e3(i)+y2(i)*e4(i)+cmb(i))/u1 + dexp(-tauc(i+1)/u0)

    
  end subroutine
  
end module


program main
  use radtran, only: two_stream
  implicit none
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: nz = 1
  real(real_kind) :: tau(nz)
  real(real_kind) :: w0(nz), u0, Rsfc, amean(nz+1), surface_radiance
  
  ! u0 = dcos(3.14159d0/3.d0)
  u0 = 0.7d0
  tau = 0.1d0
  w0  = 0.9999d0
  Rsfc = 0.25d0  
  
  call two_stream(nz, tau, w0, u0, Rsfc, amean, surface_radiance)
  
  print*,amean
  
  print*,surface_radiance

end program