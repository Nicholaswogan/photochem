
module clima_radtran_twostream
  use clima_const, only: dp
  implicit none
  private
  public :: two_stream_solar, two_stream_ir

contains
  
  pure subroutine two_stream_solar(nz, tau_in, w0_in, gt_in, u0, Rsfc, &
                                   amean, surface_radiance, fup, fdn)
    integer, intent(in) :: nz
    real(dp), intent(in) :: tau_in(nz)
    real(dp), intent(in) :: w0_in(nz)
    real(dp), intent(in) :: gt_in(nz)
    real(dp), intent(in) :: u0, Rsfc
    real(dp), intent(out) :: amean(nz+1)
    real(dp), intent(out) :: surface_radiance
    real(dp), intent(out) :: fup(nz+1), fdn(nz+1)
    
    ! local
    real(dp) :: tau(nz), w0(nz), gt(nz)
    real(dp) :: gam1(nz), gam2(nz), gam3(nz), gam4(nz)
    real(dp) :: lambda(nz), cap_gam(nz)
    real(dp) :: e1(nz), e2(nz), e3(nz), e4(nz)
    real(dp) :: tauc(nz+1), direct(nz+1)
    real(dp) :: cp0(nz), cpb(nz), cm0(nz), cmb(nz)
    real(dp) :: A(nz*2), B(nz*2), D(nz*2), E(nz*2)
    real(dp) :: y1(nz), y2(nz)
    
    integer :: i, l
    real(dp) :: wrk_real, facp, facm, et0, etb, denom, Ssfc
    
    real(dp), parameter :: u1 = 1.0_dp/sqrt(3.0_dp) ! (Quadrature)
    real(dp), parameter :: Fs_pi = 1.0_dp
    
    ! Delta-Eddington scaling (Joseph et al. 1976)
    tau = tau_in*(1.0_dp-w0_in*gt_in*gt_in)
    w0 = w0_in*(1.0_dp-gt_in*gt_in)/(1.0_dp-w0_in*gt_in*gt_in)
    gt = gt_in/(1.0_dp+gt_in)

    ! Quadrature Two-Stream coefficients (Table 1)
    gam1 = sqrt(3.0_dp)*(2.0_dp-w0*(1+gt))/2.0_dp
    gam2 = sqrt(3.0_dp)*w0*(1.0_dp-gt)/2.0_dp
    gam3 = (1.0_dp-sqrt(3.0_dp)*gt*u0)/2.0_dp
    gam4 = 1.0_dp - gam3
    ! u1 = 1.0_dp/sqrt(3.0_dp) (parameter)

    ! lambda, and capital Gamma (Equations 21, 22)
    lambda = sqrt(gam1**2.0_dp - gam2**2.0_dp)
    cap_gam = gam2 / (gam1 + lambda)
    ! cap_gam = (gam1-lambda)/gam2 ! this is the same as above
    
    ! e's (Equation 44)
    do i = 1,nz
      wrk_real = exp(-lambda(i)*tau(i))
      e1(i) = 1.0_dp + cap_gam(i)*wrk_real
      e2(i) = 1.0_dp - cap_gam(i)*wrk_real
      e3(i) = cap_gam(i) + wrk_real
      e4(i) = cap_gam(i) - wrk_real
    enddo
    
    ! tauc - cumulative optical depth at the top of each layer
    tauc(1) = 0.0_dp
    do i = 2,nz+1
      tauc(i) = tauc(i-1) + tau(i-1)
    enddo
    
    ! C+ and C- (Equation 23, 24)
    ! Fs_pi = 1.0_dp (parameter) 
    ! We take the solar flux at the top of the atmosphere to be = 1
    ! Note: Toon et al. 1989 defines Fs * pi = solar flux. So Fs = solar flux / pi
    direct(1) = u0*Fs_pi
    do i=1,nz
      
      facp = w0(i)*Fs_pi*((gam1(i)-1.0_dp/u0)*gam3(i)+gam4(i)*gam2(i))
      facm = w0(i)*Fs_pi*((gam1(i)+1.0_dp/u0)*gam4(i)+gam2(i)*gam3(i)) 
      et0 = exp(-tauc(i)/u0)
      etb = et0*exp(-tau(i)/u0)!*pi
      denom = lambda(i)**2.0_dp - 1.0_dp/(u0**2.0_dp)

      direct(i+1) = u0*Fs_pi*etb
      cp0(i) = et0*facp/denom
      cpb(i) = etb*facp/denom
      cm0(i) = et0*facm/denom
      cmb(i) = etb*facm/denom
    enddo

    Ssfc = Rsfc*direct(nz+1)
      
    ! Coefficients of tridiagonal linear system (Equations 39 - 43)
    ! Odd coeficients (Equation 41)
    A(1) = 0.0_dp
    B(1) = e1(1)
    D(1) = -e2(1)
    E(1) = 0.0_dp - cm0(1) ! assumes no downward diffuse flux at top-of-atmosphere
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
      E(l) = e2(i+1)*(cp0(i+1) - cpb(i)) - e4(i+1)*(cm0(i+1) - cmb(i))
    enddo
    l = 2*nz
    A(l) = e1(nz) - Rsfc*e3(nz)
    B(l) = e2(nz) - Rsfc*e4(nz)
    D(l) = 0.0_dp
    E(l) = Ssfc - cpb(nz) + Rsfc*cmb(nz)
    
    ! Solve tridiagonal system. e is solution
    call tridiag(nz*2,a,b,d,e)
    
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
    amean(1) = (1.0_dp/u1)*(y1(1)*e3(1)-y2(1)*e4(1)+cp0(1)) + direct(1)/u0
    do i = 1, nz
      ! J_n*4*pi = amean (Equation 49)
      amean(i+1) = (1.0_dp/u1)*(y1(i)*(e1(i)+e3(i))+y2(i)*(e2(i)+e4(i)) &
                   + cpb(i)+cmb(i)) + direct(i+1)/u0
    enddo

    ! Equations 31 and 32. But direct flux added
    fup(1) = ((y1(1)*e3(1)-y2(1)*e4(1))+cp0(1))
    fdn(1) = direct(1)
    do i = 1 , NZ
       fup(i+1) = (y1(i)*e1(i)+y2(i)*e2(i)+cpb(i))
       fdn(i+1) = (y1(i)*e3(i)+y2(i)*e4(i)+cmb(i)) + direct(i+1)
    enddo
    
    ! surface radiance (Ranjan and Sasselov 2017). photons actually hitting the ground (i think)
    i = nz
    surface_radiance = (y1(i)*e3(i)+y2(i)*e4(i)+cmb(i))/u1 + exp(-tauc(i+1)/u0)

  end subroutine
  
  subroutine two_stream_ir(nz, tau, w0, gt, emissivity, has_hard_surface, tau_min, bplanck, &
                           fup, fdn)
    use clima_const, only: pi
    integer, intent(in) :: nz
    real(dp), intent(in) :: tau(nz)
    real(dp), intent(in) :: w0(nz)
    real(dp), intent(in) :: gt(nz)
    real(dp), intent(in) :: emissivity
    logical, intent(in) :: has_hard_surface
    real(dp), intent(in) :: tau_min
    real(dp), intent(in) :: bplanck(nz+1)
    real(dp), intent(out) :: fup(nz+1), fdn(nz+1)
    
    ! local
    real(dp) :: gam1(nz), gam2(nz)
    real(dp) :: lambda(nz), cap_gam(nz)
    real(dp) :: e1(nz), e2(nz), e3(nz), e4(nz)
    real(dp) :: cp0(nz), cpb(nz), cm0(nz), cmb(nz)
    real(dp) :: A(nz*2), B(nz*2), D(nz*2), E(nz*2)
    real(dp) :: y1(nz), y2(nz)
    
    integer :: i, l
    real(dp) :: wrk_real, Ssfc, Rsfc
    real(dp) :: b0n, b1n
    real(dp) :: b_avg, b1_bot

    real(dp), parameter :: u1 = 0.5_dp ! (Hemispheric mean)
    real(dp), parameter :: norm = 2.0_dp*pi*u1
    
    ! .true. = terrestrial hard surface, .false. = gas-giant no hard surface.
    if (has_hard_surface) then
      Rsfc = 1.0_dp - emissivity
    else
      Rsfc = 0.0_dp
    endif

    ! no Delta-Eddington scaling in ir

    ! Hemispheric mean Two-Stream coefficients (Table 1)
    gam1 = 2.0_dp - w0*(1.0_dp + gt)
    gam2 = w0*(1.0_dp - gt)
    ! u1 = 0.5_dp (parameter)

    ! lambda, and capital Gamma (Equations 21, 22)
    lambda = sqrt(gam1**2.0_dp - gam2**2.0_dp)
    cap_gam = gam2 / (gam1 + lambda)
    ! cap_gam = (gam1-lambda)/gam2 ! this is the same as above
    
    ! e's (Equation 44)
    do i = 1,nz
      wrk_real = exp(-lambda(i)*tau(i))
      e1(i) = 1.0_dp + cap_gam(i)*wrk_real
      e2(i) = 1.0_dp - cap_gam(i)*wrk_real
      e3(i) = cap_gam(i) + wrk_real
      e4(i) = cap_gam(i) - wrk_real
    enddo
    
    ! C+ and C- (Equation 27)
    ! norm = 2.0_dp*pi*u1 (parameter)
    do i=1,nz
      if (tau(i) <= tau_min) then
        ! NOTE: Toon et al. effectively use a linear-in-optical-depth source
        ! function slope b1n = (B_{i+1}-B_i)/tau(i). For very optically thin layers
        ! this is numerically unstable (0/0 or noise amplification), so we instead
        ! take the tau->0 limit: constant source B ≈ mean(B_i,B_{i+1}) with zero slope.
        b_avg = 0.5_dp*(bplanck(i) + bplanck(i+1))
        b0n = b_avg
        b1n = 0.0_dp
      else
        b0n = bplanck(i)
        b1n = (bplanck(i+1) - b0n)/tau(i)
      endif
      
      cp0(i) = norm*(b0n + b1n*(1.0_dp/(gam1(i)+gam2(i))))
      cpb(i) = norm*(b0n + b1n*(tau(i) + 1.0_dp/(gam1(i)+gam2(i))))
      cm0(i) = norm*(b0n + b1n*(-1.0_dp/(gam1(i)+gam2(i))))
      cmb(i) = norm*(b0n + b1n*(tau(i) - 1.0_dp/(gam1(i)+gam2(i))))

    enddo

    if (has_hard_surface) then
      Ssfc = emissivity*pi*bplanck(nz+1) ! ground
    else
      ! PICASO-like no-hard-surface thermal lower BC:
      ! B_surface = pi * (B_bottom + mu1 * dB/dtau), with mu1=u1=0.5.
      if (tau(nz) <= tau_min) then
        b1_bot = 0.0_dp
      else
        b1_bot = (bplanck(nz+1) - bplanck(nz))/tau(nz)
      endif
      Ssfc = pi*(bplanck(nz+1) + u1*b1_bot)
    endif
      
    ! Coefficients of tridiagonal linear system (Equations 39 - 43)
    ! Odd coeficients (Equation 41)
    A(1) = 0.0_dp
    B(1) = e1(1)
    D(1) = -e2(1)
    E(1) = 0.0_dp - cm0(1) ! assumes no downward diffuse flux at top-of-atmosphere
    do i = 1, nz-1
      l = 2*i + 1
      A(l) = e2(i)*e3(i) - e4(i)*e1(i)
      B(l) = e1(i)*e1(i+1) - e3(i)*e3(i+1)
      D(l) = e3(i)*e4(i+1) - e1(i)*e2(i+1)
      E(l) = e3(i)*(cp0(i+1) - cpb(i)) + e1(i)*(cmb(i) - cm0(i+1))
    enddo
    
    ! Even coefficients (Equation 42)
    do i = 1, nz-1
      l = 2*i
      A(l) = e2(i+1)*e1(i) - e3(i)*e4(i+1)
      B(l) = e2(i)*e2(i+1) - e4(i)*e4(i+1)
      D(l) = e1(i+1)*e4(i+1) - e2(i+1)*e3(i+1)
      E(l) = e2(i+1)*(cp0(i+1) - cpb(i)) - e4(i+1)*(cm0(i+1) - cmb(i))
    enddo
    l = 2*nz
    A(l) = e1(nz) - Rsfc*e3(nz)
    B(l) = e2(nz) - Rsfc*e4(nz)
    D(l) = 0.0_dp
    E(l) = Ssfc - cpb(nz) + Rsfc*cmb(nz)
    
    ! Solve tridiagonal system. e is solution
    call tridiag(nz*2,a,b,d,e)
    
    ! unpack solution
    do i = 1, nz
      l = 2*i
      y1(i) = e(l-1)
      y2(i) = e(l)
    enddo
    
    ! Equations 31 and 32.
    fup(1) = ((y1(1)*e3(1)-y2(1)*e4(1))+cp0(1))
    fdn(1) = 0.0_dp
    do i = 1,nz
       fup(i+1) = (y1(i)*e1(i)+y2(i)*e2(i)+cpb(i))
       fdn(i+1) = (y1(i)*e3(i)+y2(i)*e4(i)+cmb(i))
    enddo
    
  end subroutine
  
  pure subroutine tridiag(n,a,b,c,d)
    integer,intent(in) :: n
    real(dp), intent(inout) :: a(n), b(n), c(n), d(n)

    integer :: i
    
    c(1) = c(1)/b(1)
    d(1) = d(1)/b(1)
    
    do i= 2, n-1
      c(i) = c(i)/(b(i) - a(i)*c(i-1))
      d(i) = (d(i) - a(i)*d(i-1)) / (b(i) - a(i)*c(i-1))
    enddo
    
    d(n) = (d(n) - a(n)*d(n-1)) / (b(n) - a(n)*c(n-1))

    do i = n-1,1,-1
      d(i) = d(i) - c(i)*d(i+1)
    enddo
  end subroutine
  
end module
