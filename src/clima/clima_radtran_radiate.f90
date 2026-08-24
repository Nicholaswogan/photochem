module clima_radtran_radiate
  use clima_const, only: dp
  implicit none

contains
  
  subroutine radiate(rtc, op, opr, rw, &
                     surface_emissivity, &
                     surface_albedo, has_hard_surface, ir_tau_min, diurnal_fac, photons_sol, &
                     zenith_u, zenith_weights, &
                     T_surface, T, &
                     fup_a, fdn_a, fup_n, fdn_n, amean, tau_band)
    use clima_radtran_types, only: OpticalProperties, OpticalPropertiesResult, &
                                   RTChannel, RadiateWork, RadiateBinWork
    use clima_radtran_types, only: SolarChannel, IRChannel
    use clima_eqns, only: planck_fcn
    use clima_radtran_twostream, only: two_stream_solar, two_stream_ir
    use clima_const, only: plank, c_light
    
    type(RTChannel), intent(in) :: rtc
    type(OpticalProperties), target, intent(in) :: op
    type(OpticalPropertiesResult), target, intent(in) :: opr
    type(RadiateWork), target, intent(inout) :: rw
    
    real(dp), intent(in) :: surface_emissivity(:) !! Surface emissivity (Needed only for IR)

    real(dp), intent(in) :: surface_albedo(:) !! Surface albedo (Needed only for solar)
    logical, intent(in) :: has_hard_surface !! If true, use hard-surface lower thermal BC.
    real(dp), intent(in) :: ir_tau_min !! Thin-layer optical-depth guard for IR two-stream source terms.
    real(dp), intent(in) :: diurnal_fac !! Diurnal averaging factor (0.5) (Needed only for solar)
    real(dp), intent(in) :: photons_sol(:) !! (nw) Average solar flux in each bin (mW/m2/Hz) (Needed only for solar)

    real(dp), intent(in) :: zenith_u(:) !! cos(zenith angles) to average stellar flux over.
                                        !! SHOULD BE SIZE 1 array for IR with any value.
    real(dp), intent(in) :: zenith_weights(:) !! Weights for integrating over zenith angle. 
                                              !! SHOULD BE SIZE 1 array for IR equal to 1.0.
    real(dp), intent(in) :: T_surface !! Surface Tempeature (K) (Needed only for solar)
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(out) :: fup_a(:,:), fdn_a(:,:) !! (nz+1,nw) mW/m2/Hz in each wavelength bin
                                                    !! at the edges of the vertical grid
    real(dp), intent(out) :: fup_n(:), fdn_n(:) !! (nz+1) mW/m2 at the edges of the vertical grid 
                                                !! (integral of fup_a and fdn_a over wavelength grid)
    real(dp), intent(out) :: amean(:,:) !! (nz+1,nw) Mean intensity in photons/cm^2/s (Only for Solar)
    real(dp), intent(out) :: tau_band(:,:) !! (nz,nw) The optical depth of each layer

    integer :: l, nz

    nz = size(opr%tau_band,1)
    
    !$omp parallel private(l)
    !$omp do
    do l = rtc%ind_start,rtc%ind_end; block
      integer :: i, ii, j, n, ll
      real(dp) :: surf_rad, avg_freq, albedo, emissivity
      type(RadiateBinWork), pointer :: rbw

      ll = l - rtc%ind_start + 1

      rbw => rw%rbw(ll)

      ! plank function, only if in the IR
      ! bplanck has units [mW sr^−1 m^−2 Hz^-1]
      if (rtc%channel_type == IRChannel) then
        avg_freq = 0.5_dp*(op%freq(l) + op%freq(l+1))
        rbw%bplanck(nz+1) = planck_fcn(avg_freq, T_surface)
        do j = 1,nz
          n = nz+1-j
          rbw%bplanck(n) = planck_fcn(avg_freq, T(j))
        enddo
        emissivity = surface_emissivity(ll)
        albedo = 0.0_dp
      elseif (rtc%channel_type == SolarChannel) then
        emissivity = 0.0_dp
        albedo = surface_albedo(ll)
      endif

      rbw%fup2 = 0.0_dp
      rbw%fdn2 = 0.0_dp
      if (rtc%channel_type == SolarChannel) then
        rbw%amean2 = 0.0_dp
      endif
      ! Loop over zenith angles and gauss points
      do ii = 1,size(zenith_u)

        rbw%fup1 = 0.0_dp
        rbw%fdn1 = 0.0_dp
        if (rtc%channel_type == SolarChannel) then
          rbw%amean1 = 0.0_dp
        endif

        do i = 1,op%kset%nbin

          if (rtc%channel_type == SolarChannel) then
            call two_stream_solar( &
              nz=nz, &
              tau_in=opr%tau(:,i,l), &
              w0_in=opr%w0(:,i,l), &
              gt_in=opr%g(:,l), &
              u0=zenith_u(ii), &
              Rsfc=albedo, &
              amean=rbw%amean0, &
              surface_radiance=surf_rad, &
              fup=rbw%fup0, &
              fdn=rbw%fdn0 &
            )
          elseif (rtc%channel_type == IRChannel) then
            call two_stream_ir( &
              nz=nz, &
              tau=opr%tau(:,i,l), &
              w0=opr%w0(:,i,l), &
              gt=opr%g(:,l), &
              emissivity=emissivity, &
              has_hard_surface=has_hard_surface, &
              tau_min=ir_tau_min, &
              bplanck=rbw%bplanck, &
              fup=rbw%fup0, &
              fdn=rbw%fdn0 &
            )
          endif

          ! Apply k weights
          rbw%fup1 = rbw%fup1 + rbw%fup0*op%kset%wbin(i)
          rbw%fdn1 = rbw%fdn1 + rbw%fdn0*op%kset%wbin(i)
          if (rtc%channel_type == SolarChannel) then
            rbw%amean1 = rbw%amean1 + rbw%amean0*op%kset%wbin(i)
          endif
        enddo

        ! Apply zenith weights
        rbw%fup2 = rbw%fup2 + rbw%fup1*zenith_weights(ii)
        rbw%fdn2 = rbw%fdn2 + rbw%fdn1*zenith_weights(ii)
        if (rtc%channel_type == SolarChannel) then
          rbw%amean2 = rbw%amean2 + rbw%amean1*zenith_weights(ii)
        endif

      enddo
      
      ! save results. Here I reverse order so that
      ! fup_a(1,l) is ground level.
      do i = 1,nz+1
        n = nz+2-i
        fup_a(i,ll) = rbw%fup2(n)
        fdn_a(i,ll) = rbw%fdn2(n)
      enddo
      if (rtc%channel_type == SolarChannel) then
        do i = 1,nz+1
          n = nz+2-i
          amean(i,ll) = rbw%amean2(n)
        enddo
      endif
      do i = 1,nz
        n = nz+1-i
        tau_band(i,ll) = opr%tau_band(n,l) ! band optical thickness
      enddo
      
    endblock; enddo
    !$omp enddo
    !$omp end parallel

    block
    integer :: i
    real(dp) :: avg_freq, avg_wavl, dfreq
  
    ! In Solar case, units for fup_a and fdn_a are unit-less.
    ! need to multiply by photons_sol (mW/m2/Hz). Then
    ! fup_a and fdn_a are in units mW/m2/Hz.
    if (rtc%channel_type == SolarChannel) then
      do l = 1,rtc%nw
        fup_a(:,l) = fup_a(:,l)*photons_sol(l)*diurnal_fac
        fdn_a(:,l) = fdn_a(:,l)*photons_sol(l)*diurnal_fac
        amean(:,l) = amean(:,l)*photons_sol(l)*diurnal_fac

        ! Convert from from mW/m^2/Hz to mW/m^2/nm
        avg_freq = 0.5_dp*(rtc%freq(l) + rtc%freq(l+1))
        avg_wavl = 1.0e9_dp*c_light/avg_freq ! nm
        amean(:,l) = amean(:,l)*(avg_freq/avg_wavl)
        ! Convert from mW/m^2/nm to photons/cm^2/s
        amean(:,l) = amean(:,l)*(avg_wavl/(plank*c_light*1.0e16_dp))*(rtc%wavl(l+1) - rtc%wavl(l))
      enddo
    endif
    
    ! Integrate fluxes over wavelength or frequency grid.
    ! Units for fup_n and fdn_n are mW/m^2.
    fup_n = 0.0_dp
    fdn_n = 0.0_dp
    do l = 1,rtc%nw
      dfreq = rtc%freq(l) - rtc%freq(l+1)
      do i = 1,nz+1
        fup_n(i) = fup_n(i) + fup_a(i,l)*dfreq
        fdn_n(i) = fdn_n(i) + fdn_a(i,l)*dfreq
      enddo
    enddo

    end block
  
  end subroutine
  
end module
  
