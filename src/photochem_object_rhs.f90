submodule(photochem_object) photochem_object_rhs
  implicit none
  
contains
  
  subroutine reaction_rates(self, nsp, nz, nrT, density, densities, rx_rates, err)
    use photochem_eqns, only: arrhenius_rate, Troe_noT2, Troe_withT2, falloff_rate
    use photochem_const, only: Rgas, k_boltz, smallest_real ! constants
    
    class(Photochem), intent(in) :: self
    integer, intent(in) :: nsp, nz, nrT
    real(real_kind), intent(in) :: density(nz)
    real(real_kind), intent(in) :: densities(nsp+1,nz)
    real(real_kind), intent(out) :: rx_rates(nz, nrT)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k, n, l, m
    real(real_kind) :: gibbs_energy(nz,self%dat%ng)
    real(real_kind) :: eff_den(nz), F(nz), k0, kinf(nz), Pr(nz)
    real(real_kind) :: gibbR_forward, gibbP_forward
    real(real_kind) :: Dg_forward
    err = ''
    
    do i = 1,self%dat%nrF
      if (self%dat%rxtypes(i) == 1) then ! elementary        
        do j = 1,nz 
          rx_rates(j,i) = arrhenius_rate(self%dat%rateparams(1,i), self%dat%rateparams(2,i), &
                                         self%dat%rateparams(3,i), self%var%temperature(j))
          
        enddo        
      elseif (self%dat%rxtypes(i) == 2) then ! three-body
        n = self%dat%num_efficient(i)
        do j = 1,self%var%nz
          eff_den(j) = density(j) * self%dat%def_eff(i)
          ! subtract the default efficiency, then add the perscribed one
          do k = 1,n ! if n is 0 then it will be skipped
            l = self%dat%eff_sp_inds(k,i)
            eff_den(j) = eff_den(j) - self%dat%def_eff(i)*densities(l,j) &
                                    + self%dat%efficiencies(k,i)*densities(l,j)
          enddo
          rx_rates(j,i) = arrhenius_rate(self%dat%rateparams(1,i), self%dat%rateparams(2,i), &
                                         self%dat%rateparams(3,i), self%var%temperature(j)) &
                                         * eff_den(j) ! we multiply density here
        enddo
        
      
      elseif (self%dat%rxtypes(i) == 3) then ! falloff
        ! compute eff_den, kinf, and Pr at all altitudes
        n = self%dat%num_efficient(i)
        do j = 1,nz
          eff_den(j) = density(j) * self%dat%def_eff(i)
          ! subtract the default efficiency, then add the perscribed one
          do k = 1,n
            l = self%dat%eff_sp_inds(k,i)
            eff_den(j) = eff_den(j) - self%dat%def_eff(i)*densities(l,j) &
                                    + self%dat%efficiencies(k,i)*densities(l,j)
          enddo
          k0  = arrhenius_rate(self%dat%rateparams(1,i), self%dat%rateparams(2,i), &
                               self%dat%rateparams(3,i), self%var%temperature(j))
          kinf(j) = arrhenius_rate(self%dat%rateparams(4,i), self%dat%rateparams(5,i), &
                                   self%dat%rateparams(6,i), self%var%temperature(j))
          kinf(j) = max(kinf(j),smallest_real)
          Pr(j) = k0*eff_den(j)/kinf(j)                        
        enddo
        
        ! compute falloff function
        if (self%dat%falloff_type(i) == 0) then ! no falloff function
          F = 1.d0
        elseif (self%dat%falloff_type(i) == 1) then ! Troe falloff without T2
          do j = 1,nz
            F(j) = Troe_noT2(self%dat%rateparams(7,i), self%dat%rateparams(8,i), &
                             self%dat%rateparams(10,i), self%var%temperature(j), Pr(j))
          enddo
        elseif (self%dat%falloff_type(i) == 2) then ! Troe falloff with T2
          do j = 1,nz
            F(j) = Troe_withT2(self%dat%rateparams(7,i), self%dat%rateparams(8,i), self%dat%rateparams(9,i), &
                               self%dat%rateparams(10,i), self%var%temperature(j), Pr(j))
          enddo
        endif
        
        ! compute rate
        do j = 1,nz
          rx_rates(j,i) = falloff_rate(kinf(j), Pr(j), F(j))
        enddo
    
      endif
    enddo ! end loop over forward reactions

    if (self%dat%reverse) then ! if there are reverse reactions
      ! compute gibbs energy at all altitudes
      call compute_gibbs_energy(self, self%var%nz, self%dat%ng, gibbs_energy, err)
      if (len_trim(err) /= 0) return
      ! compute reverse rate
      do i = self%dat%nrF+1,self%dat%nrT
        
        n = self%dat%reverse_info(i) ! Reaction number of the forward
        l = self%dat%nreactants(n) ! number of reactants for the forward reaction
        m = self%dat%nproducts(n) ! number of products for the forward reaction
        do j = 1,self%var%nz
          gibbR_forward = 0.d0
          do k = 1,l
            gibbR_forward = gibbR_forward + gibbs_energy(j,self%dat%reactants_sp_inds(k,n)-self%dat%npq)
          enddo
          gibbP_forward = 0.d0
          do k = 1,m
            gibbP_forward = gibbP_forward + gibbs_energy(j,self%dat%products_sp_inds(k,n)-self%dat%npq)
          enddo
          Dg_forward = gibbP_forward - gibbR_forward ! DG of the forward reaction (J/mol)
          ! compute the reverse rate
          rx_rates(j,i) = rx_rates(j,n) * &
                          (1.d0/exp(-Dg_forward/(Rgas * self%var%temperature(j)))) * &
                          (k_boltz*self%var%temperature(j)/1.d6)**(m-l)
        enddo
      enddo
    endif

  end subroutine
  
  subroutine compute_gibbs_energy(self, nz, ng, gibbs_energy, err)
    use photochem_eqns, only: gibbs_energy_shomate
    
    class(Photochem), intent(in) :: self
    integer, intent(in) :: nz, ng
    real(real_kind), intent(out) :: gibbs_energy(nz,ng)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, k

    err = ''
    
    do i = 1,self%dat%ng
      do j = 1,self%var%nz
        if ((self%var%temperature(j) >= self%dat%thermo_temps(1,i)) .and. &
           (self%var%temperature(j) < self%dat%thermo_temps(2,i))) then
           k = 1
        elseif ((self%var%temperature(j) >= self%dat%thermo_temps(2,i)) .and. &
              (self%var%temperature(j) <= self%dat%thermo_temps(3,i))) then
           k = 2
        else
          err = 'The temperature is not within the ranges '// &
                'given for the thermodynamic data for '//trim(self%dat%species_names(i))
          return
        endif
        ! j/mol
        call gibbs_energy_shomate(self%dat%thermo_data(1:7,k,i), self%var%temperature(j), gibbs_energy(j,i))
      enddo
    enddo
  end subroutine
  
  
end submodule