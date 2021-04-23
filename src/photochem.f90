

module photochem
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  
contains
  
  subroutine compute_reaction_rates(temperature, density, nz, nrT, reaction_rates)
    use photochem_data, only: rateparams, rxtypes, nreactants, &
                              nproducts, reactants_sp_inds, products_sp_inds, &
                              reverse_info, nrF, nrR ! all protected vars
    use photochem_const, only: Rgas, k_boltz ! constants
    use photochem_vars, only: real_nz_nsp ! pre-allocated work array
                              
    real(real_kind), intent(in) :: temperature(nz)
    real(real_kind), intent(in) :: density(nz)
    integer, intent(in) :: nz, nrT
    real(real_kind), intent(out) :: reaction_rates(nz, nrT)
    
    integer :: i, j, k, n, l, m
    real(real_kind) :: Troe, A1, A2
    real(real_kind) :: gibbR_forward, gibbP_forward
    real(real_kind) :: Dg_forward
    
    do i = 1,nrF
      if (rxtypes(i) == "falloff") then
        do j = 1,nz
          Troe = 1.d0 ! Here we must compute troe.
          A1  = rateparams(1,i) * temperature(j)**rateparams(2,i) &
                * dexp(-rateparams(3,i)/temperature(j)) * density(j)
          B1 = rateparams(4,i) * temperature(j)**rateparams(5,i) &
                * dexp(-rateparams(6,i)/temperature(j))
          reaction_rates(j,i) = A1 * Troe/(A1/B1 + 1)
        enddo
      else ! three-body or elementary are the same
        do j = 1,nz 
          reaction_rates(j,i) = rateparams(1,i) * temperature(j)**rateparams(2,i) &
                                * dexp(-rateparams(3,i)/temperature(j))
        enddo
      endif
    enddo
    
    if (nrF /= nrT) then ! if there are reverse reactions
      ! compute gibbs energy
      call compute_gibbs_energy(temperature, nz, nsp, real_nz_nsp)
      ! compute reverse rate
      do i = nrF+1,nrT
        n = reverse_info(i) ! Reaction number of the forward
        l = nreactants(n) ! number of reactants for the forward reaction
        m = nproducts(n) ! number of products for the forward reaction
        if (rxtypes(n) /= "elementary") then
          l = l - 1 ! "M" doesn't count as a reactant for three-body and fall-off
          m = m - 1
        endif
        do j = 1,nz
          gibbR_forward = 0.d0
          do k = 1,l
            gibbR_forward = gibbR_forward + real_nz_nsp(reactants_sp_inds(k,n),j)
          enddo
          gibbP_forward = 0.d0
          do k = 1,m
            gibbP_forward = gibbP_forward + real_nz_nsp(products_sp_inds(k,n),j)
          enddo
          Dg_forward = gibbsP_forward - gibbsRf_forward ! DG of the forward reaction (J/mol)
          ! compute the reverse rate
          reaction_rates(j,i) = reaction_rates(j,n) * &
                                (1.d0/dexp(-Dg_forward/(Rgas * temperature(j)))) * &
                                (k_boltz*temperature(j)/1.d6)**(m-l)
        enddo
      enddo
    endif

  end subroutine
  
  subroutine compute_gibbs_energy(temperature, nz, nsp, gibbs_energy)
    use photochem_data, only:thermo_data, thermo_temps
    use photochem_vars, only: real_nz, int_nz ! work arrays
    
    real(real_kind), intent(in) :: temperature(nz)
    integer, intent(in) :: nz, nsp
    real(real_kind), intent(in) :: gibbs_energy(nz,nsp)
    
    integer :: i, j, k
    real(real_kind) :: enthalpy, entropy, TT
    
    ! Shomate parameterization
    do j = 1,nz
      real_nz(j) = temperature(j)/1000.d0 
      if ((temperature(j) >= thermo_temps(1,i)) .and. &
         (temperature(j) < thermo_temps(2,i))) then
         int_nz(j) = 1
      elseif ((temperature(j) >= thermo_temps(2,i)) .and. &
            (temperature(j) <= thermo_temps(3,i))) then
         int_nz(j) = 2
      else
        print*,'PhotoError: The temperature is not within the ranges ', &
               'given for the thermodynamic data.'
        stop
      endif 
    enddo

    do i = 1,nsp
      do j = 1,nz
        k = int_nz(j)
        TT = real_nz(j)
        enthalpy = thermo_data(1,k,i)*TT + (thermo_data(2,k,i)*TT**2)/2.d0 &
                 + (thermo_data(3,k,i)*TT**3)/3.d0  + (thermo_data(4,k,i)*TT**4)/4.d0 &
                 - thermo_data(5,k,i)/TT + thermo_data(6,k,i)
        entropy = thermo_data(1,k,i)*dlog(TT) + thermo_data(2,k,i)*TT &
                + (thermo_data(3,k,i)*TT**2)/2.d0 + (thermo_data(4,k,i)*TT**3)/3.d0 &
                - thermo_data(5,k,i)/(2.d0 * TT**2) + thermo_data(7,k,i)
        gibbs_energy(j,i) = enthalpy*1000.d0 - temperature(j)*entropy
        ! j/mol
      enddo
    enddo
  end subroutine
  
end module