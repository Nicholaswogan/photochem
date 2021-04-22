

module photochem
  implicit none
  private
  integer,parameter :: real_kind = kind(1.0d0)
  
contains
  
  subroutine chem_reaction_rates(temperature, density, nz, nrT, reaction_rates)
    use photochem_data, only: rateparams, rxtypes, nreactants, &
                              nproducts, reactants_sp_inds, products_sp_inds, &
                              reverse2forward, nrF, nrR ! all protected vars
    use photochem_const, only: Rgas, k_boltz ! protected
    use photochem_vars, only: ! things that change between iterations
    use photochem_wrk, only: gibbs_energy_wrk ! module of preallocated work arrays
                              
    real(real_kind), intent(in) :: temperature(nz)
    real(real_kind), intent(in) :: density(nz)
    integer, intent(in) :: nz, nrT
    real(real_kind), intent(out) :: reaction_rates(nz, nrT)
    
    do i = 1,nrF
      if (rxtypes(i) == "falloff") then
        do j = 1,nz
          Troe = 1.d0 ! Here we must compute troe!
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
    
    if (nrF /= nrT) then
      ! compute gibbs energy
      call compute_gibbs_energy(temperature, nz, nsp, gibbs_energy_wrk)
      ! compute reverse rate
      do i = nrF+1,nrT
        l = nreactants(i)
        m = nproducts(i)
        n = reverse_info(i)
        if (rxtypes(n) /= "elementary") then
          l = l - 1 ! "M" doesn't count as a reactant for three-body and fall-off
          m = m - 1
        endif
        do j = 1,nz
          gibbR = 0.d0
          do k = 1,l
            gibbR = gibbR + gibbs_energy_wrk(reactants_sp_inds(k,i),j)
          enddo
          gibbP = 0.d0
          do k = 1,m
            gibbP = gibbP + gibbs_energy_wrk(products_sp_inds(k,i),j)
          enddo
          Dgibbs = gibbsR - gibbsP
          ! compute rate
          reaction_rates(j,i) = reaction_rates(j,n) * dexp(-Dgibbs/(Rgas * temperature(j))) * &
                                (k_boltz*temperature(j))**(m-l)
        enddo
      enddo
    endif
    
  end subroutine
  
  subroutine compute_gibbs_energy(temperature, nz, nsp, gibbs_energy)
    
  end subroutine
  

end module