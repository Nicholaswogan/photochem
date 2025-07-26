
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine deallocate_productionloss(ptr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    deallocate(pl)
  end subroutine

  subroutine deallocate_conservationfluxes(ptr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    deallocate(f)
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine productionloss_production_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    dim1 = size(pl%production,1)
    dim2 = size(pl%production,2)
  end subroutine
  subroutine productionloss_production_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    arr = pl%production
  end subroutine

  subroutine productionloss_loss_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    dim1 = size(pl%loss,1)
    dim2 = size(pl%loss,2)
  end subroutine
  subroutine productionloss_loss_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1, dim2)
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    arr = pl%loss
  end subroutine

  subroutine productionloss_integrated_production_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    dim1 = size(pl%integrated_production,1)
  end subroutine
  subroutine productionloss_integrated_production_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    arr = pl%integrated_production
  end subroutine
  
  subroutine productionloss_integrated_loss_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    dim1 = size(pl%integrated_loss,1)
  end subroutine
  subroutine productionloss_integrated_loss_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    arr = pl%integrated_loss
  end subroutine
  
  subroutine productionloss_production_rx_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    dim1 = size(pl%production_rx)
  end subroutine
  subroutine productionloss_production_rx_get(ptr, dim1, names) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: names(dim1*m_str_len+1)
    type(ProductionLoss), pointer :: pl
    integer :: i, j, k
    
    call c_f_pointer(ptr, pl)
    do i = 1,dim1
      do j = 1,m_str_len
        k = j + (i - 1) * m_str_len
        names(k) = pl%production_rx(i)(j:j)
      enddo
    enddo
    names(dim1*m_str_len+1) = c_null_char
  end subroutine
  
  subroutine productionloss_loss_rx_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    dim1 = size(pl%loss_rx)
  end subroutine
  subroutine productionloss_loss_rx_get(ptr, dim1, names) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: names(dim1*m_str_len+1)
    type(ProductionLoss), pointer :: pl
    integer :: i, j, k
    
    call c_f_pointer(ptr, pl)
    do i = 1,dim1
      do j = 1,m_str_len
        k = j + (i - 1) * m_str_len
        names(k) = pl%loss_rx(i)(j:j)
      enddo
    enddo
    names(dim1*m_str_len+1) = c_null_char
  end subroutine
  


  subroutine conservationfluxes_chemical_production_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%chemical_production,1)
  end subroutine
  subroutine conservationfluxes_chemical_production_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%chemical_production
  end subroutine

  subroutine conservationfluxes_chemical_loss_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%chemical_loss,1)
  end subroutine
  subroutine conservationfluxes_chemical_loss_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%chemical_loss
  end subroutine

  subroutine conservationfluxes_rainout_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%rainout,1)
  end subroutine
  subroutine conservationfluxes_rainout_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%rainout
  end subroutine

  subroutine conservationfluxes_special_h2o_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%special_H2O,1)
  end subroutine
  subroutine conservationfluxes_special_h2o_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%special_H2O
  end subroutine

  subroutine conservationfluxes_condensation_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%condensation,1)
  end subroutine
  subroutine conservationfluxes_condensation_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%condensation
  end subroutine

  subroutine conservationfluxes_evaporation_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%evaporation,1)
  end subroutine
  subroutine conservationfluxes_evaporation_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%evaporation
  end subroutine

  subroutine conservationfluxes_custom_rates_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%custom_rates,1)
  end subroutine
  subroutine conservationfluxes_custom_rates_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%custom_rates
  end subroutine

  subroutine conservationfluxes_lower_boundary_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%lower_boundary,1)
  end subroutine
  subroutine conservationfluxes_lower_boundary_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%lower_boundary
  end subroutine

  subroutine conservationfluxes_upper_boundary_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%upper_boundary,1)
  end subroutine
  subroutine conservationfluxes_upper_boundary_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%upper_boundary
  end subroutine

  subroutine conservationfluxes_net_flux_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%net_flux,1)
  end subroutine
  subroutine conservationfluxes_net_flux_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%net_flux
  end subroutine

  subroutine conservationfluxes_columns_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%columns,1)
  end subroutine
  subroutine conservationfluxes_columns_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%columns
  end subroutine

  subroutine conservationfluxes_timescale_of_change_get_size(ptr, dim1) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    dim1 = size(f%timescale_of_change,1)
  end subroutine
  subroutine conservationfluxes_timescale_of_change_get(ptr, dim1, arr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ConservationFluxes), pointer :: f
    call c_f_pointer(ptr, f)
    arr = f%timescale_of_change
  end subroutine