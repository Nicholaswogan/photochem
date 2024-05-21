
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine deallocate_productionloss(ptr) bind(c)
    type(c_ptr), value, intent(in) :: ptr
    type(ProductionLoss), pointer :: pl
    call c_f_pointer(ptr, pl)
    deallocate(pl)
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
  