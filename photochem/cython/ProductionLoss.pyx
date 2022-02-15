cimport ProductionLoss_pxd as pl_pxd

cdef class ProductionLoss:
  cdef void *_ptr

  def __cinit__(self):
    # never allocate. Only allow creation by functions.
    pass

  def __dealloc__(self):
    pl_pxd.deallocate_productionloss(&self._ptr)
    self._ptr = NULL 
      
  property production:
    def __get__(self):
      cdef int dim1, dim2
      pl_pxd.productionloss_production_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      pl_pxd.productionloss_production_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property loss:
    def __get__(self):
      cdef int dim1, dim2
      pl_pxd.productionloss_loss_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      pl_pxd.productionloss_loss_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property integrated_production:
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_integrated_production_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.productionloss_integrated_production_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property integrated_loss:
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_integrated_loss_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.productionloss_integrated_loss_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property production_rx:
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_production_rx_get_size(&self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      pl_pxd.productionloss_production_rx_get(&self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
    
  property loss_rx:
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_loss_rx_get_size(&self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      pl_pxd.productionloss_loss_rx_get(&self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
  
