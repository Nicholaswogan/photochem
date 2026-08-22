cimport ProductionLoss_pxd as pl_pxd

cdef class ProductionLoss:
  """Reaction-resolved production and loss rates for one chemical species.

  Instances are returned by :meth:`EvoAtmosphere.production_and_loss`.
  Reaction columns are ordered from largest to smallest vertically integrated
  rate.
  """

  cdef pl_pxd.ProductionLoss *_ptr

  def __cinit__(self):
    self._ptr = NULL

  def __dealloc__(self):
    pl_pxd.deallocate_productionloss(self._ptr)
    self._ptr = NULL 
      
  property production:
    """ndarray, shape (nz, nproduction)

    Production rate from each reaction at every layer, in molecules/cm^3/s.
    """
    def __get__(self):
      cdef int dim1, dim2
      pl_pxd.productionloss_production_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      pl_pxd.productionloss_production_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property loss:
    """ndarray, shape (nz, nloss)

    Loss rate from each reaction at every layer, in molecules/cm^3/s. Gas
    rainout is included as an additional loss process.
    """
    def __get__(self):
      cdef int dim1, dim2
      pl_pxd.productionloss_loss_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      pl_pxd.productionloss_loss_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property integrated_production:
    """ndarray, shape (nproduction,)

    Vertically integrated production rate from each reaction, in
    molecules/cm^2/s.
    """
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_integrated_production_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.productionloss_integrated_production_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property integrated_loss:
    """ndarray, shape (nloss,)

    Vertically integrated loss rate from each reaction, in molecules/cm^2/s.
    """
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_integrated_loss_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.productionloss_integrated_loss_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property production_rx:
    """list[str], shape (nproduction,)

    Reaction equations corresponding to the columns of :attr:`production`.
    """
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_production_rx_get_size(self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      pl_pxd.productionloss_production_rx_get(self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
    
  property loss_rx:
    """list[str], shape (nloss,)

    Reaction equations or process labels corresponding to the columns of
    :attr:`loss`. The additional gas-rainout loss is labeled ``"rainout"``.
    """
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_loss_rx_get_size(self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      pl_pxd.productionloss_loss_rx_get(self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
  
