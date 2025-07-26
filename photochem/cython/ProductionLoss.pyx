cimport ProductionLoss_pxd as pl_pxd

cdef class ProductionLoss:
  """A class containing data describing the reactions that produce and destroy 
  a species. This class is produced when calling the `production_and_loss` routine.
  """

  cdef pl_pxd.ProductionLoss *_ptr

  def __cinit__(self):
    self._ptr = NULL

  def __dealloc__(self):
    pl_pxd.deallocate_productionloss(self._ptr)
    self._ptr = NULL 
      
  property production:
    """ndarray[double,dim=2], shape (nz,nproduction). The rate the molecule
    is produced in molecules/cm^3 at each atmospheric layer from each reaction
    that produces the molecule
    """
    def __get__(self):
      cdef int dim1, dim2
      pl_pxd.productionloss_production_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      pl_pxd.productionloss_production_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property loss:
    """ndarray[double,dim=2], shape (nz,nloss). The rate the molecule
    is destroyed in molecules/cm^3 at each atmospheric layer from each reaction
    that destroys the molecule
    """
    def __get__(self):
      cdef int dim1, dim2
      pl_pxd.productionloss_loss_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      pl_pxd.productionloss_loss_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property integrated_production:
    """ndarray[double,dim=1], shape (nproduction). The vertically-integrated production
    rate of the molecule in molecules/cm^2 for each reaction that produces the molecule.
    """
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_integrated_production_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.productionloss_integrated_production_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property integrated_loss:
    """ndarray[double,dim=1], shape (nloss). The vertically-integrated production
    rate of the molecule in molecules/cm^2 for each reaction that produces the molecule.
    """
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_integrated_loss_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.productionloss_integrated_loss_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property production_rx:
    "List, shape (nproduction). The reaction equations that produce the molecule."
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_production_rx_get_size(self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      pl_pxd.productionloss_production_rx_get(self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
    
  property loss_rx:
    "List, shape (nloss). The reaction equations that destroy the molecule."
    def __get__(self):
      cdef int dim1
      pl_pxd.productionloss_loss_rx_get_size(self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      pl_pxd.productionloss_loss_rx_get(self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
  
cdef class ConservationFluxes:

  cdef pl_pxd.ConservationFluxes *_ptr

  def __cinit__(self):
    self._ptr = NULL

  def __dealloc__(self):
    pl_pxd.deallocate_conservationfluxes(self._ptr)
    self._ptr = NULL 

  property chemical_production:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_chemical_production_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_chemical_production_get(self._ptr, &dim1, <double *>arr.data)
      return arr
  
  property chemical_loss:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_chemical_loss_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_chemical_loss_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property rainout:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_rainout_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_rainout_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property special_H2O:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_special_h2o_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_special_h2o_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property condensation:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_condensation_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_condensation_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property evaporation:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_evaporation_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_evaporation_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property custom_rates:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_custom_rates_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_custom_rates_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property lower_boundary:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_lower_boundary_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_lower_boundary_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property upper_boundary:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_upper_boundary_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_upper_boundary_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property net_flux:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_net_flux_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_net_flux_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property columns:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_columns_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_columns_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property timescale_of_change:
    def __get__(self):
      cdef int dim1
      pl_pxd.conservationfluxes_timescale_of_change_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double)
      pl_pxd.conservationfluxes_timescale_of_change_get(self._ptr, &dim1, <double *>arr.data)
      return arr
