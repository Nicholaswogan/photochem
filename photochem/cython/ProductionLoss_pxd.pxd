from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct ProductionLoss:
    pass
  struct ConservationFluxes:
    pass

cdef extern void deallocate_productionloss(ProductionLoss *ptr)
cdef extern void deallocate_conservationfluxes(ConservationFluxes *ptr)

cdef extern void productionloss_production_get_size(ProductionLoss *ptr, int *dim1, int *dim2)
cdef extern void productionloss_production_get(ProductionLoss *ptr, int *dim1, int *dim2, double *arr)

cdef extern void productionloss_loss_get_size(ProductionLoss *ptr, int *dim1, int *dim2)
cdef extern void productionloss_loss_get(ProductionLoss *ptr, int *dim1, int *dim2, double *arr)

cdef extern void productionloss_integrated_production_get_size(ProductionLoss *ptr, int *dim1)
cdef extern void productionloss_integrated_production_get(ProductionLoss *ptr, int *dim1, double *arr)

cdef extern void productionloss_integrated_loss_get_size(ProductionLoss *ptr, int *dim1)
cdef extern void productionloss_integrated_loss_get(ProductionLoss *ptr, int *dim1, double *arr)

cdef extern void productionloss_production_rx_get_size(ProductionLoss *ptr, int *dim1)
cdef extern void productionloss_production_rx_get(ProductionLoss *ptr, int *dim1, char *names)

cdef extern void productionloss_loss_rx_get_size(ProductionLoss *ptr, int *dim1)
cdef extern void productionloss_loss_rx_get(ProductionLoss *ptr, int *dim1, char *names)


cdef extern void conservationfluxes_chemical_production_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_chemical_production_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_chemical_loss_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_chemical_loss_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_rainout_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_rainout_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_special_h2o_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_special_h2o_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_condensation_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_condensation_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_evaporation_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_evaporation_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_custom_rates_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_custom_rates_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_lower_boundary_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_lower_boundary_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_upper_boundary_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_upper_boundary_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_net_flux_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_net_flux_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_columns_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_columns_get(ConservationFluxes *ptr, int *dim1, double *arr)

cdef extern void conservationfluxes_timescale_of_change_get_size(ConservationFluxes *ptr, int *dim1)
cdef extern void conservationfluxes_timescale_of_change_get(ConservationFluxes *ptr, int *dim1, double *arr)