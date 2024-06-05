from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

cdef extern from *:
  struct ProductionLoss:
    pass

cdef extern void deallocate_productionloss(ProductionLoss *ptr)

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