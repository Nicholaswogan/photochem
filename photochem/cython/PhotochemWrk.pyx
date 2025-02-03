cimport PhotochemWrk_pxd as wrk_pxd

cdef class PhotochemWrkEvo(PhotochemWrk):

  property pressure_hydro:
    """ndarray[double,dim=1], shape (nz). The hydrostatic pressure at the center of each 
    atmospheric layer (dynes/cm^2).
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrkevo_pressure_hydro_get_size(<wrk_pxd.PhotochemWrkEvo *>self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrkevo_pressure_hydro_get(<wrk_pxd.PhotochemWrkEvo *>self._ptr, &dim1, <double *>arr.data)
      return arr

cdef class PhotochemWrk:
  """This class contains data that changes during each step 
  when integrating the photochemical model
  """

  cdef wrk_pxd.PhotochemWrk *_ptr

  def __cinit__(self):
    self._ptr = NULL

  property nsteps_total:
    "int. Total number of steps in a robust integration."
    def __get__(self):
      cdef int val
      wrk_pxd.photochemwrk_nsteps_total_get(self._ptr, &val)
      return val

  property nsteps:
    "int. Number of integration steps excuted. Updated after every successful step."
    def __get__(self):
      cdef int val
      wrk_pxd.photochemwrk_nsteps_get(self._ptr, &val)
      return val

  property t_history:
    """ndarray[double,dim=1], shape (500). History of times at previous integration steps. 
    Index 1 is current, while index 2, 3, 4 are previous steps. Updated after every successful step.
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_t_history_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double, order="F")
      wrk_pxd.photochemwrk_t_history_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property mix_history:
    """ndarray[double,dim=3], shape (nq,nz,500). History of mixing ratios at previous integration steps. 
    Index 1 is current, while index 2, 3, 4 are previous steps. Updated after every successful step.
    """
    def __get__(self):
      cdef int dim1,dim2,dim3
      wrk_pxd.photochemwrk_mix_history_get_size(self._ptr, &dim1, &dim2, &dim3)
      cdef ndarray arr = np.empty((dim1,dim2,dim3), np.double, order="F")
      wrk_pxd.photochemwrk_mix_history_get(self._ptr, &dim1, &dim2, &dim3, <double *>arr.data)
      return arr

  property longdy:
    "double. Normalized change in mixing ratios over some number of integrations steps."
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_longdy_get(self._ptr, &val)
      return val

  property longdydt:
    """double. Normalized change in mixing ratios divided by change in time over some 
    number of integrations steps.
    """
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_longdydt_get(self._ptr, &val)
      return val

  property tn:
    "double. The current time of integration."
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_tn_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wrk_pxd.photochemwrk_tn_set(self._ptr, &val)

  property usol:
    """ndarray[double,dim=2], shape (nq,nz). Current gas concentrations in the atmosphere
    in units and molecules/cm^3.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_usol_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_usol_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=2] usol_new_):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_usol_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray usol_new = np.asfortranarray(usol_new_)
      if usol_new.shape[0] != dim1 or usol_new.shape[1] != dim2:
        raise PhotoException("Input usol is the wrong size.")
      wrk_pxd.photochemwrk_usol_set(self._ptr, &dim1, &dim2, <double *>usol_new.data)  
  
  property pressure:
    """ndarray[double,dim=1], shape (nz). The pressure at the center of each 
    atmospheric layer (dynes/cm^2).
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_pressure_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_pressure_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property density:
    """ndarray[double,dim=1], shape (nz). The total number density at the 
    center of each atmospheric layer (molecules/cm^3).
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_density_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_density_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property densities:
    """ndarray[double,dim=2], shape (nsp+1,nz). The number density (molecules/cm^3)
    or particle density (particles/cm^3) of each molecules or particle at each atmospheric
    layer.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_densities_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_densities_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property rx_rates:
    """ndarray[double,dim=2], shape (nz,nrT). Reaction rate constants in various units
    involving molecules cm^3 and s. These rates include 3rd body contributions.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_rx_rates_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_rx_rates_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property mubar:
    """ndarray[double,dim=1], shape (nz). The mean molar mass of each atmospheric layer
    (g/mol)
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_mubar_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_mubar_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property prates:
    "ndarray[double,dim=2], shape (nz,kj). The rates of each photolysis reaction (1/s)"
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_prates_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_prates_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property amean_grd:
    """ndarray[double,dim=2], shape (nz,nw). The mean irradiance at each 
    atmospheric layer in each wavelength bin. Assumes the total flux is 1.
    So, `amean_grd*photon_flux` is [photons/cm^2/s]
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_amean_grd_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_amean_grd_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property optical_depth:
    """ndarray[double,dim=2], shape (nz,nw). The optical depth at each atmospheric
    layer in each wavelength bin.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_optical_depth_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_optical_depth_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property surf_radiance:
    """ndarray[double,dim=1], shape (nw). The light hitting the ground.
    Assumes the total flux is 1. So, `surf_radiance*photon_flux` is [photons/cm^2/s]
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_surf_radiance_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_surf_radiance_get(self._ptr, &dim1, <double *>arr.data)
      return arr
  
