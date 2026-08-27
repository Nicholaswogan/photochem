cimport PhotochemWrk_pxd as wrk_pxd

cdef class PhotochemWrk:
  """Prepared atmosphere, solver state, and runtime diagnostics.

  This object is a live borrowed view available as [EvoAtmosphere.wrk][photochem.EvoAtmosphere.wrk];
  it should not be constructed directly. Most diagnostics describe the most
  recently prepared atmospheric state. Array-valued getters return copies;
  assigning ``usol`` or ``tn`` updates the live model state.
  """

  cdef wrk_pxd.PhotochemWrk *_ptr

  def __cinit__(self):
    self._ptr = NULL

  property surface_pressure:
    "double. Surface pressure derived from the current atmospheric column, in bars."
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_surface_pressure_get(self._ptr, &val)
      return val

  property n_toa_pressure_updates:
    """int. Number of successful automatic TOA pressure updates."""
    def __get__(self):
      cdef int val
      wrk_pxd.photochemwrk_n_toa_pressure_updates_get(self._ptr, &val)
      return val

  property n_toa_pressure_failures:
    """int. Number of failed automatic TOA pressure maintenance attempts."""
    def __get__(self):
      cdef int val
      wrk_pxd.photochemwrk_n_toa_pressure_failures_get(self._ptr, &val)
      return val

  property nsteps_since_toa_pressure_update:
    """int. Accepted steps since the last successful TOA pressure update."""
    def __get__(self):
      cdef int val
      wrk_pxd.photochemwrk_nsteps_since_toa_pressure_update_get(self._ptr, &val)
      return val

  property pressure_hydro:
    """ndarray[double,dim=1], shape (nz). The hydrostatic pressure at the center of each
    atmospheric layer (dyn/cm^2).
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_pressure_hydro_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_pressure_hydro_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property scale_height:
    "ndarray[double,dim=1], shape (nz). Atmospheric scale height (cm)."
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_scale_height_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_scale_height_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property wfall:
    """ndarray[double,dim=2], shape (np,nz). Particle settling velocities
    (cm/s).
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_wfall_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_wfall_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property gas_sat_den:
    """ndarray[double,dim=2], shape (np,nz). Saturation number densities
    (molecules/cm^3) for condensing particle gas phases. Entries for
    non-condensing particles are zero.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_gas_sat_den_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_gas_sat_den_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property molecules_per_particle:
    """ndarray[double,dim=2], shape (np,nz). Condensate molecules represented
    by each particle.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_molecules_per_particle_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_molecules_per_particle_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property rainout_rates:
    """ndarray[double,dim=2], shape (nq,nz). First-order gas rainout loss
    rates (1/s). Values are zero where rainout is inactive.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_rainout_rates_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_rainout_rates_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property nsteps_total:
    "int. Total number of accepted steps in a robust integration."
    def __get__(self):
      cdef int val
      wrk_pxd.photochemwrk_nsteps_total_get(self._ptr, &val)
      return val

  property nerrors_total:
    "int. Total number of failed steps in a robust integration."
    def __get__(self):
      cdef int val
      wrk_pxd.photochemwrk_nerrors_total_get(self._ptr, &val)
      return val

  property robust_stepper_initialized:
    "bool. True while the active CVODE stepper belongs to a robust integration."
    def __get__(self):
      cdef bool val
      wrk_pxd.photochemwrk_robust_stepper_initialized_get(self._ptr, &val)
      return val

  property nsteps:
    """int. Number of accepted steps in the current integration segment."""
    def __get__(self):
      cdef int val
      wrk_pxd.photochemwrk_nsteps_get(self._ptr, &val)
      return val

  property t_history:
    """ndarray[float], shape (500,). Recent accepted integration times.

    Element 0 is current; subsequent elements contain progressively older
    entries, in seconds. Updated after each accepted step.
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_t_history_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty((dim1), np.double, order="F")
      wrk_pxd.photochemwrk_t_history_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property mix_history:
    """ndarray[float], shape (nq, nz, 500). Recent mixing-ratio history.

    ``mix_history[:, :, 0]`` is current; subsequent slices contain
    progressively older entries. Updated after each accepted step.
    """
    def __get__(self):
      cdef int dim1,dim2,dim3
      wrk_pxd.photochemwrk_mix_history_get_size(self._ptr, &dim1, &dim2, &dim3)
      cdef ndarray arr = np.empty((dim1,dim2,dim3), np.double, order="F")
      wrk_pxd.photochemwrk_mix_history_get(self._ptr, &dim1, &dim2, &dim3, <double *>arr.data)
      return arr

  property longdy:
    """float. Normalized mixing-ratio change over the convergence interval."""
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_longdy_get(self._ptr, &val)
      return val

  property longdydt:
    """float. ``longdy`` divided by elapsed history time, in s^-1."""
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_longdydt_get(self._ptr, &val)
      return val

  property tn:
    """float. Current integration time in seconds.

    Right-hand-side and Jacobian evaluations update this value and pass it to
    the time-dependent photon-flux callback.
    """
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_tn_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      wrk_pxd.photochemwrk_tn_set(self._ptr, &val)

  property usol:
    """ndarray[float], shape (nq, nz). Current evolved number densities.

    Gas and condensed-material entries are in molecules/cm^3.
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
        raise PhotoException("usol must have shape (nq, nz).")
      wrk_pxd.photochemwrk_usol_set(self._ptr, &dim1, &dim2, <double *>usol_new.data)  
  
  property pressure:
    """ndarray[double,dim=1], shape (nz). The pressure at the center of each 
    atmospheric layer (dyn/cm^2).
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
    """ndarray[float], shape (nsp + 1, nz). Chemistry number densities.

    The first ``nsp`` rows follow ``dat.species_names[:nsp]``. Particle rows
    are in particles/cm^3, gas rows are in molecules/cm^3, and the final row
    is the unit density used for ``hv``.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_densities_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_densities_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property rx_rates:
    """ndarray[float], shape (nz, nrT). Effective reaction rate coefficients.

    Units depend on reaction order. Falloff and third-body contributions are
    already included.
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
    """ndarray[float], shape (nz, kj). Photolysis frequencies in s^-1."""
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_prates_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_prates_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
  
  property amean_grd:
    """ndarray[float], shape (nz, nw). Dimensionless mean-radiance factors.

    Multiplying by ``var.photon_flux`` gives photons/cm^2/s in each bin.
    """
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_amean_grd_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_amean_grd_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property optical_depth:
    """ndarray[float], shape (nz, nw). Dimensionless layer optical depth."""
    def __get__(self):
      cdef int dim1, dim2
      wrk_pxd.photochemwrk_optical_depth_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wrk_pxd.photochemwrk_optical_depth_get(self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr
      
  property surf_radiance:
    """ndarray[float], shape (nw,). Dimensionless surface-radiance factors.

    Multiplying by ``var.photon_flux`` gives photons/cm^2/s at the surface.
    """
    def __get__(self):
      cdef int dim1
      wrk_pxd.photochemwrk_surf_radiance_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wrk_pxd.photochemwrk_surf_radiance_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property VH2_esc:
    """float. Effective H2 effusion velocity representing escape, in cm/s."""
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_vh2_esc_get(self._ptr, &val)
      return val

  property VH_esc:
    """float. Effective H effusion velocity representing escape, in cm/s."""
    def __get__(self):
      cdef double val
      wrk_pxd.photochemwrk_vh_esc_get(self._ptr, &val)
      return val
