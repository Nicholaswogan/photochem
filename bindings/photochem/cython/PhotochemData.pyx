cimport PhotochemData_pxd as dat_pxd

cdef class PhotochemData:
  """Read-mostly mechanism and planetary data owned by a model.

  This object is a live borrowed view available as :attr:`EvoAtmosphere.dat`;
  it should not be constructed directly. Array-valued properties return copies.
  """

  cdef dat_pxd.PhotochemData *_ptr

  def __cinit__(self):
    self._ptr = NULL

  property nq:
    """int. Number of particle and gas species evolved by transport."""
    def __get__(self):
      cdef int nq
      dat_pxd.photochemdata_nq_get(self._ptr, &nq)
      return nq
      
  property np:
    """int. Number of particle species."""
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_np_get(self._ptr, &val)
      return val
      
  property ng:
    """int. Number of gas species, including short-lived gases."""
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_ng_get(self._ptr, &val)
      return val

  property nsl:
    """int. Number of short-lived gas species."""
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_nsl_get(self._ptr, &val)
      return val
      
  property nll:
    """int. Number of long-lived gas species."""
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_nll_get(self._ptr, &val)
      return val
      
  property nsp:
    """int. Total number of chemical species, excluding ``hv`` and ``M``."""
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_nsp_get(self._ptr, &val)
      return val
      
  property nw:
    """int. Number of wavelength bins."""
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_nw_get(self._ptr, &val)
      return val

  property gas_rainout:
    "bool. Whether tropospheric gas rainout is enabled."
    def __get__(self):
      cdef bool val
      dat_pxd.photochemdata_gas_rainout_get(self._ptr, &val)
      return val

  property planet_mass:
    "Planet mass (g)."
    def __get__(self):
      cdef double val
      dat_pxd.photochemdata_planet_mass_get(self._ptr, &val)
      return val

  property planet_radius:
    "Planet radius (cm)."
    def __get__(self):
      cdef double val
      dat_pxd.photochemdata_planet_radius_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      dat_pxd.photochemdata_planet_radius_set(self._ptr, &val)

  property species_names:
    """list[str], length ``nsp + 2``. Ordered model species names.

    Particle species precede gases; the final two entries are ``"hv"`` and
    ``"M"``.
    """
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_species_names_get_size(self._ptr, &dim1)
      cdef ndarray species_names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      dat_pxd.photochemdata_species_names_get(self._ptr, &dim1, <char *>species_names_c.data)
      return c2stringarr(species_names_c, S_STR_LEN, dim1)
      
  property species_composition:
    """ndarray[int], shape (natoms, nsp + 2). Species elemental composition."""
    def __get__(self):
      cdef int dim1, dim2
      dat_pxd.photochemdata_species_composition_get_size(self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.intc, order="F")
      dat_pxd.photochemdata_species_composition_get(self._ptr, &dim1, &dim2, <int *>arr.data)
      return arr

  property atoms_names:
    """list[str], length ``natoms``. Element names in composition-array order."""
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_atoms_names_get_size(self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      dat_pxd.photochemdata_atoms_names_get(self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, S_STR_LEN, dim1)
      
  property reaction_equations:
    """list[str], length ``nrT``. Reaction equations in reaction-rate order."""
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_reaction_equations_get_size(self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      dat_pxd.photochemdata_reaction_equations_get(self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
      
  property photonums:
    """ndarray[int], shape (kj,). One-based reaction numbers for photolysis."""
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_photonums_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.int32)
      dat_pxd.photochemdata_photonums_get(self._ptr, &dim1, <int *>arr.data)
      return arr
      
  property wavl:
    """ndarray[float], shape (nw + 1,). Wavelength-bin edges in nm."""
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_wavl_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      dat_pxd.photochemdata_wavl_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property species_mass:
    """ndarray[float], shape (nsp,). Species molar masses in g/mol."""
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_species_mass_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      dat_pxd.photochemdata_species_mass_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property species_redox:
    """ndarray[float], shape (nsp,). Redox coefficient of each species."""
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_species_redox_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      dat_pxd.photochemdata_species_redox_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property particle_sat:
    """list[SaturationData], length ``np``. Particle saturation models.

    Entries are live borrowed views owned by the model.
    """
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_particle_sat_get_size(self._ptr, &dim1)
      cdef atom_pxd.SaturationData **arrp = <atom_pxd.SaturationData **> malloc(dim1 * sizeof(atom_pxd.SaturationData *))
      dat_pxd.photochemdata_particle_sat_get(self._ptr, &dim1, arrp)
      arr1 = []
      for i in range(dim1):
        tmp = SaturationData()
        tmp._ptr = arrp[i]
        arr1.append(tmp)
      free(arrp)
      return arr1
    
  property H_escape_coeff:
    """float. Zahnle hydrogen-escape coefficient in molecules/cm^2/s."""
    def __get__(self):
      cdef double val
      dat_pxd.photochemdata_h_escape_coeff_get(self._ptr, &val)
      return val
