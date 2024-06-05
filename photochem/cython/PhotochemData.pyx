cimport PhotochemData_pxd as dat_pxd

cdef class PhotochemData:
  """This class contains data which rarely changes once a photochemical model has
  has been initialized.
  """

  cdef dat_pxd.PhotochemData *_ptr

  def __cinit__(self):
    self._ptr = NULL

  property nq:
    "The number of atmospheric species which evolve according to the PDEs"
    def __get__(self):
      cdef int nq
      dat_pxd.photochemdata_nq_get(self._ptr, &nq)
      return nq
      
  property np:
    "The number of particles"
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_np_get(self._ptr, &val)
      return val
      
  property ng:
    "The number of gases in the model"
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_ng_get(self._ptr, &val)
      return val

  property nsl:
    "The number of short-lived gases"
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_nsl_get(self._ptr, &val)
      return val
      
  property nll:
    "The number of long-lived gases"
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_nll_get(self._ptr, &val)
      return val
      
  property nsp:
    "The total number of species"
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_nsp_get(self._ptr, &val)
      return val
      
  property nw:
    "The number of wavelength bins"
    def __get__(self):
      cdef int val
      dat_pxd.photochemdata_nw_get(self._ptr, &val)
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
    """List, shape (nsp+2). A list of the species in the model (particles and gases). 
    The last two elements are 'hv' and 'M'.
    """
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_species_names_get_size(self._ptr, &dim1)
      cdef ndarray species_names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      dat_pxd.photochemdata_species_names_get(self._ptr, &dim1, <char *>species_names_c.data)
      return c2stringarr(species_names_c, S_STR_LEN, dim1)
      
  property atoms_names:
    "List, shape (natoms). The atoms in the model"
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_atoms_names_get_size(self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      dat_pxd.photochemdata_atoms_names_get(self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, S_STR_LEN, dim1)
      
  property reaction_equations:
    "List, shape (nRT). A list of all reaction equations"
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_reaction_equations_get_size(self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      dat_pxd.photochemdata_reaction_equations_get(self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
      
  property photonums:
    "ndarray[int,dim=1], shape (kj). The reaction number of each photolysis reaction"
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_photonums_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.int32)
      dat_pxd.photochemdata_photonums_get(self._ptr, &dim1, <int *>arr.data)
      return arr
      
  property wavl:
    "ndarray[int,dim=1], shape (nw). The wavelength bins (nm)"
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_wavl_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      dat_pxd.photochemdata_wavl_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property species_mass:
    "ndarray[double,dim=1], shape (nsp). The molar mass of each species (g/mol)"
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_species_mass_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      dat_pxd.photochemdata_species_mass_get(self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property species_redox:
    "ndarray[double,dim=1], shape (nsp). The redox state of each molecule"
    def __get__(self):
      cdef int dim1
      dat_pxd.photochemdata_species_redox_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      dat_pxd.photochemdata_species_redox_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property particle_sat:
    "list, shape(np)."
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
    