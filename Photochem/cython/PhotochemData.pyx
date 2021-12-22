
cdef extern void allocate_photochemdata(void *ptr)
cdef extern void deallocate_photochemdata(void *ptr)

cdef extern void photochemdata_nq_get(void *ptr, int *nq)
cdef extern void photochemdata_np_get(void *ptr, int *nq)
cdef extern void photochemdata_nsp_get(void *ptr, int *nq)
cdef extern void photochemdata_ng_get(void *ptr, int *nq)
cdef extern void photochemdata_nsl_get(void *ptr, int *nq)
cdef extern void photochemdata_nll_get(void *ptr, int *nq)
cdef extern void photochemdata_nw_get(void *ptr, int *nq)

cdef extern void photochemdata_species_names_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_species_names_get(void *ptr, int *dim1, char* species_names)

cdef extern void photochemdata_atoms_names_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_atoms_names_get(void *ptr, int *dim1, char* names)

cdef extern void photochemdata_reaction_equations_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_reaction_equations_get(void *ptr, int *dim1, char* names)

cdef extern void photochemdata_photonums_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_photonums_get(void *ptr, int *dim1, int *arr)

cdef extern void photochemdata_wavl_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_wavl_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemdata_species_mass_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_species_mass_get(void *ptr, int *dim1, double *arr)

cdef extern void photochemdata_species_redox_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_species_redox_get(void *ptr, int *dim1, double *arr)

cdef class PhotochemData:
  cdef void *_ptr
  cdef bint _destroy

  def __cinit__(self, bint alloc = True):
    if alloc:
      allocate_photochemdata(&self._ptr)
      self._destroy = True
    else:
      self._destroy = False

  def __dealloc__(self):
    if self._destroy:
      deallocate_photochemdata(&self._ptr)
      self._ptr = NULL
  
  property nq:
    def __get__(self):
      cdef int nq
      photochemdata_nq_get(&self._ptr, &nq)
      return nq
      
  property np:
    def __get__(self):
      cdef int val
      photochemdata_np_get(&self._ptr, &val)
      return val
      
  property ng:
    def __get__(self):
      cdef int val
      photochemdata_ng_get(&self._ptr, &val)
      return val

  property nsl:
    def __get__(self):
      cdef int val
      photochemdata_nsl_get(&self._ptr, &val)
      return val
      
  property nll:
    def __get__(self):
      cdef int val
      photochemdata_nll_get(&self._ptr, &val)
      return val
      
  property nsp:
    def __get__(self):
      cdef int val
      photochemdata_nsp_get(&self._ptr, &val)
      return val
      
  property nw:
    def __get__(self):
      cdef int val
      photochemdata_nw_get(&self._ptr, &val)
      return val

  property species_names:
    def __get__(self):
      cdef int dim1
      photochemdata_species_names_get_size(&self._ptr, &dim1)
      cdef ndarray species_names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      photochemdata_species_names_get(&self._ptr, &dim1, <char *>species_names_c.data)
      return c2stringarr(species_names_c, S_STR_LEN, dim1)
      
  property atoms_names:
    def __get__(self):
      cdef int dim1
      photochemdata_atoms_names_get_size(&self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      photochemdata_atoms_names_get(&self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, S_STR_LEN, dim1)
      
  property reaction_equations:
    def __get__(self):
      cdef int dim1
      photochemdata_reaction_equations_get_size(&self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*M_STR_LEN + 1, 'S1')
      photochemdata_reaction_equations_get(&self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, M_STR_LEN, dim1)
      
  property photonums:
    def __get__(self):
      cdef int dim1
      photochemdata_photonums_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.int32)
      photochemdata_photonums_get(&self._ptr, &dim1, <int *>arr.data)
      return arr
      
  property wavl:
    def __get__(self):
      cdef int dim1
      photochemdata_wavl_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemdata_wavl_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property species_mass:
    def __get__(self):
      cdef int dim1
      photochemdata_species_mass_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemdata_species_mass_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
      
  property species_redox:
    def __get__(self):
      cdef int dim1
      photochemdata_species_redox_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      photochemdata_species_redox_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

cdef c2stringarr(ndarray c_str_arr, int str_len, int arr_len):  
  bs = c_str_arr[:-1].tobytes()
  return [bs[i:i+str_len].decode().strip() for i in range(0, str_len*arr_len, str_len)]
    
    