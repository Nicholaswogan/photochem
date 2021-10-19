
cdef extern void allocate_photochemdata(void *ptr)
cdef extern void deallocate_photochemdata(void *ptr)

cdef extern void photochemdata_nq_get(void *ptr, int *nq)

cdef extern void photochemdata_species_names_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_species_names_get(void *ptr, int *dim1, char* species_names)

cdef extern void photochemdata_atoms_names_get_size(void *ptr, int *dim1)
cdef extern void photochemdata_atoms_names_get(void *ptr, int *dim1, char* names)

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
      
  property species_names:
    def __get__(self):
      cdef int dim1
      photochemdata_species_names_get_size(&self._ptr, &dim1)
      cdef ndarray species_names_c = np.empty(dim1*S_STR_LEN + 1, 'c')
      photochemdata_species_names_get(&self._ptr, &dim1, <char *>species_names_c.data)
      return c2stringarr(species_names_c, S_STR_LEN, dim1)
      
  property atoms_names:
    def __get__(self):
      cdef int dim1
      photochemdata_atoms_names_get_size(&self._ptr, &dim1)
      cdef ndarray names_c = np.empty(dim1*S_STR_LEN + 1, 'c')
      photochemdata_atoms_names_get(&self._ptr, &dim1, <char *>names_c.data)
      return c2stringarr(names_c, S_STR_LEN, dim1)


cdef c2stringarr(ndarray c_str_arr, int str_len, int arr_len):
  cdef int i, j, k
  tmp = [' ' for i in range(str_len)]
  str_arr = []
  for i in range(arr_len):
      for j in range(str_len):
          k = j + i * str_len
          tmp[j] = c_str_arr[k].decode('utf-8')
      str_arr.append(''.join(tmp).strip())
  return str_arr
    
    
    
    