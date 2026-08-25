from numpy cimport ndarray, uint8_t
from libc.stdint cimport uintptr_t
from libcpp cimport bool as cbool
from cpython.object cimport PyObject_GenericSetAttr
cimport ChemEquiAnalysis_pxd as cea_pxd
import numpy as np
import ctypes as ct
import os

DEF S_STR_LEN = 20;
DEF ERR_LEN = 1024;

cdef class ChemEquiAnalysis:
  """A chemical equilibrium solver.
  """

  cdef cea_pxd.ChemEquiAnalysis *_ptr
  cdef cbool _init_called

  def __cinit__(self, *args, **kwargs):
    self._init_called = False
    self._ptr = cea_pxd.allocate_chemequianalysis()

  def __dealloc__(self):
    cea_pxd.deallocate_chemequianalysis(self._ptr)

  def __getattribute__(self, name):
    if not self._init_called:
      raise EquilibrateException('The "__init__" method of ChemEquiAnalysis has not been called.')
    return super().__getattribute__(name)

  def __setattr__(self, name, value):
    if not self._init_called:
      raise EquilibrateException('The "__init__" method of ChemEquiAnalysis has not been called.')
    PyObject_GenericSetAttr(self, name, value)

  def __init__(self, str thermofile, atoms = None, species = None):           
    """Initializes the chemical equilibrium solver given an input thermodynamic file.
    The file can only have ".yaml" format. `atoms` and `species` are optional inputs 
    with the following effects:
    - If `atoms` are specified but not `species`, then the code will consider all 
      species with the input `atoms`.
    - If `species` are specified but not `atoms`, then the code will consider all
      atoms corresponding to the input `species`.
    - If neither `atoms` or `species` are specified, then the code will consider all
      atoms and species in the input file.

    Parameters
    ----------
    thermofile : str
        Path to file describing the thermodynamic data of all species.
    atoms : list, optional
        Names of atoms to include.
    species : list, optional
        Names of species to include.
    """

    self._init_called = True

    # convert strings to char
    cdef bytes thermofile_b = pystring2cstring(thermofile)
    cdef char *thermofile_c = thermofile_b
    cdef char err[ERR_LEN+1]

    cdef cbool atoms_present = False
    cdef int atoms_dim = 1
    cdef ndarray atoms_c = np.zeros(atoms_dim*S_STR_LEN + 1, 'S1')
    if atoms is not None:
      atoms_present = True
      atoms_dim = len(atoms)
      atoms_c = list2cstring(atoms, S_STR_LEN)

    cdef cbool species_present = False
    cdef int species_dim = 1
    cdef ndarray species_c = np.zeros(species_dim*S_STR_LEN + 1, 'S1')
    if species is not None:
      species_present = True
      species_dim = len(species)
      species_c = list2cstring(species, S_STR_LEN)

    if atoms_present and species_present:
      raise EquilibrateException('atoms and species cannot both be provided.')
    
    # Initialize
    cea_pxd.chemequianalysis_create_wrapper(self._ptr, thermofile_c,
                                            &atoms_present, &atoms_dim, <char *>atoms_c.data,
                                            &species_present, &species_dim, <char *>species_c.data,
                                            err)
    if len(err.strip()) > 0:
      raise EquilibrateException(err.decode("utf-8").strip())

  def solve(self, double P, double T, molfracs_atoms = None, molfracs_species = None):
    """Computes chemical equilibrium given input atom or species mole fractions.
    If successful, then the equilibrium composition will be stored in a number
    of attributes (e.g., self%molfracs_species).

    Parameters
    ----------
    P : double
        Pressure in dynes/cm^2
    T : double
        Temperature in Kelvin
    molfracs_atoms : ndarray[double,ndim=1], optional
        Atom mole fractions in the same order and length as self.atoms_names.
    molfracs_species : ndarray[double,ndim=1], optional
        Species mole fractions in the same order and length as self.species_names.

    Results
    -------
    converged : bool
        If true, then the calculation successfully achieved chemical equilibrium
        to within the specified tolerances.
    """

    cdef ndarray[double, ndim=1] molfracs_atoms_ = np.empty(1,dtype=np.double)
    cdef cbool molfracs_atoms_present = False
    if molfracs_atoms is not None:
      molfracs_atoms_present = True
      molfracs_atoms_ = molfracs_atoms
    cdef int molfracs_atoms_dim = molfracs_atoms_.shape[0]

    cdef ndarray[double, ndim=1] molfracs_species_ = np.empty(1,dtype=np.double)
    cdef cbool molfracs_species_present = False
    if molfracs_species is not None:
      molfracs_species_present = True
      molfracs_species_ = molfracs_species
    cdef int molfracs_species_dim = molfracs_species_.shape[0]

    cdef cbool converged
    cdef char err[ERR_LEN+1]

    cea_pxd.chemequianalysis_solve_wrapper(self._ptr, &P, &T,
                         &molfracs_atoms_present, &molfracs_atoms_dim, <double *>molfracs_atoms_.data, 
                         &molfracs_species_present, &molfracs_species_dim, <double *>molfracs_species_.data, 
                         &converged, err)

    if len(err.strip()) > 0:
      raise EquilibrateException(err.decode("utf-8").strip())

    return converged

  def solve_metallicity(self, double P, double T, double metallicity, CtoO = None):
    """Computes chemical equilibrium given an input metallicity and, optionally,
    a C/O ratio. If successful, then the equilibrium composition will be stored in a number
    of attributes (e.g., self%molfracs_species).

    Parameters
    ----------
    P : double
        Pressure in dynes/cm^2
    T : double
        Temperature in Kelvin
    metallicity : double
        Metallicity relative to the Sun
    CtoO : double, optional
        The C/O ratio relative to solar. CtoO = 1 would be the same composition as solar.

    Results
    -------
    converged : bool
        If true, then the calculation successfully achieved chemical equilibrium
        to within the specified tolerances.
    """

    cdef double CtoO_ = 1.0
    cdef cbool CtoO_present = False
    if CtoO is not None:
      CtoO_present = True
      CtoO_ = CtoO

    cdef cbool converged
    cdef char err[ERR_LEN+1]

    cea_pxd.chemequianalysis_solve_metallicity_wrapper(self._ptr, &P, &T,
                         &metallicity, &CtoO_present, &CtoO_, 
                         &converged, err)

    if len(err.strip()) > 0:
      raise EquilibrateException(err.decode("utf-8").strip())

    return converged

  property atoms_names:
    "List. Names of atoms"
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_atoms_names_get_size(self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_atoms_names_get(self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property species_names:
    "List. Names of species"
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_species_names_get_size(self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_species_names_get(self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property gas_names:
    "List. Names of gases"
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_gas_names_get_size(self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_gas_names_get(self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property condensate_names:
    "List. Names of condensates"
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_condensate_names_get_size(self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_condensate_names_get(self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property molfracs_atoms_sun:
    "ndarray[double,ndim=1]. Assumed mole fractions of each atom in the Sun."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_sun_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_sun_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self,ndarray[double, ndim=1] arr):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_sun_get_size(self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise EquilibrateException(
          'molfracs_atoms_sun must have shape (number_of_atoms,).'
        )
      cea_pxd.chemequianalysis_molfracs_atoms_sun_set(self._ptr, &dim1, <double *>arr.data)

  property molfracs_atoms:
    "ndarray[double,ndim=1]. Mole fractions of each atom."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species:
    "ndarray[double,ndim=1]. Mole fractions of each species."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property massfracs_species:
    "ndarray[double,ndim=1]. Mass fractions of each species."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_massfracs_species_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_massfracs_species_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_atoms_gas:
    "ndarray[double,ndim=1]. Mole fractions of atoms in gas phase."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_gas_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_gas_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species_gas:
    "ndarray[double,ndim=1]. Mole fractions of species in gas phase."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_gas_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_gas_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_atoms_condensate:
    "ndarray[double,ndim=1]. Mole fractions of atoms in condensed phase."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_condensate_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_condensate_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species_condensate:
    "ndarray[double,ndim=1]. Mole fractions of species in condensed phase."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_condensate_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_condensate_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property mubar:
    "float. Mean molecular weight (g/mol)."
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_mubar_get(self._ptr, &val)
      return val

  property verbose:
    "bool. Determines amount of printing."
    def __get__(self):
      cdef cbool val
      cea_pxd.chemequianalysis_verbose_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      cea_pxd.chemequianalysis_verbose_set(self._ptr, &val)

  property use_prev_guess:
    "bool. Use previous converged solution as initial guess."
    def __get__(self):
      cdef cbool val
      cea_pxd.chemequianalysis_use_prev_guess_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      cea_pxd.chemequianalysis_use_prev_guess_set(self._ptr, &val)

  property mass_tol:
    """float. Degree to which mass will be balanced. Gordon & McBride's default 
    is 1.0e-6, but it seems like 1.0e-2 is OK.
    """
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_mass_tol_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      cea_pxd.chemequianalysis_mass_tol_set(self._ptr, &val)

# utils
cdef pystring2cstring(str pystring):
  # add a null c char, and convert to byes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring

cdef c2stringarr(ndarray c_str_arr, int str_len, int arr_len):  
  bs = c_str_arr[:-1].tobytes()
  return [bs[i:i+str_len].decode().strip() for i in range(0, str_len*arr_len, str_len)]

cdef list2cstring(list arr, int str_len):
  arr_c = np.zeros(len(arr)*str_len + 1, 'S1')
  for i in range(len(arr)):
    if len(arr[i]) > str_len:
          raise EquilibrateException(
            'A list entry exceeds the supported string length.'
          )
    arr_c[i*str_len:(i+1)*str_len] = b' '
    arr_c[i*str_len:i*str_len+len(arr[i])] = np.array([elem.encode('utf-8') for elem in arr[i]])
  return arr_c

class EquilibrateException(Exception):
    """Error reported by the compiled Equilibrate model or its Python wrapper."""
