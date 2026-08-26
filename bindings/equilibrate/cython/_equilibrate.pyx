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
  """Chemical-equilibrium solver and state from its most recent solve.

  Name and mass properties are established during construction. Composition
  and thermodynamic result properties are updated by ``solve`` and
  ``solve_metallicity``; check the returned convergence flag before using
  them. Array-valued properties return copies.
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
    """Initialize a chemical-equilibrium solver from thermodynamic data.

    Python accepts YAML thermodynamic files. Supplying only ``atoms`` selects
    every compatible species; supplying only ``species`` selects every atom in
    those species; and supplying neither uses the entire file. The two
    selections are mutually exclusive in Python.

    Parameters
    ----------
    thermofile : str
        Path to a YAML thermodynamic file.
    atoms : sequence of str, optional
        Nonempty atom-name selection.
    species : sequence of str, optional
        Nonempty species-name selection.

    Raises
    ------
    EquilibrateException
        If the file or selections are invalid, or an entry exceeds the
        supported 20-character length.
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
    """Compute chemical equilibrium from atom or species mole fractions.

    Supply exactly one composition array. It is normalized internally. On a
    completed solve, the result properties contain the latest equilibrium
    state and the return value reports whether the iteration converged.

    Parameters
    ----------
    P : float
        Pressure (dyn/cm^2).
    T : float
        Temperature (K).
    molfracs_atoms : ndarray, shape (na,), optional
        Nonnegative atom mole fractions in ``atoms_names`` order.
    molfracs_species : ndarray, shape (ns,), optional
        Nonnegative species mole fractions in ``species_names`` order.

    Returns
    -------
    bool
        Whether the equilibrium iteration met its tolerances.

    Raises
    ------
    EquilibrateException
        If the inputs are invalid or the compiled solver reports an error.
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
    """Compute chemical equilibrium from metallicity and an optional C/O ratio.

    Heavy-element abundances in ``molfracs_atoms_sun`` are scaled relative to
    H and He. If ``CtoO`` is supplied, carbon and oxygen are redistributed while
    preserving their combined abundance. Result properties are updated as in
    ``solve``.

    Parameters
    ----------
    P : float
        Pressure (dyn/cm^2).
    T : float
        Temperature (K).
    metallicity : float
        Positive metallicity relative to the reference solar composition.
    CtoO : float, optional
        Positive C/O ratio relative to the reference solar value. A value of 1
        preserves solar C/O.

    Returns
    -------
    bool
        Whether the equilibrium iteration met its tolerances.

    Raises
    ------
    EquilibrateException
        If metallicity or C/O is nonpositive, required C or O atoms are absent,
        or the compiled solver reports an error.
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
    """list[str]: Atom names in the ordering used by atom-indexed arrays."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_atoms_names_get_size(self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_atoms_names_get(self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property species_names:
    """list[str]: All species names in the ordering used by full-species arrays."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_species_names_get_size(self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_species_names_get(self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property species_mass:
    """ndarray, shape (ns,): Read-only copy of species molar masses (g/mol)."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_species_mass_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_species_mass_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property gas_names:
    """list[str]: Gas names in the ordering used by gas-phase arrays."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_gas_names_get_size(self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_gas_names_get(self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property gas_mass:
    """ndarray, shape (ng,): Read-only copy of gas molar masses (g/mol)."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_gas_mass_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_gas_mass_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property condensate_names:
    """list[str]: Condensate names in condensed-phase array ordering."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_condensate_names_get_size(self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_condensate_names_get(self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property condensate_mass:
    """ndarray, shape (nc,): Read-only copy of condensate molar masses (g/mol)."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_condensate_mass_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_condensate_mass_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_atoms_sun:
    """ndarray, shape (na,): Writable reference solar atom mole fractions.

    Values follow ``atoms_names`` order. Getting this property returns a copy;
    assign the complete array to update the reference used by
    ``solve_metallicity``.
    """
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
    """ndarray, shape (na,): Read-only copy of atom mole fractions across all phases."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species:
    """ndarray, shape (ns,): Read-only copy of species mole fractions across all phases."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property massfracs_species:
    """ndarray, shape (ns,): Read-only copy of species mass fractions across all phases."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_massfracs_species_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_massfracs_species_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_atoms_gas:
    """ndarray, shape (na,): Read-only copy of gas-phase atom mole fractions."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_gas_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_gas_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species_gas:
    """ndarray, shape (ng,): Read-only copy of gas-phase species mole fractions."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_gas_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_gas_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_atoms_condensate:
    """ndarray, shape (na,): Read-only copy of condensed-phase atom mole fractions."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_condensate_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_condensate_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species_condensate:
    """ndarray, shape (nc,): Read-only copy of condensed species mole fractions."""
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_condensate_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_condensate_get(self._ptr, &dim1, <double *>arr.data)
      return arr

  property mubar:
    """float: Mean molecular weight of the equilibrium gas phase (g/mol)."""
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_mubar_get(self._ptr, &val)
      return val

  property nabla_ad:
    """float: Adiabatic logarithmic temperature gradient, ``dln(T)/dln(P)``."""
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_nabla_ad_get(self._ptr, &val)
      return val

  property gamma2:
    """float: Second adiabatic exponent, ``1 / (1 - nabla_ad)``."""
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_gamma2_get(self._ptr, &val)
      return val

  property rho:
    """float: Equilibrium gas mass density (g/cm^3)."""
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_rho_get(self._ptr, &val)
      return val

  property c_pe:
    """float: Equilibrium specific heat at constant pressure (erg/(g K))."""
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_c_pe_get(self._ptr, &val)
      return val

  property verbose:
    """bool: Whether the equilibrium solver prints iteration details."""
    def __get__(self):
      cdef cbool val
      cea_pxd.chemequianalysis_verbose_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      cea_pxd.chemequianalysis_verbose_set(self._ptr, &val)

  property use_prev_guess:
    """bool: Whether to initialize from the previous converged solution."""
    def __get__(self):
      cdef cbool val
      cea_pxd.chemequianalysis_use_prev_guess_get(self._ptr, &val)
      return val
    def __set__(self, cbool val):
      cea_pxd.chemequianalysis_use_prev_guess_set(self._ptr, &val)

  property mass_tol:
    """float: Dimensionless mass-balance convergence tolerance."""
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_mass_tol_get(self._ptr, &val)
      return val
    def __set__(self, double val):
      cea_pxd.chemequianalysis_mass_tol_set(self._ptr, &val)

# utils
cdef pystring2cstring(str pystring):
  # Add a null C character and convert to bytes.
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
