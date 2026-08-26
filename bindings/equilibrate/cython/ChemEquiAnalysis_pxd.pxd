from libcpp cimport bool as cbool
cdef extern from "<stdbool.h>":
  pass

cdef extern from "_equilibrate.h":
  struct ChemEquiAnalysis:
    pass

# Allocate and destroy
cdef extern ChemEquiAnalysis *allocate_chemequianalysis();
cdef extern void deallocate_chemequianalysis(ChemEquiAnalysis *ptr);

# Wrappers for functions

cdef extern void chemequianalysis_create_wrapper(ChemEquiAnalysis *ptr, char *thermofile, 
                         cbool *atoms_present, int *atoms_dim, char *atoms,
                         cbool *species_present, int *species_dim, char *species,
                         char *err);
cdef extern void chemequianalysis_solve_wrapper(ChemEquiAnalysis *ptr, double *P, double *T,
                         cbool *molfracs_atoms_present, int *molfracs_atoms_dim, double *molfracs_atoms, 
                         cbool *molfracs_species_present, int *molfracs_species_dim, double *molfracs_species, 
                         cbool *converged, char *err)
cdef extern void chemequianalysis_solve_metallicity_wrapper(ChemEquiAnalysis *ptr, double *P, double *T,
                         double *metallicity, cbool *CtoO_present, double *CtoO, 
                         cbool *converged, char *err)

# Getters and setters

cdef extern void chemequianalysis_atoms_names_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_atoms_names_get(ChemEquiAnalysis *ptr, int *dim1, char* arr)

cdef extern void chemequianalysis_species_names_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_species_names_get(ChemEquiAnalysis *ptr, int *dim1, char* arr)
cdef extern void chemequianalysis_species_mass_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_species_mass_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_gas_names_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_gas_names_get(ChemEquiAnalysis *ptr, int *dim1, char* arr)
cdef extern void chemequianalysis_gas_mass_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_gas_mass_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_condensate_names_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_condensate_names_get(ChemEquiAnalysis *ptr, int *dim1, char* arr)
cdef extern void chemequianalysis_condensate_mass_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_condensate_mass_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_atoms_sun_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_atoms_sun_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)
cdef extern void chemequianalysis_molfracs_atoms_sun_set(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_atoms_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_atoms_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_species_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_species_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_massfracs_species_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_massfracs_species_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_atoms_gas_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_atoms_gas_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_species_gas_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_species_gas_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_atoms_condensate_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_atoms_condensate_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_species_condensate_get_size(ChemEquiAnalysis *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_species_condensate_get(ChemEquiAnalysis *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_mubar_get(ChemEquiAnalysis *ptr, double *val)
cdef extern void chemequianalysis_nabla_ad_get(ChemEquiAnalysis *ptr, double *val)
cdef extern void chemequianalysis_gamma2_get(ChemEquiAnalysis *ptr, double *val)
cdef extern void chemequianalysis_rho_get(ChemEquiAnalysis *ptr, double *val)
cdef extern void chemequianalysis_c_pe_get(ChemEquiAnalysis *ptr, double *val)

cdef extern void chemequianalysis_verbose_get(ChemEquiAnalysis *ptr, cbool *val)
cdef extern void chemequianalysis_verbose_set(ChemEquiAnalysis *ptr, cbool *val)

cdef extern void chemequianalysis_use_prev_guess_get(ChemEquiAnalysis *ptr, cbool *val)
cdef extern void chemequianalysis_use_prev_guess_set(ChemEquiAnalysis *ptr, cbool *val)

cdef extern void chemequianalysis_mass_tol_get(ChemEquiAnalysis *ptr, double *val)
cdef extern void chemequianalysis_mass_tol_set(ChemEquiAnalysis *ptr, double *val)
