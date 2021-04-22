
.SUFFIXES:
FC=gfortran

FCFLAGS = -g -fbounds-check

IDIR=dependencies/modules
LDIR=dependencies/lib
FFLAGS=-I$(IDIR) -L$(LDIR)
LIBS=-lyaml -lstringifor



hellomake: 
	$(FC) photochem_types.f90 photochem_io.f90 main.f90 $(FFLAGS) $(LIBS) -o test.run