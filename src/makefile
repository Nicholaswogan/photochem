

FCC=gfortran
IDIR=dependencies/modules
LDIR=dependencies/lib
LIBS=-lyaml -lstringifor
FFLAGS=-I$(IDIR) -L$(LDIR)

hellomake:
	$(FCC) test.f90 $(FFLAGS) $(LIBS) -o test.run