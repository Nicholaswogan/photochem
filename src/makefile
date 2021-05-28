
FC=gfortran
FCFLAGS=-g -fbounds-check -std=f2018#-Wunused-variable
SRC=photochem_const.f90 photochem_types.f90 photochem_io.f90 main.f90
OBJ=$(SRC:.f90=.o)

IDIR=dependencies/modules
LDIR=dependencies/lib
FFLAGS=-I$(IDIR) -L$(LDIR)
LIBS=-lyaml -lstringifor
OBJLIBS=dependencies/lib/binning.o

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS) 

photo.run: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $(OBJ) $(OBJLIBS) $(FFLAGS) $(LIBS) 
	
clean:
	@rm -f *.mod *.o photo.run