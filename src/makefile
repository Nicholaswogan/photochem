
FC=gfortran
FCFLAGS= -g -std=f2008 -O3#-Wunused-variable#  -fcheck=all -fopenmp -O3
SRC=photochem_const.f90 \
		photochem_types.f90 \
		photochem_vars.f90 \
		photochem_input.f90 \
		photochem_wrk.f90 \
		photochem_data.f90 \
		photochem_setup.f90 \
		photochem_radtran.f90 \
		photochem.f90 \
		main.f90
OBJ=$(SRC:.f90=.o)

IDIR=dependencies/modules
LDIR=dependencies/lib
FFLAGS=-I$(IDIR) -L$(LDIR)
LIBS=-lyaml -lbinning -lsundials_fcvode_mod -lsundials_cvode -llapack

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FFLAGS) $(LIBS)

photo.run: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $(OBJ) $(OBJLIBS) $(FFLAGS) $(LIBS)

clean:
	@rm -f *.mod *.o photo.run
