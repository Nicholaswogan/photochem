gfortran \
photochem_const.f90 \
photochem_types.f90 \
photochem_vars.f90 \
photochem_input.f90 \
photochem_wrk.f90 \
photochem_data.f90 \
photochem_setup.f90 \
photochem_radtran.f90 \
photochem.f90 \
-shared -fPIC -o photochem_jl.so -O3 \
-Idependencies/modules -Ldependencies/lib \
-lyaml -lbinning -lsundials_fcvode_mod -lsundials_cvode
    
