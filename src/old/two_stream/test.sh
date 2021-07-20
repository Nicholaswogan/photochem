gfortran two_stream.f90 -O3 -llapack -lblas -o photo
gfortran two_stream.f90 -O3 -ffast-math -llapack -lblas -o photof
./photo
./photo
./photof
./photof

time ./photo
time ./photof

# f2py -c two_stream.f90 -m two_stream --opt="-O3" -llapack -lblas
# f2py -c two_stream.f90 -m two_streamf --opt="-O3 -ffast-math" -llapack -lblas

