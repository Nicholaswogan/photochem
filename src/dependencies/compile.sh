

cd ../cvode-5.7.0
mkdir build_dir
cd build_dir
cmake \
-DCMAKE_INSTALL_PREFIX=../../dependencies \
-DCMAKE_BUILD_TYPE=Release \
-DEXAMPLES_ENABLE_C=OFF \
-DEXAMPLES_ENABLE_F77=OFF \
-DEXAMPLES_ENABLE_F90=OFF \
-DEXAMPLES_ENABLE_F2003=OFF \
-DBUILD_FORTRAN77_INTERFACE=OFF \
-DBUILD_FORTRAN_MODULE_INTERFACE=ON \
-DBUILD_SHARED_LIBS=OFF \
-DCMAKE_Fortran_COMPILER=gfortran \
../
make install

cd ../../fortran-yaml
mkdir build_dir
cd build_dir
cmake \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_INSTALL_PREFIX=../../dependencies \
../
make install

cp -r modules ../../dependencies
mv ../../dependencies/fortran/* ../../dependencies/modules

# interpolation
cd ../../binning
gfortran -c binning.f -O3
gfortran -c interp.f90 -O3
ar rcs libbinning.a *.o
rm *.o
mv libbinning.a ../dependencies/lib

cd ../minpack
gfortran -c *.f -O3
ar rcs libminpack.a *.o
rm *.o
mv libminpack.a ../dependencies/lib



