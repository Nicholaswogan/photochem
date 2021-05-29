

cd ../cvode-5.7.0
mkdir build_dir
cd build_dir
cmake \
-DBUILD_EXAMPLES=OFF \
-DCMAKE_INSTALL_PREFIX=../../dependencies \
-DEXAMPLES_ENABLE_F77=ON \
-DBUILD_FORTRAN77_INTERFACE=ON \
-DBUILD_SHARED_LIBS=OFF \
-DCMAKE_Fortran_COMPILER=gfortran \
../
make install

cd ../../fortran-yaml
mkdir build_dir
cd build_dir
cmake \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_INSTALL_PREFIX=../../dependencies \
../
make install

cp -r modules ../../dependencies

# Download and build Stringifor
cd ../../
python -m pip install FoBiS.py
git clone https://github.com/szaghi/StringiFor
cd StringiFor
git submodule update --init
FoBiS.py build -mode stringifor-static-gnu

cp lib/libstringifor.a ../dependencies/lib
cp lib/mod/* ../dependencies/modules

# interpolation
cd ../binning
gfortran -c binning.f -O3
gfortran -c interp.f90 -O3
ar rcs libbinning.a *.o
rm *.o
mv libbinning.a ../dependencies/lib



