# Photochem

`Photochem` is a photochemical model of planet's atmospheres. It is currently under construction. 

# Compiling and testing

You need `CMake`, a Fortran compiler, and a C compiler. You can then compile `libphotochem` with the following commands.

```sh
mkdir build
cd build
cmake ..
make
```

After building, you can run the test, which is a model run of the modern Earth.

```sh
ctest # runs the test model with ctest
./test.run # runs the test case
```