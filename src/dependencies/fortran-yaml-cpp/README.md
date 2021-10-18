# fortran-yaml-cpp

This is [YAML](http://yaml.org) parser for Fortran matching the [YAML 1.2 spec](https://yaml.org/spec/1.2.2/).

This package uses the C++ package [yaml-cpp](https://github.com/jbeder/yaml-cpp) to parse yaml documents, then stores the data in Fortran derived types created by [fortran-yaml](https://github.com/BoldingBruggeman/fortran-yaml). Hence the name `fortran-yaml-cpp`.

## Building

First, clone this repository

```
git clone --recursive https://github.com/Nicholaswogan/fortran-yaml-cpp.git
```

Next, navigate to the root directory of `fortran-yaml-cpp`. Finally, run the following commands to build with `cmake` and run the test and example.

```
mkdir build
cd build
cmake ..
make
./test_yaml
./example
```