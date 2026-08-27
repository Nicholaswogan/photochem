# Building from Source

Most users should install Photochem from conda-forge as described on the [Installation](../installation.md) page. A source build is useful when you want to modify the Python, Cython, C, or Fortran implementation.

There are two supported source-build workflows. Installing through Python builds the compiled extensions and installs the Python package into the active environment. Building directly with CMake exposes all build options and also builds the Fortran tests and executables.

## Get the source code

Clone the repository and enter its root directory:

```sh
git clone https://github.com/Nicholaswogan/photochem.git
cd photochem
```

These instructions build the current development version. To build a released version instead, check out its Git tag before continuing.

## Create the development environment

The repository includes `environment-dev.yml` with the runtime libraries, build tools, compilers, and model data needed by both workflows:

```sh
mamba env create -f environment-dev.yml
conda activate photochem-dev
```

The environment uses the platform-appropriate C and C++ compilers supplied by conda-forge and explicitly installs GNU Fortran. Point the build system to `gfortran` after activating the environment:

```sh
export FC="$CONDA_PREFIX/bin/gfortran"
```

Photochem requires `gfortran` 14 or newer.

## Option 1: Install the Python package from source

This is the preferred workflow if you plan to use Photochem through Python. From the repository root, tell CMake to search the active Conda environment and install the package without creating a separate build environment:

```sh
export CMAKE_ARGS="-DCMAKE_PREFIX_PATH=$CONDA_PREFIX"
python -m pip install --no-deps --no-build-isolation . -v
```

The command invokes the CMake build through `scikit-build-core`, compiles Photochem with OpenMP support, and installs the resulting Python package into `photochem-dev`.

Verify the installation from outside the repository so that Python does not import the source tree in place of the installed package:

```sh
cd ..
python -c "import photochem; print(photochem.__version__)"
cd photochem
```

Run the installation command again after changing compiled Cython, C, or Fortran code. Add `--force-reinstall` if `pip` reports that the current version is already installed.

## Option 2: Build everything directly with CMake

The direct CMake workflow builds the Fortran libraries and tests as well as the Python extensions. Configure an out-of-source build from the repository root:

```sh
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DPython_EXECUTABLE="$CONDA_PREFIX/bin/python" \
  -DPython3_EXECUTABLE="$CONDA_PREFIX/bin/python" \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DBUILD_EXECUTABLES=ON \
  -DBUILD_PYTHON_PHOTOCHEM=ON \
  -DBUILD_WITH_OPENMP=ON
```

Compile all configured targets, then copy the compiled Python extensions into the source-tree `photochem/` package:

```sh
cmake --build build -j
cmake --install build
```

This CMake installation step does not create a conventional Python package installation in the Conda environment. It enables in-place Python use from the repository root, which is convenient while editing the source.

Run representative checks with:

```sh
./build/tests/clima/test_adiabat
./build/tests/clima/test_radtran
./build/tests/equilibrate/test_equilibrate
./build/tests/photochem/test_api
python -c "import photochem; print(photochem.__version__)"
```

## Troubleshooting

If CMake cannot find a compiler or library, first confirm that `photochem-dev` is active, that `FC` points to `gfortran` in that environment, and that `CMAKE_PREFIX_PATH` points to `$CONDA_PREFIX`.

CMake stops with an explicit error if GNU Fortran is older than version 14. Recreate or update the development environment rather than trying to build with the older compiler.

HDF5 is supplied through the Conda environment. If HDF5 is not found, confirm that the environment is active and configure a fresh build directory with `-DCMAKE_PREFIX_PATH="$CONDA_PREFIX"`.

CMake caches compiler and dependency locations. Use a fresh build directory after changing environments, compilers, or important configuration options.
