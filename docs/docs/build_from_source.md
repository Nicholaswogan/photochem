# Building from source

Building Photochem from source can be challenging, which is why `conda` is the [preferred installation route](installation.md). But in order to change the underlying Fortran, you need to be able to build the code from source.

First, lets create a new conda environment called `test` with all the photochem dependencies

```sh
conda create -n test -c conda-forge photochem_clima_data python numpy scipy pyyaml numba astropy h5py threadpoolctl scikit-build cmake ninja cython fypp pip c-compiler cxx-compiler fortran-compiler git jupyter matplotlib

conda activate test
```

Here is a breakdown of the dependencies:

- `photochem_clima_data` - All model data (e.g., photolysis cross sections etc.) 
- `python numpy scipy pyyaml numba astropy h5py threadpoolctl` - Packages needed to run Photochem
- `scikit-build cmake ninja cython fypp pip c-compiler cxx-compiler fortran-compiler` - Packages needed to build Photochem. Critically, this includes C and Fortran compilers.
- `git` - Needed to download the Photochem source from the internet
- `jupyter matplotlib` - Only needed to run notebooks/plotting

Next, we need to download Photochem from the internet with `git`.

```sh
git clone https://github.com/Nicholaswogan/photochem.git
cd photochem
```

We need to set some environment variables to help the build find libraries.

```sh
export CMAKE_ARGS="-DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DCMAKE_POSITION_INDEPENDENT_CODE=ON"
export CONDA_PREFIX_SAVE=$CONDA_PREFIX
unset CONDA_PREFIX
```

Finally, we can install the code, then reset your conda prefix.

```sh
python -m pip install --no-deps --no-build-isolation . -v
export CONDA_PREFIX=$CONDA_PREFIX_SAVE
```

## Building in place without `pip`

It is also useful to build the code without installing it. Run the following commands from the root directory of the Photochem repository.

```sh
mkdir build
cd build

export CONDA_PREFIX_SAVE=$CONDA_PREFIX
unset CONDA_PREFIX

# Configure
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX_SAVE -DBUILD_PYTHON_PHOTOCHEM=ON -DBUILD_WITH_OPENMP=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_BUILD_TYPE=Release -DSKBUILD_CMAKE_MODULE_DIR=$(python -c "from skbuild import __file__; print(__file__.strip('__init__.py')+'resources/cmake')")

# Build and install
cmake --build . -j && cmake --install .
```

New binary filles will have now appeared in the `photochem/` directory. So now, you can run python code at the root of the `photochem` directory, and it will use the binaries you just built.