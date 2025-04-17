# Photochem

`Photochem` is a photochemical and climate model of planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere, or evolve atmospheres through time. The code also contains 1-D climate models and a chemical equilibrium solver.

## Installation

### Option 1: Conda

`conda install -c conda-forge photochem`

### Option 2: From source

You need a Fortran compiler (`gfortran>=9.30`, [install instructions here](https://fortran-lang.org/learn/os_setup/install_gfortran)) and C compiler (e.g. install with `conda install -c conda-forge clang`)

Create a `conda` environment with all dependencies

```sh
conda create -n photochem -c conda-forge python numpy scipy pyyaml numba scikit-build cython cmake ninja pip hdf5 fypp
```

Clone this Gitub repository: 

```sh
git clone --depth=1 --recursive https://github.com/Nicholaswogan/photochem.git
```

Navigate to the root directory with a terminal, activate your new `conda` environment, then install with pip:

```sh
conda activate photochem
python -m pip install --no-deps --no-build-isolation .
```

## Examples/Tutorial

Check out the `examples` directory.

## History

In the 1980s Kevin Zahnle and Jim Kasting wrote the `Atmos` photochemical model in Fortran 77. An updated version of this code is maintained at [this link](https://github.com/VirtualPlanetaryLaboratory/atmos) by some excellent people at NASA Goddard. In December 2020, I reworked the `Atmos` photochemical model in Fortran 90, and made a Python wrapper to it using `numpy.f2py`. This resulted in [PhotochemPy](https://github.com/Nicholaswogan/PhotochemPy). PhotochemPy has several fundamental limitations that makes it challenging to build upon. So, starting Spring 2021, I began re-writing the model in Modern Fortran (2008) which resulted in this package.

## Contact

If you have questions email me: nicholaswogan at gmail

## Funding and Acknowledgements

Funding for the development of Photochem comes from
- [The Virtual Planetary Laboratory](https://depts.washington.edu/naivpl/content/welcome-virtual-planetary-laboratory)
- [Simons Collaboration on the Origins of Life](https://www.simonsfoundation.org/life-sciences/origins-of-life/simons-collaboration-on-the-origins-of-life/)
- [AEThER](https://planets.carnegiescience.edu/)

This model was built in collaboration with
- David Catling
- Kevin Zahnle
- Mark Claire
- Sandra Bastelberger
- Shawn Domagal-Goldman
- Josh Krissansen-Totton
