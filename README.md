# Photochem [![codecov](https://codecov.io/gh/Nicholaswogan/Photochem/branch/main/graph/badge.svg?token=ZTCXVTG371)](https://codecov.io/gh/Nicholaswogan/Photochem) [![Build Status](https://app.travis-ci.com/Nicholaswogan/Photochem.svg?branch=main)](https://app.travis-ci.com/Nicholaswogan/Photochem)

**NOTE: Photochem has not been published in a peer-reviewed journal yet. Please do not use this model in a paper that will be published. Wait for me to publish the model first, so you are able to cite it. Email me if you have questions (wogan@uw.edu)**

`Photochem` is a photochemical model of planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere, or evolve atmospheres through time.

## Installation

You need a Fortran compiler (`gfortran>=9.30`, [install instructions here](https://fortran-lang.org/learn/os_setup/install_gfortran)) and C compiler (e.g. install with `conda install -c conda-forge clang`)

Create a `conda` environment with all dependencies

```sh
conda create -n photochem -c conda-forge python numpy scipy pyyaml numba scikit-build cython
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

**Fortran library only:** 

You can build `libphotochem` with CMake (install CMake with `conda install -c anaconda cmake`). Download or clone this repository, then from the root directory of the repository run

```sh
mkdir build
cd build
cmake ..
cmake --build . -j
```

## Examples/Tutorial

Check out the `examples` directory.

## History

In the 1980s Kevin Zahnle and Jim Kasting wrote the `Atmos` photochemical model in Fortran 77. An updated version of this code is maintained at [this link](https://github.com/VirtualPlanetaryLaboratory/atmos) by some excellent people at NASA Goddard. In December 2020, I reworked the `Atmos` photochemical model in Fortran 90, and made a Python wrapper to it using `numpy.f2py`. This resulted in [PhotochemPy](https://github.com/Nicholaswogan/PhotochemPy). PhotochemPy has several fundamental limitations that makes it challenging to build upon. So, starting Spring 2021, I began re-writing the model in Modern Fortran (2008) which resulted in this package.

## Contact

If you have questions email me: wogan@uw.edu

## Funding and Acknowledgements

Funding for the development of Photochem comes from
- [The Virtual Planetary Laboratory](https://depts.washington.edu/naivpl/content/welcome-virtual-planetary-laboratory)
- [Simons Collaboration on the Origins of Life](https://www.simonsfoundation.org/life-sciences/origins-of-life/simons-collaboration-on-the-origins-of-life/)
- [AEThER](https://planets.carnegiescience.edu/)

This model was build in collaboration with
- David Catling
- Kevin Zahnle
- Sandra Bastelberger
- Shawn Domagal-Goldman
- Josh Krissansen-Totton
