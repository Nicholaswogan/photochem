# Photochem [![codecov](https://codecov.io/gh/Nicholaswogan/Photochem/branch/main/graph/badge.svg?token=ZTCXVTG371)](https://codecov.io/gh/Nicholaswogan/Photochem) [![Build Status](https://app.travis-ci.com/Nicholaswogan/Photochem.svg?branch=main)](https://app.travis-ci.com/Nicholaswogan/Photochem)

`Photochem` is a photochemical model of planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere, or evolve atmospheres through time.

## Installation

**Requirements:**
- Python >= 3.7
- Fortran and C compiler. I suggest the GNU compiler collection (includes `gfortran`, `gcc`, etc.). If you are using a Mac, install it with Homebrew: `brew install gcc.`

**Python Module:** 
- Clone or download the Github repository.
- Navigate to the root directory, then install with `python -m pip install .`

**Fortran library:** 

If you prefer to use the code exclusively in Fortran, that is OK too. You can build `libphotochem` with CMake. Download or clone this repository, then from the root directory of the repository run

```sh
mkdir build
cd build
cmake ..
make -j
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

This model was build in collaboration with
- David Catling
- Kevin Zahnle
- Sandra Bastelberger
- Giada Arney
- Shawn Domagal-Goldman
- Josh Krissansen-Totton
