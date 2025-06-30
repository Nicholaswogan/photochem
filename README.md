# Photochem

`Photochem` is a photochemical and climate model of planet's atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the steady-state chemical composition of an atmosphere, or evolve atmospheres through time. The code also contains 1-D climate models and a chemical equilibrium solver.

## Installation

`Photochem` can be installed with conda:

```sh
conda install -c conda-forge photochem
```

For more installation instruction see the "Documentation" section below.

## Documentation

`Photochem` does not have a formal documentation website. A website will soon be built and linked here. In the mean time, the best available documentation is the following tutorials from a `Photochem` workshop at the University of Arizona: https://github.com/Nicholaswogan/UofA_Photochem_Workshop. Furthermore, there are a few tutorials in the `examples` directory of this repository.

## History

In the 1980s Kevin Zahnle and Jim Kasting wrote the `Atmos` photochemical model in Fortran 77. An updated version of this code is maintained at [this link](https://github.com/VirtualPlanetaryLaboratory/atmos) by some excellent people at NASA Goddard. In December 2020, I reworked the `Atmos` photochemical model in Fortran 90, and made a Python wrapper to it using `numpy.f2py`. This resulted in [PhotochemPy](https://github.com/Nicholaswogan/PhotochemPy). PhotochemPy has several fundamental limitations that makes it challenging to build upon. So, starting Spring 2021, I began re-writing the model in Modern Fortran (2008) which resulted in this package.
