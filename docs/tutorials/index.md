# Tutorials

## Summary

- [Rocky Planet Photochemistry](rocky_planet_photochemistry.py) demonstrates how to simulate the steady-state photochemistry of a rocky planet, using the modern Earth as an example.
- [Rocky Planet Climate](rocky_planet_climate.py) shows how to model the climate of a rocky planet like Earth. The tutorial covers constructing simple P-T profiles, computing their radiative properties, and computing full radiative-convective equilibrium.
- [Gas Giant Photochemistry](gas_giant_photochemistry.py) covers how to simulate the photochemistry of a gas-rich planet that has no surface, using WASP-39b as an example.

## How to run these tutorials

The tutorials on this website have already been executed. To modify parameters and run the calculations yourself, install Photochem and the tutorial dependencies described in the [installation guide](../installation.md).

Clone the Photochem repository, enter the tutorials directory, and start Jupyter:

```sh
git clone https://github.com/Nicholaswogan/photochem.git
cd photochem/docs/tutorials
jupyter notebook
```

The tutorials are stored as Jupytext `.py` files. To open one as a notebook, right-click the tutorial file in the Jupyter file browser and select **Open With → Notebook**.
