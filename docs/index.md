# Photochem

**One-dimensional photochemistry, climate, and chemical equilibrium for planetary atmospheres.**

Photochem is an open-source modeling package for exploring how radiation, chemistry, transport, and climate shape planetary atmospheres. Its Python interfaces can evolve atmospheric composition through time, find chemical steady states, calculate radiative-convective climates, and determine thermochemical-equilibrium compositions.

[Install Photochem](installation.md){ .md-button .md-button--primary } [Start with a tutorial](tutorials/index.md){ .md-button } [Browse the API](reference/index.md){ .md-button }

## Choose a modeling workflow

<div class="grid cards" markdown>

-   **Rocky-planet photochemistry**

    Use `photochem.EvoAtmosphere` to model chemistry and vertical transport on an altitude grid, with surface and upper-atmosphere boundary conditions. Integrate an evolving atmosphere or seek a chemical steady state.

-   **Gas-giant photochemistry**

    Use `EvoAtmosphereGasGiant` from `photochem.extensions` for pressure-grid models of gas-rich planets and brown dwarfs, including deep atmospheres initialized from thermochemical equilibrium.

-   **Climate**

    Use `photochem.clima.AdiabatClimate` to calculate one-dimensional radiative-convective temperature structures, including atmospheres with multiple condensible species.

-   **Chemical equilibrium**

    Use `photochem.equilibrate.ChemEquiAnalysis` to calculate equilibrium compositions from elemental abundances at specified temperatures and pressures.

</div>

## Get started

- [Installation](installation.md) covers the supported ways to install Photochem and verify the installation.
- [Tutorials](tutorials/index.md) develop complete atmospheric-modeling workflows in executable notebooks.
- [Input Files](input-files/index.md) describes reaction mechanisms, model settings, stellar spectra, atmosphere profiles, and equilibrium species.
- [API Reference](reference/index.md) documents the public Python classes and functions.

## Capabilities

- Time-dependent integration and steady-state solutions for atmospheric photochemistry.
- Flexible chemical mechanisms and thermodynamic data stored in readable input files.
- Photolysis, gas and particle transport, condensation and evaporation, rainout, and configurable atmospheric boundary conditions.
- Separate workflows suited to rocky planets and gas-rich planets.
- Radiative-convective climate calculations and thermochemical-equilibrium chemistry in the same Python package.
- Modern Fortran and C solvers exposed through a Python interface.

## Scientific applications

Photochem is designed for comparative planetary-atmosphere studies spanning Solar System worlds, early terrestrial environments, and exoplanets. Example applications include:

- modern and ancient rocky-planet atmospheres;
- origin-of-life chemistry following large impacts on the early Earth;
- atmospheric biosignatures and habitability;
- temperate mini-Neptunes, hot Jupiters, and brown dwarfs; and
- preparing pressure-temperature-composition profiles for spectrum models.

Published applications include modeling [post-impact chemistry on the early Earth](https://doi.org/10.3847/PSJ/aced83) and interpreting [JWST observations of K2-18b](https://doi.org/10.3847/2041-8213/ad2616).

## Project and citation

Photochem is developed openly on [GitHub](https://github.com/Nicholaswogan/photochem). Bug reports, feature requests, and contributions are welcome through the repository.

If you use Photochem in published research, please cite the scientific work most relevant to your application. Preliminary citation guidance is available in the [repository README](https://github.com/Nicholaswogan/photochem#citation); a complete **How to Cite Photochem** page will be added before the v1.0 documentation release.

!!! note "Documentation status"

    The v1.0 documentation is under active development. Pages marked as under construction establish stable navigation destinations for material added in later documentation passes.
