# Photochem

**One-dimensional photochemistry, climate, and chemical equilibrium for planetary atmospheres.**

Photochem is an open-source code for simulating how chemistry, radiation, transport, and climate shape planetary atmospheres. Its Python interfaces can evolve atmospheric composition through time, find photochemical steady states, calculate radiative-convective-equilibrium climates, and determine thermochemical-equilibrium compositions.

[Install Photochem](installation.md){ .md-button .md-button--primary } [Start with a tutorial](tutorials/index.md){ .md-button } [Browse the API](reference/index.md){ .md-button }

## Choose a modeling workflow

<div class="grid cards" markdown>

-   **Rocky-planet photochemistry**

    Use `photochem.EvoAtmosphere` to model the atmospheric photochemistry of a rocky planet like Earth or Venus.

-   **Gas-giant photochemistry**

    Use `EvoAtmosphereGasGiant` from `photochem.extensions` for models of gas-rich planets including deep atmospheres initialized from thermochemical equilibrium.

-   **Climate**

    Use `photochem.clima.AdiabatClimate` to calculate one-dimensional radiative-convective temperature structures, including atmospheres with multiple condensible species.

-   **Chemical equilibrium**

    Use `photochem.equilibrate.ChemEquiAnalysis` to calculate equilibrium compositions from elemental abundances at specified temperatures and pressures.

</div>

## Get started

- [Installation](installation.md) covers the supported ways to install Photochem and verify the installation.
- [Tutorials](tutorials/index.md) develop complete atmospheric-modeling workflows in executable notebooks.
- [API Reference](reference/index.md) documents the public Python classes and functions.
