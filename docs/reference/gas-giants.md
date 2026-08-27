# Gas-giant extension API

`EvoAtmosphereGasGiant` extends `EvoAtmosphere` with initialization and domain-management tools for gas-rich planets. It initializes composition from chemical equilibrium along a supplied climate profile, manages the deeper photochemical domain, and provides helpers for saving and restoring model state.

Pressure-based and altitude-based atmospheric initialization are general Photochem capabilities and are not, by themselves, distinctions between rocky-planet and gas-giant workflows.

::: photochem.extensions.EvoAtmosphereGasGiant
    options:
      inherited_members: false
      members:
        - initialize_to_climate_equilibrium_PT
        - reinitialize_to_new_climate_PT
        - return_atmosphere_climate_grid
        - return_atmosphere
        - robust_step
        - model_state_to_dict
        - initialize_from_dict

The inherited photochemical methods and `dat`, `var`, and `wrk` state views are documented in the [Photochem API](photochem.md). Gas-giant-specific configuration and saved climate profiles are stored in `pc.gdat`.

::: photochem.extensions.gasgiants.GasGiantData
    options:
      members: false
