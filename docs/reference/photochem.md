# Photochem API

`EvoAtmosphere` is the primary interface for one-dimensional photochemical models. Constructing the object loads the reaction mechanism, settings, and stellar flux; the initial atmosphere can be supplied at construction or initialized later from arrays.

## Construction

```python
photochem.EvoAtmosphere(
    mechanism_file,
    settings_file,
    flux_file,
    atmosphere_txt=None,
    data_dir=None,
)
```

The constructor is shown explicitly because Python's runtime inspection does not expose constructor parameters for compiled extension classes. The generated class and member documentation below comes from the installed extension.

::: photochem.EvoAtmosphere
    options:
      docstring_options:
        warn_unknown_params: false
      members: null

## Model state

An `EvoAtmosphere` instance exposes three non-owning state views. Use `pc.dat` for static mechanism data, `pc.var` for prepared configuration and atmospheric profiles, and `pc.wrk` for mutable integration state. These objects are created and owned by `EvoAtmosphere`; users should not construct them directly.

### `pc.dat`

::: photochem._photochem.PhotochemData
    options:
      members: null

Particle saturation entries in `pc.dat.particle_sat` are borrowed state views.

::: photochem._photochem.SaturationData
    options:
      members: null

### `pc.var`

::: photochem._photochem.PhotochemVars
    options:
      members: null

Condensation parameters in `pc.var.cond_params` are borrowed state views.

::: photochem._photochem.CondensationParameters
    options:
      members: null

The top-of-atmosphere pressure controller is available as `pc.var.toa_pressure_maintenance`.

::: photochem._photochem.TOAPressureMaintenance
    options:
      members: null

### `pc.wrk`

::: photochem._photochem.PhotochemWrk
    options:
      members: null

## Diagnostic results

`pc.production_and_loss(...)` returns a `ProductionLoss` object containing the vertically resolved production and loss terms for a selected species.

::: photochem._photochem.ProductionLoss
    options:
      members: null
