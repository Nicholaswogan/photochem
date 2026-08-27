# Clima API

`AdiabatClimate` is the primary interface for one-dimensional radiative-convective climate calculations. It constructs atmospheric profiles, computes top-of-atmosphere fluxes, solves for surface temperature, and can perform a radiative-convective-equilibrium calculation.

## Construction

```python
photochem.clima.AdiabatClimate(
    species_file,
    settings_file,
    flux_file,
    data_dir=None,
    double_radiative_grid=True,
)
```

The constructor is shown explicitly because Python's runtime inspection does not expose constructor parameters for compiled extension classes. The generated class and member documentation below comes from the installed extension.

::: photochem.clima.AdiabatClimate
    options:
      docstring_options:
        warn_unknown_params: false
      members: null

## Radiative-transfer state

The radiative-transfer object is available as `clima.rad`. It is created and owned by `AdiabatClimate`; users should not construct it directly.

::: photochem._clima.Radtran
    options:
      members: null

Workspace arrays and longwave and shortwave channel metadata are available through `clima.rad.wrk`, `clima.rad.lw`, and `clima.rad.sw`.

::: photochem._clima.ClimaRadtranWrk
    options:
      members: null

::: photochem._clima.RTChannel
    options:
      members: null

## Rebinning utilities

::: photochem.clima.rebin

::: photochem.clima.rebin_with_errors
