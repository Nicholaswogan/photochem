# Equilibrate API

`ChemEquiAnalysis` computes thermochemical equilibrium for a specified pressure, temperature, and elemental or species composition. Results are exposed as properties on the same object after a successful solve.

## Construction

```python
photochem.equilibrate.ChemEquiAnalysis(
    thermofile,
    atoms=None,
    species=None,
)
```

The Python interface accepts YAML thermodynamic data. Legacy CEA `.inp` files are not supported through Python. The constructor is shown explicitly because Python's runtime inspection does not expose constructor parameters for compiled extension classes.

::: photochem.equilibrate.ChemEquiAnalysis
    options:
      docstring_options:
        warn_unknown_params: false
      members: null
