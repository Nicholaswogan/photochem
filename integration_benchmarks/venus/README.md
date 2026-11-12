# Venus integration benchmark

Port of the nominal legacy SolarSystem Venus setup. All case-specific inputs
are copied here; no script imports from or reads the original project. The
usual installed `photochem_clima_data` database is still required, as for other
Photochem runs. This is a solver benchmark, not a new validation of Venus
atmospheric physics.

## Starts and inputs

- Default **saved-atmosphere start**: `atmosphere.txt`, the legacy nominal Venus
  output. It is a text warm start, not a current-version checkpoint. Do not
  assume it is already converged under the repaired/current model.
- `--fresh`: construct the atmosphere from `venus_seiff1985.txt` (VIRA), extend
  its temperature constantly from 100 to 110 km, and interpolate Kzz from
  `venus_rimmer2021.dat`. Surface pressure is 92.1 bar. The original prescribed
  gas fractions are normalized: CO2 0.965, N2 0.035, SO2 100 ppm, H2O 35 ppm,
  OCS 5 ppm, CO 10 ppm, HCl 500 ppb, and H2 4.5 ppb. Other evolved species start
  at an absolute density of `1e-40 cm^-3`, matching the legacy initializer.
- `flux.txt` is the original `SunNow.txt`; `settings.yaml` retains its 1.93 photon
  scale factor, Venus mass/radius, 100 layers, and closed H2 boundaries.
- `mechanism_legacy.yaml` is the unmodified original nominal mechanism;
  `mechanism.yaml` is the thermodynamically repaired input used for all methods.

The VIRA file identifies Seiff et al. (1985),
[doi:10.1016/0273-1177(85)90197-8](https://doi.org/10.1016/0273-1177(85)90197-8).
The Kzz input identifies the Rimmer et al. (2021) dataset,
[doi:10.7910/DVN/AKUTME](https://doi.org/10.7910/DVN/AKUTME).
`legacy_source_hashes.json` records the original source/input file identities.
Those source-project paths are provenance labels, not runtime dependencies.

## Physical and numerical controls

`initialize.py` owns the controls reapplied by the native comparison driver:

- All particle condensation parameters: smooth factor 5, condensation rate
  coefficient 100, evaporation rate coefficient 100.
- H2SO4 particle radius `2e-4 cm` (2 microns); water particle radius `1e-3 cm`
  (10 microns). Other radii remain at their defaults.
- Built-in sulfuric-acid optical properties (`palmer1975`).
- The legacy nominal Rimmer UV absorber: pure absorption, cutoff at 400 nm.
  Crisp prescribed clouds, Haus absorber alternatives, and the optional extra
  Rimmer chlorine reaction set are not included in this nominal case.
- Upwind molecular diffusion as in the old wrapper, with the current analytical
  chemistry Jacobian (`jacobian_method=3`) for all compared backends. The earlier
  45-second experiments retained legacy autodiff; those generated results are
  archived under `results/legacy-autodiff/`. The current result tables and plots
  use analytical derivatives only. Use `--jacobian-method autodiff` to reproduce
  the earlier choice.
- CVODE `atol=1e-17`, `rtol=1e-3`; the legacy time target is `3.15e13 s`.
  The comparison disables history-based early convergence (`conv_longdy=0`).
- Default PTC uses exactly the same step policy, initial timestep, and `1e-20`
  correction-scale abundance floor as the Archean case. That floor is distinct
  from the CVODE absolute tolerance. No Venus-specific PTC tuning is applied.

## Necessary port differences

1. **Thermodynamic continuity.** The current loader rejects the legacy joins.
   `repair_mechanism.py` uses `make_thermo_continuous`, retaining the segment at
   298.15 K and matching Cp/H/S at adjacent joins. It adjusts 181 segments across
   96 species. `thermo_corrections.json` lists every correction. Reactions,
   particles, and non-thermodynamic species metadata are unchanged. This changes
   some Gibbs energies and therefore derived reverse rates; it is not merely a
   file-format edit. This port is not a bitwise reproduction of the old model.
2. **Explicit initialization.** The fresh start uses the supported
   `initialize_atmosphere_z` API instead of the old custom wrapper and dummy
   atmosphere file. Hydrostatic interpolation may differ slightly numerically.
   Obsolete `bottom`/`top` settings keys were removed; the initializer/input
   atmosphere determines the same 0–110 km domain.
3. **UV opacity coordinate.** The removed `var.tauc/w0c/g0c` setters are replaced
   with `set_custom_optical_properties(..., tau_mode=1)` against initial
   hydrostatic pressure. The imposed extinction times layer thickness matches
   the legacy altitude formula exactly at initialization. During evolution the
   prescription follows pressure, whereas the old prescription stayed at fixed
   altitude layers. The migration is documented rather than claiming exact
   equivalence after the pressure distribution changes.
4. **Jacobian selection.** The legacy wrapper explicitly selected autodiff. The
   current benchmark uses the analytical Jacobian at the user's request. This is
   a numerical-method change, not a change to the specified chemical reactions.
5. **Robust integration.** Use current native robust stepping instead of the old
   Python wrapper's resets and negative-abundance restart heuristics. Model-state
   export/import does not preserve these controls; the shared driver reapplies
   the case namelist. The current installed photochemical data are used rather
   than requiring an old private data directory.

## Running

From the repository root in the build's Python/compiler environment:

```sh
PYTHONPATH=. python integration_benchmarks/venus/check_setup.py
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case venus --modes cvode cvode1 ptc ptc-residual ptc-floor12 --jacobian-method analytical --steps 100000 --seconds 45
PYTHONPATH=. python integration_benchmarks/run_ptc.py --case venus --fresh --modes cvode cvode1 ptc ptc-residual ptc-floor12 --jacobian-method analytical --steps 100000 --seconds 300
PYTHONPATH=. python integration_benchmarks/compare_ptc.py --case venus
```

The common driver must have a configured/built Photochem test target; see the
[parent instructions](../README.md). `--reuse-inputs` avoids rebuilding a prepared
start. Longer CPU/step budgets can be selected explicitly. Results, native
executables, and logs are ignored under `results/`.

To reproduce the mechanism repair:

```sh
PYTHONPATH=. python integration_benchmarks/venus/repair_mechanism.py
```

`check_setup.py` checks that reactions and particle data were preserved, both
starts are valid, physical controls/radii are applied, and the installed initial
UV optical depth matches the legacy formula. It performs no long integration
and remains a local port check, outside CI.
