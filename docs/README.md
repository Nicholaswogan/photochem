# Building the documentation

Activate a development environment of your choice, then install the documentation dependencies:

```sh
mamba install -c conda-forge \
  mkdocs mkdocs-material mkdocs-jupyter mkdocstrings-python black \
  jupytext matplotlib
```

From the repository root, start the local preview server with:

```sh
mkdocs serve
```

Create a production build and treat every warning as an error with:

```sh
mkdocs build --strict
```

The generated site is written to `site/`, which is ignored by Git. The API reference is extracted from the locally installed Photochem package, so rebuild and install the compiled extensions before checking API changes. The compiled `EvoAtmosphere`, `AdiabatClimate`, and `ChemEquiAnalysis` constructors are declared manually on their API pages because their Cython constructor slots do not expose signatures to runtime inspection; their methods and properties remain generated from the package. Tutorials are stored as Jupytext `py:percent` files under `docs/tutorials/`. MkDocs executes them from the repository root during the build, stops on notebook errors, and includes the Jupytext source with each rendered tutorial. Tutorial paths should support both the repository root used by MkDocs and `docs/tutorials` used for manual execution.

Tutorials can also be converted or run manually from `docs/tutorials/`. For example:

```sh
cd docs/tutorials
jupytext --to notebook rocky_planet_photochemistry.py
```

Write Markdown prose with one physical line per paragraph; do not hard-wrap paragraphs or list items to a fixed column width. This keeps the source easier to edit directly. Code blocks, tables, and other structures should retain the line breaks required by their syntax.

## Automation and deployment

The `documentation` GitHub Actions workflow builds the complete site on pushes to `dev` and `main`, on pull requests targeting either branch, and when manually requested. Only a push to `main` uploads and deploys the GitHub Pages artifact. The root website contains the latest documentation from `main`; historical versioned documentation is deferred until the project has another published documentation version.

After generated outputs were removed, the initial v0.9 production documentation builds took about 3.5 to 4 minutes on the development machine. The workflow records the MkDocs build duration in its job summary. The gas-giant tutorial intentionally downloads the TOI-193 spectrum from the MUSCLES archive, so the automated build also tests that supported network workflow. This external dependency is acceptable for v0.9, but failures should be monitored and the decision revisited if the archive proves unreliable in CI.
