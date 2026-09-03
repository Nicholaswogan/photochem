# Building the documentation

## Prepare the development environment and compiled extensions

Follow the [building-from-source guide](development/building-from-source.md) to clone the repository and create the `photochem-dev` environment. From the repository root, the essential environment setup is:

```sh
mamba env create -f environment-dev.yml
conda activate photochem-dev
export FC="$CONDA_PREFIX/bin/gfortran"
```

Next, follow Option 2 in the source-build guide to configure a direct CMake build with `BUILD_PYTHON_PHOTOCHEM=ON`. Build and install the configured targets with:

```sh
cmake --build build -j
cmake --install build
```

For a direct CMake build, the install step copies `_photochem`, `_clima`, and `_equilibrate` into the source-tree `photochem/` directory. The documentation commands below place the repository root on `PYTHONPATH`, so these in-place extensions are required. Repeat the CMake build and install steps after changing Fortran, C, or Cython source. Documentation-only and pure-Python changes do not require recompilation.

## Install the documentation dependencies

With `photochem-dev` still active, install the documentation dependencies:

```sh
mamba install -c conda-forge \
  mkdocs mkdocs-material mkdocs-jupyter mkdocstrings-python black \
  jupytext matplotlib
python -m pip install --upgrade "mkdocs-jupyter>=0.26.3,<1"
```

The Mamba command installs MkDocs-Jupyter and its dependencies from conda-forge. The pip command then upgrades MkDocs-Jupyter itself because version 0.26.3 or newer is required for persistent tutorial caching and a sufficiently recent version is not currently available from conda-forge.

## Preview and build the site

MkDocs executes every tutorial from `docs/tutorials`, matching the directory used for interactive tutorial work. From the repository root, start the local preview server with:

```sh
cd docs/tutorials
PYTHONPATH=../.. mkdocs serve --config-file ../../mkdocs.yml
```

This is most useful for developing or writing documentation, as you can actively see the site changing. Create a production build with:

```sh
cd docs/tutorials
PYTHONPATH=../.. mkdocs build --strict --config-file ../../mkdocs.yml
```

The generated site is written to the repository's `site/` directory, which is ignored by Git. The API reference is extracted from the in-place source-tree package, so rebuild and install the compiled extensions before checking Cython API changes. Tutorials are stored as Jupytext `py:percent` files under `docs/tutorials/`. MkDocs stops on notebook errors and includes the Jupytext source with each rendered tutorial.

Tutorials can also be converted or run manually from `docs/tutorials/`. For example, from the repository root:

```sh
cd docs/tutorials
jupytext --to notebook rocky_planet_photochemistry.py
```

Write Markdown prose with one physical line per paragraph; do not hard-wrap paragraphs or list items to a fixed column width. This keeps the source easier to edit directly. Code blocks, tables, and other structures should retain the line breaks required by their syntax.

## Automation and deployment

The `documentation` GitHub Actions workflow builds the complete site on pushes to `dev` and `main`, on pull requests targeting either branch, and when manually requested. Only a push to `main` uploads and deploys the GitHub Pages artifact. The root website contains the latest documentation from `main`; historical versioned documentation is deferred until the project has another published documentation version.
