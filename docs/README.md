# Building the documentation

Install the documentation dependencies into the `photochem` conda environment:

```sh
mamba install -n photochem -c conda-forge \
  mkdocs mkdocs-material mkdocs-jupyter jupytext matplotlib
```

Activate the environment and start the local preview server from the repository
root:

```sh
conda activate photochem
mkdocs serve
```

Create a production build and treat every warning as an error with:

```sh
mkdocs build --strict
```

The generated site is written to `site/`, which is ignored by Git. Tutorials
are stored as Jupytext `py:percent` files under `docs/tutorials/`. MkDocs
executes them during the build, stops on notebook errors, and includes the
Jupytext source with each rendered tutorial.
