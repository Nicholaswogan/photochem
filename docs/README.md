# Building the documentation

Activate a development environment of your choice, then install the documentation dependencies:

```sh
mamba install -c conda-forge \
  mkdocs mkdocs-material mkdocs-jupyter jupytext matplotlib
```

From the repository root, start the local preview server with:

```sh
mkdocs serve
```

Create a production build and treat every warning as an error with:

```sh
mkdocs build --strict
```

The generated site is written to `site/`, which is ignored by Git. Tutorials are stored as Jupytext `py:percent` files under `docs/tutorials/`. MkDocs executes them during the build, stops on notebook errors, and includes the Jupytext source with each rendered tutorial.

Write Markdown prose with one physical line per paragraph; do not hard-wrap paragraphs or list items to a fixed column width. This keeps the source easier to edit directly. Code blocks, tables, and other structures should retain the line breaks required by their syntax.
