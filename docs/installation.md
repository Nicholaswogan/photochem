# Installation

Photochem is distributed through conda-forge for Linux and macOS. Windows users should run Photochem through the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install).

## Install Photochem

If you already have Conda, install Photochem into the active environment with:

```sh
conda install -c conda-forge photochem
```

An isolated environment is optional but recommended when starting a new project:

```sh
conda create -n my_env_name -c conda-forge photochem
conda activate my_env_name
```

Verify that Python can import the package:

```sh
python -c "import photochem; print(photochem.__version__)"
```

The installation automatically includes `photochem_clima_data`, which contains the bundled reaction mechanisms, photochemical cross sections, climate opacity data, and other model data used by Photochem. These data are located automatically, so a separate model-data download is not required.

## Tutorial dependencies

Photochem itself does not require JupyterLab, Jupytext, or Matplotlib. However, to run the documentation tutorials, you must install them in the same environment:

```sh
conda install -c conda-forge jupyterlab jupytext matplotlib
```

## What if I don't have Conda?

There are several ways to obtain Conda. If you do not already have it, we recommend installing it via [Miniforge](https://github.com/conda-forge/miniforge), a minimal Conda distribution configured to use conda-forge.

## Building from source

The conda-forge package is the recommended installation method. Source builds are intended primarily for developers and are covered in [Building from Source](development/building-from-source.md).
