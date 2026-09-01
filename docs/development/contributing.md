# Contributing

Contributions to Photochem are welcome!

## Before starting

For a substantial change, consider opening a [GitHub issue](https://github.com/Nicholaswogan/photochem/issues) first. An early discussion can confirm that the proposed approach fits the project and may identify relevant code or existing work.

Set up a development environment and verify that you can compile Photochem by following [Building from Source](building-from-source.md).

## Contribution workflow

Fork Photochem on GitHub, then clone your fork and configure the main Photochem repository as the `upstream` remote:

```sh
git clone https://github.com/YOUR-USERNAME/photochem.git
cd photochem
git remote add upstream https://github.com/Nicholaswogan/photochem.git
```

Bring your fork's `main` branch up to date before beginning a change:

```sh
git switch main
git fetch upstream
git merge --ff-only upstream/main
git push origin main
```

Create a new branch from the updated `main` branch using a short, descriptive name:

```sh
git switch -c improve-climate-validation
```

Make and test a focused change. Commit the relevant files with a clear message, then push the branch to your fork:

```sh
git add path/to/changed-file
git commit -m "Improve climate input validation"
git push -u origin improve-climate-validation
```

Open a pull request from that branch to the Photochem `main` branch. Continue pushing review changes to the same branch until the pull request is complete.

## Pull requests

A pull request should contain one core change. Avoid combining multiple features, unrelated bug fixes, broad reformatting, or cleanup in one contribution. Update tests and documentation when behavior or a public interface changes.

Include:

- a clear title;
- the motivation for the change;
- a concise summary of the implementation;
- links to related issues; and
- the tests or other checks used to verify the result.

## General style

- Follow the style of the surrounding Python, Cython, C, or Fortran code.
- Use descriptive variable names. Comments should explain reasoning or non-obvious behavior rather than restating the code.
- Use NumPy-style docstrings for public Python interfaces. Cython docstrings should reflect the corresponding Fortran documentation, which is the source of truth for compiled behavior.
- Use FORD-style docstrings for public Fortran interfaces.
- Cite the relevant paper and table/equation when implementing a published scientific method.
- Generally follow CGS units, but regardless, always document the units of inputs and outputs.
- Write Markdown prose with one physical line per paragraph.
