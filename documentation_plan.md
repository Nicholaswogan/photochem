# Photochem v1.0 documentation plan

## Goals

The v1.0 documentation should make Photochem approachable through runnable, scientifically meaningful tutorials without creating a second implementation of the package that must be maintained separately. It should cover Photochem, its gas-giant extension, Clima, and Equilibrate as parts of one package.

The documentation will:

- use MkDocs and be published as a versioned website;
- use Markdown for ordinary pages;
- use Jupytext `py:percent` files for notebook tutorials;
- emphasize complete scientific workflows over long theoretical derivations;
- generate the Python API reference from source docstrings as much as possible;
- document all public input-file formats;
- explain the essential physics and numerical behavior, while linking to papers and other authoritative sources for detailed derivations; and
- execute its tutorials when the documentation is built so that examples do not silently drift from the code.

The University of Arizona Photochem workshop, former repository examples preserved in Git history, the documentation proposal, and selected PICASO documentation are source material. They are not specifications: all borrowed material must be updated for the v1.0 API and rewritten where clarity can be improved. The old `docs` branch may be consulted for useful prose, but it is too old to merge as the foundation of the new site.

## Site structure

- **Home**
- **Installation**
- **Tutorials**
  - Rocky Planet Photochemistry
  - Gas Giant Photochemistry
  - Rocky Planet Climate
  - Chemical Equilibrium
  - The Habitable Zone
  - Changing Model Data
  - Climate Then Photochemistry: TRAPPIST-1e
  - Preparing Photochem Results for Spectrum Calculations
- **Input Files**
  - Reaction Mechanisms
  - Photochem Settings
  - Stellar Flux Files
  - Atmosphere Profiles
  - Clima Species and Settings
  - Equilibrate Thermodynamic Data
- **Explanation**, with sections for the Photochemical Model, Climate Model, and Chemical Equilibrium
- **API Reference**
  - Photochem
  - Gas-Giant Extension
  - Clima
  - Equilibrate
  - Utilities
- **Development**
  - Building from Source
  - Contributing
  - Future Directions
- **How to Cite Photochem**

The tutorial navigation will remain flat for v1.0. Individual tutorials can link to related tutorials without introducing more navigation levels.

## Content principles

### Tutorials

Tutorials are the center of the documentation. Each tutorial should:

- begin with the scientific task and the result that will be produced;
- introduce only the concepts needed for that task;
- use current public APIs rather than internal implementation details;
- keep inputs with the tutorial or use clearly identified package data;
- state relevant units and array ordering when each quantity is introduced;
- use cgs units unless a documented interface explicitly uses other units;
- state that atmospheric profiles are ordered from the bottom to the top when profile arrays are first introduced;
- include plots or diagnostics that help the user assess whether the model ran correctly;
- finish with a short summary and links to relevant input and API reference pages; and
- run from a clean environment without relying on state produced by another tutorial.

Only the Jupytext `.py` source will be committed. Generated `.ipynb` files and notebook outputs will not be committed. Tutorial execution products should be written to temporary or ignored build locations.

The principal tutorials should construct useful models from transparent starting ingredients instead of loading preconverged answers. In particular, the Rocky Planet Photochemistry tutorial should use the packaged `zahnle_earth` mechanism (or create a clearly motivated variant), a committed Modern Earth settings file, and a simple initial pressure-temperature-Kzz and composition profile. It should initialize the atmosphere through the public array-based API and integrate the model itself. A final steady-state `atmosphere.txt` file should not be a required tutorial input.

Pressure-based and altitude-based initialization are general coordinate choices available to planetary workflows. The documentation must not present pressure grids as specific to gas giants or altitude grids as specific to rocky planets; the gas-giant extension is distinguished by its equilibrium initialization and deep-atmosphere workflow.

The same tutorial should demonstrate creating the stellar spectrum with the `photochem.utils.stars` API rather than treating a finished Photochem stellar flux file as an unexplained input. The modern-Sun reference spectrum is packaged in `photochem_clima_data`, so `solar_spectrum` can construct this input reproducibly without network access.

The Atmosphere Profiles tutorial section or input-reference page should still show how to write and reload `atmosphere.txt`. It should present that format as a convenient atmospheric-profile initialization and output format, not as a complete model restart. It does not preserve all mutable model configuration, including persistent pressure-temperature-Kzz profile configuration.

### Input-file reference

Input pages are reference material rather than extended tutorials. Each page should document:

- the purpose of the file;
- its overall structure;
- required and optional entries;
- accepted values, defaults, shapes, and units;
- important validation rules and interactions between options;
- a minimal valid example; and
- links to complete files used by tutorials.

The Fortran readers and validation routines are the source of truth for these pages. Full input files should not be copied into several Markdown pages. A canonical runnable file should live with a tutorial, while the reference page explains its fields and links to it.

The Atmosphere Profiles page should distinguish among legacy atmosphere text files, direct pressure- and altitude-based initialization, saved model state, and `evolve` trajectory output. It must state prominently that an atmosphere text file is not a complete model reinitialization or restart format and does not preserve settings such as persistent pressure-temperature-Kzz profiles.

### Explanation

The Explanation content should be one Markdown page with sections for the photochemical model, climate model, and chemical-equilibrium solver. The methods and appendices in `codex_reference/Wogan_2025_Photochem-main/` should be its primary starting point, updated wherever the v1.0 implementation has changed since the paper was drafted.

Each model section should answer:

1. What problem does this component solve?
2. What are its governing equations or optimization conditions?
3. What are the major physical assumptions?
4. What numerical approach is used?
5. Which papers or derivations provide more detail?
6. Which tutorials and API objects demonstrate the implementation?

These pages should be sufficient to interpret model inputs and outputs, but should not reproduce long derivations already maintained in papers or theses.

### API reference

API pages should be generated from the installed package rather than copying docstrings into Markdown. The Cython docstrings reflect the corresponding Fortran documentation, with the Fortran interfaces remaining the source of truth for compiled behavior.

The initial implementation should use `mkdocstrings` with runtime inspection of the compiled extension modules. Cython method signatures and class, method, and property docstrings are already available at runtime. Constructor signatures currently require special attention. The first attempt should enable Cython's `embedsignature` compiler directive and test the rendered result. Small manually curated constructor declarations are preferable to maintaining a large parallel stub API if embedded signatures are insufficient.

API pages should explicitly select public classes and members. They should not dump every name found in compiled modules. Nested state objects such as `PhotochemData`, `PhotochemVars`, and `PhotochemWrk` should be presented in the way users encounter them through `pc.dat`, `pc.var`, and `pc.wrk`.

### External packages and experimental extensions

The spectrum tutorial will prepare a generic pressure-temperature-composition profile suitable for a spectrum model. It may mention PICASO as one possible consumer, but PICASO will not be a dependency of the main documentation build. The tutorial will not claim that Photochem itself computes spectra.

PICASO's gas-giant climate/photochemistry workflow may be linked as a related external example rather than duplicated. The experimental `hotrocks.py` extension is outside the primary v1.0 documentation because of its larger optional dependency stack and support requirements.

## Technical design

The initial stack is expected to contain:

- MkDocs;
- Material for MkDocs, with minimal or no custom theme overrides;
- `mkdocs-jupyter` for Jupytext tutorial rendering and execution;
- Jupytext using the `py:percent` format;
- `mkdocstrings` with its Python handler for API generation; and
- MathJax or KaTeX support for equations.

The preferred repository layout is:

```text
mkdocs.yml
docs/
    README.md
    index.md
    installation.md
    tutorials/
    input-files/
    explanation/
    reference/
    development/
    assets/
    requirements.txt
```

The MkDocs configuration should remain at the repository root. Documentation builds should execute from `docs/tutorials` so tutorial input paths match interactive use, while the repository root is placed on `PYTHONPATH` for the in-place Photochem package. Documentation dependencies should be isolated from Photochem's runtime dependencies and constrained sufficiently to make builds reproducible without requiring frequent manual upgrades.

First-party Markdown prose should use one physical line per paragraph, without hard-wrapping at a fixed column width. List items should likewise remain on one physical line where Markdown syntax permits; code blocks, tables, and other structures retain the line breaks required by their syntax.

## Implementation sequence

Work should proceed in bounded passes. Each pass should leave the repository in a buildable, reviewable state and should avoid mixing broad content writing with unrelated code changes.

### Pass 1: Establish the MkDocs site

- [x] Create the root `mkdocs.yml` and `docs/` directory structure.
- [x] Add the minimal documentation dependency file.
- [x] Configure Material, navigation, search, code highlighting, admonitions, and mathematical rendering.
- [x] Configure Jupytext notebook discovery without adding real tutorials yet.
- [x] Add placeholder pages only where needed to validate navigation.
- [x] Establish `docs/README.md` as the concise maintainer guide for building and previewing the documentation. Keep it current as dependencies, tutorial execution, validation, and deployment evolve.
- [x] Verify that `mkdocs build --strict` succeeds from a development environment containing the documented dependencies. Clean-environment automation is added in Pass 9.

The goal of this pass is a plain but functional website. Custom branding, versioning, deployment, and extensive styling should not delay it.

### Pass 2: Write the home page

- [x] Explain what Photochem is and the planetary problems it addresses.
- [x] Introduce the four supported surfaces: rocky-planet Photochem, gas-giant Photochem, Clima, and Equilibrate.
- [x] Provide clear routes for a new user: install, run the first tutorial, read the input reference, or open the API reference.
- [x] Add a concise feature list and representative scientific applications.
- [x] Add preliminary citation and repository links without duplicating the full How to Cite Photochem page.
- [x] Render and inspect the page at desktop and narrow viewport widths.

### Pass 3: Establish generated API documentation

- [x] Enable embedded Cython signatures and rebuild the installed package.
- [x] Configure `mkdocstrings` to inspect the compiled extension modules.
- [x] Prototype pages for `EvoAtmosphere`, `AdiabatClimate`, and `ChemEquiAnalysis`.
- [x] Verify constructor and method signatures, property documentation, NumPy docstring sections, cross-links, and nested state objects.
- [x] Use small manual constructor declarations for the three compiled extension classes whose constructor slots cannot be inspected at runtime.
- [x] Add curated pages for Photochem, the gas-giant extension, Clima, Equilibrate, and utilities.
- [x] Verify that the strict documentation build fails if an API page cannot be generated.

This pass should first solve the extraction mechanism, then populate the reference. It should not rewrite API prose already present in source docstrings.

### Pass 4: Write the installation page

- [x] Make conda-forge installation the preferred user path.
- [x] Include environment creation, activation, installation, and a minimal import/version verification command.
- [x] Explain that `photochem_clima_data` is installed automatically and provides the bundled model data.
- [x] Explain the additional packages needed to run and edit tutorials.
- [x] Include concise Windows/WSL guidance and state that native Windows binaries are not supported.
- [x] Direct advanced users to the Building from Source page rather than placing the full compiler toolchain on the main installation page.
- [x] Test the installation and verification commands in a clean environment using the current conda-forge package.

Completion of passes 1 through 4 forms the first documentation milestone: a working site with a useful landing page, installed-package API reference, and tested installation instructions.

### Pass 5: Document input files (deferred until after the initial release)

Input-file documentation is outside the scope of the initial documentation release. Complete this pass after the initial site is published.

Write the input reference one format at a time, checking every statement against the corresponding Fortran reader and validation code:

- [ ] Reaction Mechanisms *(deferred)*
- [ ] Photochem Settings *(deferred)*
- [ ] Stellar Flux Files *(deferred)*
- [ ] Atmosphere Profiles *(deferred)*
- [ ] Clima Species and Settings *(deferred)*
- [ ] Equilibrate Thermodynamic Data *(deferred)*

Use small test fixtures to confirm examples where practical. Link each page to the tutorials that use the format.

### Pass 6: Build the tutorials

For the initial documentation release, port and modernize these four core tutorials:

- [x] Rocky Planet Photochemistry
- [x] Rocky Planet Climate
- [x] Gas Giant Photochemistry
- [ ] The Habitable Zone *(deferred)*

The remaining tutorials are deferred until after the initial release:

- [ ] Chemical Equilibrium *(deferred)*
- [ ] Changing Model Data *(deferred)*
- [ ] Climate Then Photochemistry: TRAPPIST-1e *(deferred)*
- [ ] Preparing Photochem Results for Spectrum Calculations *(deferred)*

Likely source material is:

| New tutorial | Primary source material |
| --- | --- |
| Rocky Planet Photochemistry | Workshop `ModernEarth_COMPLETED.ipynb` and the former v0.8.4 `examples/Tutorial.ipynb` |
| Gas Giant Photochemistry | Workshop `GasGiants.ipynb` and the former v0.8.4 `examples/GasGiants.ipynb` |
| Rocky Planet Climate | Workshop `ClimateTutorial_COMPLETED.ipynb` |
| Chemical Equilibrium | New, based on the Equilibrate API and tests |
| The Habitable Zone | Workshop `HabitableZone_COMPLETED.ipynb` |
| Changing Model Data | Workshop `ModelData.ipynb` and the current utility APIs |
| Climate Then Photochemistry: TRAPPIST-1e | Workshop `TRAPPIST1e_COMPLETED.ipynb` |
| Spectrum preparation | New, based on public Photochem output dictionaries |

For every tutorial:

- [ ] Audit the old notebook against the current API.
- [ ] Identify and copy only inputs still scientifically appropriate.
- [ ] Convert the content into a canonical Jupytext `py:percent` file.
- [ ] Rewrite the narrative and code for clarity.
- [ ] Execute it from a clean state with an explicit timeout.
- [ ] Check numerical output for finite and scientifically reasonable results.
- [ ] Render the page and inspect plots, equations, navigation, and downloads.
- [ ] Add links to related input, explanation, and API pages.

### Pass 7: Write the Explanation page (deferred until after the initial release)

The Explanation page is outside the scope of the initial documentation release. After release, write one concise page with Photochemical Model, Climate Model, and Chemical Equilibrium sections. Use `codex_reference/Wogan_2025_Photochem-main/` as the primary source, supplemented by the proposal, workshop background, relevant papers, and current source code. Prefer citations over duplicating full published derivations. Check equations and scientific claims carefully and connect them to the exact public interfaces used by tutorials.

### Pass 8: Complete development and citation material

- [x] Write current source-build instructions and verify them in a clean build.
- [x] Write contribution expectations, test commands, documentation workflow, and pull-request guidance.
- [ ] Curate valuable unfinished roadmap items into Future Directions. *(deferred)*
- [ ] Write the How to Cite Photochem page, including software, model-data, component, and tutorial-specific scientific citations. *(deferred)*

## v0.9 initial documentation and package release

Photochem v0.9 is an intentional early documentation release rather than completion of the full v1.0 documentation scope. It publishes the useful material that is ready now: installation and source-build instructions, generated API documentation, and complete tutorials for rocky-planet photochemistry, rocky-planet climate, and gas-giant photochemistry. The v0.9 checklist below is authoritative for this milestone; the full v1.0 requirements in Passes 9 and 10 remain below.

### v0.9 Pass 1: Clean the published documentation scope

- [x] Define the v0.9 tutorial set as Rocky Planet Photochemistry, Rocky Planet Climate, and Gas Giant Photochemistry; defer The Habitable Zone.
- [x] Delete everything under `examples/`; tracked examples remain recoverable from Git history.
- [x] Remove the Input Files placeholder page, its navigation entry, and links that imply the input reference is available.
- [x] Remove the Future Directions placeholder page and its navigation entry.
- [x] Remove the home-page development-status notice and other language that makes the deliberately scoped v0.9 site look unfinished.
- [x] Remove the Habitable Zone “in preparation” notice from the tutorial overview.
- [x] Remove notebook checkpoints and generated documentation artifacts from the working tree.
- [x] Search the repository for stale references to `examples/`, deleted placeholders, and material promised only for v1.0.

### v0.9 Pass 2: Prepare release-facing material

- [x] Rewrite the repository `README.md` to describe the current package and link to installation, tutorials, API documentation, source builds, contribution guidance, and the published site.
- [x] Remove README references to the workshop as the primary documentation and to the deleted `examples/` directory.
- [x] Audit release-facing installation commands, documentation links, version references, and package claims for consistency.
- [x] Draft GitHub release notes and a final release checklist for v0.9.0.
- [x] Retain a concise citation section in the README that directs users to cite Wogan et al. (2025) until the full How to Cite Photochem page is written.

### v0.9 Pass 3: Automate and deploy the site

- [x] Add a separate documentation workflow that creates a clean environment, builds and installs Photochem, installs the documentation dependencies, and runs `mkdocs build --strict`.
- [x] Execute all three v0.9 tutorials in the automated build, including the live TOI-193 MUSCLES download.
- [x] Run documentation builds without deployment on development branches and pull requests.
- [x] Upload and deploy the completed site to GitHub Pages only from `main`, with appropriate permissions, concurrency control, and a bounded timeout.
- [x] Enable GitHub Pages with GitHub Actions as its deployment source and verify the expected `https://nicholaswogan.github.io/photochem/` URL.
- [x] Document the v0.9 policy that the root website contains the latest documentation from `main`; historical tagged documentation is deferred until another published documentation version exists.
- [x] Record clean-build duration and evaluate whether the external MUSCLES dependency is acceptably reliable in automation.

### v0.9 Pass 4: Perform the release audit

- [x] Read every published page in navigation order and remove unnecessary duplication.
- [x] Inspect all three rendered tutorials and confirm their numerical results remain finite and scientifically reasonable.
- [x] Verify documentation commands, internal links, external links, downloads, and API cross-references.
- [x] Build the source distribution and wheel in a clean environment, install the wheel, and verify `photochem.__version__ == "0.9.0"`.
- [x] Run the existing release-oriented test suite and confirm the GitHub Actions test workflow passes.
- [x] Confirm a clean checkout builds the complete documentation and no code, tests, or documentation depend on the deleted `examples/` directory.
- [x] Review the final diff for generated files, placeholders, stale version references, and unintended changes.

### v0.9 Pass 5: Publish v0.9.0

- [ ] Commit the audited release changes on `dev` and merge them into `main`.
- [ ] Wait for the test and documentation workflows on `main` and inspect any failures.
- [ ] Verify the deployed GitHub Pages site from the `main` build.
- [ ] Tag `v0.9.0` and create the GitHub release from the prepared release notes.
- [ ] Update the conda-forge package and verify installation from the released channel.

The following documentation work is explicitly deferred beyond v0.9: the input-file reference, Explanation page, Chemical Equilibrium tutorial, Habitable Zone tutorial, Changing Model Data tutorial, climate-then-photochemistry tutorial, spectrum-preparation tutorial, Future Directions, full citation page, historical versioned documentation, and automated external-link checking.

The v0.9 milestone is complete when the site contains no placeholder destinations, the repository README points to the deployed documentation, all three tutorials and the generated API reference build in automation, the release artifacts and installation have been verified, the site is live on GitHub Pages, and the `v0.9.0` release has been published with its deferrals stated clearly.

### Pass 9: Automate, deploy, and version the site

- [ ] Add a clean documentation-build job that installs Photochem and the docs dependencies.
- [ ] Execute all first-party tutorials as part of the documentation build.
- [ ] Publish successful builds to GitHub Pages.
- [ ] Add versioned documentation for tagged releases and a clear stable/latest policy.
- [ ] Add link checking and other inexpensive validation once the site structure is stable.
- [ ] Measure build time and separate optional external examples if they make the primary documentation build unreliable.

### Pass 10: v1.0 documentation and release audit

- [ ] Read the complete site in navigation order and remove duplication.
- [ ] Verify all commands, downloads, cross-references, and external links.
- [ ] Confirm every public component has an installation path, tutorial, input reference where applicable, explanation, and API reference.
- [ ] Refresh the repository README to point users to the site and first tutorial.
- [ ] Complete the migration guide, changelog, release notes, and final tutorial validation required by roadmap item 9.8.

## Working method for Codex passes

Codex will perform much of this implementation. To keep the work accurate and reviewable, each pass should follow this pattern:

1. Inspect the relevant current source, tests, and legacy material.
2. State the intended scope and any scientific assumptions before editing.
3. Implement one bounded part of the plan.
4. Build and execute the affected documentation locally.
5. Inspect rendered output when layout, plots, or notebook presentation change.
6. Run source and link checks appropriate to the pass.
7. Report what changed, what was verified, and any decisions requiring review.
8. Leave unrelated user changes and reference material untouched.

Major checkpoints for user review are the first working site, the API-reference prototype, each completed tutorial, the Explanation page, and the first deployed version. Scientific interpretation and tutorial scope should be reviewed at those checkpoints rather than postponed until the final audit.

## Completion criteria

The documentation effort is complete for v1.0 when:

- the MkDocs site builds reproducibly and is deployed with versioning;
- installation instructions have been tested;
- the curated public Python API is generated from current package docstrings;
- every listed input format is documented against its source parser;
- all eight first-party tutorials execute successfully during a documentation build;
- the Explanation page accurately describes all three implemented models and cites deeper sources;
- source-build, contribution, future-direction, and citation pages are present;
- the repository README points to the published documentation; and
- the release documentation, migration notes, and final validation required by the v1.0 roadmap are complete.
