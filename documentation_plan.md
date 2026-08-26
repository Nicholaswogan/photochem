# Photochem v1.0 documentation plan

## Goals

The v1.0 documentation should make Photochem approachable through runnable,
scientifically meaningful tutorials without creating a second implementation
of the package that must be maintained separately. It should cover Photochem,
its gas-giant extension, Clima, and Equilibrate as parts of one package.

The documentation will:

- use MkDocs and be published as a versioned website;
- use Markdown for ordinary pages;
- use Jupytext `py:percent` files for notebook tutorials;
- emphasize complete scientific workflows over long theoretical derivations;
- generate the Python API reference from source docstrings as much as possible;
- document all public input-file formats;
- explain the essential physics and numerical behavior, while linking to papers
  and other authoritative sources for detailed derivations; and
- execute its tutorials when the documentation is built so that examples do
  not silently drift from the code.

The University of Arizona Photochem workshop, the existing repository
examples, the documentation proposal, and selected PICASO documentation are
source material. They are not specifications: all borrowed material must be
updated for the v1.0 API and rewritten where clarity can be improved. The old
`docs` branch may be consulted for useful prose, but it is too old to merge as
the foundation of the new site.

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
- **Explanation**, with sections for the Photochemical Model, Climate Model,
  and Chemical Equilibrium
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

The tutorial navigation will remain flat for v1.0. Individual tutorials can
link to related tutorials without introducing more navigation levels.

## Content principles

### Tutorials

Tutorials are the center of the documentation. Each tutorial should:

- begin with the scientific task and the result that will be produced;
- introduce only the concepts needed for that task;
- use current public APIs rather than internal implementation details;
- keep inputs with the tutorial or use clearly identified package data;
- state relevant units and array ordering when each quantity is introduced;
- use cgs units unless a documented interface explicitly uses other units;
- state that atmospheric profiles are ordered from the bottom to the top when
  profile arrays are first introduced;
- include plots or diagnostics that help the user assess whether the model ran
  correctly;
- finish with a short summary and links to relevant input and API reference
  pages; and
- run from a clean environment without relying on state produced by another
  tutorial.

Only the Jupytext `.py` source will be committed. Generated `.ipynb` files and
notebook outputs will not be committed. Tutorial execution products should be
written to temporary or ignored build locations.

The principal tutorials should construct useful models from transparent
starting ingredients instead of loading preconverged answers. In particular,
the Rocky Planet Photochemistry tutorial should use the packaged
`zahnle_earth` mechanism (or create a clearly motivated variant), a committed
Modern Earth settings file, and a simple initial pressure-temperature-Kzz and
composition profile. It should initialize the atmosphere through the public
array-based API and integrate the model itself. A final steady-state
`atmosphere.txt` file should not be a required tutorial input.

The same tutorial should demonstrate creating the stellar spectrum with the
`photochem.utils.stars` API rather than treating a finished Photochem stellar
flux file as an unexplained input. This must not make documentation builds
depend on a live network service. The current `solar_spectrum` implementation
downloads its modern-Sun reference spectrum on every call, so the tutorial
must wait until that utility has a reproducible offline source or cache with a
validated checksum. A committed source spectrum is an acceptable fallback,
but a live download during every build is not.

The Atmosphere Profiles tutorial section or input-reference page should still
show how to write and reload `atmosphere.txt`. It should present that format as
a convenient atmospheric-profile initialization and output format, not as a
complete model restart. It does not preserve all mutable model configuration,
including persistent pressure-temperature-Kzz profile configuration.

### Input-file reference

Input pages are reference material rather than extended tutorials. Each page
should document:

- the purpose of the file;
- its overall structure;
- required and optional entries;
- accepted values, defaults, shapes, and units;
- important validation rules and interactions between options;
- a minimal valid example; and
- links to complete files used by tutorials.

The Fortran readers and validation routines are the source of truth for these
pages. Full input files should not be copied into several Markdown pages. A
canonical runnable file should live with a tutorial, while the reference page
explains its fields and links to it.

The Atmosphere Profiles page should distinguish among legacy atmosphere text
files, direct pressure- and altitude-based initialization, saved model state,
and `evolve` trajectory output. It must state prominently that an atmosphere
text file is not a complete model reinitialization or restart format and does
not preserve settings such as persistent pressure-temperature-Kzz profiles.

### Explanation

The Explanation content should be one Markdown page with sections for the
photochemical model, climate model, and chemical-equilibrium solver. The
methods and appendices in
`codex_reference/Wogan_2025_Photochem-main/` should be its primary starting
point, updated wherever the v1.0 implementation has changed since the paper
was drafted.

Each model section should answer:

1. What problem does this component solve?
2. What are its governing equations or optimization conditions?
3. What are the major physical assumptions?
4. What numerical approach is used?
5. Which papers or derivations provide more detail?
6. Which tutorials and API objects demonstrate the implementation?

These pages should be sufficient to interpret model inputs and outputs, but
should not reproduce long derivations already maintained in papers or theses.

### API reference

API pages should be generated from the installed package rather than copying
docstrings into Markdown. The Cython docstrings reflect the corresponding
Fortran documentation, with the Fortran interfaces remaining the source of
truth for compiled behavior.

The initial implementation should use `mkdocstrings` with runtime inspection
of the compiled extension modules. Cython method signatures and class, method,
and property docstrings are already available at runtime. Constructor
signatures currently require special attention. The first attempt should
enable Cython's `embedsignature` compiler directive and test the rendered
result. Small manually curated constructor declarations are preferable to
maintaining a large parallel stub API if embedded signatures are insufficient.

API pages should explicitly select public classes and members. They should not
dump every name found in compiled modules. Nested state objects such as
`PhotochemData`, `PhotochemVars`, and `PhotochemWrk` should be presented in the
way users encounter them through `pc.dat`, `pc.var`, and `pc.wrk`.

### External packages and experimental extensions

The spectrum tutorial will prepare a generic pressure-temperature-composition
profile suitable for a spectrum model. It may mention PICASO as one possible
consumer, but PICASO will not be a dependency of the main documentation build.
The tutorial will not claim that Photochem itself computes spectra.

PICASO's gas-giant climate/photochemistry workflow may be linked as a related
external example rather than duplicated. The experimental `hotrocks.py`
extension is outside the primary v1.0 documentation because of its larger
optional dependency stack and support requirements.

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

The MkDocs configuration should remain at the repository root so local build
commands are conventional. Documentation dependencies should be isolated from
Photochem's runtime dependencies and constrained sufficiently to make builds
reproducible without requiring frequent manual upgrades.

## Implementation sequence

Work should proceed in bounded passes. Each pass should leave the repository in
a buildable, reviewable state and should avoid mixing broad content writing
with unrelated code changes.

### Pass 1: Establish the MkDocs site

- [x] Create the root `mkdocs.yml` and `docs/` directory structure.
- [x] Add the minimal documentation dependency file.
- [x] Configure Material, navigation, search, code highlighting, admonitions,
  and mathematical rendering.
- [x] Configure Jupytext notebook discovery without adding real tutorials yet.
- [x] Add placeholder pages only where needed to validate navigation.
- [x] Establish `docs/README.md` as the concise maintainer guide for building
  and previewing the documentation. Keep it current as dependencies, tutorial
  execution, validation, and deployment evolve.
- [x] Verify that `mkdocs build --strict` succeeds from a development
  environment containing the documented dependencies. Clean-environment
  automation is added in Pass 9.

The goal of this pass is a plain but functional website. Custom branding,
versioning, deployment, and extensive styling should not delay it.

### Pass 2: Write the home page

1. Explain what Photochem is and the planetary problems it addresses.
2. Introduce the four supported surfaces: rocky-planet Photochem, gas-giant
   Photochem, Clima, and Equilibrate.
3. Provide clear routes for a new user: install, run the first tutorial, read
   the input reference, or open the API reference.
4. Add a concise feature list and representative scientific applications.
5. Add preliminary citation and repository links without duplicating the full
   How to Cite Photochem page.
6. Render and inspect the page at desktop and narrow viewport widths.

### Pass 3: Establish generated API documentation

1. Enable embedded Cython signatures and rebuild the installed package.
2. Configure `mkdocstrings` to inspect the compiled extension modules.
3. Prototype pages for `EvoAtmosphere`, `AdiabatClimate`, and
   `ChemEquiAnalysis`.
4. Verify constructor and method signatures, property documentation, NumPy
   docstring sections, cross-links, and nested state objects.
5. Decide whether any small manual constructor declarations are necessary.
6. Add curated pages for Photochem, the gas-giant extension, Clima,
   Equilibrate, and utilities.
7. Add a strict documentation build check that fails when API pages cannot be
   generated.

This pass should first solve the extraction mechanism, then populate the
reference. It should not rewrite API prose already present in source
docstrings.

### Pass 4: Write the installation page

1. Make conda-forge installation the preferred user path.
2. Include environment creation, activation, installation, and a minimal
   import/version verification command.
3. Explain the additional packages needed to run and edit tutorials.
4. Include concise Windows/WSL guidance where necessary.
5. Link advanced users to Building from Source rather than placing the full
   compiler toolchain on the main installation page.
6. Test every command in a clean environment or against the release candidate
   package channel available at the time.

Completion of passes 1 through 4 forms the first documentation milestone: a
working site with a useful landing page, installed-package API reference, and
tested installation instructions.

### Pass 5: Document input files

Write the input reference one format at a time, checking every statement
against the corresponding Fortran reader and validation code:

1. Reaction Mechanisms
2. Photochem Settings
3. Stellar Flux Files
4. Atmosphere Profiles
5. Clima Species and Settings
6. Equilibrate Thermodynamic Data

Use small test fixtures to confirm examples where practical. Link each page to
the tutorials that use the format.

### Pass 6: Build the tutorials

Port and modernize tutorials individually in this order unless dependencies or
scientific review suggest a better sequence:

1. Rocky Planet Photochemistry
2. Rocky Planet Climate
3. Chemical Equilibrium
4. Gas Giant Photochemistry
5. The Habitable Zone
6. Changing Model Data
7. Climate Then Photochemistry: TRAPPIST-1e
8. Preparing Photochem Results for Spectrum Calculations

Likely source material is:

| New tutorial | Primary source material |
| --- | --- |
| Rocky Planet Photochemistry | Workshop `ModernEarth_COMPLETED.ipynb` and `examples/Tutorial.ipynb` |
| Gas Giant Photochemistry | Workshop `GasGiants.ipynb` and `examples/GasGiants.ipynb` |
| Rocky Planet Climate | Workshop `ClimateTutorial_COMPLETED.ipynb` |
| Chemical Equilibrium | New, based on the Equilibrate API and tests |
| The Habitable Zone | Workshop `HabitableZone_COMPLETED.ipynb` |
| Changing Model Data | Workshop `ModelData.ipynb` and the current utility APIs |
| Climate Then Photochemistry: TRAPPIST-1e | Workshop `TRAPPIST1e_COMPLETED.ipynb` |
| Spectrum preparation | New, based on public Photochem output dictionaries |

For every tutorial:

1. Audit the old notebook against the current API.
2. Identify and copy only inputs still scientifically appropriate.
3. Convert the content into a canonical Jupytext `py:percent` file.
4. Rewrite the narrative and code for clarity.
5. Execute it from a clean state with an explicit timeout.
6. Check numerical output for finite and scientifically reasonable results.
7. Render the page and inspect plots, equations, navigation, and downloads.
8. Add links to related input, explanation, and API pages.

### Pass 7: Write the Explanation page

Write one concise page with Photochemical Model, Climate Model, and Chemical
Equilibrium sections. Use `codex_reference/Wogan_2025_Photochem-main/` as the
primary source, supplemented by the proposal, workshop background, relevant
papers, and current source code. Prefer citations over duplicating full
published derivations. Check equations and scientific claims carefully and
connect them to the exact public interfaces used by tutorials.

### Pass 8: Complete development and citation material

1. Write current source-build instructions and verify them in a clean build.
2. Write contribution expectations, test commands, documentation workflow, and
   pull-request guidance.
3. Curate valuable unfinished roadmap items into Future Directions.
4. Write the How to Cite Photochem page, including software, model-data,
   component, and tutorial-specific scientific citations.

### Pass 9: Automate, deploy, and version the site

1. Add a clean documentation-build job that installs Photochem and the docs
   dependencies.
2. Execute all first-party tutorials as part of the documentation build.
3. Publish successful builds to GitHub Pages.
4. Add versioned documentation for tagged releases and a clear stable/latest
   policy.
5. Add link checking and other inexpensive validation once the site structure
   is stable.
6. Measure build time and separate optional external examples if they make the
   primary documentation build unreliable.

### Pass 10: v1.0 documentation and release audit

1. Read the complete site in navigation order and remove duplication.
2. Verify all commands, downloads, cross-references, and external links.
3. Confirm every public component has an installation path, tutorial, input
   reference where applicable, explanation, and API reference.
4. Refresh the repository README to point users to the site and first tutorial.
5. Complete the migration guide, changelog, release notes, and final tutorial
   validation required by roadmap item 9.8.

## Working method for Codex passes

Codex will perform much of this implementation. To keep the work accurate and
reviewable, each pass should follow this pattern:

1. Inspect the relevant current source, tests, and legacy material.
2. State the intended scope and any scientific assumptions before editing.
3. Implement one bounded part of the plan.
4. Build and execute the affected documentation locally.
5. Inspect rendered output when layout, plots, or notebook presentation change.
6. Run source and link checks appropriate to the pass.
7. Report what changed, what was verified, and any decisions requiring review.
8. Leave unrelated user changes and reference material untouched.

Major checkpoints for user review are the first working site, the API-reference
prototype, each completed tutorial, the Explanation page, and the first
deployed version. Scientific interpretation and tutorial scope should be
reviewed at those checkpoints rather than postponed until the final audit.

## Completion criteria

The documentation effort is complete for v1.0 when:

- the MkDocs site builds reproducibly and is deployed with versioning;
- installation instructions have been tested;
- the curated public Python API is generated from current package docstrings;
- every listed input format is documented against its source parser;
- all eight first-party tutorials execute successfully during a documentation
  build;
- the Explanation page accurately describes all three implemented models and
  cites deeper sources;
- source-build, contribution, future-direction, and citation pages are present;
- the repository README points to the published documentation; and
- the release documentation, migration notes, and final validation required by
  the v1.0 roadmap are complete.
