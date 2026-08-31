# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Rocky planet photochemistry
#
# This tutorial constructs a one-dimensional model of modern Earth's atmosphere, evolves its chemistry to steady state, and diagnoses the processes controlling methane. The setup is adapted from the modern-Earth model of [Wogan et al. (2025)](https://doi.org/10.3847/PSJ/ae0e1c).
#
# Photochem's photochemical model is contained in the `EvoAtmosphere` class. Constructing an `EvoAtmosphere` requires three files: a chemical mechanism, a settings file, and a stellar-flux file. We will generate the mechanism and stellar spectrum, while the settings file is kept with this tutorial so it can be inspected and modified.

# %%
from pathlib import Path
import tempfile
import time

import matplotlib.pyplot as plt
import numpy as np

from photochem import EvoAtmosphere
from photochem.utils import stars, zahnle_rx_and_thermo_files


# %% [markdown]
# ## Construct the photochemical model
#
# We first locate the two version-controlled inputs that are specific to this tutorial. The fallback makes the same source work when it is executed either from the repository root or from the `docs/tutorials` directory used by the documentation builder.

# %%
tutorial_inputs = Path("docs/tutorials/rocky_planet_photochemistry")
if not tutorial_inputs.is_dir():
    tutorial_inputs = Path("rocky_planet_photochemistry")

settings_file = tutorial_inputs / "settings.yaml"
profile_file = tutorial_inputs / "modern_earth_profile.txt"


# %% [markdown]
# The settings file specifies the planet, water and rainout behavior, particle properties, and lower and upper boundary conditions. This particular file is based on the modern-Earth calculation of Wogan et al. (2025). The [Input Files](../input-files/) section will describe the format in detail.
#
# Photochem includes a general Zahnle reaction network. Here, `zahnle_rx_and_thermo_files` creates a mechanism containing H, He, N, O, C, and S chemistry while excluding chlorine. The generated mechanism and stellar spectrum are temporary products rather than tutorial inputs that must be maintained by hand.

# %%
temporary_directory = tempfile.TemporaryDirectory()
work_directory = Path(temporary_directory.name)

mechanism_file = work_directory / "zahnle_earth_HNOCHeS.yaml"
zahnle_rx_and_thermo_files(
    atoms_names=["H", "He", "N", "O", "C", "S"],
    rxns_filename=str(mechanism_file),
    thermo_filename=None,
    remove_reaction_particles=True,
)


# %% [markdown]
# Photochemical rates depend on the wavelength-dependent stellar irradiance. `solar_spectrum` starts from the packaged Thuillier et al. (2004) ATLAS-1 solar reference spectrum, appends the infrared spectrum, and rebins it to a resolution appropriate for Photochem. No network connection is required.

# %%
stellar_flux_file = work_directory / "modern_sun.txt"
stars.solar_spectrum(outputfile=str(stellar_flux_file))


# %% [markdown]
# We can now construct `EvoAtmosphere`. We deliberately omit an atmosphere file: the constructor reads the chemistry, settings, and stellar spectrum, but atmospheric pressure, temperature, eddy diffusion, and composition will be initialized separately.

# %%
pc = EvoAtmosphere(
    str(mechanism_file),
    str(settings_file),
    str(stellar_flux_file),
)
pc.var.verbose = 0

print(f"Number of species: {pc.dat.nq}")
print(f"Atmosphere initialized: {pc.atmosphere_initialized}")


# %% [markdown]
# ## Initialize the atmospheric state
#
# The atmospheric profile combines the equatorial January CIRA-86 pressure-temperature structure with an eddy-diffusion profile based on Massie and Hunten (1981). The equatorial profile is warmer at the surface than a global-mean Earth profile, but it maintains consistency with the modern-Earth benchmark on which this setup is based.

# %%
altitude_km, pressure_bar, temperature, eddy = np.loadtxt(profile_file).T
pressure = pressure_bar * 1.0e6  # bar to dyn/cm^2

fig, temperature_axis = plt.subplots(figsize=(5.5, 5.0))
eddy_axis = temperature_axis.twiny()

temperature_axis.plot(temperature, pressure_bar, color="C3", linewidth=2)
eddy_axis.plot(eddy, pressure_bar, color="C0", linewidth=2)

temperature_axis.set_xlabel("Temperature (K)", color="C3")
eddy_axis.set_xlabel(r"Eddy diffusion (cm$^2$ s$^{-1}$)", color="C0")
temperature_axis.set_ylabel("Pressure (bar)")
temperature_axis.set_yscale("log")
temperature_axis.invert_yaxis()
eddy_axis.set_xscale("log")
temperature_axis.tick_params(axis="x", colors="C3")
eddy_axis.tick_params(axis="x", colors="C0")
temperature_axis.grid(alpha=0.25)
plt.show()


# %% [markdown]
# `initialize_atmosphere_p` maps this pressure-based profile onto the model grid. Setting `persistent=True` retains temperature and eddy diffusion as functions of hydrostatic pressure while the composition evolves. Rainout is enabled, so persistent initialization also needs the tropopause pressure; we interpolate it at the 11 km tropopause specified in `settings.yaml`.
#
# We do not supply `mix`. Photochem therefore infers the initial composition from gases with fixed-partial-pressure lower boundary conditions. N2, O2, and CO2 fill the dry atmosphere in their prescribed proportions. H2O participates in the mixture until it reaches the configured 40% relative-humidity limit, after which it is capped and cold trapped upward. Other gases and all particles begin at negligible abundance. We provide a 10 micron radius for water particles, following the benchmark setup.

# %%
tropopause_altitude_km = 11.0
tropopause_pressure = 10.0 ** np.interp(
    tropopause_altitude_km,
    altitude_km,
    np.log10(pressure),
)
particle_radius = {"H2Oaer": np.full(pressure.size, 1.0e-3)}

pc.initialize_atmosphere_p(
    pressure,
    temperature,
    eddy,
    mix=None,
    particle_radius=particle_radius,
    persistent=True,
    tropopause_pressure=tropopause_pressure,
    target_pressure=pressure[-1],
)
pc.var.atol = 1.0e-23

print(f"Atmosphere initialized: {pc.atmosphere_initialized}")
print(f"Surface pressure: {pc.wrk.surface_pressure:.4f} bar")
print(f"Tropopause pressure: {tropopause_pressure / 1.0e6:.3e} bar")


# %% [markdown]
# Most information in an `EvoAtmosphere` is organized under three attributes. `pc.dat` contains model definitions that generally do not change, such as species and reactions. `pc.var` contains configurable quantities such as temperature, eddy diffusion, and solver tolerances. `pc.wrk` contains the evolving numerical state and prepared diagnostics, including pressure and `usol`.
#
# `pc.wrk.usol` contains number densities rather than mixing ratios. For routine inspection, `mole_fraction_dict` is more convenient: it returns the atmospheric structure and each species' mixing-ratio profile.

# %%
initial_atmosphere = pc.mole_fraction_dict()

print("Model-grid temperature shape:", pc.var.temperature.shape)
print("ODE-state shape:", pc.wrk.usol.shape)
print("Available atmosphere fields:", list(initial_atmosphere)[:8], "...")


# %%
def plot_composition(atmosphere, species, title, axis=None):
    """Plot selected mixing ratios against pressure."""
    if axis is None:
        _, axis = plt.subplots(figsize=(5.5, 5.0))
    for name in species:
        axis.plot(atmosphere[name], atmosphere["pressure"] / 1.0e6, label=name)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.invert_yaxis()
    axis.set_xlim(1.0e-15, 1.2)
    axis.set_xlabel("Mixing ratio")
    axis.set_ylabel("Pressure (bar)")
    axis.set_title(title)
    axis.grid(alpha=0.25)
    axis.legend()
    return axis


plot_composition(
    initial_atmosphere,
    ["H2O", "N2", "O2", "CO2"],
    "Initial atmosphere",
)
plt.show()


# %% [markdown]
# ## Evolve the atmosphere to steady state
#
# The simplest way to integrate is `pc.find_steady_state()`. It manages solver restarts and recovery internally, then returns whether the configured convergence criteria were met. Here we use the underlying `initialize_robust_stepper` and `robust_step` methods instead so that we can retain intermediate states and watch the atmosphere develop.

# %%
evolving_species = ["H2O", "O3", "CH4", "CO", "N2O"]


def save_snapshot(model):
    """Copy the small subset of model state needed for the evolution plot."""
    atmosphere = model.mole_fraction_dict()
    snapshot = {
        "time": model.wrk.tn,
        "pressure": atmosphere["pressure"].copy(),
    }
    for name in evolving_species:
        snapshot[name] = atmosphere[name].copy()
    return snapshot


snapshots = [save_snapshot(pc)]
pc.initialize_robust_stepper(pc.wrk.usol)

converged = False
give_up = False
robust_steps = 0
start_time = time.perf_counter()
while not converged and not give_up:
    give_up, converged = pc.robust_step()
    robust_steps += 1
    snapshots.append(save_snapshot(pc))
    if robust_steps % 100 == 0 or converged or give_up:
        print(
            f"Robust step {robust_steps:4d}: "
            f"t = {pc.wrk.tn:.3e} s, converged = {converged}"
        )

elapsed = time.perf_counter() - start_time
final_time = pc.wrk.tn
pc.destroy_stepper()

if not converged:
    raise RuntimeError("The Modern Earth calculation did not reach steady state.")

print(f"Reached steady state at t = {final_time:.3e} s in {elapsed:.1f} seconds.")


# %% [markdown]
# The panels below retain four stages of the integration. The robust stepper advances through increasingly long chemical timescales, so evenly spaced saved states provide a useful view of the atmosphere assembling from its simple initial composition.

# %%
snapshot_indices = np.unique(
    np.linspace(0, len(snapshots) - 1, 4, dtype=int)
)
fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex=True, sharey=True)

for axis, snapshot_index in zip(axes.flat, snapshot_indices):
    snapshot = snapshots[snapshot_index]
    atmosphere = {
        "pressure": snapshot["pressure"],
        **{name: snapshot[name] for name in evolving_species},
    }
    if snapshot["time"] == 0.0:
        title = "Initial state"
    else:
        title = f"t = {snapshot['time']:.1e} s"
    plot_composition(atmosphere, evolving_species, title, axis=axis)

for axis in axes.flat:
    axis.legend(fontsize=8)
plt.tight_layout()
plt.show()


# %% [markdown]
# The final atmosphere contains the chemically produced trace gases that were absent from the inferred initial state. The dominant background gases remain close to their prescribed lower-boundary abundances, while water follows the imposed relative humidity and cold trap.

# %%
steady_atmosphere = pc.mole_fraction_dict()
final_species = ["H2O", "N2", "O2", "CO2", "O3", "CH4", "CO", "N2O"]

plot_composition(steady_atmosphere, final_species, "Steady-state atmosphere")
plt.show()


# %% [markdown]
# ## Diagnose methane production and loss
#
# `production_and_loss` decomposes the tendency of one species into resolved production and loss contributions. These can include chemical reactions, vertical transport, boundary exchange, rainout, condensation and evaporation, distributed sources, custom rates, and hydrogen escape where applicable. Contributions are ordered by their vertically integrated importance.

# %%
methane_budget = pc.production_and_loss("CH4", pc.wrk.usol)

number_of_processes = 3
pressure_plot = steady_atmosphere["pressure"] / 1.0e6
fig, axis = plt.subplots(figsize=(7.0, 5.0))

for index in range(min(number_of_processes, methane_budget.production.shape[1])):
    axis.plot(
        methane_budget.production[:, index],
        pressure_plot,
        label=f"Production: {methane_budget.production_rx[index]}",
    )
for index in range(min(number_of_processes, methane_budget.loss.shape[1])):
    axis.plot(
        methane_budget.loss[:, index],
        pressure_plot,
        linestyle="--",
        label=f"Loss: {methane_budget.loss_rx[index]}",
    )

axis.set_xscale("log")
axis.set_yscale("log")
axis.invert_yaxis()
axis.set_xlabel(r"Rate (molecules cm$^{-3}$ s$^{-1}$)")
axis.set_ylabel("Pressure (bar)")
axis.grid(alpha=0.25)
axis.legend(fontsize=8)
plt.show()


# %%
print("Leading column-integrated CH4 production terms:")
for label, rate in zip(
    methane_budget.production_rx[:number_of_processes],
    methane_budget.integrated_production[:number_of_processes],
):
    print(f"  {label}: {rate:.3e} molecules cm^-2 s^-1")

print("Leading column-integrated CH4 loss terms:")
for label, rate in zip(
    methane_budget.loss_rx[:number_of_processes],
    methane_budget.integrated_loss[:number_of_processes],
):
    print(f"  {label}: {rate:.3e} molecules cm^-2 s^-1")


# %% [markdown]
# ## Save and reload the atmosphere
#
# `out2atmosphere_txt` writes composition, temperature, eddy diffusion, and atmospheric structure to a human-readable profile. This is a convenient initialization file, but it is not a complete checkpoint: the reaction mechanism, stellar spectrum, boundary conditions, and other settings still come from the three files supplied to the constructor, and persistent pressure-profile configuration is not stored in the atmosphere file.

# %%
restart_file = work_directory / "modern_earth_steady.txt"
pc.out2atmosphere_txt(str(restart_file), overwrite=True)

pc_restart = EvoAtmosphere(
    str(mechanism_file),
    str(settings_file),
    str(stellar_flux_file),
    str(restart_file),
)
pc_restart.var.verbose = 0
pc_restart.var.atol = 1.0e-23

restart_start = time.perf_counter()
restart_converged = pc_restart.find_steady_state()
restart_elapsed = time.perf_counter() - restart_start

print(f"Restart converged: {restart_converged}")
print(f"Restart integration time: {restart_elapsed:.2f} seconds")
print(f"Accepted steps in the final solver session: {pc_restart.wrk.nsteps}")


# %% [markdown]
# Starting from the saved steady-state profile should require much less work than constructing the atmosphere from the simple inferred composition. The result nevertheless depends on the separately supplied mechanism, settings, and stellar spectrum, which is why an atmosphere file should be treated as a reusable profile rather than a self-contained model restart.

# %%
temporary_directory.cleanup()
