# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Rocky planet photochemistry: Modern Earth
#
# This tutorial constructs a one-dimensional photochemical model of modern Earth's atmosphere, evolves its chemistry to steady state, and analyzes the results. The setup is adapted from the modern-Earth model of [Wogan et al. (2025)](https://doi.org/10.3847/PSJ/ae0e1c).

# %%
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np

from photochem import EvoAtmosphere
from photochem.utils import stars, zahnle_rx_and_thermo_files

# Support execution by MkDocs from the repository root and interactive
# execution from docs/tutorials.
tutorial_directory = Path("docs/tutorials/rocky_planet_photochemistry")
if not tutorial_directory.is_dir():
    tutorial_directory = Path("rocky_planet_photochemistry")


# %% [markdown]
# ## Initializing the photochemical model
#
# Photochem's photochemical model is contained in the `EvoAtmosphere` class. Constructing an `EvoAtmosphere` requires three files: a chemical mechanism, a settings file, and a stellar-flux file.
#
# **Chemical mechanism:** Photochem includes the Zahnle reaction network. Here, we use `zahnle_rx_and_thermo_files` to create a subset of that network that contains H, He, N, O, C, and S chemistry.

# %%
zahnle_rx_and_thermo_files(
    atoms_names=["H", "He", "N", "O", "C", "S"],
    rxns_filename=str(tutorial_directory / "zahnle_earth_HNOCHeS.yaml"),
    thermo_filename=None,
    remove_reaction_particles=True,
)


# %% [markdown]
# **Settings file:** We use the settings file `rocky_planet_photochemistry/settings.yaml`, which is based on the modern-Earth calculation of [Wogan et al. (2025)](https://doi.org/10.3847/PSJ/ae0e1c). It specifies the planet's mass and radius, boundary conditions, and other parameters.
#
# **Stellar flux file:** To generate a stellar-flux file, we use the `solar_spectrum` routine. It starts from the packaged Thuillier et al. (2004) ATLAS-1 solar reference spectrum, extends it to long wavelengths with a solar blackbody, and rebins it to a resolution appropriate for Photochem.

# %%
_ = stars.solar_spectrum(outputfile=str(tutorial_directory / "modern_sun.txt"))


# %% [markdown]
# We can now initialize an instance of `EvoAtmosphere`:

# %%
pc = EvoAtmosphere(
    str(tutorial_directory / "zahnle_earth_HNOCHeS.yaml"),
    str(tutorial_directory / "settings.yaml"),
    str(tutorial_directory / "modern_sun.txt"),
)
pc.var.verbose = 0  # Turn off printing


# %% [markdown]
# At this point the `pc` object has been constructed, which means the chemical kinetics, thermodynamics, stellar spectrum, etc. are set. However, the initial atmospheric state (temperature, eddy diffusion, and composition) has not yet been established:

# %%
print(f"Atmosphere initialized: {pc.atmosphere_initialized}")

# %% [markdown]
# ## Setting the atmospheric state
#
# Now we must set the initial atmospheric state. This includes the temperature, eddy diffusion, and mixing ratio profiles. For the temperature and eddy diffusion, we will use `rocky_planet_photochemistry/modern_earth_profile.txt`, which contains the January CIRA-86 pressure-temperature structure and the eddy-diffusion profile based on Massie and Hunten (1981). The equatorial profile is warmer at the surface than a global-mean Earth profile, but it maintains consistency with the modern-Earth benchmark in [Wogan et al. (2025)](https://doi.org/10.3847/PSJ/ae0e1c) on which this setup is based.

# %%
altitude_km, pressure_bar, temperature, eddy = np.loadtxt(
    tutorial_directory / "modern_earth_profile.txt"
).T
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
# To set this P-T-$K_{zz}$ profile in the code, we will use `initialize_atmosphere_p`. This function maps the pressure-based profile onto the model grid. Setting `persistent=True` retains temperature and eddy diffusion as functions of hydrostatic pressure while the composition evolves. Rainout is enabled, so persistent initialization also needs the tropopause pressure; we interpolate it at the 11 km tropopause specified in `settings.yaml`.
#
# `initialize_atmosphere_p` optionally accepts a `mix` argument, which is the initial volume mixing ratios of the atmosphere. Here, we do not supply `mix`. When this is the case, the code guesses a sensible initial composition from gases with fixed-partial-pressure lower boundary conditions.

# %%
tropopause_altitude_km = 11.0
tropopause_pressure = 10.0 ** np.interp(
    tropopause_altitude_km,
    altitude_km,
    np.log10(pressure),
)
# Use a 10 micron (1.0e-3 cm) radius for water particles.
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

print(f"Atmosphere initialized: {pc.atmosphere_initialized}")
print(f"Surface pressure: {pc.wrk.surface_pressure:.4f} bar")
print(f"Tropopause pressure: {tropopause_pressure / 1.0e6:.3e} bar")


# %%
def plot_composition(atmosphere, species, title, axis=None, legend=True):
    """Plot selected mixing ratios against pressure."""
    if axis is None:
        _, axis = plt.subplots(figsize=(5.5, 5.0))
    for name in species:
        axis.plot(atmosphere[name], atmosphere["pressure"] / 1.0e6, label=name)
    axis.set_xscale("log")
    axis.set_yscale("log")
    if not axis.yaxis_inverted():
        axis.invert_yaxis()
    axis.set_xlim(1.0e-15, 1.2)
    axis.set_xlabel("Mixing ratio")
    axis.set_ylabel("Pressure (bar)")
    axis.set_title(title)
    axis.grid(alpha=0.25)
    if legend:
        axis.legend()
    return axis

initial_atmosphere = pc.mole_fraction_dict()
plot_composition(
    initial_atmosphere,
    ["H2O", "N2", "O2", "CO2"],
    "Initial atmosphere",
)
plt.show()

# %% [markdown]
# ## Organization of `EvoAtmosphere`
#
# Most information in an `EvoAtmosphere` is organized under three attributes. `pc.dat` contains model definitions that generally do not change, such as species and reactions. `pc.var` contains configurable quantities such as temperature, eddy diffusion, and solver tolerances. `pc.wrk` contains the evolving numerical state, including the number densities in `pc.wrk.usol`.
#
# The function `mole_fraction_dict` is convenient because it returns the most important model state as a dictionary.

# %%
print("Model-grid temperature shape:", pc.var.temperature.shape)
print("ODE-state shape:", pc.wrk.usol.shape)
print("Available atmosphere fields:", list(initial_atmosphere)[:8], "...")

# %% [markdown]
# ## Evolve the atmosphere to steady state
#
# The simplest way to integrate to a steady state is `pc.find_steady_state()`. It manages solver restarts and recovery internally, then returns whether the configured convergence criteria are met. Here we use the underlying `initialize_robust_stepper` and `robust_step` methods instead so that we can retain intermediate states and watch the atmosphere develop.

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
fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True, sharey=True)

show_legend = True
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
    plot_composition(
        atmosphere,
        evolving_species,
        title,
        axis=axis,
        legend=show_legend,
    )
    show_legend = False

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
axis.legend(ncol=2, bbox_to_anchor=(0.5, 1.01), loc="lower center", fontsize=8)
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
pc.out2atmosphere_txt(str(tutorial_directory / "modern_earth_steady.txt"), overwrite=True)

pc_restart = EvoAtmosphere(
    str(tutorial_directory / "zahnle_earth_HNOCHeS.yaml"),
    str(tutorial_directory / "settings.yaml"),
    str(tutorial_directory / "modern_sun.txt"),
    str(tutorial_directory / "modern_earth_steady.txt"),
)
pc_restart.var.verbose = 0


# %% [markdown]
# The atmosphere file does not preserve the persistent P-T-$K_{zz}$ configuration. We therefore restore the original pressure-based profiles and tropopause pressure before integrating the restarted model.

# %%
pc_restart.set_press_temp_edd_profile(
    pressure,
    temperature,
    eddy,
    trop_p=tropopause_pressure,
    hydro_pressure=True,
    maintain_toa_pressure=True,
    target_pressure=pressure[-1],
)

restart_start = time.perf_counter()
restart_converged = pc_restart.find_steady_state()
restart_elapsed = time.perf_counter() - restart_start

print(f"Restart converged: {restart_converged}")
print(f"Restart integration time: {restart_elapsed:.2f} seconds")
print(f"Accepted steps in the final solver session: {pc_restart.wrk.nsteps}")


# %% [markdown]
# Starting from the saved steady-state profile should require much less work than constructing the atmosphere from the simple inferred composition. The result nevertheless depends on the separately supplied mechanism, settings, and stellar spectrum, which is why an atmosphere file should be treated as a reusable profile rather than a self-contained model restart.
