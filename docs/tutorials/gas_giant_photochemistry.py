# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:percent
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
# # Gas giant photochemistry
#
# This tutorial uses `EvoAtmosphereGasGiant` to simulate the photochemistry of the hot Saturn WASP-39b. We focus on the production of sulfur dioxide motivated by the analysis of [Tsai et al. (2023)](https://doi.org/10.1038/s41586-023-05902-2). The calculation illustrates the gas-giant workflow but does not reproduce that study exactly.

# %%
from pathlib import Path
import time

from astropy import constants
import matplotlib.pyplot as plt
import numpy as np

from photochem.extensions import EvoAtmosphereGasGiant
from photochem.utils import stars, zahnle_rx_and_thermo_files

# Support execution by MkDocs from the repository root and interactive
# execution from docs/tutorials.
tutorial_directory = Path("docs/tutorials/gas_giant_photochemistry")
if not tutorial_directory.is_dir():
    tutorial_directory = Path("gas_giant_photochemistry")


# %% [markdown]
# ## Initializing the photochemical model
#
# `EvoAtmosphereGasGiant` adapts Photochem's photochemical model for hydrogen-rich planets. It joins a photochemical atmosphere to a deeper atmosphere assumed to be in chemical equilibrium and provides tools for initializing the composition from metallicity and C/O ratio.
#
# The first required input is a chemical mechanism. Here, `zahnle_rx_and_thermo_files` creates reaction and thermodynamic files containing H, He, N, O, C, and S chemistry. Reaction particles are removed because they do not work well for planets with hot and deep atmospheres.

# %%
mechanism_file = str(tutorial_directory / "photochem_rxns.yaml")
thermo_file = str(tutorial_directory / "photochem_thermo.yaml")

zahnle_rx_and_thermo_files(
    atoms_names=["H", "He", "N", "O", "C", "S"],
    rxns_filename=mechanism_file,
    thermo_filename=thermo_file,
    remove_reaction_particles=True,
)


# %% [markdown]
# We also need a stellar-flux file that includes ultraviolet wavelengths important for photochemistry. The [MUSCLES Treasury Survey](https://archive.stsci.edu/prepds/muscles/) provides panchromatic spectra for nearby stars. WASP-39 is not in the survey, so we use `closest_muscles_to_Teff` to select the available star with the closest effective temperature. `muscles_spectrum` downloads that spectrum, extends it to long wavelengths with a blackbody, and scales its total irradiance to the equilibrium temperature specified by `Teq`.

# %%
equilibrium_temperature = 1166.0  # K
host_star_temperature = 5400.0  # K
comparison_star = stars.closest_muscles_to_Teff(host_star_temperature)
stellar_flux_file = str(tutorial_directory / "toi193_spectrum.txt")

print(
    f"Closest MUSCLES star: {comparison_star['name']} "
    f"({comparison_star['st_teff']:.0f} K)"
)

wavelength, stellar_flux = stars.muscles_spectrum(
    comparison_star["name"],
    outputfile=stellar_flux_file,
    Teq=equilibrium_temperature,
)


# %% [markdown]
# The downloaded spectrum contains the observed ultraviolet emission that would be absent from a photospheric model alone.

# %%
fig, axis = plt.subplots(figsize=(6.5, 4.0))
axis.plot(wavelength, stellar_flux, color="black", linewidth=1)
axis.set_xscale("log")
axis.set_yscale("log")
axis.set_xlim(1.0, 1000.0)
axis.set_xlabel("Wavelength (nm)")
axis.set_ylabel(r"Flux at WASP-39b (mW m$^{-2}$ nm$^{-1}$)")
axis.set_title(f"{comparison_star['name']} MUSCLES spectrum")
axis.grid(alpha=0.25)
plt.show()


# %% [markdown]
# We initialize `EvoAtmosphereGasGiant` with the mechanism, stellar flux, planetary mass and radius, and the separate thermodynamic file. Unlike `EvoAtmosphere`, this extension does not require a user-supplied settings file because its constructor configures gas-giant defaults. We use the 83° solar zenith angle adopted by Tsai et al. (2023) and set the diurnal factor to one so the supplied beam is not reduced by an additional day-night average.

# %%
planet_mass = 0.28 * constants.M_jup.cgs.value
planet_radius = 1.279 * constants.R_jup.cgs.value

pc = EvoAtmosphereGasGiant(
    mechanism_file,
    stellar_flux_file,
    planet_mass,
    planet_radius,
    solar_zenith_angle=83.0,
    thermo_file=thermo_file,
)
pc.gdat.verbose = False
pc.var.diurnal_fac = 1.0


# %% [markdown]
# `pc.gdat.gas` is a `ChemEquiAnalysis` object used to calculate equilibrium chemistry in the deep atmosphere. It contains a reference solar elemental composition that defines what `metallicity` and `CtoO` mean during initialization. To follow Tsai et al. (2023) more closely, we replace that reference with the elemental abundances used in their calculation.

# %%
tsai_elemental_abundances = {
    "H": 1.0,
    "He": 0.0838,
    "C": 2.95e-4,
    "N": 7.08e-5,
    "O": 5.37e-4,
    "S": 1.41e-5,
}
abundance_total = sum(tsai_elemental_abundances.values())
pc.gdat.gas.molfracs_atoms_sun = np.array(
    [
        tsai_elemental_abundances[atom] / abundance_total
        for atom in pc.gdat.gas.atoms_names
    ]
)


# %% [markdown]
# ## Pressure, temperature, and eddy diffusion
#
# We use the 10×-solar-metallicity morning-terminator profile distributed with the [VULCAN model](https://github.com/shami-EEG/VULCAN/blob/3ed92cf0222316d238efe9059959d29621962b17/atm/atm_W39b_10Xsolar_Twhole_morning_TP_20deg.txt). It averages ±10° around the western terminator in an exo-FMS general circulation model. We also follow Tsai et al. (2023) for the eddy-diffusion profile: $K_{zz}$ is fixed at $5\times10^7$ cm² s⁻¹ at pressures above 5 bar and increases in proportion to $P^{-1/2}$ at lower pressures.

# %%
pressure, temperature = np.loadtxt(
    tutorial_directory / "wasp39b_morning_pt.txt"
).T
pressure_bar = pressure / 1.0e6

# Eddy diffusion profile used by Tsai et al. (2023), in cm^2/s.
eddy = np.where(
    pressure_bar > 5.0,
    5.0e7,
    5.0e7 * np.sqrt(5.0 / pressure_bar),
)

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
# `initialize_to_climate_equilibrium_PT` calculates equilibrium chemistry along the supplied profile, estimates important chemical quench levels, and builds the photochemical grid deep enough to include them. The initial state includes an approximate treatment of quenched species. `metallicity=10` means ten times the solar elemental metallicity, whereas `CtoO=1` means the solar C/O ratio rather than an absolute C/O ratio of one.

# %%
pc.initialize_to_climate_equilibrium_PT(
    pressure,
    temperature,
    eddy,
    metallicity=10.0,
    CtoO=1.0,
)

print(
    "Prescribed pressure range: "
    f"{pressure_bar[0]:.1f} to {pressure_bar[-1]:.1e} bar"
)
print(
    "Explicit photochemical range: "
    f"{pc.wrk.pressure[0] / 1.0e6:.2f} to "
    f"{pc.wrk.pressure[-1] / 1.0e6:.1e} bar"
)


# %% [markdown]
# The `pc.gdat` attribute stores gas-giant-specific information, including the metallicity, C/O ratio, and the `ChemEquiAnalysis` object used for equilibrium calculations. The usual `pc.dat`, `pc.var`, and `pc.wrk` attributes remain available for the reaction mechanism, configurable variables, and evolving numerical state.
#
# `return_atmosphere` combines the explicit photochemical domain with the deeper equilibrium atmosphere by default.

# %%
initial_atmosphere = pc.return_atmosphere()

fig, axis = plt.subplots(figsize=(5.5, 5.0))
for species in ["H2", "He", "H2O", "CO", "CH4"]:
    axis.plot(
        initial_atmosphere[species],
        initial_atmosphere["pressure"] / 1.0e6,
        linewidth=2,
        label=species,
    )

axis.set_xscale("log")
axis.set_yscale("log")
axis.set_xlim(1.0e-10, 1.2)
axis.set_ylim(2.0e2, 1.0e-7)
axis.set_xlabel("Mixing ratio")
axis.set_ylabel("Pressure (bar)")
axis.set_title("Initialized atmosphere")
axis.grid(alpha=0.25)
axis.legend()
plt.show()


# %% [markdown]
# ## Evolve the atmosphere to steady state
#
# The simplest way to integrate the atmosphere is `pc.find_steady_state()`. As in the rocky planet tutorial, we instead use the robust stepper directly so that we can save intermediate states and see the photochemistry develop.

# %%
sulfur_species = ["H2S", "S", "S2", "SO", "SO2"]


def save_snapshot(model):
    """Copy the subset of atmospheric state used in the evolution plot."""
    atmosphere = model.return_atmosphere()
    snapshot = {
        "time": model.wrk.tn,
        "pressure": atmosphere["pressure"].copy(),
    }
    for species in sulfur_species:
        snapshot[species] = atmosphere[species].copy()
    return snapshot


snapshots = [save_snapshot(pc)]
pc.initialize_robust_stepper(pc.wrk.usol)

converged = False
give_up = False
robust_steps = 0
start_time = time.perf_counter()
try:
    while not converged and not give_up:
        give_up, converged = pc.robust_step()
        robust_steps += 1
        if robust_steps % 300 == 0 or converged or give_up:
            snapshots.append(save_snapshot(pc))
            print(
                f"Robust step {robust_steps:4d}: "
                f"t = {pc.wrk.tn:.3e} s, converged = {converged}"
            )
finally:
    pc.destroy_stepper()

elapsed = time.perf_counter() - start_time
if not converged:
    raise RuntimeError("The WASP-39b calculation did not reach steady state.")

print(f"Reached steady state in {elapsed:.1f} seconds.")


# %%
def plot_sulfur(atmosphere, title, axis=None, legend=True):
    """Plot selected sulfur-bearing gases against pressure."""
    if axis is None:
        _, axis = plt.subplots(figsize=(5.5, 5.0))
    for species in sulfur_species:
        axis.plot(
            atmosphere[species],
            atmosphere["pressure"] / 1.0e6,
            linewidth=2,
            label=species,
        )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(1.0e-12, 1.0e-3)
    axis.set_ylim(10.0, 1.0e-7)
    axis.set_xlabel("Mixing ratio")
    axis.set_ylabel("Pressure (bar)")
    axis.set_title(title)
    axis.grid(alpha=0.25)
    if legend:
        axis.legend()
    return axis


snapshot_indices = np.unique(
    np.linspace(0, len(snapshots) - 1, 4, dtype=int)
)
fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)

for panel_index, (axis, snapshot_index) in enumerate(
    zip(axes.flat, snapshot_indices)
):
    snapshot = snapshots[snapshot_index]
    title = (
        "Initial state"
        if snapshot["time"] == 0.0
        else f"t = {snapshot['time']:.1e} s"
    )
    plot_sulfur(
        snapshot,
        title,
        axis=axis,
        legend=panel_index == 0,
    )

plt.tight_layout()
plt.show()


# %% [markdown]
# Hydrogen sulfide supplied from the deep atmosphere is photolyzed and converted into several sulfur-bearing products. Passing `equilibrium=True` to `return_atmosphere` evaluates equilibrium composition across the same returned pressure grid, allowing a direct comparison. The solid curves below are the photochemical solution, while the dotted curves show the equilibrium abundances.

# %%
steady_atmosphere = pc.return_atmosphere()
equilibrium_atmosphere = pc.return_atmosphere(equilibrium=True)

fig, axis = plt.subplots(figsize=(6.5, 5.0))
colors = dict(zip(sulfur_species, ["C1", "gold", "0.5", "C0", "C4"]))

for species in sulfur_species:
    axis.plot(
        steady_atmosphere[species],
        steady_atmosphere["pressure"] / 1.0e6,
        color=colors[species],
        linewidth=2,
        label=species,
    )
    axis.plot(
        equilibrium_atmosphere[species],
        equilibrium_atmosphere["pressure"] / 1.0e6,
        color=colors[species],
        linestyle=":",
        linewidth=1.5,
    )

axis.set_xscale("log")
axis.set_yscale("log")
axis.set_xlim(1.0e-12, 1.0e-3)
axis.set_ylim(10.0, 1.0e-7)
axis.set_xlabel("Mixing ratio")
axis.set_ylabel("Pressure (bar)")
axis.set_title("Photochemical (solid) and equilibrium (dotted)")
axis.grid(alpha=0.25)
axis.legend()
plt.show()


# %% [markdown]
# The enhanced SO₂ between roughly 10⁻² and 10⁻⁶ bar broadly reproduces the photochemical behavior found by Tsai et al. (2023). Differences are expected because TOI-193 is a proxy for the WASP-39 stellar spectrum and the chemical networks are not identical.
