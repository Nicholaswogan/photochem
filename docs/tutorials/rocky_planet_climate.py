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
# # Rocky planet climate
#
# Photochem's climate model is contained in the `AdiabatClimate` class. It constructs one-dimensional atmospheric profiles, calculates their radiative properties, and solves for approximate or full radiative-convective equilibrium. The radiative calculations are multithreaded; here we limit them to four threads.

# %%
import numpy as np
from astropy import constants
from matplotlib import pyplot as plt
from threadpoolctl import threadpool_limits

from photochem.clima import AdiabatClimate
from photochem.utils import (
    settings_file_for_climate,
    species_file_for_climate,
    stars,
)

_ = threadpool_limits(limits=4)

# %% [markdown]
# ## Initialization
#
# Initializing the model requires three files: a species file, a settings file, and a stellar-flux file. Below, we generate each file on the fly.
#
# **Species file:** This file defines the gases included in the climate model, their elemental compositions and thermodynamic properties, and which gases may condense. Here, H₂O and CO₂ are condensable while N₂ is not.

# %%
species_file_for_climate(
    filename="rocky_planet_climate/species.yaml",
    species=["H2O", "CO2", "N2"],
    condensates=["H2O", "CO2"],
)

# %% [markdown]
# **Settings file:** This file defines the atmospheric grid, planetary mass and radius, surface albedo, and radiative-transfer configuration. We use Earth's mass and radius with a surface albedo of 0.2.

# %%
settings_file_for_climate(
    filename="rocky_planet_climate/settings.yaml",
    planet_mass=float(constants.M_earth.cgs.value),
    planet_radius=float(constants.R_earth.cgs.value),
    surface_albedo=0.2,
)

# %% [markdown]
# **Stellar flux file:** This file supplies the wavelength-dependent stellar irradiance used for shortwave radiative transfer. `solar_spectrum` constructs the modern solar spectrum from the packaged ATLAS-1 reference, extends it to long wavelengths with a solar blackbody, and rebins it for Photochem's photochemical and climate models.

# %%
_ = stars.solar_spectrum(outputfile="rocky_planet_climate/modern_sun.txt")

# %% [markdown]
# Now we can initialize `AdiabatClimate`:

# %%
c = AdiabatClimate(
    species_file="rocky_planet_climate/species.yaml",
    settings_file="rocky_planet_climate/settings.yaml",
    flux_file="rocky_planet_climate/modern_sun.txt",
)

# %% [markdown]
# Once initialized, the `AdiabatClimate` object contains all the information needed to do climate calculations.

# %%
print(f"Climate species: {c.species_names}")

# %% [markdown]
# ## Constructing an atmosphere
#
# Given a surface temperature and surface partial pressures, `AdiabatClimate.make_profile` constructs a multispecies pseudoadiabat following Equation 1 of [Graham et al. (2021)](https://doi.org/10.3847/PSJ/ac214c), connected to an isothermal stratosphere. Condensable material in excess of its saturation vapor pressure is placed in a surface reservoir.

# %%
# Surface partial pressures in species_names order
surface_partial_pressures = np.full(len(c.species_names), 1.0e-15)
surface_partial_pressures[c.species_names.index("H2O")] = 270.0  # bar
surface_partial_pressures[c.species_names.index("N2")] = 1.0
surface_partial_pressures[c.species_names.index("CO2")] = 10.0
surface_partial_pressures *= 1.0e6  # bar to dyn/cm^2

surface_temperature = 280.0  # K
c.T_trop = 150.0  # K

c.make_profile(surface_temperature, surface_partial_pressures)

# %% [markdown]
# The object `c` now has attributes that describe the atmosphere:
#
# - `c.P`, pressures in atmospheric layers (dyn cm⁻²)
# - `c.P_surf`, surface pressure (dyn cm⁻²)
# - `c.T`, temperature in atmospheric layers (K)
# - `c.f_i`, mixing ratios of all species in atmospheric layers (2-D array)
# - `c.z`, altitude at the center of each atmospheric layer (cm)
# - `c.dz`, thickness of each atmospheric layer (cm)
# - `c.N_atmos`, reservoir of each gas in the atmosphere (mol cm⁻²)
# - `c.N_surface`, reservoir of each gas on the surface (mol cm⁻²)
# - `c.N_ocean`, matrix of gases dissolved in oceans formed by each condensate (mol cm⁻²)
#
# Here is an illustrative plot:

# %%
fig, ax = plt.subplots(figsize=(5, 4))

for index, species in enumerate(c.species_names):
    ax.plot(c.f_i[:, index], c.P / 1.0e6, linewidth=2, label=species)
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(c.P_surf / 1.0e6, c.P_top / 1.0e6)
ax.legend(bbox_to_anchor=(1.1, 1.02), loc="upper left")
ax.grid()
ax.set_ylabel("Pressure (bar)")
ax.set_xlabel("Mixing ratio")

temperature_axis = ax.twiny()

temperature_axis.plot(c.T, c.P / 1.0e6, "k--", linewidth=2, label="Temperature")
temperature_axis.set_xlabel("Temperature (K)")
temperature_axis.legend(bbox_to_anchor=(1.1, 0.2), loc="upper left")

plt.show()

# %% [markdown]
# In the figure above, the tropopause is near 10⁻² bar. Below that, the slope of the tropospheric temperature profile changes near 3 bar when CO₂ begins to condense.
#
# `c.N_atmos` and `c.N_surface` contain one reservoir per species. `c.N_ocean[:, j]` contains the gases dissolved in the ocean formed by species `j`. All reservoirs have units of mol cm⁻². Only species that are saturated at the surface have a surface reservoir. In this calculation, water is the only species with a surface reservoir, which can be interpreted as an ocean. We have not configured gas dissolution, so the dissolved reservoirs are zero.

# %%
print(f"Atmospheric reservoirs (mol/cm^2): {c.N_atmos}")
print(f"Surface reservoirs (mol/cm^2): {c.N_surface}")
print(f"Ocean-dissolved reservoir matrix (mol/cm^2): {c.N_ocean}")

# %% [markdown]
# The surface pressure is about 11.01 bar from 10 bar of CO₂, 1 bar of N₂, and a little water vapor.

# %%
print(f"Surface pressure: {c.P_surf / 1.0e6:.3f} bar")

# %% [markdown]
# Two related routines construct atmospheres from different inputs:
#
# - `AdiabatClimate.make_column` constructs an atmosphere from total volatile reservoirs in mol cm⁻².
# - `AdiabatClimate.make_profile_bg_gas` constructs an atmosphere in which one molecule is treated as the background gas.
#

# %% [markdown]
# ## Radiative transfer
#
# Atmospheric-construction routines have corresponding methods that also calculate radiative properties. For example, `AdiabatClimate.TOA_fluxes` calls `AdiabatClimate.make_profile`, performs shortwave and longwave radiative transfer, and returns the net incoming stellar radiation and outgoing longwave radiation at the top of the atmosphere.

# %%
surface_partial_pressures = np.full(len(c.species_names), 1.0e-10)
surface_partial_pressures[c.species_names.index("H2O")] = 270.0  # bar
surface_partial_pressures[c.species_names.index("N2")] = 1.0
surface_partial_pressures[c.species_names.index("CO2")] = 400.0e-6
surface_partial_pressures *= 1.0e6  # bar to dyn/cm^2

surface_temperature = 280.0  # K
c.T_trop = 215.0  # K

absorbed_stellar_radiation, outgoing_longwave_radiation = c.TOA_fluxes(
    surface_temperature,
    surface_partial_pressures,
)
print(f"Absorbed stellar radiation: {absorbed_stellar_radiation / 1.0e3:.2f} W/m^2")
print(f"Outgoing longwave radiation: {outgoing_longwave_radiation / 1.0e3:.2f} W/m^2")

# %% [markdown]
# We can look at the emission of the planet as a function of wavelength.

# %%
frequency_edges = c.rad.ir.freq  # Hz
frequency_centers = 0.5 * (frequency_edges[1:] + frequency_edges[:-1])
wavelength_centers = 1.0e6 * constants.c.value / frequency_centers  # microns

# Convert the top-of-atmosphere flux from mW/m^2/Hz to W/m^2/micron.
spectral_flux = (
    1.0e-3
    * c.rad.wrk_ir.fup_a[-1, :]
    * frequency_centers
    / wavelength_centers
)
wavelength_edges = c.rad.ir.wavl / 1.0e3  # nm to microns

fig, ax = plt.subplots(figsize=(8, 5))

ax.stairs(spectral_flux, wavelength_edges, linewidth=2, color="black")

ax.set_xscale("log")
ax.set_xlim(2, 100)
ax.set_ylim(0, ax.get_ylim()[1])
ax.grid(alpha=0.4)

ax.set_ylabel("Top-of-atmosphere thermal flux\n" + r"(W m$^{-2}$ $\mu$m$^{-1}$)")
ax.set_xlabel(r"Wavelength ($\mu$m)")

plt.show()

# %% [markdown]
# ## Steady-state climate modeling (approximate climates)
#
# Up to this point we have constructed atmospheres and calculated their radiative properties. We can obtain an approximate steady-state climate by retaining the assumed pseudoadiabat and isothermal stratosphere, then solving for the surface temperature at which absorbed stellar radiation balances outgoing longwave radiation. `AdiabatClimate.surface_temperature` performs this nonlinear solve.
#
# Below, we compute a steady-state climate for a generic habitable planet.

# %%
c.T_trop = 215.0  # Prescribed stratospheric temperature (K)
c.RH = np.full(len(c.species_names), 0.5)

surface_temperature = c.surface_temperature(
    surface_partial_pressures,
    T_guess=280.0,
)
print(f"Approximate equilibrium surface temperature: {surface_temperature:.2f} K")

# %% [markdown]
# ## Steady-state climate modeling (full radiative-convective equilibrium)
#
# The calculation above restricts the atmosphere to an adiabat connected to an isothermal stratosphere. `AdiabatClimate.RCE` instead adjusts the full temperature profile to solve the more computationally expensive radiative-convective equilibrium problem. RCE is sensitive to its initial temperature profile, so we initialize it with the approximate equilibrium solution calculated above.

# %%
c.verbose = False
converged = c.RCE(
    surface_partial_pressures,
    c.T_surf,
    c.T,
    c.convecting_with_below,
)
if not converged:
    raise RuntimeError("The radiative-convective equilibrium solve did not converge.")

print(f"RCE converged: {converged}")
print(f"RCE surface temperature: {c.T_surf:.2f} K")

# %%
fig, ax = plt.subplots(figsize=(5, 4))

for index, species in enumerate(c.species_names):
    ax.plot(c.f_i[:, index], c.P / 1.0e6, linewidth=2, label=species)
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(c.P_surf / 1.0e6, c.P_top / 1.0e6)
ax.legend(bbox_to_anchor=(1.1, 1.02), loc="upper left")
ax.grid()
ax.set_ylabel("Pressure (bar)")
ax.set_xlabel("Mixing ratio")

temperature_axis = ax.twiny()

temperature_axis.plot(c.T, c.P / 1.0e6, "k--", linewidth=2, label="Temperature")
temperature_axis.set_xlabel("Temperature (K)")
temperature_axis.legend(bbox_to_anchor=(1.1, 0.2), loc="upper left")
temperature_axis.set_xlim(130, 300)

plt.show()

# %%
fig, ax = plt.subplots(figsize=(5, 4))

total_radiative_flux = (
    c.rad.wrk_sol.fdn_n[1:-2:2]
    - c.rad.wrk_sol.fup_n[1:-2:2]
    + c.rad.wrk_ir.fdn_n[1:-2:2]
    - c.rad.wrk_ir.fup_n[1:-2:2]
)
total_radiative_flux *= 1.0e-3  # mW/m^2 to W/m^2

ax.plot(total_radiative_flux, c.P / 1.0e6)
ax.axvline(0.0, color="black", linewidth=1, linestyle="--")

ax.set_yscale("log")
ax.set_ylim(c.P_surf / 1.0e6, c.P_top / 1.0e6)
ax.set_xlabel(r"Net downward radiative flux (W/m$^2$)")
ax.set_ylabel("Pressure (bar)")
ax.grid(alpha=0.4)

plt.show()

# %% [markdown]
# The net radiative flux approaches zero in regions that are in radiative equilibrium. In convective regions, departures from zero indicate the radiative imbalance that is offset by convective energy transport.
#
# This tutorial demonstrated the main `AdiabatClimate` workflow: generating its input files, constructing a condensable atmosphere, calculating top-of-atmosphere and wavelength-dependent fluxes, and solving for both approximate and full radiative-convective equilibrium climates.
