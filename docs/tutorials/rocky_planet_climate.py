# # Rocky planet climate
#
# The climate model is contained within the `AdiabatClimate` class, which we import below along with several other standard modules. The climate model uses multithreading so we can set the threadcount to a good value for your given machine.

# +
from pathlib import Path
import numpy as np
from astropy import constants
from matplotlib import pyplot as plt

# Climate model
from photochem.clima import AdiabatClimate
from photochem.utils import species_file_for_climate, settings_file_for_climate
from photochem.utils import stars

# Also set the thread count
from threadpoolctl import threadpool_limits
_ = threadpool_limits(limits=4)

# Support execution by MkDocs from the repository root and interactive
# execution from docs/tutorials.
tutorial_directory = Path("docs/tutorials/rocky_planet_climate")
if not tutorial_directory.is_dir():
    tutorial_directory = Path("rocky_planet_climate")
# -

# ## Initialization
#
# Initializing the code requires three files: A species file, a settings file and a stellar spectrum file. Below, we generate each file on-the-fly.
#
# **Species file:** The species file specifies the different molecules in the code...

species_file = str(tutorial_directory / "species.yaml")
species_file_for_climate(
    filename=species_file,
    species=['H2O', 'CO2', 'N2'],
    condensates=['H2O', 'CO2']
)

# **Settings file:** The settings file...

settings_file = str(tutorial_directory / "settings.yaml")
settings_file_for_climate(
    filename=settings_file,
    planet_mass=float(constants.M_earth.cgs.value),
    planet_radius=float(constants.R_earth.cgs.value),
    surface_albedo=0.2
)

# **Stellar flux file:** We can generate a solar spectrum with...

flux_file = str(tutorial_directory / "modern_sun.txt")
_ = stars.solar_spectrum(outputfile=flux_file)

# Now we can initialize `AdiabatClimate`:

c = AdiabatClimate(
    species_file=species_file,
    settings_file=settings_file,
    flux_file=flux_file,
)

# Once initialized, the `AdiabatClimate` object contains all the information needed to do climate calculations.

print(c.species_names) # species in the model

# ## Constructing an atmosphere
#
# Given the surface temperature and volatile inventories `AdiabatClimate` is able to draw a temperature profile and distribute gases in the atmosphere. One function that does this is `AdiabatClimate.make_profile`, which is illustrated below.

# +
# The surface volatile inventories
P_i = np.ones(len(c.species_names))*1e-15
P_i[c.species_names.index('H2O')] = 270 # bar
P_i[c.species_names.index('N2')] = 1
P_i[c.species_names.index('CO2')] = 10
P_i *= 1e6 # Convert to dynes/cm^2

T_surf = 280 # surface temperature (K)
c.T_trop = 150 # tropopause temperature (K). 

# Integrates a multispecies pseudoadiabat upward following
# Equation (1) in Graham et al. 2021 (doi.org/10.3847/PSJ/ac214c)
c.make_profile(T_surf, P_i)
# -

# The object `c` now has attributes that describe the atmosphere
# - `c.P`, pressures in atmospheric layers (dynes/cm$^2$)
# - `c.P_surf` surface pressure (dynes/cm$^2$)
# - `c.T`, temperature in atmospheric layers (K)
# - `c.f_i`, mixing ratios of all species in atmospheric layers (2-D array)
# - `c.z`, altitude of center of atmospheric layers (cm)
# - `c.dz`, thickness of each atmospheric layers (cm)
# - `c.N_atmos`, reservoir of each gas in the atmosphere (mol/cm$^2$)
# - `c.N_surface`, reservoir of each gas on the surface (mol/cm$^2$)
# - `c.N_ocean`, reservoir of each gas dissolved in oceans (mol/cm$^2$)
#
# Here is an illustrative plot:

# +
fig,ax = plt.subplots(1,1,figsize=[5,4])

for i in range(len(c.species_names)):
    ax.plot(c.f_i[:,i], c.P/1e6, lw=2, label=c.species_names[i])
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(c.P[0]/1e6,c.P_top/1e6)
ax.legend(ncol=1,bbox_to_anchor=(1.1, 1.02), loc='upper left')
ax.grid()
ax.set_ylabel('Pressure (bar)')
ax.set_xlabel('Mixing ratio')

ax1 = ax.twiny()

ax1.plot(c.T, c.P/1e6, 'k--', lw=2, label='Temperature')
ax1.set_xlabel('Temperature (K)')
ax1.legend(ncol=1,bbox_to_anchor=(1.1, .2), loc='upper left')

plt.show()
# -

# In the figure above, the tropopause is at $10^{-2}$ bars. Below that, in the troposphere, the slope of the temperature profile changes at about 3 bars pressure. This is caused by the start of CO$_2$ condensation.
#
# `c.N_atmos`, `c.N_surface` and `c.N_ocean` contains the reservoir of volatiles in mol/cm$^{2}$ in the atmosphere, surface, and dissolved in oceans. Only gases that are saturated at surface have a surface reservoir. In this case, only water has a surface reservoir, which can be thought of as an ocean. We have not included gas dissolution in oceans, so the model reports no dissolved gases. Later in this notebook we implement gas dissolution in oceans.  

print(c.N_atmos)

print(c.N_surface)

print(c.N_ocean)

# The surface pressure is about 11.01 bars from 10 bars of CO$_2$, 1 bar of N$_2$ and a little water vapor.

print(c.P_surf/1e6) # bars

# The code contains two other routines for constructing atmosphere, which you can check out if you are interested:
#
# - `AdiabatClimate.make_column` - Constructs an atmosphere given surface volatile reservoirs in moles/cm$^2$.
# - `AdiabatClimate.make_profile_bg_gas` - Constructs an atmosphere where one molecule is treated as a background gas.
#

# ## Radiative transfer
#
# Routines for constructing atmospheres have corresponding methods that both construct and atomsphere and subsequently compute that atmosphere's radiative properties. For example, the routine `AdiabatClimate.TOA_fluxes` calls `AdiabatClimate.make_profile`, then does solar and infrared radiative transfer, returning the fluxes at the top of the atmosphere.

# +
P_i = np.ones(len(c.species_names))*1e-10
P_i[c.species_names.index('H2O')] = 270 # bar
P_i[c.species_names.index('N2')] = 1
P_i[c.species_names.index('CO2')] = 400e-6
P_i *= 1e6 # Convert to dynes/cm^2

T_surf = 280
c.T_trop = 215

# Call TOA fluxes
ASR, OLR = c.TOA_fluxes(T_surf, P_i)
print(ASR) # mW/m^2
print(OLR) # mW/m^2
# -

# We can look at the emission of the planet as a function of wavelength.

# +
freq = c.rad.ir.freq # Hz
freq_av = (freq[1:]+freq[:-1])/2 # Hz
wv_av = 1e6*constants.c.value/freq_av # microns
# c.rad.wrk_ir.fup_a is mW/m^2/Hz. Here I convert to W/m^2/um
F = 1e-3*c.rad.wrk_ir.fup_a[-1,:]*(freq_av/(wv_av))

fig,ax = plt.subplots(1,1,figsize=[8,5])

ax.plot(c.rad.ir.wavl[1:]/1e3, F, drawstyle='steps-pre',lw=2,c='k')

ax.set_xscale("log")
ax.set_xlim(2,100)
ax.set_ylim(0,ax.get_ylim()[1])
ax.grid(alpha=0.4)

ax.set_ylabel('Top of atmosphere thermal radiance\n'+r'(W m$^{-2}$ $\mu$m$^{-1}$)')
ax.set_xlabel(r'Wavelength ($\mu$m)')

plt.show()
# -

# ## Steady-state climate modeling (approximate climates)
#
# Up to this point we have merely constructed atmospheres and computed their radiative properties. In this section we calculate approximate steady-state climates, where we assume an adiabat connected to an isothermal stratosphere. To find steady states, we can use the `AdiabatClimate.surface_temperature` routine which numerically solves (via a variation of Newton's method) for the surface temperature which results in an atmosphere that balances incoming solar and outgoing long-wave radiation.
#
# Below, we compute a steady-state climate for a generic habitable planet.

# +
P_i = np.ones(len(c.species_names))*1e-10
P_i[c.species_names.index('H2O')] = 270 # bar
P_i[c.species_names.index('N2')] = 1
P_i[c.species_names.index('CO2')] = 400e-6
P_i *= 1e6 # Convert to dynes/cm^2

c.T_trop = 215 # Have to set the stratosphere temp.
c.RH = np.ones(len(c.species_names))*0.5 # Change relative humidity to 0.5

# Call TOA fluxes
T_surf = c.surface_temperature(P_i, T_guess=280)
print(T_surf)
# -

# ## Steady-state climate modeling (full radiative-convective equilibrium)
#
# In the calculation above, we made climate a much easier problem by making assumptions about the temperature profile. Here, we use the `AdiabatClimate.RCE` routine to do the more challenging (and computationally expensive) problem of computing full radiative-convective equilibrium (RCE). The `AdiabatClimate.RCE` routine can fail to converge if given a bad starting guess for the temperature profile.

# +
P_i = np.ones(len(c.species_names))*1e-10
P_i[c.species_names.index('H2O')] = 270 # bar
P_i[c.species_names.index('N2')] = 1
P_i[c.species_names.index('CO2')] = 400e-6
P_i *= 1e6 # Convert to dynes/cm^2

c.RH = np.ones(len(c.species_names))*0.5 # Change relative humidity to 0.5

# Get a starting guess for the temperature profile
T_surf = c.surface_temperature(P_i, T_guess=280)

# Solve for RCE with that starting guess 
converged = c.RCE(P_i, c.T_surf, c.T, c.convecting_with_below)
print(c.T_surf)

# +
fig,ax = plt.subplots(1,1,figsize=[5,4])

for i in range(len(c.species_names)):
    ax.plot(c.f_i[:,i], c.P/1e6, lw=2, label=c.species_names[i])
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(c.P_surf/1e6,c.P_top/1e6)
ax.legend(ncol=1,bbox_to_anchor=(1.1, 1.02), loc='upper left')
ax.grid()
ax.set_ylabel('Pressure (bar)')
ax.set_xlabel('Mixing ratio')

ax1 = ax.twiny()

ax1.plot(c.T, c.P/1e6, 'k--', lw=2, label='Temperature')
ax1.set_xlabel('Temperature (K)')
ax1.legend(ncol=1,bbox_to_anchor=(1.1, .2), loc='upper left')
ax1.set_xlim(130,300)

plt.show()

# +
fig,ax = plt.subplots(1,1,figsize=[5,4])

f_total = (c.rad.wrk_sol.fdn_n[1:-2:2] - c.rad.wrk_sol.fup_n[1:-2:2]) + (c.rad.wrk_ir.fdn_n[1:-2:2] - c.rad.wrk_ir.fup_n[1:-2:2])
f_total *= 1e-3 # W/m^2

ax.plot(f_total, c.P/1e6)

ax.set_yscale('log')
ax.set_ylim(c.P_surf/1e6,1e-4)
ax.set_xlabel('Total radiative flux (W/m$^2$)')
ax.set_ylabel('Pressure (bar)')

plt.show()
# -


