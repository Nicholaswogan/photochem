"""Port of the nominal legacy Venus chemistry case, with all inputs local.

No climate solve, equilibrium initialization, Crisp clouds, or extra Rimmer
chlorine reactions are used. Built-in H2SO4 optics and the Rimmer UV absorber
are retained. All dimensional values below use CGS except wavelength (nm).
"""
from pathlib import Path
import shutil
import numpy as np
from photochem import EvoAtmosphere

HERE = Path(__file__).resolve().parent
MECHANISM = 'mechanism.yaml'
SETTINGS = 'settings.yaml'
CONTROLS = dict(atol=1e-17, rtol=1e-3, equilibrium_time=3.15e13,
                smooth_factor=5., k_cond=100., k_evap=100.,
                upwind_molec_diff=True, jacobian_method=3)


def initialize(destination, fresh=False):
    shutil.copyfile(HERE / 'flux.txt', destination / 'flux.txt')
    pc = EvoAtmosphere(str(HERE / 'mechanism.yaml'), str(HERE / 'settings.yaml'),
                       str(HERE / 'flux.txt'), None if fresh else str(HERE / 'atmosphere.txt'))
    if fresh:
        z, temperature, pressure = np.loadtxt(HERE / 'venus_seiff1985.txt', skiprows=2).T
        z = np.append(z, 110.) * 1e5
        temperature = np.append(temperature, temperature[-1])
        _, _, _, kzz_input, z_input = np.loadtxt(HERE / 'venus_rimmer2021.dat', skiprows=2).T
        edd = 10**np.interp(z, z_input, np.log10(kzz_input))
        surface_mix = {'CO2': .965, 'N2': .035, 'SO2': 100e-6, 'H2O': 35e-6,
                       'OCS': 5e-6, 'CO': 10e-6, 'HCl': 500e-9, 'H2': 4.5e-9}
        mix = {name: np.full(z.size, value) for name, value in surface_mix.items()}
        pc.initialize_atmosphere_z(z, temperature, edd, pressure[0]*1e6, mix=mix)
        # Legacy setup used an absolute density floor for unspecified species.
        densities = pc.wrk.usol
        for i, name in enumerate(pc.dat.species_names[:densities.shape[0]]):
            if name not in mix:
                densities[i, :] = 1e-40
        pc.wrk.usol = densities
    # Retain legacy diffusion; use the current analytical Jacobian for benchmarks.
    pc.var.upwind_molec_diff = True
    pc.var.jacobian_method = 3
    for params in pc.var.cond_params:
        params.smooth_factor = 5.
        params.k_cond = 100.
        params.k_evap = 100.
    radii = pc.var.particle_radius
    for name, radius in {'H2SO4aer': 2e-4, 'H2Oaer': 1e-3}.items():
        radii[pc.dat.species_names.index(name), :] = radius
    pc.var.particle_radius = radii
    pc.prep_atmosphere(pc.wrk.usol)
    install_uv_absorber(pc)
    pc.var.atol = 1e-17
    pc.var.equilibrium_time = 3.15e13
    pc.var.verbose = 0
    return pc


def install_uv_absorber(pc):
    """Install nominal Rimmer extinction/cm against initial hydrostatic P.

    The formula was applied on fixed altitude layers by the legacy code. The
    current API preserves its initial profile on a pressure coordinate instead.
    Pressure extrapolation follows the current constant-endpoint policy.
    Wavelength is in nm.
    """
    z = pc.var.z
    wavelength = .5*(pc.dat.wavl[1:]+pc.dat.wavl[:-1])
    coefficient = np.zeros_like(z)
    coefficient[(z > 58e5) & (z <= 67e5)] = .056e-5
    above = z > 67e5
    coefficient[above] = .056e-5*np.exp(-(z[above]-67e5)/3e5)
    extinction = coefficient[:, None]*np.exp(-(wavelength[None, :]-360.)/100.)
    extinction[:, wavelength > 400.] = 0.
    pc.set_custom_optical_properties(wavelength, pc.wrk.pressure_hydro,
                                    extinction, np.zeros_like(extinction),
                                    np.zeros_like(extinction), tau_mode=1)
