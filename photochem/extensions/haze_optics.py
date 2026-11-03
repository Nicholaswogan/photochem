"""Mie optical properties of compact, monodisperse haze particles."""

import numpy as np
import numba as nb
from miepython.mie_jit import _single_sphere_nb


# Khare et al. (1984), Icarus 60, 127-137.
# Wavelengths in cm and refractive indices in the convention m = n - i*k.
KHARE_WAVELENGTH = np.array([
    2.0700e-02*1e-4,
    3.1200e-02*1e-4,
    4.1500e-02*1e-4,
    4.4900e-02*1e-4,
    5.8800e-02*1e-4,
    7.8500e-02*1e-4,
    8.0000e-02*1e-4,
    9.2500e-02*1e-4,
    1.0160e-01*1e-4,
    1.0970e-01*1e-4,
    1.1590e-01*1e-4,
    1.1810e-01*1e-4,
    1.2150e-01*1e-4,
    1.3620e-01*1e-4,
    1.6310e-01*1e-4,
    1.9680e-01*1e-4,
    2.3840e-01*1e-4,
    2.6380e-01*1e-4,
    2.9520e-01*1e-4,
    3.2630e-01*1e-4,
    3.5420e-01*1e-4,
    3.8740e-01*1e-4,
    4.1330e-01*1e-4,
    4.4280e-01*1e-4,
    5.6350e-01*1e-4,
    6.8880e-01*1e-4,
    8.7310e-01*1e-4,
    1.0160e+00*1e-4,
    1.1480e+00*1e-4,
    1.1920e+00*1e-4,
    1.3570e+00*1e-4,
    1.3810e+00*1e-4,
    1.8020e+00*1e-4,
    1.8130e+00*1e-4,
    1.8230e+00*1e-4,
    1.8730e+00*1e-4,
    2.0190e+00*1e-4,
    2.2140e+00*1e-4,
    2.3930e+00*1e-4,
    2.4030e+00*1e-4,
    2.4220e+00*1e-4,
    2.6950e+00*1e-4,
    2.7430e+00*1e-4,
    2.8180e+00*1e-4,
    2.9380e+00*1e-4,
    3.0090e+00*1e-4,
    3.2460e+00*1e-4,
    3.4630e+00*1e-4,
    3.6680e+00*1e-4,
    3.9480e+00*1e-4,
    4.2170e+00*1e-4,
    4.4280e+00*1e-4,
    4.4920e+00*1e-4,
    4.5920e+00*1e-4,
    4.6260e+00*1e-4,
    4.8810e+00*1e-4,
    5.1660e+00*1e-4,
    5.2530e+00*1e-4,
    5.4380e+00*1e-4,
    5.7400e+00*1e-4,
    5.8480e+00*1e-4,
    5.9610e+00*1e-4,
    6.3260e+00*1e-4,
    6.4570e+00*1e-4,
    6.6660e+00*1e-4,
    7.0440e+00*1e-4,
    7.6530e+00*1e-4,
    8.7310e+00*1e-4,
    1.0510e+01*1e-4,
    1.1070e+01*1e-4,
    1.1700e+01*1e-4,
    1.2910e+01*1e-4,
    1.4420e+01*1e-4,
    1.7710e+01*1e-4,
    1.8230e+01*1e-4,
    2.2140e+01*1e-4,
    3.1000e+01*1e-4,
    3.6460e+01*1e-4,
    5.6350e+01*1e-4,
    8.1570e+01*1e-4,
    1.2150e+02*1e-4,
    1.4090e+02*1e-4,
    1.7220e+02*1e-4,
    2.2960e+02*1e-4,
    3.8740e+02*1e-4,
    5.6350e+02*1e-4,
    6.8880e+02*1e-4,
    7.7490e+02*1e-4,
    8.5000e+02*1e-4,
    9.2000e+02*1e-4,
])

KHARE_REFRACTIVE_INDEX = np.array([
    0.9200 - 0.0490j,
    0.8500 - 0.1400j,
    0.8020 - 0.3100j,
    0.8120 - 0.3800j,
    0.9630 - 0.6200j,
    1.3300 - 0.7700j,
    1.3700 - 0.7700j,
    1.5800 - 0.6700j,
    1.6700 - 0.5800j,
    1.7200 - 0.5000j,
    1.7500 - 0.4300j,
    1.7500 - 0.4000j,
    1.7400 - 0.3700j,
    1.7000 - 0.2700j,
    1.6500 - 0.2400j,
    1.6600 - 0.2200j,
    1.6800 - 0.2100j,
    1.6800 - 0.1800j,
    1.6600 - 0.1500j,
    1.6400 - 0.1300j,
    1.6300 - 0.1100j,
    1.6600 - 0.0910j,
    1.6900 - 0.0760j,
    1.7200 - 0.0600j,
    1.7000 - 0.0230j,
    1.6800 - 0.0088j,
    1.6600 - 0.0024j,
    1.6500 - 0.0010j,
    1.6500 - 0.0006j,
    1.6500 - 0.0005j,
    1.6400 - 0.0004j,
    1.6400 - 0.0004j,
    1.6400 - 0.0004j,
    1.6400 - 0.0004j,
    1.6400 - 0.0004j,
    1.6300 - 0.0005j,
    1.6300 - 0.0007j,
    1.6300 - 0.0009j,
    1.6200 - 0.0008j,
    1.6200 - 0.0009j,
    1.6200 - 0.0012j,
    1.6000 - 0.0041j,
    1.5900 - 0.0110j,
    1.5800 - 0.0240j,
    1.5900 - 0.0600j,
    1.6100 - 0.0750j,
    1.6500 - 0.0560j,
    1.6500 - 0.0350j,
    1.6500 - 0.0210j,
    1.6400 - 0.0130j,
    1.6300 - 0.0100j,
    1.6100 - 0.0110j,
    1.6100 - 0.0140j,
    1.6100 - 0.0280j,
    1.6200 - 0.0270j,
    1.6100 - 0.0100j,
    1.5800 - 0.0076j,
    1.5800 - 0.0087j,
    1.5500 - 0.0110j,
    1.4800 - 0.0290j,
    1.4400 - 0.0700j,
    1.4300 - 0.1500j,
    1.5800 - 0.2800j,
    1.6400 - 0.2600j,
    1.6900 - 0.2200j,
    1.6900 - 0.1700j,
    1.7200 - 0.1300j,
    1.7100 - 0.0920j,
    1.6700 - 0.0750j,
    1.6600 - 0.0790j,
    1.6400 - 0.0970j,
    1.6700 - 0.1400j,
    1.6700 - 0.1700j,
    1.7400 - 0.2100j,
    1.7600 - 0.2100j,
    1.8000 - 0.1800j,
    1.8100 - 0.1400j,
    1.8100 - 0.1500j,
    1.8600 - 0.1900j,
    1.9300 - 0.2100j,
    2.0200 - 0.1900j,
    2.0300 - 0.1600j,
    2.0400 - 0.1600j,
    2.0700 - 0.1400j,
    2.1200 - 0.1000j,
    2.1500 - 0.0700j,
    2.1600 - 0.0470j,
    2.1600 - 0.0290j,
    2.1700 - 0.0100j,
    2.1700 - 0.0030j,
])
KHARE_WAVELENGTH.setflags(write=False)
KHARE_REFRACTIVE_INDEX.setflags(write=False)


@nb.njit
def haze_optics(wavelength, particle_radius, particle_number_density, refractive_index=None):
    """Compute Mie optical properties.

    Parameters
    ----------
    wavelength : ndarray, shape (n_wavelength,)
        Positive, strictly increasing wavelengths (cm).
    particle_radius : ndarray, shape (n_points,)
        Positive particle radii (cm), ordered from bottom to top.
    particle_number_density : ndarray, shape (n_points,)
        Nonnegative particle number densities (particles/cm^3).
    refractive_index : ndarray of complex, shape (n_wavelength,), optional
        Index at each wavelength, using m = n - i*k with n > 0 and k >= 0.
        Defaults to Khare et al. (1984) tholins, interpolating real and
        imaginary parts linearly in log wavelength with constant extrapolation
        beyond the tabulated wavelength range.

    Returns
    -------
    dtau_dz : ndarray, shape (n_points, n_wavelength)
        Extinction coefficient (1/cm); not integrated layer optical depth.
    w0 : ndarray, shape (n_points, n_wavelength)
        Single-scattering albedo (dimensionless).
    g0 : ndarray, shape (n_points, n_wavelength)
        Scattering asymmetry parameter (dimensionless).

    """
    wavelength = np.asarray(wavelength, dtype=np.float64)
    particle_radius = np.asarray(particle_radius, dtype=np.float64)
    particle_number_density = np.asarray(particle_number_density, dtype=np.float64)
    if wavelength.ndim != 1 or wavelength.size == 0 or not np.all(np.isfinite(wavelength)):
        raise ValueError('wavelength must be a nonempty, finite, one-dimensional array.')
    if particle_radius.ndim != 1 or particle_radius.size == 0 or not np.all(np.isfinite(particle_radius)):
        raise ValueError('particle_radius must be a nonempty, finite, one-dimensional array.')
    if (particle_number_density.ndim != 1 or particle_number_density.size == 0
            or not np.all(np.isfinite(particle_number_density))):
        raise ValueError('particle_number_density must be a nonempty, finite, one-dimensional array.')
    if np.any(wavelength <= 0.0) or np.any(np.diff(wavelength) <= 0.0):
        raise ValueError('wavelength must be positive and strictly increasing.')
    if np.any(particle_radius <= 0.0):
        raise ValueError('particle_radius must be positive.')
    if particle_number_density.shape != particle_radius.shape or np.any(particle_number_density < 0.0):
        raise ValueError('particle_number_density must be nonnegative and match particle_radius.')
    if refractive_index is None:
        log_wavelength = np.log(wavelength)
        log_table = np.log(KHARE_WAVELENGTH)
        refractive_index = (np.interp(log_wavelength, log_table, KHARE_REFRACTIVE_INDEX.real)
                            + 1j*np.interp(log_wavelength, log_table, KHARE_REFRACTIVE_INDEX.imag))
    else:
        refractive_index = np.asarray(refractive_index, dtype=np.complex128)
        if (refractive_index.shape != wavelength.shape
                or not np.all(np.isfinite(refractive_index))
                or np.any(refractive_index.real <= 0.0)
                or np.any(refractive_index.imag > 0.0)):
            raise ValueError('refractive_index must match wavelength and have finite n > 0, k >= 0 in n - i*k.')
    shape = (particle_radius.size, wavelength.size)
    dtau_dz = np.zeros(shape)
    w0 = np.zeros(shape)
    g0 = np.zeros(shape)
    for i in range(particle_radius.size):
        if particle_number_density[i] == 0.0:
            continue
        area_density = np.pi*particle_radius[i]**2*particle_number_density[i]
        for j in range(wavelength.size):
            x = 2.0*np.pi*particle_radius[i]/wavelength[j]
            # Explicitly use miepython's compiled kernel, independent of its
            # process-wide MIEPYTHON_USE_JIT setting.
            qext, qsca, _, asymmetry = _single_sphere_nb(refractive_index[j], x, 0, True)
            dtau_dz[i, j] = area_density*qext
            if qext > 0.0:
                w0[i, j] = qsca/qext
            if qsca > 0.0:
                g0[i, j] = asymmetry
    if (not np.all(np.isfinite(dtau_dz)) or not np.all(np.isfinite(w0))
            or not np.all(np.isfinite(g0))):
        raise RuntimeError('Mie calculation produced nonfinite optical properties.')
    return dtau_dz, w0, g0
