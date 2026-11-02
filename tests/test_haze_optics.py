"""Smoke checks for compact-particle haze optical properties."""

import numpy as np
from miepython import efficiencies_mx
from photochem.extensions.haze_optics import (
    haze_optics, KHARE_WAVELENGTH, KHARE_REFRACTIVE_INDEX, _haze_optics,
)


def test_haze_optics():
    wavelength = np.array([2e-5, 5e-5, 1e-4])
    radius = np.array([1e-7, 1e-5, 1e-4])
    number = np.array([0.0, 10.0, 100.0])
    result = haze_optics(wavelength, radius, number)
    doubled = haze_optics(wavelength, radius, 2*number)
    for array in result:
        assert array.shape == (3, 3)
        assert np.all(np.isfinite(array))
        assert np.all(array[0] == 0.0)
    np.testing.assert_allclose(doubled[0], 2*result[0])
    np.testing.assert_array_equal(doubled[1], result[1])
    np.testing.assert_array_equal(doubled[2], result[2])
    assert np.all((result[1] >= 0.0) & (result[1] <= 1.0))
    assert np.all(np.abs(result[2]) <= 1.0)
    assert _haze_optics.nopython_signatures

    # Direct efficiencies and a Rayleigh-limit check independently constrain
    # radius/wavelength units, geometric area, and refractive-index convention.
    index = np.full(3, 1.5 - 0.1j)
    ext, albedo, asymmetry = haze_optics(wavelength, radius, number, index)
    qext, qsca, _, g = efficiencies_mx(index, 2*np.pi*radius[1]/wavelength)
    np.testing.assert_allclose(ext[1], number[1]*np.pi*radius[1]**2*qext)
    np.testing.assert_allclose(albedo[1], qsca/qext)
    np.testing.assert_allclose(asymmetry[1], g)
    ext, albedo, _ = haze_optics(wavelength, np.array([1e-8]), np.array([1.0]),
                                np.full(3, 1.5+0j))
    x = 2*np.pi*1e-8/wavelength
    rayleigh = (8/3)*x**4*((1.5**2-1)/(1.5**2+2))**2
    np.testing.assert_allclose(ext[0]/(np.pi*1e-16), rayleigh, rtol=1e-5)
    np.testing.assert_allclose(albedo, 1.0)
    for array in haze_optics(wavelength, radius, number, np.ones(3, dtype=complex)):
        np.testing.assert_array_equal(array, 0.0)

    # Default interpolation reproduces table values at the original knots.
    defaults = haze_optics(KHARE_WAVELENGTH, radius[1:2], number[1:2])
    explicit = haze_optics(KHARE_WAVELENGTH, radius[1:2], number[1:2],
                           KHARE_REFRACTIVE_INDEX)
    for actual, expected in zip(defaults, explicit):
        np.testing.assert_allclose(actual, expected)

    # Outside the table, hold the refractive index at the nearest endpoint.
    outside = np.array([KHARE_WAVELENGTH[0]/2, KHARE_WAVELENGTH[-1]*2])
    defaults = haze_optics(outside, radius[1:2], number[1:2])
    explicit = haze_optics(outside, radius[1:2], number[1:2],
                           KHARE_REFRACTIVE_INDEX[[0, -1]])
    for actual, expected in zip(defaults, explicit):
        np.testing.assert_allclose(actual, expected)


def test_invalid_inputs():
    args = dict(wavelength=np.array([5e-5]), particle_radius=np.array([1e-5]),
                particle_number_density=np.array([1.0]))
    for update in (
        dict(wavelength=[0.0]), dict(wavelength=[5e-5, 4e-5]),
        dict(particle_radius=[-1.0]),
        dict(particle_number_density=[-1.0]), dict(particle_number_density=[1., 2.]),
        dict(refractive_index=[1.5+0.1j]), dict(refractive_index=[complex('nan')]),
        dict(refractive_index=[1.5+0j, 1.5+0j]),
    ):
        try:
            haze_optics(**dict(args, **update))
        except ValueError:
            pass
        else:
            raise AssertionError(f'Invalid inputs accepted: {update}')


if __name__ == '__main__':
    test_haze_optics()
    test_invalid_inputs()
    print('Haze optics smoke checks passed.')
