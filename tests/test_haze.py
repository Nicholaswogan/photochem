"""Fast standalone checks for primary-particle haze production and transport."""

import numpy as np
from photochem.extensions import haze
from photochem.extensions.haze import fall_velocity, mckay_microphysics


def atmosphere(altitude):
    ones = np.ones_like(altitude)
    return dict(altitude=altitude, temperature=250.0*ones,
                air_density=1e-6*ones, air_number_density=2e16*ones,
                gravity=980.0*ones, production_rate=1e-20*ones)


def fluxes(inputs, result, seed_radius=1e-7):
    radius, number, mass_flux = result
    speed = fall_velocity(inputs['gravity'], radius, 1.0, inputs['air_density'],
                          inputs['air_number_density'], inputs['temperature'])
    np.testing.assert_allclose(number*(4*np.pi/3)*radius**3*speed, mass_flux,
                               rtol=1e-12, atol=0)
    assert mass_flux[-1] == number[-1] == 0.0
    assert np.all(radius >= seed_radius)
    return number*speed


def test_empty_and_mass_budget():
    inputs = atmosphere(np.array([0., 1e5, 3e5, 6e5, 1e6]))
    inputs['production_rate'][:] = 0
    r, n, C = mckay_microphysics(**inputs)
    np.testing.assert_array_equal(r, np.full(5, 1e-7))
    assert np.all(n == 0) and np.all(C == 0)
    inputs['production_rate'][:3] = [1e-20, 2e-20, 1e-20]
    result = mckay_microphysics(**inputs)
    np.testing.assert_allclose(result[2], [6e-15, 4.5e-15, 1.5e-15, 0, 0], rtol=1e-14)
    assert np.all(result[1][3:] == 0)
    fluxes(inputs, result)
    original = inputs['production_rate'].copy()
    scaled = mckay_microphysics(**inputs, column_production=3e-14)
    np.testing.assert_allclose(scaled[2][0], 3e-14, rtol=1e-14)
    np.testing.assert_array_equal(original, inputs['production_rate'])
    r, n, C = mckay_microphysics(**inputs, column_production=0.)
    assert np.all(n == 0) and np.all(C == 0)


def test_production_without_coagulation():
    inputs = atmosphere(np.linspace(0, 1e6, 9))
    for method in ('RK45', 'Radau', 'BDF', 'LSODA'):
        result = mckay_microphysics(**inputs, charging_radius=1e-30, method=method)
        # Each mass increment introduces exactly q/m0 new primary particles.
        np.testing.assert_allclose(result[0], 1e-7, rtol=1e-12)
        F = fluxes(inputs, result)
        np.testing.assert_allclose(F, result[2]/((4*np.pi/3)*1e-21), rtol=1e-12)


def test_coagulation_and_refinement():
    inputs = atmosphere(np.linspace(0, 1e6, 9))
    result = mckay_microphysics(**inputs)
    F = fluxes(inputs, result)
    assert np.all(result[0][:-1] > 1e-7)
    assert np.all(F[:-1] < result[2][:-1]/((4*np.pi/3)*1e-21))
    fine = atmosphere(np.linspace(0, 1e6, 17))
    refined = mckay_microphysics(**fine, rtol=1e-8)
    for a, b in zip(result, refined):
        np.testing.assert_allclose(a, b[::2], rtol=2e-3, atol=0)
    # A source-free lower region can only lose number flux to coagulation.
    inputs['production_rate'][:5] = 0
    result = mckay_microphysics(**inputs)
    F = fluxes(inputs, result)
    assert np.all(np.diff(F[:4]) >= 0)


def test_cap_and_sweep():
    inputs = atmosphere(np.linspace(0, 1e6, 9))
    baseline = mckay_microphysics(**inputs)
    cap = 1.1e-7
    capped = mckay_microphysics(**inputs, max_radius=cap)
    assert np.all(capped[0] <= cap)
    assert np.any(capped[0] == cap)
    fluxes(inputs, capped)
    np.testing.assert_array_equal(capped[2], baseline[2])
    primary = mckay_microphysics(**inputs, max_radius=1e-7)
    np.testing.assert_allclose(primary[0], 1e-7, rtol=1e-12)
    fluxes(inputs, primary)
    # A source taper followed by a thin source pulse exercises both reaching
    # the cap and leaving it when freshly produced primary particles arrive.
    pulse = atmosphere(np.r_[0., .01, np.linspace(1e5, 1e6, 10)])
    pulse['production_rate'][:] = 0
    pulse['production_rate'][-4:] = 1e-21
    pulse['production_rate'][0] = 1e-12
    released = mckay_microphysics(**pulse, max_radius=5e-7)
    assert released[0][1] == 5e-7
    assert released[0][0] < released[0][1]
    fluxes(pulse, released)
    for threshold in (0., 3.5e5, 1e6, 1.1e6):
        result = mckay_microphysics(**inputs, sweep_below=threshold)
        removed = inputs['altitude'] < threshold
        assert np.all(result[1][removed] == 0) and np.all(result[2][removed] == 0)
        for a, b in zip(result, baseline):
            np.testing.assert_allclose(a[~removed], b[~removed], rtol=1e-12, atol=0)
        expected = result[0][np.flatnonzero(~removed)[0]] if np.any(~removed) else 1e-7
        assert np.all(result[0][removed] == expected)


def test_integration_failures():
    from types import SimpleNamespace
    from unittest.mock import patch
    inputs = atmosphere(np.linspace(0, 1e6, 5))
    with patch.object(haze, 'solve_ivp', return_value=SimpleNamespace(
            success=False, message='forced integration failure')):
        try:
            mckay_microphysics(**inputs)
        except RuntimeError as error:
            assert 'forced integration failure' in str(error)
        else:
            raise AssertionError('Solver failure was ignored.')


def test_invalid_inputs():
    inputs = atmosphere(np.linspace(0, 1e6, 5))
    for updates in ({'altitude': inputs['altitude'][::-1]}, {'temperature': np.ones(4)},
                    {'temperature': np.ones((5, 1))}, {'air_number_density': np.zeros(5)},
                    {'production_rate': -np.ones(5)}, {'production_rate': np.full(5, np.nan)},
                    {'seed_radius': 0.}, {'charging_radius': 0.}, {'max_radius': 1e-8},
                    {'max_radius': np.inf}, {'sweep_below': np.nan}, {'particle_density': 1e-7},
                    {'rtol': 0.}, {'atol': np.inf}, {'method': 'unknown'}, {'column_production': -1.},
                    {'column_production': 1., 'production_rate': np.zeros(5)}):
        try:
            mckay_microphysics(**dict(inputs, **updates))
        except ValueError:
            pass
        else:
            raise AssertionError(f'Invalid inputs accepted: {updates}')


if __name__ == '__main__':
    test_empty_and_mass_budget()
    test_production_without_coagulation()
    test_coagulation_and_refinement()
    test_cap_and_sweep()
    test_integration_failures()
    test_invalid_inputs()
    print('Haze microphysics smoke checks passed.')
