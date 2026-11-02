"""Standalone smoke checks for the supplied-grid haze microphysics solver."""

from unittest.mock import patch
from contextlib import contextmanager

import numpy as np

from photochem.extensions import haze

from photochem.extensions.haze import fall_velocity, mckay_microphysics


@contextmanager
def track_integration():
    """Inspect interval calls using the same driver body without JIT binding."""
    with patch.object(haze, 'mckay_microphysics', haze.mckay_microphysics.py_func):
        with patch.object(haze, '_integrate_log_radius', wraps=haze._integrate_log_radius) as integration:
            yield integration


def atmosphere(altitude):
    ones = np.ones_like(altitude)
    return dict(
        altitude=altitude,
        temperature=250.0*ones,
        air_density=1.0e-6*ones,
        air_number_density=2.0e16*ones,
        gravity=980.0*ones,
        production_rate=1.0e-20*ones,
    )


def test_empty_and_seed_regions():
    inputs = atmosphere(np.array([0.0, 1.0e5, 3.0e5, 6.0e5, 1.0e6]))
    inputs['production_rate'][:] = 0.0
    radius, number, flux = mckay_microphysics(**inputs)
    np.testing.assert_array_equal(radius, np.full(5, 1.0e-7))
    np.testing.assert_array_equal(number, np.zeros(5))
    np.testing.assert_array_equal(flux, np.zeros(5))

    inputs['production_rate'][:3] = [1.0e-20, 2.0e-20, 1.0e-20]
    radius, number, flux = mckay_microphysics(**inputs)
    np.testing.assert_allclose(flux, [6.0e-15, 4.5e-15, 1.5e-15, 0.0, 0.0])
    np.testing.assert_array_equal(radius[2:], np.full(3, 1.0e-7))
    np.testing.assert_array_equal(number[3:], np.zeros(2))
    assert np.all(number[:3] > 0.0)
    assert np.all(radius[:-1] >= radius[1:])

    speed = fall_velocity(
        inputs['gravity'], radius, 1.0, inputs['air_density'],
        inputs['air_number_density'], inputs['temperature']
    )
    mass = (4.0*np.pi/3.0)*radius**3
    np.testing.assert_allclose(number*mass*speed, flux, rtol=1.0e-13)


def test_production_without_coagulation():
    inputs = atmosphere(np.linspace(0.0, 1.0e6, 13))
    # Exponential charge suppression underflows to zero, isolating mass addition.
    radius, number, flux = mckay_microphysics(
        **inputs, charging_radius=1.0e-30, rtol=1.0e-9, atol=1.0e-11
    )
    expected = 1.0e-7*(flux[:-1]/flux[-2])**(1.0/3.0)
    np.testing.assert_allclose(radius[:-1], expected, rtol=2.0e-8)
    speed = fall_velocity(
        inputs['gravity'], radius, 1.0, inputs['air_density'],
        inputs['air_number_density'], inputs['temperature']
    )
    # Without nucleation or coagulation, downward number flux stays constant.
    np.testing.assert_allclose(number[:-1]*speed[:-1], number[-2]*speed[-2], rtol=1.0e-7)


def test_column_normalization():
    inputs = atmosphere(np.linspace(0.0, 1.0e6, 9))
    original = inputs['production_rate'].copy()
    target = 3.0e-14
    result = mckay_microphysics(**inputs, column_production=target)
    np.testing.assert_allclose(result[2][0], target, rtol=1.0e-14)
    np.testing.assert_array_equal(inputs['production_rate'], original)

    inputs['production_rate'] = 7.0*original
    rescaled = mckay_microphysics(**inputs, column_production=target)
    for left, right in zip(result, rescaled):
        np.testing.assert_allclose(left, right, rtol=1.0e-8)
    radius, number, flux = mckay_microphysics(**inputs, column_production=0.0)
    np.testing.assert_array_equal(radius, np.full(9, 1.0e-7))
    np.testing.assert_array_equal(number, np.zeros(9))
    np.testing.assert_array_equal(flux, np.zeros(9))


def test_integration_convergence():
    # Keep the same seed location and the same continuous atmosphere/source
    # while refining the sampling below it; otherwise the boundary also changes.
    zseed, ztop = 9.0e6, 1.0e7
    coarse = atmosphere(np.r_[np.linspace(0.0, zseed, 17), ztop])
    fine = atmosphere(np.r_[np.linspace(0.0, zseed, 33), ztop])
    default = mckay_microphysics(**coarse)
    tight = mckay_microphysics(**coarse, rtol=1.0e-11, atol=1.0e-13)
    refined = mckay_microphysics(**fine, rtol=1.0e-11, atol=1.0e-13)
    # Number density amplifies radius error; allow 0.05% with the default rtol.
    for actual, reference in zip(default, tight):
        np.testing.assert_allclose(actual, reference, rtol=5.0e-4)
    for actual, reference in zip(tight, refined):
        np.testing.assert_allclose(actual[:-1], reference[:-1:2], rtol=1.0e-7)
        np.testing.assert_allclose(actual[-1], reference[-1])


def test_max_radius():
    inputs = atmosphere(np.linspace(0.0, 1.0e6, 13))
    baseline = mckay_microphysics(**inputs)
    cap = 1.0e-6
    with track_integration() as integration:
        radius, number, flux = haze.mckay_microphysics(**inputs, max_radius=cap)
        assert 0 < integration.call_count < inputs['altitude'].size - 2
    assert np.all(radius <= cap)
    assert np.any(radius == cap)
    np.testing.assert_allclose(radius, np.minimum(baseline[0], cap), rtol=1.0e-8)
    np.testing.assert_array_equal(flux, baseline[2])
    speed = fall_velocity(
        inputs['gravity'], radius, 1.0, inputs['air_density'],
        inputs['air_number_density'], inputs['temperature']
    )
    np.testing.assert_allclose(number*(4*np.pi/3)*radius**3*speed, flux, rtol=1.0e-13)

    unchanged = mckay_microphysics(**inputs, max_radius=2*baseline[0][0])
    for actual, expected in zip(unchanged, baseline):
        np.testing.assert_array_equal(actual, expected)
    with track_integration() as integration:
        radius, number, flux = haze.mckay_microphysics(**inputs, max_radius=1.0e-7)
        assert integration.call_count == 0
    np.testing.assert_array_equal(radius, np.full(13, 1.0e-7))
    assert np.all(number[:-1] > 0.0)


def test_sweep_below():
    inputs = atmosphere(np.linspace(0.0, 1.0e6, 9))
    altitude = inputs['altitude']
    for cap in (None, 1.0e-6):
        baseline = mckay_microphysics(**inputs, max_radius=cap, column_production=3.0e-14)
        for threshold in (altitude[3], 0.5*(altitude[3] + altitude[4]),
                          altitude[0], altitude[0] - 1, altitude[-1], altitude[-1] + 1):
            with track_integration() as integration:
                radius, number, flux = haze.mckay_microphysics(
                    **inputs, max_radius=cap, sweep_below=threshold,
                    column_production=3.0e-14
                )
                expected_intervals = max(
                    altitude.size - 2 - np.count_nonzero(altitude < threshold), 0
                )
                if cap is None:
                    assert integration.call_count == expected_intervals
                else:
                    assert integration.call_count <= expected_intervals
            removed = altitude < threshold
            assert np.all(number[removed] == 0.0)
            assert np.all(flux[removed] == 0.0)
            np.testing.assert_allclose(number[~removed], baseline[1][~removed], rtol=1.0e-13)
            np.testing.assert_allclose(flux[~removed], baseline[2][~removed], rtol=1.0e-13)
            np.testing.assert_allclose(radius[~removed], baseline[0][~removed], rtol=1.0e-13)
            if np.any(~removed):
                expected_radius = baseline[0][np.flatnonzero(~removed)[0]]
            else:
                expected_radius = 1.0e-7
            assert np.all(radius[removed] == expected_radius)
            speed = fall_velocity(
                inputs['gravity'], radius, 1.0, inputs['air_density'],
                inputs['air_number_density'], inputs['temperature']
            )
            np.testing.assert_allclose(number*(4*np.pi/3)*radius**3*speed, flux, rtol=1.0e-13)

    # Production can be nonzero, yet lie entirely below the retained region.
    inputs['production_rate'][3:] = 0.0
    with track_integration() as integration:
        radius, number, flux = haze.mckay_microphysics(**inputs, sweep_below=altitude[5])
        assert integration.call_count == 0
    np.testing.assert_array_equal(radius, np.full(9, 1.0e-7))
    np.testing.assert_array_equal(number, np.zeros(9))
    np.testing.assert_array_equal(flux, np.zeros(9))

    inputs['production_rate'][:] = 0.0
    radius, number, flux = mckay_microphysics(
        **inputs, max_radius=1.0e-6, sweep_below=altitude[3]
    )
    np.testing.assert_array_equal(radius, np.full(9, 1.0e-7))
    np.testing.assert_array_equal(number, np.zeros(9))
    np.testing.assert_array_equal(flux, np.zeros(9))


def test_integrator_failures():
    state = np.array([250.0, 1.0e-6, 2.0e16, 980.0, 1.0e-20])
    kwargs = dict(
        log_radius=np.log(1.0e-7), thickness=1.0e5,
        state_above=state, state_below=state, flux_above=1.0e-15,
        particle_density=1.0, charging_radius=9.0e-6, rtol=1.0e-6, atol=1.0e-10,
    )
    _, capped, status = haze._integrate_log_radius(**kwargs, max_steps=1)
    assert not capped and status == 3
    _, capped, status = haze._integrate_log_radius(**dict(kwargs, flux_above=0.0))
    assert not capped and status == 1
    _, capped, status = haze._integrate_log_radius(
        **dict(kwargs, thickness=np.nextafter(0.0, 1.0))
    )
    assert not capped and status == 2


def test_nonfinite_trial_recovery():
    # A coarse interval on a steep source tail generates nonfinite RK trial
    # states. The integrator must reject those trials rather than abort.
    above = np.array([278.177727, 7.41031595e-12, 1.54381582e11, 945.968033, 3.99761877e-34])
    below = np.array([227.568339, 2.75189088e-11, 5.73310600e11, 948.300025, 1.05086645e-28])
    initial = -17.12841732
    thickness = 797812.5744476914
    args = (thickness, above, below, 9.796764923311649e-29, 1.2876324544508917, 9.0e-6)
    result, capped, status = haze._integrate_log_radius(initial, *args, 1.0e-6, 1.0e-10)
    assert status == 0 and not capped

    assert np.isfinite(result) and result >= initial


def test_invalid_inputs():
    inputs = atmosphere(np.linspace(0.0, 1.0e6, 5))
    cases = [
        {'altitude': inputs['altitude'][::-1]},
        {'temperature': np.ones(4)},
        {'temperature': np.ones((5, 1))},
        {'air_number_density': np.zeros(5)},
        {'production_rate': -np.ones(5)},
        {'production_rate': np.full(5, np.nan)},
        {'seed_radius': 0.0},
        {'charging_radius': 0.0},
        {'max_radius': 1.0e-8},
        {'max_radius': np.inf},
        {'max_radius': np.nan},
        {'sweep_below': np.nan},
        {'sweep_below': np.inf},
        {'particle_density': 1.0e-7},
        {'rtol': 0.0},
        {'atol': np.inf},
        {'column_production': -1.0},
        {'column_production': 1.0, 'production_rate': np.zeros(5)},
    ]
    for updates in cases:
        try:
            mckay_microphysics(**dict(inputs, **updates))
        except ValueError:
            pass
        else:
            raise AssertionError(f'Invalid inputs accepted: {updates}')


if __name__ == '__main__':
    test_empty_and_seed_regions()
    test_production_without_coagulation()
    test_column_normalization()
    test_integration_convergence()
    test_max_radius()
    test_sweep_below()
    test_integrator_failures()
    test_nonfinite_trial_recovery()
    test_invalid_inputs()
    print('Haze microphysics smoke checks passed.')
