"""Smoke test for the installed Clima Python wrapper. See README.md."""

from pathlib import Path

import numpy as np

from photochem import clima


TEST_DIR = Path(__file__).resolve().parent


def test_clima_python_surface():
    """Clima's public arrays, nested views, and rebin utilities are usable."""
    c = clima.AdiabatClimate(
        str(TEST_DIR / "species.yaml"),
        str(TEST_DIR / "settings_adiabat.yaml"),
        str(TEST_DIR / "sun.txt"),
        double_radiative_grid=False,
    )

    partial_pressures = 1.0e6 * np.array(
        [270.0, 400.0e-6, 1.0, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10]
    )
    isr, olr = c.TOA_fluxes(280.0, partial_pressures)
    assert np.isfinite(isr) and isr > 0.0
    assert np.isfinite(olr) and olr > 0.0

    ng = len(c.species_names)
    npart = len(c.particle_names)
    nz = c.P.size
    assert c.f_i_surf.shape == (ng,)
    assert c.P.shape == c.T.shape == c.z.shape == c.dz.shape == (nz,)
    assert c.f_i.shape == c.densities.shape == (nz, ng)
    assert c.pdensities.shape == c.pradii.shape == (nz, npart)
    assert c.N_atmos.shape == c.N_surface.shape == (ng,)
    assert c.N_ocean.shape == (ng, ng)

    pressure_copy = c.P
    pressure_copy[0] = -1.0
    assert c.P[0] > 0.0

    relative_humidity = c.RH
    relative_humidity[0] = 0.75
    assert c.RH[0] != 0.75
    c.RH = relative_humidity
    assert c.RH[0] == 0.75

    c.convective_hysteresis_frac_on = 0.03
    c.convective_hysteresis_frac_off = 0.04
    c.convective_hysteresis_min = 0.002
    c.compute_solar_in_jac = True
    assert c.convective_hysteresis_frac_on == 0.03
    assert c.convective_hysteresis_frac_off == 0.04
    assert c.convective_hysteresis_min == 0.002
    assert c.compute_solar_in_jac
    c.rce_solve_strategy = clima.RCE_SOLVE_HYBRJ_ONLY
    assert c.rce_solve_strategy == clima.RCE_SOLVE_HYBRJ_ONLY
    assert clima.RCE_SOLVE_PTC_THEN_HYBRJ == 2
    assert clima.RCE_SOLVE_HYBRJ_THEN_PTC_THEN_HYBRJ == 3

    rad = c.rad
    assert rad.ir.wavl.shape == rad.ir.freq.shape
    assert rad.sol.wavl.shape == rad.sol.freq.shape
    assert rad.wrk_ir.fup_n.shape == (nz + 1,)
    assert rad.wrk_sol.fdn_n.shape == (nz + 1,)
    assert rad.wrk_ir.fup_a.shape[0] == nz + 1
    assert rad.wrk_ir.tau_band.shape[0] == nz

    original_flux = rad.bolometric_flux()
    rad.set_bolometric_flux(0.9 * original_flux)
    np.testing.assert_allclose(rad.bolometric_flux(), 0.9 * original_flux)
    rad.set_bolometric_flux(original_flux)

    values, errors = clima.rebin_with_errors(
        np.array([0.0, 1.0, 2.0]),
        np.array([2.0, 4.0]),
        np.array([1.0, 2.0]),
        np.array([0.0, 2.0]),
    )
    np.testing.assert_allclose(values, [3.0])
    np.testing.assert_allclose(errors, [np.sqrt(1.25)])


def test_custom_optical_properties():
    """The integer opacity mode reaches Fortran through the Python interface."""
    c = clima.AdiabatClimate(
        str(TEST_DIR / "species.yaml"),
        str(TEST_DIR / "settings_adiabat.yaml"),
        str(TEST_DIR / "sun.txt"),
        double_radiative_grid=False,
    )
    partial_pressures = 1.0e6 * np.array(
        [270.0, 400.0e-6, 1.0, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10]
    )
    baseline = c.TOA_fluxes(280.0, partial_pressures)
    wavelength = np.array([1.0, 1.0e9])
    pressure = np.array([1.0e9, 1.0e-9])
    dtau = np.zeros((2, 2))
    w0 = np.full((2, 2), 0.5)
    g0 = np.full((2, 2), 0.3)
    for mode in (1, 2):
        c.rad.set_custom_optical_properties(wavelength, pressure, dtau, w0, g0, tau_mode=mode)
        np.testing.assert_allclose(c.TOA_fluxes(280.0, partial_pressures), baseline)

    dtau[:] = 1.0e-6
    c.rad.set_custom_optical_properties(wavelength, pressure, dtau, w0, g0, tau_mode=1)
    c.TOA_fluxes(280.0, partial_pressures)
    per_cm_tau = c.rad.wrk_sol.tau_band
    c.rad.set_custom_optical_properties(wavelength, pressure, dtau, w0, g0, tau_mode=2)
    c.TOA_fluxes(280.0, partial_pressures)
    assert np.max(per_cm_tau - c.rad.wrk_sol.tau_band) > 1.0e-3

    # Sliced coordinate arrays must reach Fortran in the same order as copies.
    dtau[:] = [[1.0e-6, 2.0e-6], [3.0e-6, 4.0e-6]]
    c.rad.set_custom_optical_properties(wavelength, pressure, dtau, w0, g0, tau_mode=2)
    contiguous_flux = c.TOA_fluxes(280.0, partial_pressures)
    sliced_wavelength = np.array([1.0, 99.0, 1.0e9])[::2]
    sliced_pressure = np.array([1.0e-9, 99.0, 1.0e9])[::-2]
    c.rad.set_custom_optical_properties(sliced_wavelength, sliced_pressure, dtau, w0, g0, tau_mode=2)
    np.testing.assert_allclose(c.TOA_fluxes(280.0, partial_pressures), contiguous_flux)

    try:
        c.rad.set_custom_optical_properties(wavelength, pressure, dtau, w0, g0, tau_mode=-1)
    except clima.ClimaException:
        pass
    else:
        raise AssertionError("Invalid tau_mode was accepted")
    c.rad.unset_custom_optical_properties()
    np.testing.assert_allclose(c.TOA_fluxes(280.0, partial_pressures), baseline)


def main():
    test_clima_python_surface()
    test_custom_optical_properties()


if __name__ == "__main__":
    main()
