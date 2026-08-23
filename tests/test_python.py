"""Smoke test for the installed Python wrappers. See README.md."""

import tempfile
from pathlib import Path

import numpy as np

from photochem import EvoAtmosphere, PhotoException, zahnle_earth


def _make_file_atmosphere():
    pc = EvoAtmosphere(
        zahnle_earth,
        "../examples/ModernEarth/settings.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        "../examples/ModernEarth/atmosphere.txt",
        data_dir="../data",
    )
    pc.var.verbose = 0
    return pc


def test_evolve_uses_fixed_grid():
    """Evolution does not move the grid and rejects mismatched restarts."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        filename = Path(tmp_dir) / "evolution.dat"
        pc = _make_file_atmosphere()
        initial_usol = pc.wrk.usol.copy()

        assert pc.evolve(
            str(filename),
            0.0,
            initial_usol,
            np.array([1.0, 2.0]),
            overwrite=True,
        )
        assert pc.var.top_atmos == 1.0e7

        mismatched = _make_file_atmosphere()
        mismatched.update_vertical_grid(TOA_alt=0.9 * mismatched.var.top_atmos)
        try:
            mismatched.evolve(
                str(filename),
                0.0,
                mismatched.wrk.usol.copy(),
                np.array([3.0]),
                restart_from_file=True,
            )
        except PhotoException as exc:
            assert "fixed grid" in str(exc)
        else:
            raise AssertionError("restart accepted a mismatched fixed grid")


def test_static_construction():
    pc = EvoAtmosphere(
        "no_particle_test.yaml",
        "test_settings_minimal.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        data_dir="../data",
    )

    assert not pc.atmosphere_initialized
    assert pc.dat.nq > 0
    assert pc.var.nz > 0

    try:
        pc.gas_fluxes()
    except PhotoException as exc:
        assert "atmosphere is not initialized" in str(exc)
    else:
        raise AssertionError("gas_fluxes accepted an uninitialized atmosphere")

    z = np.array([0.0, 5.0e6, 1.0e7])
    temperature = np.array([300.0, 240.0, 180.0])
    edd = np.array([1.0e5, 1.0e6, 1.0e7])
    mix = {"H2": np.ones(z.size)}
    pc.initialize_atmosphere_z(z, temperature, edd, 1.0e6, mix)

    assert pc.atmosphere_initialized


def test_toa_pressure_maintenance_api():
    """The optional TOA-maintenance settings and counters are Python-visible."""
    pc = EvoAtmosphere(
        "no_particle_test.yaml",
        "test_settings_minimal.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        data_dir="../data",
    )

    maintenance = pc.var.toa_pressure_maintenance
    assert maintenance.enabled is False
    assert maintenance.pressure_factor == 3.0

    maintenance.enabled = True
    maintenance.target_pressure = 2.5e-4
    maintenance.pressure_factor = 4.0
    maintenance.nsteps_between_updates = 7
    maintenance.max_failures = 2

    assert maintenance.enabled is True
    assert maintenance.target_pressure == 2.5e-4
    assert maintenance.pressure_factor == 4.0
    assert maintenance.nsteps_between_updates == 7
    assert maintenance.max_failures == 2

    # Counters are available before a robust integration starts and have no
    # side effects from merely configuring the mode.
    assert pc.wrk.n_toa_pressure_updates == 0
    assert pc.wrk.n_toa_pressure_failures == 0
    assert pc.wrk.nsteps_since_toa_pressure_update == 0


def test_persistent_profile_controls_toa_maintenance():
    """Persistent-profile options reach the shared Fortran maintenance state."""
    pc = EvoAtmosphere(
        "no_particle_test.yaml",
        "test_settings_minimal.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        data_dir="../data",
    )
    z = np.array([0.0, 5.0e6, 1.0e7])
    pc.initialize_atmosphere_z(
        z,
        np.array([300.0, 240.0, 180.0]),
        np.array([1.0e5, 1.0e6, 1.0e7]),
        1.0e6,
        {"H2": np.ones(z.size)},
    )

    pressure = np.array([
        2.0 * pc.wrk.surface_pressure * 1.0e6,
        0.5 * pc.wrk.pressure_hydro[-1],
    ])
    temperature = np.array([300.0, 180.0])
    edd = np.array([3.0e7, 4.0e5])

    pc.set_press_temp_edd_profile(pressure, temperature, edd)
    assert pc.var.toa_pressure_maintenance.enabled is True
    assert pc.var.toa_pressure_maintenance.target_pressure == 0.1

    explicit_target = 2.5e-4
    pc.set_press_temp_edd_profile(
        pressure, temperature, edd, target_pressure=explicit_target
    )
    assert pc.var.toa_pressure_maintenance.enabled is True
    assert pc.var.toa_pressure_maintenance.target_pressure == explicit_target

    pc.set_press_temp_edd_profile(
        pressure, temperature, edd, maintain_toa_pressure=False
    )
    assert pc.var.toa_pressure_maintenance.enabled is False
    pc.clear_press_temp_edd_profile()


def test_robust_initial_toa_pressure_preflight():
    """Robust initialization corrects an out-of-band TOA before CVODE starts."""
    pc = EvoAtmosphere(
        "no_particle_test.yaml",
        "test_settings_minimal.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        data_dir="../data",
    )
    z = np.array([0.0, 5.0e6, 1.0e7])
    pc.initialize_atmosphere_z(
        z,
        np.array([300.0, 240.0, 180.0]),
        np.array([1.0e5, 1.0e6, 1.0e7]),
        1.0e6,
        {"H2": np.ones(z.size)},
    )
    pressure = np.array([
        2.0 * pc.wrk.surface_pressure * 1.0e6,
        0.5 * pc.wrk.pressure_hydro[-1],
    ])
    pc.set_press_temp_edd_profile(
        pressure,
        np.array([300.0, 180.0]),
        np.array([3.0e7, 4.0e5]),
    )

    maintenance = pc.var.toa_pressure_maintenance
    maintenance.target_pressure = 0.95 * pc.wrk.pressure[-1]
    maintenance.pressure_factor = 1.01
    maintenance.nsteps_between_updates = 100
    pc.initialize_robust_stepper(pc.wrk.usol)
    assert pc.wrk.nsteps_total == 0
    assert pc.wrk.n_toa_pressure_updates == 1
    assert pc.wrk.nsteps_since_toa_pressure_update == 0
    assert np.isclose(
        pc.wrk.pressure[-1] / maintenance.target_pressure, 1.0, rtol=0.0, atol=2.0e-5
    )
    pc.destroy_stepper()


def test_gas_giant_static_construction():
    from photochem.extensions.gasgiants import EvoAtmosphereGasGiant

    pc = EvoAtmosphereGasGiant(
        "no_particle_test.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        5.972e27,
        6.371e8,
        thermo_file="no_particle_test.yaml",
        data_dir="../data",
    )

    assert not pc.atmosphere_initialized

    pressure = np.array([1.0e6, 1.0e5, 1.0e4, 1.0e2])
    temperature = np.array([300.0, 260.0, 220.0, 180.0])
    edd = np.array([1.0e5, 3.0e5, 1.0e6, 1.0e7])
    mix = {"H2": np.ones(pressure.size)}
    pc.initialize_atmosphere_p(
        pressure, temperature, edd, mix, persistent=True
    )

    assert pc.atmosphere_initialized


def _make_initialized_gas_giant():
    from photochem.extensions.gasgiants import EvoAtmosphereGasGiant

    pc = EvoAtmosphereGasGiant(
        "no_particle_test.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        5.972e27,
        6.371e8,
        nz=20,
        thermo_file="no_particle_test.yaml",
        data_dir="../data",
    )
    pc.gdat.verbose = False

    # The climate grid extends above the 0.1 dyn/cm^2 photochemical target.
    pressure = np.array([
        1.0e7, 1.0e6, 1.0e4, 1.0e2, 1.0, 1.0e-2, 1.0e-3
    ])
    temperature = np.array([
        320.0, 300.0, 260.0, 220.0, 180.0, 150.0, 140.0
    ])
    edd = np.array([
        5.0e4, 1.0e5, 3.0e5, 1.0e6, 1.0e7, 3.0e7, 5.0e7
    ])
    mix = {"H2": np.ones(pressure.size)}

    pc.gdat.P_clima_grid = pressure.copy()
    pc.gdat.T_clima_grid = temperature.copy()
    pc.gdat.Kzz_clima_grid = edd.copy()
    pc.gdat.P_desired = pressure.copy()
    pc.gdat.T_desired = temperature.copy()
    pc.gdat.Kzz_desired = edd.copy()
    pc.gdat.metallicity = 1.0
    pc.gdat.CtoO = 1.0
    pc.gdat.ind_b = 0
    pc._initialize_atmosphere(
        pressure, temperature, edd, np.zeros(pressure.size), mix
    )
    return pc


def test_gas_giant_returns_climate_levels_above_photochemical_grid():
    pc = _make_initialized_gas_giant()
    sol = pc.return_atmosphere_climate_grid()
    above = pc.gdat.P_clima_grid < pc.wrk.pressure_hydro[-1]

    assert np.count_nonzero(above) >= 2
    assert np.array_equal(sol['temperature'], pc.gdat.T_clima_grid)
    assert np.array_equal(sol['Kzz'], pc.gdat.Kzz_clima_grid)

    species_names = pc.dat.species_names[:(-2-pc.dat.nsl)]
    h2_index = species_names.index('H2')
    h2_at_photochemical_top = (
        pc.wrk.usol[h2_index, -1] / pc.wrk.density[-1]
    )
    assert np.allclose(sol['H2'][above], h2_at_photochemical_top)


def test_gas_giant_reinitializes_from_public_climate_grid():
    pc = _make_initialized_gas_giant()
    pressure = pc.gdat.P_clima_grid.copy()
    initial = pc.return_atmosphere_climate_grid()
    gas_names = pc.dat.species_names[pc.dat.np:(-2-pc.dat.nsl)]
    mix = {sp: initial[sp].copy() for sp in gas_names}
    temperature = pc.gdat.T_clima_grid + np.linspace(
        15.0, 5.0, pressure.size
    )
    edd = 2.0 * pc.gdat.Kzz_clima_grid

    pc.reinitialize_to_new_climate_PT(
        pressure, temperature, edd, mix
    )
    sol = pc.return_atmosphere_climate_grid()
    above = pressure < pc.wrk.pressure_hydro[-1]

    assert np.count_nonzero(above) >= 2
    assert np.array_equal(pc.gdat.T_clima_grid, temperature)
    assert np.array_equal(pc.gdat.Kzz_clima_grid, edd)
    assert np.array_equal(sol['temperature'], temperature)
    assert np.array_equal(sol['Kzz'], edd)
    for sp in gas_names:
        assert np.all(np.isfinite(sol[sp]))
        assert np.all(sol[sp] >= 0.0)
        assert np.allclose(sol[sp][above], sol[sp][above][0])


def test_gas_giant_uses_shared_robust_stepper():
    from photochem.extensions.gasgiants import EvoAtmosphereGasGiant

    # Initialization and the steady-state loop are inherited directly. The
    # remaining robust_step override only adds progress presentation.
    assert "initialize_robust_stepper" not in EvoAtmosphereGasGiant.__dict__
    assert "find_steady_state" not in EvoAtmosphereGasGiant.__dict__

    pc = _make_initialized_gas_giant()
    maintenance = pc.var.toa_pressure_maintenance
    assert pc.var.nerrors_before_giveup == 10
    assert pc.var.nsteps_before_conv_check == 300
    assert pc.var.nsteps_before_reinit == 1000
    assert pc.var.nsteps_before_giveup == 100_000
    assert maintenance.enabled
    assert np.isclose(maintenance.target_pressure, 0.1)
    assert maintenance.pressure_factor == 3.0
    assert maintenance.nsteps_between_updates == 1000

    pc.initialize_robust_stepper(pc.wrk.usol)
    assert pc.wrk.robust_stepper_initialized
    assert pc.wrk.nsteps_total == 0
    assert pc.wrk.nerrors_total == 0

    # Exercise the progress wrapper over a shared step with an immediate,
    # deterministic convergence condition.
    pc.var.equilibrium_time = 0.0
    give_up, converged = pc.robust_step()
    assert not give_up
    assert converged
    assert pc.wrk.nsteps_total == 1

    pc.destroy_stepper()
    assert not pc.wrk.robust_stepper_initialized

    # The inherited steady-state loop dispatches through the shared robust
    # initializer and the presentation-only step wrapper.
    assert pc.find_steady_state()
    assert pc.wrk.robust_stepper_initialized
    pc.destroy_stepper()


def test_gas_giant_shared_limits_and_state_restore():
    pc = _make_initialized_gas_giant()
    custom_target = 0.2
    pc.destroy_stepper()
    pc.set_press_temp_edd_profile(
        pc.gdat.P_desired, pc.gdat.T_desired, pc.gdat.Kzz_desired,
        hydro_pressure=True, target_pressure=custom_target
    )

    # The exact shared accepted-step ceiling replaces the legacy Python
    # counter, which allowed one extra step.
    pc.var.equilibrium_time = 1.0e100
    pc.var.nsteps_before_conv_check = 0
    pc.var.nsteps_before_reinit = 2
    pc.var.nsteps_before_giveup = 1
    pc.var.conv_longdy = -1.0
    pc.var.conv_longdydt = -1.0
    pc.initialize_robust_stepper(pc.wrk.usol)
    give_up, converged = pc.robust_step()
    assert give_up
    assert not converged
    assert pc.wrk.nsteps_total == 1

    state = pc.model_state_to_dict()
    pc.initialize_from_dict(state)
    assert not pc.wrk.robust_stepper_initialized
    assert np.array_equal(pc.gdat.T_clima_grid, state['T_clima_grid'])
    assert np.array_equal(pc.gdat.Kzz_clima_grid, state['Kzz_clima_grid'])

    legacy_state = state.copy()
    legacy_state.pop('T_clima_grid')
    legacy_state.pop('Kzz_clima_grid')
    pc.initialize_from_dict(legacy_state)
    nclima = pc.gdat.P_clima_grid.size
    assert np.array_equal(
        pc.gdat.T_clima_grid, state['T_desired'][:nclima]
    )
    assert np.array_equal(
        pc.gdat.Kzz_clima_grid, state['Kzz_desired'][:nclima]
    )
    maintenance = pc.var.toa_pressure_maintenance
    assert maintenance.enabled
    assert np.isclose(maintenance.target_pressure, custom_target)
    assert maintenance.pressure_factor == 3.0
    assert maintenance.nsteps_between_updates == 1000


def test_initialize_atmosphere_z_no_particles():
    pc = EvoAtmosphere(
        "no_particle_test.yaml",
        "test_settings_minimal.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        "../examples/ModernEarth/atmosphere.txt",
        data_dir="../data",
    )
    z = np.array([0.0, 5.0e6, 1.0e7])
    temperature = np.array([300.0, 240.0, 180.0])
    edd = np.array([1.0e5, 1.0e6, 1.0e7])
    names = pc.dat.species_names[:pc.dat.nq]
    usol = pc.wrk.usol
    dominant_gas = np.argmax(usol[:, 0])
    mix = {names[dominant_gas]: np.ones(z.size)}
    pc.initialize_atmosphere_z(z, temperature, edd, 1.0e6, mix)
    assert pc.dat.np == 0
    assert pc.var.top_atmos == z[-1]
    assert np.all(pc.wrk.pressure_hydro > 0.0)


def test_initialize_atmosphere_z_particles():
    pc = EvoAtmosphere(
        zahnle_earth,
        "../examples/ModernEarth/settings.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        "../examples/ModernEarth/atmosphere.txt",
        data_dir="../data",
    )
    z = np.array([0.0, 5.0e6, 1.0e7])
    temperature = np.array([300.0, 240.0, 180.0])
    edd = np.array([1.0e5, 1.0e6, 1.0e7])
    species_names = pc.dat.species_names[:pc.dat.nq]
    particle_names = species_names[:pc.dat.np]
    dominant_gas = pc.dat.np + np.argmax(
        pc.wrk.usol[pc.dat.np:pc.dat.nq, 0]
    )
    mix = {species_names[dominant_gas]: np.ones(z.size)}
    radii = {particle_names[0]: np.full(z.size, 2.0e-5)}

    pc.initialize_atmosphere_z(
        z, temperature, edd, 1.0e6, mix, particle_radius=radii
    )

    assert pc.var.top_atmos == z[-1]
    assert np.all(np.isfinite(pc.wrk.usol))
    assert np.all(pc.wrk.pressure_hydro > 0.0)
    assert np.all(np.isfinite(pc.wrk.density))
    assert np.all(pc.wrk.density > 0.0)
    assert set(mix).isdisjoint(particle_names)
    assert np.allclose(pc.var.particle_radius[0], 2.0e-5)
    assert np.allclose(pc.var.particle_radius[1:], 1.0e-5)


def test_initialize_atmosphere_p_no_particles():
    pc = EvoAtmosphere(
        "no_particle_test.yaml",
        "test_settings_minimal.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        "../examples/ModernEarth/atmosphere.txt",
        data_dir="../data",
    )
    pressure = np.array([1.0e6, 1.0e5, 1.0e4, 1.0e2])
    temperature = np.array([300.0, 260.0, 220.0, 180.0])
    edd = np.array([1.0e5, 3.0e5, 1.0e6, 1.0e7])
    names = pc.dat.species_names[:pc.dat.nq]
    dominant_gas = np.argmax(pc.wrk.usol[:, 0])
    mix = {names[dominant_gas]: np.ones(pressure.size)}

    pc.initialize_atmosphere_p(
        pressure, temperature, edd, mix, persistent=True
    )

    assert pc.var.bottom_atmos == 0.0
    assert pc.var.top_atmos > 0.0
    assert np.all(np.diff(pc.wrk.pressure_hydro) < 0.0)
    assert pc.var.toa_pressure_maintenance.enabled is True
    assert pc.var.toa_pressure_maintenance.target_pressure == 0.1

    pc.initialize_atmosphere_p(
        pressure, temperature, edd, mix,
        persistent=True, maintain_toa_pressure=False
    )
    assert pc.var.toa_pressure_maintenance.enabled is False

    target_pressure = 0.25
    pc.initialize_atmosphere_p(
        pressure, temperature, edd, mix,
        persistent=True, target_pressure=target_pressure
    )
    assert pc.var.toa_pressure_maintenance.enabled is True
    assert pc.var.toa_pressure_maintenance.target_pressure == target_pressure

    try:
        pc.initialize_atmosphere_p(
            pressure, temperature, edd, mix, target_pressure=target_pressure
        )
    except PhotoException as exc:
        assert "persistent" in str(exc)
    else:
        raise AssertionError("nonpersistent pressure initialization accepted a TOA target")

    pc.clear_press_temp_edd_profile()


def test_initialize_atmosphere_p_particles():
    pc = EvoAtmosphere(
        zahnle_earth,
        "../examples/ModernEarth/settings.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        "../examples/ModernEarth/atmosphere.txt",
        data_dir="../data",
    )
    pressure = np.array([1.0e6, 1.0e5, 1.0e4, 1.0e2])
    temperature = np.array([300.0, 260.0, 220.0, 180.0])
    edd = np.array([1.0e5, 3.0e5, 1.0e6, 1.0e7])
    species_names = pc.dat.species_names[:pc.dat.nq]
    particle_names = species_names[:pc.dat.np]
    dominant_gas = pc.dat.np + np.argmax(
        pc.wrk.usol[pc.dat.np:pc.dat.nq, 0]
    )
    mix = {species_names[dominant_gas]: np.ones(pressure.size)}
    radii = {particle_names[0]: np.full(pressure.size, 2.0e-5)}

    pc.initialize_atmosphere_p(
        pressure, temperature, edd, mix, particle_radius=radii,
        persistent=True, tropopause_pressure=1.0e5
    )

    assert pc.var.top_atmos > 0.0
    assert np.all(np.isfinite(pc.wrk.usol))
    assert np.allclose(pc.var.particle_radius[0], 2.0e-5)
    assert np.allclose(pc.var.particle_radius[1:], 1.0e-5)
    pc.clear_press_temp_edd_profile()


def test_wrapper():

    pc = EvoAtmosphere(
        zahnle_earth,
        "../examples/ModernEarth/settings.yaml",
        "../examples/ModernEarth/Sun_now.txt",
        "../examples/ModernEarth/atmosphere.txt",
        data_dir="../data"
    )
    
    print(pc.dat.nq)
    print(pc.dat.np)
    print(pc.dat.ng)
    print(pc.dat.nsl)
    print(pc.dat.nll)
    print(pc.dat.nsp)
    print(pc.dat.nw)
    print(pc.dat.planet_mass)
    print(pc.dat.planet_radius)
    print(pc.dat.species_names[0])
    print(pc.dat.atoms_names[0])
    print(pc.dat.photonums[0])
    print(pc.dat.wavl[0])
    print(pc.dat.species_mass[0])
    print(pc.dat.species_redox[0])
    print()
    print(pc.var.nz)
    print(pc.var.top_atmos)
    print(pc.var.bottom_atmos)
    print(pc.var.trop_alt)
    print(pc.var.trop_ind)
    print(pc.var.z[0])
    # print(pc.var.photon_flux_fcn)
    pc.var.photon_flux_fcn = None
    print(pc.var.cond_params[0].k_cond)
    print(pc.var.temperature[0])
    print(pc.var.edd[0])
    pc.var.custom_binary_diffusion_fcn = None
    print(pc.var.photon_flux[0])
    print(pc.var.grav[0])
    print(pc.wrk.surface_pressure)
    print(pc.var.max_error_reinit_attempts)
    print(pc.var.rtol)
    print(pc.var.atol)
    print(pc.var.mxsteps)
    print(pc.var.equilibrium_time)
    print(pc.var.conv_hist_factor)
    print(pc.var.conv_min_mix)
    print(pc.var.conv_longdy)
    print(pc.var.conv_longdydt)
    assert pc.var.jacobian_method == 1
    pc.var.jacobian_method = 2
    assert pc.var.jacobian_method == 2
    pc.var.jacobian_method = 1
    print(pc.var.jacobian_method)
    print(pc.var.epsj)
    print(pc.var.verbose)
    print()
    print(pc.wrk.nsteps)
    print(pc.wrk.longdy)
    print(pc.wrk.longdydt)
    print(pc.wrk.tn)
    print(pc.wrk.usol[0,0])
    print(pc.wrk.pressure[0])
    print(pc.wrk.density[0])
    print(pc.wrk.densities[0,0])
    print(pc.wrk.rx_rates[0,0])
    print(pc.wrk.mubar[0])
    print(pc.wrk.prates[0,0])
    print(pc.wrk.amean_grd[0,0])
    print(pc.wrk.optical_depth[0,0])
    print(pc.wrk.surf_radiance[0])

    pc.prep_atmosphere(pc.wrk.usol)
    with tempfile.TemporaryDirectory() as tmp_dir:
        atmosphere_path = Path(tmp_dir) / "atmosphere.txt"
        pc.out2atmosphere_txt(str(atmosphere_path), overwrite=True)
    pc.gas_fluxes()
    pc.set_lower_bc('O',bc_type='vdep',vdep=0)
    pc.set_upper_bc('O',bc_type='veff',veff=0)
    pc.set_rate_fcn('CH4',None)
    pc.set_temperature(pc.var.temperature)
    pc.set_press_temp_edd(pc.wrk.pressure,pc.var.temperature,pc.var.edd,1e-1*1e6)
    P = np.array([2.0*pc.wrk.surface_pressure*1.0e6,
                  0.5*pc.wrk.pressure_hydro[-1]])
    T = np.array([300.0, 180.0])
    edd = np.array([3.0e7, 4.0e5])
    pc.set_press_temp_edd_profile(
        P, T, edd, pc.wrk.pressure_hydro[pc.var.nz//2],
        hydro_pressure=True
    )
    pc.update_vertical_grid(TOA_pressure=1e-7*1e6)
    pc.clear_press_temp_edd_profile()
    pc.initialize_stepper(pc.wrk.usol)
    print(pc.wrk.t_history[0])
    print(pc.wrk.mix_history[0,0,0])
    pc.step()
    converged = pc.check_for_convergence()
    pc.destroy_stepper()
    pl = pc.production_and_loss('CH4',pc.wrk.usol)
    print(pl.production_rx[0])
    print(pl.loss_rx[0])
    print(pl.production[0,0])
    print(pl.loss[0,0])
    print(pl.integrated_production[0])
    print(pl.integrated_loss[0])

def main():

    test_static_construction()
    test_evolve_uses_fixed_grid()
    test_gas_giant_static_construction()
    test_gas_giant_returns_climate_levels_above_photochemical_grid()
    test_gas_giant_reinitializes_from_public_climate_grid()
    test_gas_giant_uses_shared_robust_stepper()
    test_gas_giant_shared_limits_and_state_restore()
    test_wrapper()
    test_initialize_atmosphere_z_particles()
    test_initialize_atmosphere_z_no_particles()
    test_initialize_atmosphere_p_particles()
    test_initialize_atmosphere_p_no_particles()

if __name__ == "__main__":
    main()
