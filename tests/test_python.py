"""Smoke test for the installed Python wrappers. See README.md."""

import numpy as np

from photochem import EvoAtmosphere, zahnle_earth


def test_initialize_atmosphere_z_no_particles():
    pc = EvoAtmosphere(
        "no_particle_test.yaml",
        "test_settings_top_atmospherefile.yaml",
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
    print(pc.var.usol_init[0,0])
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
    print(pc.var.at_photo_equilibrium)
    print(pc.var.surface_pressure)
    print(pc.var.max_error_reinit_attempts)
    print(pc.var.rtol)
    print(pc.var.atol)
    print(pc.var.mxsteps)
    print(pc.var.equilibrium_time)
    print(pc.var.conv_hist_factor)
    print(pc.var.conv_min_mix)
    print(pc.var.conv_longdy)
    print(pc.var.conv_longdydt)
    print(pc.var.autodiff)
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
    pc.out2atmosphere_txt('tmp.txt',overwrite=True)
    pc.gas_fluxes()
    pc.set_lower_bc('O',bc_type='vdep',vdep=0)
    pc.set_upper_bc('O',bc_type='veff',veff=0)
    pc.set_rate_fcn('CH4',None)
    pc.set_temperature(pc.var.temperature)
    pc.set_press_temp_edd(pc.wrk.pressure,pc.var.temperature,pc.var.edd,1e-1*1e6)
    P = np.array([2.0*pc.var.surface_pressure*1.0e6,
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

    test_wrapper()
    test_initialize_atmosphere_z_particles()
    test_initialize_atmosphere_z_no_particles()

if __name__ == "__main__":
    main()
