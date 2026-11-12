"""Local setup checks for the ported Venus case; no long integration."""
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
import yaml
from initialize import initialize

HERE = Path(__file__).resolve().parent


def main():
    original = yaml.safe_load((HERE / 'mechanism_legacy.yaml').read_text())
    repaired = yaml.safe_load((HERE / 'mechanism.yaml').read_text())
    assert original['reactions'] == repaired['reactions']
    assert original['particles'] == repaired['particles']
    for old, new in zip(original['species'], repaired['species']):
        assert {k: v for k, v in old.items() if k != 'thermo'} == {
            k: v for k, v in new.items() if k != 'thermo'}
    with TemporaryDirectory() as tmp:
        for fresh in (False, True):
            pc = initialize(Path(tmp), fresh)
            assert pc.var.jacobian_method == 3
            assert pc.wrk.usol.shape == (114, 100)
            assert np.all(np.isfinite(pc.wrk.usol))
            assert np.all(pc.wrk.usol >= 0)
            assert np.all(np.diff(pc.var.z) > 0)
            assert abs(pc.wrk.surface_pressure - 92.1) < .05
            assert np.isclose(pc.var.temperature[-1], 172.2)
            for name, radius in {'H2SO4aer': 2e-4, 'H2Oaer': 1e-3}.items():
                i = pc.dat.species_names.index(name)
                np.testing.assert_allclose(pc.var.particle_radius[i], radius)
            for params in pc.var.cond_params:
                assert (params.smooth_factor, params.k_cond, params.k_evap) == (5., 100., 100.)
            pc.prep_atmosphere(pc.wrk.usol)
            z = pc.var.z
            wavelength = .5*(pc.dat.wavl[:-1] + pc.dat.wavl[1:])
            coefficient = np.where(z > 67e5, .056e-5*np.exp(-(z-67e5)/3e5),
                                   np.where(z > 58e5, .056e-5, 0.))
            extinction = coefficient[:, None]*np.exp(-(wavelength[None, :]-360.)/100.)
            extinction[:, wavelength > 400.] = 0.
            # At the reference atmosphere, the pressure-aware API reproduces
            # exactly the legacy imposed extinction times layer thickness.
            np.testing.assert_allclose(pc.wrk.tauc, extinction*pc.var.dz[:, None],
                                       rtol=1e-9, atol=1e-25)
            np.testing.assert_array_equal(pc.wrk.w0c, 0)
            np.testing.assert_array_equal(pc.wrk.g0c, 0)
    print('Venus fresh/restart setup and initial UV opacity checks passed')


if __name__ == '__main__':
    main()
