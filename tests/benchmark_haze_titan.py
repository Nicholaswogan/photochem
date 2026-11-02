"""Compare Python haze microphysics with exported native McKay results.

Run from the repository root::

    PYTHONPATH=. python tests/benchmark_haze_titan.py

Only mckay_haze_reference.txt is needed; no private Fortran or compiler is
required. Atmospheric grid, source samples, material density, charging
radius and primary-particle radius come from that file.
The source is neither regenerated nor renormalized, and no sweeping is added.

This compares the two physical/numerical implementations, not identical
algorithms. Our model uses Stokes-Cunningham settling, adaptive integration,
and linear source interpolation. McKay uses its native settling blend,
stepwise radius update, and AVERGE source quadrature. Its top number density
is nonzero while its top mass flux is zero. Our revised source instead
creates primary particles (J=q/m0) with zero incoming mass and number flux.
The source closure is therefore deliberately different from McKay's J=0.
The resulting PDF is ignored by Git; this manual comparison is outside CI.
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from photochem.extensions.haze import mckay_microphysics


HERE = Path(__file__).resolve().parent


def load_reference():
    """Read the self-describing space-separated native CGS reference file."""
    path = HERE/'mckay_haze_reference.txt'
    metadata = {}
    for line in path.read_text().splitlines():
        if line.startswith('# ') and ': ' in line:
            key, value = line[2:].split(': ', 1)
            metadata[key] = value
    columns = metadata['columns'].split()
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] != len(columns) or not np.all(np.isfinite(raw)):
        raise ValueError('Invalid McKay reference table.')
    reference = dict(zip(columns, raw.T))
    if not np.all(np.diff(reference['altitude_cm']) > 0):
        raise ValueError('Reference must be ordered bottom to top.')
    return reference, metadata


def run_benchmark():
    """Solve on the reference grid with native mass production and J=q/m0."""
    reference, metadata = load_reference()
    density = float(metadata['particle_density_g_cm3'])
    monomer_radius = float(metadata['monomer_radius_cm'])
    charging_radius = float(metadata['charging_radius_cm'])
    # Preserve the exported source samples and their original normalization.
    radius, number, flux = mckay_microphysics(
        reference['altitude_cm'], reference['temperature_K'],
        reference['air_density_g_cm3'], reference['air_number_density_cm3'],
        reference['gravity_cm_s2'], reference['production_g_cm3_s'],
        particle_density=density, monomer_radius=monomer_radius,
        charging_radius=charging_radius, rtol=1e-6,
    )
    internal_column = np.sum(
        .5*(reference['production_g_cm3_s'][1:] + reference['production_g_cm3_s'][:-1])
        *np.diff(reference['altitude_cm'])
    )
    assert number[-1] == flux[-1] == 0.0
    np.testing.assert_allclose(flux[0], internal_column, rtol=1e-12)
    assert np.all(np.isfinite(radius)) and np.all(np.isfinite(number))
    return reference, dict(radius=radius, number=number, flux=flux,
                          internal_column=internal_column)


def main():
    reference, result = run_benchmark()
    z = reference['altitude_cm']/1e5
    fig, axes = plt.subplots(1, 4, figsize=(13, 6.5), sharey=True)
    ours, native = '#2166ac', '#b2182b'
    axes[0].semilogx(result['radius']*1e4, z, color=ours, label='Python (J=q/m₀)')
    axes[0].semilogx(reference['particle_radius_cm']*1e4, z, color=native,
                     ls='--', label='Native McKay Fortran')
    axes[0].set_xlabel('Radius (µm)')
    axes[1].semilogx(np.where(result['number'] > 0, result['number'], np.nan), z, color=ours)
    axes[1].semilogx(reference['particle_number_density_cm3'], z, color=native, ls='--')
    axes[1].set_xlabel('Number density (cm$^{-3}$)')
    source = reference['production_g_cm3_s']
    axes[2].plot(source/np.max(source), z, color='0.3')
    axes[2].set_xlabel('Source / maximum\n(same samples for both)')
    axes[3].plot(result['flux']/1e-14, z, color=ours)
    axes[3].plot(reference['mass_flux_g_cm2_s']/1e-14, z, color=native, ls='--')
    axes[3].set_xlabel('Downward mass flux\n($10^{-14}$ g cm$^{-2}$ s$^{-1}$)')
    axes[0].set_ylabel('Altitude (km)')
    for ax in axes:
        ax.set_ylim(0, z[-1])
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.grid(alpha=.2)
    fig.suptitle('Titan haze: identical atmospheric inputs and native production samples')
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center',
               bbox_to_anchor=(.5, .94), ncol=2, frameon=False)
    fig.text(.5, .025,
             'Native initial atmosphere; no rainout or climate iteration.\n'
             'Python creates primary particles (J=q/m₀); McKay adds mass to existing particles (J=0).',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .09, 1, .89))
    output = HERE/'benchmark_haze_titan.pdf'
    fig.savefig(output)
    plt.close(fig)
    print(f'Saved {output}')
    print(f'Native internal column: {reference["mass_flux_g_cm2_s"][0]:.6e} g/cm^2/s')
    print(f'Linear-source integral: {result["internal_column"]:.6e} g/cm^2/s')
    print('Altitude   Radius (µm): Python / native    Density (cm^-3): Python / native')
    for height in (100, 200, 350, 500, 600):
        r = np.interp(height, z, result['radius'])*1e4
        r_native = np.interp(height, z, reference['particle_radius_cm'])*1e4
        n = np.interp(height, z, result['number'])
        n_native = np.interp(height, z, reference['particle_number_density_cm3'])
        print(f'{height:4d} km          {r:.4f} / {r_native:.4f}'
              f'                     {n:.3g} / {n_native:.3g}')


if __name__ == '__main__':
    main()
