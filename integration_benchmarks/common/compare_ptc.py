"""Compare generated snapshots and plot local solver progress (not a CI test)."""
import argparse
import json
from pathlib import Path
import numpy as np
from photochem import EvoAtmosphereState

BENCHMARKS = Path(__file__).resolve().parent.parent


def compare(directory):
    reference = EvoAtmosphereState(str(directory / 'cvode1.h5')).to_dict()
    result = {}
    def reached_target(mode):
        summary = directory / (mode + '.summary')
        if not summary.exists():
            return False
        return any(line.split() == ['finished', 'T'] for line in summary.read_text().splitlines())
    reference_converged = reached_target('cvode1')
    for path in sorted(directory.glob('*.h5')):
        if path.name in ('initial.h5', 'cvode1.h5'):
            continue
        state = EvoAtmosphereState(str(path)).to_dict()
        flags = {'reference_reached_target': reference_converged,
                 'candidate_reached_target': reached_target(path.stem),
                 'modified_physics': path.stem in ('ptc-periodic', 'ptc-toa')}
        # Forced TOA experiments use different grids and are not direct comparisons.
        if not np.array_equal(state['z'], reference['z']):
            result[path.stem] = dict(flags, different_grid=True)
            continue
        a, b = reference['usol'], state['usol']
        gas_start = reference['nq'] - reference['ng'] + reference['nsl']
        abundance = a / np.maximum(a[gas_start:].sum(axis=0), 1e-100)
        ca, cb = a @ reference['dz'], b @ state['dz']
        columns = ca > 1e-10 * ca.sum()
        item = {**flags, 'max_column_relative_above_1e-10': float(np.max(np.abs(cb-ca)[columns]/ca[columns])),
                'max_temperature_difference_K': float(np.max(np.abs(state['temperature']-reference['temperature'])))}
        for floor in (1e-20, 1e-10, 1e-6):
            mask = abundance > floor
            item[f'max_profile_relative_above_{floor}'] = float(np.max(np.abs(b-a)[mask]/a[mask]))
        result[path.stem] = item
    (directory / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def plot(case):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, name in zip(axes, ('restart', 'fresh')):
        directory = BENCHMARKS / case / 'results' / name
        for mode in ('cvode', 'cvode1', 'ptc', 'ptc-residual', 'ptc-floor12'):
            path = directory / (mode + '.csv')
            if not path.exists():
                continue
            data = np.genfromtxt(path, delimiter=',', names=True, usecols=(0, 1))
            data = np.atleast_1d(data)
            ax.semilogy(data['step'], data['time'], label=mode)
        controls_path = directory / 'controls.json'
        target = json.loads(controls_path.read_text())['equilibrium_time'] if controls_path.exists() else 1e17
        ax.axhline(target, color='k', linestyle=':', linewidth=1)
        ax.set(xlabel='Accepted step', ylabel='Time / pseudo-time (s)', title=name)
        ax.legend(fontsize=8)
    fig.savefig(BENCHMARKS / case / 'results' / 'progress.pdf')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['archean', 'venus'], default='archean')
    args = parser.parse_args()
    for name in ('restart', 'fresh'):
        directory = BENCHMARKS / args.case / 'results' / name
        if (directory / 'cvode1.h5').exists():
            print(args.case, name, json.dumps(compare(directory), indent=2))
    plot(args.case)


if __name__ == '__main__':
    main()
