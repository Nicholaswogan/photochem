"""Local, opt-in PTC comparisons; run with the Python used to build Photochem.

Prepare inputs without modifying state1.pkl, compile the Fortran driver, and
write generated files beneath results/ (ignored by Git). No Python API changes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import importlib.util
import shlex
import shutil
import subprocess
import sys
import time

import photochem_clima_data

HERE = Path(__file__).resolve().parent
BENCHMARKS = HERE.parent
ROOT = BENCHMARKS.parent


def prepare(destination, case, fresh=False):
    destination.mkdir(parents=True, exist_ok=True)
    case_dir = BENCHMARKS / case
    spec = importlib.util.spec_from_file_location('benchmark_case', case_dir / 'initialize.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    pc = module.initialize(destination, fresh)
    controls = module.CONTROLS
    shutil.copyfile(case_dir / module.MECHANISM, destination / 'mechanism.yaml')
    shutil.copyfile(case_dir / module.SETTINGS, destination / 'settings.yaml')
    pc.model_state().save_to_h5(str(destination / 'initial.h5'))
    save_controls(destination, controls)


def save_controls(destination, controls):
    (destination / 'controls.json').write_text(json.dumps(controls, indent=2) + '\n')
    values = []
    for name, value in controls.items():
        formatted = ('.true.' if value else '.false.') if isinstance(value, bool) else str(value)
        values.append(f'  {name} = {formatted}')
    (destination / 'controls.nml').write_text('&benchmark_controls\n' + '\n'.join(values) + '\n/\n')


def compile_driver(build, output):
    """Reuse CMake's platform-specific compiler/link configuration."""
    args = shlex.split((build / 'tests/photochem/CMakeFiles/photochem_test_api.dir/link.txt').read_text())
    command = []
    skip = False
    for arg in args:
        if skip:
            skip = False
            continue
        if arg == '-o':
            command += ['-o', str(output)]
            skip = True
        elif arg.endswith('/test_api.f90.o'):
            command += ['-I' + str(build / 'modules'), '-I' + str(build / 'fortran_STATIC'),
                        '-ffree-line-length-none', str(HERE / 'ptc_driver.f90')]
        elif arg.startswith('../../lib/'):
            command.append(str(build / 'lib' / Path(arg).name))
        else:
            command.append(arg)
    env = dict(os.environ)
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    with (output.parent / "compile.log").open("w") as stream:
        subprocess.run(command, check=True, env=env, stdout=stream, stderr=subprocess.STDOUT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['archean', 'venus'], default='archean')
    parser.add_argument('--fresh', action='store_true')
    parser.add_argument('--jacobian-method', choices=['analytical', 'autodiff', 'finite-difference'],
                        help='Override the case Jacobian consistently for every backend')
    parser.add_argument('--reuse-inputs', action='store_true', help='Use previously prepared HDF5 inputs')
    parser.add_argument('--build', type=Path, default=ROOT / 'build')
    parser.add_argument('--modes', nargs='+', choices=['cvode', 'cvode1', 'cvode1-tight',
                        'ptc', 'ptc-residual', 'ptc-periodic', 'ptc-toa', 'ptc-floor12'], default=['cvode', 'cvode1', 'ptc', 'ptc-residual'])
    parser.add_argument('--steps', type=int, default=7000)
    parser.add_argument('--seconds', type=float, default=120)
    args = parser.parse_args()
    if args.steps < 1 or not (0 < args.seconds < float('inf')):
        parser.error('--steps and --seconds must be positive and finite')
    if args.case == 'venus' and any(mode in ('ptc-periodic', 'ptc-toa') for mode in args.modes):
        parser.error('Venus uses fixed altitude profiles; profile/TOA perturbation modes are Archean-only')
    output = BENCHMARKS / args.case / 'results' / ('fresh' if args.fresh else 'restart')
    if not args.reuse_inputs:
        prepare(output, args.case, args.fresh)
    if args.reuse_inputs and not (output / 'controls.nml').exists():
        parser.error('No prepared case controls found; rerun without --reuse-inputs')
    if args.jacobian_method:
        controls = json.loads((output / 'controls.json').read_text())
        controls['jacobian_method'] = {'analytical': 3, 'autodiff': 1, 'finite-difference': 2}[args.jacobian_method]
        save_controls(output, controls)
    print('Case controls:', (output / 'controls.json').read_text(), flush=True)
    hashes = {name: hashlib.sha256((output / name).read_bytes()).hexdigest()
              for name in ['mechanism.yaml', 'settings.yaml', 'flux.txt', 'initial.h5', 'controls.nml']}
    (output / 'input_hashes.json').write_text(json.dumps(hashes, indent=2) + '\n')
    executable = output / 'ptc_driver'
    compile_driver(args.build.resolve(), executable)
    data = Path(photochem_clima_data.__file__).parent / 'data'
    env = dict(os.environ, OMP_NUM_THREADS='1')
    for mode in args.modes:
        begin = time.perf_counter()
        with (output / (mode + '.log')).open('w') as stream:
            subprocess.run([str(executable), str(output), str(data), mode, str(output / mode),
                            str(args.steps), str(args.seconds)], stdout=stream,
                           stderr=subprocess.STDOUT, check=True, env=env,
                           timeout=args.seconds * 3 + 60)
        print(mode, 'wall_seconds=', round(time.perf_counter() - begin, 3), flush=True)
        print((output / (mode + '.summary')).read_text(), flush=True)


if __name__ == '__main__':
    main()
