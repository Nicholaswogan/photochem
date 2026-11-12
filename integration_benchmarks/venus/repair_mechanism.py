"""Repair legacy thermodynamic joins for the current strict mechanism loader."""
from pathlib import Path
import json
import yaml
from photochem.utils import make_thermo_continuous

HERE = Path(__file__).resolve().parent


def main():
    network = yaml.load((HERE / 'mechanism_legacy.yaml').read_text(), Loader=yaml.SafeLoader)
    repaired, corrections = make_thermo_continuous(network)
    (HERE / 'mechanism.yaml').write_text(yaml.safe_dump(repaired, sort_keys=False))
    (HERE / 'thermo_corrections.json').write_text(json.dumps(corrections, indent=2) + '\n')
    print('Adjusted thermodynamic segments:', len(corrections))


if __name__ == '__main__':
    main()
