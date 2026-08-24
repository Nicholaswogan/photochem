import numpy as np
from pathlib import Path

from photochem.equilibrate import ChemEquiAnalysis

def main():
    atoms = [
        'H ',
        'He',
        'C ',
        'N ',
        'O ',
        'Na',
        'Mg',
        'Al',
        'Si',
        'P ',
        'S ',
        'Cl',
        'K ',
        'Ca',
        'Ti',
        'V ',
        'Fe',
        'Ni'
    ]
    X = np.array([
        9.207539305000000e-01,
        7.836886940000000e-02,
        2.478241000000000e-04,
        6.225060569498810e-05,
        4.509658000000000e-04,
        1.600086943532050e-06,
        3.665587420553620e-05,
        2.595000000000000e-06,
        2.979500000000000e-05,
        2.366702019976680e-07,
        1.213790073460400e-05,
        2.911679584995890e-07,
        9.866056119256769e-08,
        2.014390114292550e-06,
        8.206228043663590e-08,
        7.836886940899920e-09,
        2.911679584995890e-05,
        1.528071168062810e-06
    ])

    # 
    thermofile = Path(__file__).with_name('thermo_easy_chem_simp_own.yaml')
    cea = ChemEquiAnalysis(str(thermofile), atoms=atoms)
    converged = cea.solve(1.0e6, 1000.0, molfracs_atoms=X)
    
    # Check the solution
    ind = cea.species_names.index('H2')
    assert np.isclose(cea.molfracs_species[ind],8.531731613887943e-01,rtol=1e-4)

if __name__ == '__main__':
    main()
