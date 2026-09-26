"""Tests for continuity repair of gas-phase thermodynamic fits."""

from copy import deepcopy
import unittest

import yaml

from photochem.utils import check_thermo_continuity, make_thermo_continuous
from photochem.utils.thermo_continuity import _enthalpy_entropy, _heat_capacity


class ThermoContinuityTest(unittest.TestCase):
    def test_repair_matches_all_properties_outward_from_reference(self):
        network = {"species": []}
        adjustable = {"Shomate": (0, 5, 6), "NASA7": (0, 5, 6),
                      "NASA9": (2, 7, 8)}
        for model, length in (("Shomate", 7), ("NASA7", 7), ("NASA9", 9)):
            rows = [[0.0] * length for _ in range(3)]
            for index, row in enumerate(rows):
                row[adjustable[model][0]] = float(index + 2)
            network["species"].append({
                "name": model,
                "thermo": {"model": model,
                           "temperature-ranges": [100.0, 200.0, 400.0, 800.0],
                           "data": rows},
            })
        original = deepcopy(network)

        repaired, changes = make_thermo_continuous(network, reference_temperature=298.15)
        self.assertEqual(network, original)
        self.assertEqual(len(changes), 6)
        self.assertEqual(check_thermo_continuity(repaired), [])
        for before, after in zip(original["species"], repaired["species"]):
            model = after["thermo"]["model"]
            old_rows = before["thermo"]["data"]
            rows = after["thermo"]["data"]
            self.assertEqual(rows[1], old_rows[1])  # 298.15 K reference segment
            for row, old_row in zip(rows, old_rows):
                self.assertEqual([v for i, v in enumerate(row) if i not in adjustable[model]],
                                 [v for i, v in enumerate(old_row) if i not in adjustable[model]])
            for temperature, left, right in ((200.0, 0, 1), (400.0, 1, 2)):
                left_h, left_s = _enthalpy_entropy(model, rows[left], temperature)
                right_h, right_s = _enthalpy_entropy(model, rows[right], temperature)
                self.assertAlmostEqual(left_h, right_h, delta=1e-6)
                self.assertAlmostEqual(left_s, right_s, delta=1e-8)
                self.assertAlmostEqual(left_h - temperature * left_s,
                                       right_h - temperature * right_s, delta=1e-6)
                self.assertAlmostEqual(_heat_capacity(model, rows[left], temperature),
                                       _heat_capacity(model, rows[right], temperature), delta=1e-8)

        serialized = yaml.safe_dump(repaired)
        self.assertEqual(check_thermo_continuity(yaml.safe_load(serialized)), [])
        repaired_twice, further_changes = make_thermo_continuous(repaired)
        self.assertEqual(further_changes, [])
        self.assertEqual(repaired_twice, repaired)

    def test_checker_reports_heat_capacity_jump_alone(self):
        network = {"species": [{"name": "X", "thermo": {
            "model": "Shomate", "temperature-ranges": [200, 500, 1000],
            "data": [[1, 0, 0, 0, 0, 0, 0],
                     [2, 0, 0, 0, 0, -0.5, 0.6931471805599453]],
        }}]}
        jumps = check_thermo_continuity(network)
        self.assertEqual(len(jumps), 1)
        self.assertAlmostEqual(jumps[0]["delta_heat_capacity"], 1)
        self.assertAlmostEqual(jumps[0]["delta_enthalpy"], 0, delta=1e-9)
        self.assertAlmostEqual(jumps[0]["delta_entropy"], 0, delta=1e-9)

    def test_equilibrate_condensates_are_skipped(self):
        network = {"species": [
            {"name": "X", "thermo": {
                "model": "Shomate", "temperature-ranges": [200, 500, 1000],
                "data": [[1, 0, 0, 0, 0, 0, 0],
                         [2, 0, 0, 0, 0, 0, 0]],
            }},
            {"name": "Xaer", "condensate": True, "thermo": {
                "model": "Shomate", "temperature-ranges": [100, 150, 200],
                "data": [[1, 0, 0, 0, 0, 0, 0],
                         [2, 0, 0, 0, 0, 0, 0]],
            }},
        ]}
        original = deepcopy(network)
        self.assertEqual([jump["species"] for jump in check_thermo_continuity(network)],
                         ["X"])

        repaired, corrections = make_thermo_continuous(network)
        self.assertEqual(network, original)
        self.assertEqual(repaired["species"][1], original["species"][1])
        self.assertEqual([change["species"] for change in corrections], ["X"])
        self.assertEqual(check_thermo_continuity(repaired), [])


if __name__ == "__main__":
    unittest.main()
