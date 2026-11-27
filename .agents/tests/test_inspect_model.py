# SPDX-License-Identifier: LGPL-3.0-only
"""Inspect small native meshes without a compiled ThinCurr runtime."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy


SCRIPT = Path(__file__).resolve().parents[1] / 'skills/thincurr-model/scripts/inspect_model.py'


class InspectModelTests(unittest.TestCase):
    def inspect(self, points, cells, mesh_sets, *args):
        with tempfile.TemporaryDirectory() as directory:
            mesh = Path(directory) / 'mesh.h5'
            with h5py.File(mesh, 'w') as mesh_file:
                mesh_file['mesh/R'] = points
                mesh_file['mesh/LC'] = cells
                mesh_file['mesh/REG'] = numpy.ones(len(cells), dtype='i4')
                for name, values in mesh_sets.items():
                    mesh_file['mesh/' + name] = values
            return subprocess.run(
                [sys.executable, str(SCRIPT), str(mesh), '--no-setup', *args],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False,
            )

    def test_legacy_sideset_supplies_closure(self):
        points = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        cells = [(1, 3, 2), (1, 2, 4), (2, 3, 4), (3, 1, 4)]
        sets = {'NUM_SIDESETS': [1], 'SIDESET0001': [1]}
        result = self.inspect(points, cells, sets)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn('SIDESET0001', result.stdout)

        disk = self.inspect(points[:3], cells[:1], sets)
        self.assertEqual(disk.returncode, 1, disk.stdout)
        self.assertIn('File defines 1 closure(s), but the topology needs 0', disk.stdout)

    def test_legacy_jumpers_do_not_count_as_holes(self):
        # A square annulus has one hole. Its only nodeset is reserved for a jumper.
        points = [(-2, -2, 0), (2, -2, 0), (2, 2, 0), (-2, 2, 0),
                  (-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0)]
        cells = [(1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
                 (3, 4, 8), (3, 8, 7), (4, 1, 5), (4, 5, 8)]
        sets = {'NUM_NODESETS': [1], 'NODESET0001': [5, 6, 7, 8]}
        valid = self.inspect(points, cells, sets)
        self.assertEqual(valid.returncode, 0, valid.stdout)
        for start in ('1', '-1'):
            result = self.inspect(points, cells, sets, '--jumper-start', start)
            self.assertEqual(result.returncode, 1, result.stdout)
            self.assertIn('Legacy nodesets define 0 hole(s), but the topology needs 1', result.stdout)
        for start in ('2', '-2'):
            result = self.inspect(points, cells, sets, '--jumper-start', start)
            self.assertEqual(result.returncode, 1, result.stdout)
            self.assertIn('exceeds number of nodesets', result.stdout)
