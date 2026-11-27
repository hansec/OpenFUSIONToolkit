# SPDX-License-Identifier: LGPL-3.0-only
"""Exercise uninstall safety from the command line in an isolated project."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


INSTALL = Path(__file__).resolve().parents[1] / 'install.py'


class UninstallSafetyTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.project = Path(directory.name)

    def install(self, *args, expected_exit=0):
        result = subprocess.run(
            [sys.executable, str(INSTALL), '--scope', 'project', '--project-dir',
             str(self.project), *args],
            stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, expected_exit, result.stdout + result.stderr)
        return result.stdout + result.stderr

    def test_uninstall_requires_force(self):
        marker = self.project / '.claude/skills/thincurr-run/user-file'
        marker.parent.mkdir(parents=True)
        marker.write_text('user content')
        target = ('--target', 'anthropic', '--skills', 'thincurr-run')
        self.assertIn('exists', self.install(*target))
        self.assertIn('--force', self.install(*target, '--uninstall', expected_exit=2))
        self.assertIn(str(marker.parent), self.install(*target, '--uninstall', '--dry-run'))
        self.assertEqual(marker.read_text(), 'user content')
        self.install(*target, '--uninstall', '--force')
        self.assertFalse(marker.exists())
