# SPDX-License-Identifier: LGPL-3.0-only
"""Exercise shared-skill retention during uninstall in an isolated project."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


INSTALL = Path(__file__).resolve().parents[1] / 'install.py'


class SharedSkillsUninstallTests(unittest.TestCase):
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

    def test_uninstall_one_tool_keeps_shared_skills(self):
        self.install('--target', 'openai', 'google')
        skills = self.project / '.agents/skills/thincurr-run'
        codex = self.project / '.codex/agents'
        gemini = self.project / '.gemini/agents'
        self.assertTrue(list(codex.iterdir()) and list(gemini.iterdir()))
        self.assertIn('Keeping shared skills', self.install('--target', 'openai', '--uninstall', '--force'))
        self.assertTrue(skills.is_dir())
        self.assertFalse(list(codex.iterdir()))
        self.assertTrue(list(gemini.iterdir()))
        self.install('--target', 'openai', 'google', '--uninstall', '--force')
        self.assertFalse(skills.exists())
        self.assertFalse(list(gemini.iterdir()))
