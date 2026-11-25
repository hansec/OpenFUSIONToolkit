# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../../README.md).
# Not part of the regression suite.
"""Run a command and print its exit code and peak RSS (children)."""
import resource
import subprocess
import sys
rc = subprocess.call(sys.argv[1:])
print('EXITCODE=%d' % rc)
print('Maximum resident set size (kbytes): %d' % resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
