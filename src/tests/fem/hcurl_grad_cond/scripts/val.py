# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import load, restrict_free
from basis import assemble
d, p = sys.argv[1], int(sys.argv[2])
Mo, info, colors, cells = load(d)
Mp, lay = assemble(d, p, 'oft')
A, fi, nid = restrict_free(Mo, info)
B, _, _ = restrict_free(Mp.tocsr(), info)
D = abs(A - B)
print('p', p, 'free n', A.shape[0], 'max|diff| =', D.max(), 'max|M| =', abs(A).max(), 'rel', D.max() / abs(A).max())
