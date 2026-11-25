# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse.linalg as spla
from analyze import load, restrict_free
from basis import assemble
d, p = sys.argv[1], int(sys.argv[2])
Mo, info, colors, cells = load(d)
F = lambda x: np.stack([np.sin(2*x[:,1]+x[:,2]), np.cos(x[:,0]*x[:,2]), np.exp(0.3*x[:,0])*x[:,1]], axis=1)
for kind in ('oft', 'sz'):
    M, lay, b = assemble(d, p, kind, field=F)
    A, fi, nid = restrict_free(M.tocsr(), info)
    bb = b[nid >= 0]
    x = spla.spsolve(A.tocsc(), bb)
    print(kind, 'n', A.shape[0], 'b^T M^-1 b = %.15e' % (bb @ x))
