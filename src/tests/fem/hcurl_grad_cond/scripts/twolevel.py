# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse.linalg as spla
from analyze import load, restrict_free, Schwarz, blocks_from_color
d = sys.argv[1]
lag = '--lag' in sys.argv
M, info, colors, cells = load(d, 'lag_' if lag else '')
if lag:
    Mf, finfo, newid = M.tocsr(), info, np.arange(M.shape[0])
    # coarse: vertex DoFs (P1) of each component
    coarse = np.nonzero(finfo[:, 1] == 0)[0]
else:
    Mf, finfo, newid = restrict_free(M, info)
    # coarse p=1 space: Whitney edge + first H1 edge DoF (hierarchical injection)
    coarse = np.nonzero((finfo[:, 1] == 1) & (finfo[:, 3] == 0))[0]
keep = newid >= 0
n = Mf.shape[0]
Mc = Mf[coarse][:, coarse].tocsc()
luc = spla.splu(Mc)
def Cc(r):
    z = np.zeros(n); z[coarse] = luc.solve(r[coarse]); return z
key = finfo[:, 1] * 10**7 + finfo[:, 2]
_, inv = np.unique(key, return_inverse=True)
o = np.argsort(inv, kind='stable'); eblks = np.split(o, np.cumsum(np.bincount(inv))[:-1])
smoothers = {'entity8': Schwarz(Mf, blocks_from_color(colors['entity'], keep)),
             'metis8': Schwarz(Mf, blocks_from_color(colors['metis'], keep)),
             'perentity': Schwarz(Mf, eblks)}
blks, asgs = cells['ov1p']
bl = [newid[b][newid[b] >= 0] for b in blks]; asg = [newid[a][newid[a] >= 0] for a in asgs]
smoothers['ras_ov1p'] = Schwarz(Mf, bl, assigned=asg)
def specrad(fn):
    E = spla.LinearOperator((n, n), matvec=fn, dtype=float)
    return np.max(np.abs(spla.eigs(E, k=2, which='LM', return_eigenvectors=False, tol=1e-6)))
def extremal(fn):
    A = spla.LinearOperator((n, n), matvec=fn, dtype=float)
    lmax = np.max(spla.eigs(A, k=1, which='LR', return_eigenvectors=False, tol=1e-6).real)
    s = 1.05 * lmax
    B = spla.LinearOperator((n, n), matvec=lambda x: s * x - fn(x), dtype=float)
    lmin = s - np.max(spla.eigs(B, k=1, which='LR', return_eigenvectors=False, tol=1e-6).real)
    return lmin, lmax
print('==', d, 'lag' if lag else 'hcurl', 'n', n, 'ncoarse', len(coarse))
for name, S in smoothers.items():
    out = []
    if name != 'ras_ov1p':
        lmin, lmax = extremal(lambda x: S.apply(Mf @ x) + Cc(Mf @ x))
        out.append('2L-additive kappa=%.3g' % (lmax / lmin))
    # two-grid "up" cycle: coarse correction then one smoothing step (undamped)
    def E(x):
        y = x - Cc(Mf @ x)
        return y - S.apply(Mf @ y)
    out.append('rho(TG coarse+1 smooth)=%.3g' % specrad(E))
    print('  %-10s' % name, ' '.join(out), flush=True)
