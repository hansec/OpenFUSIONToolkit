# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from analyze import load, restrict_free
d = sys.argv[1]
M, info, colors, cells = load(d)
Mf, fi, nid = restrict_free(M, info)
key = fi[:, 1] * 10**7 + fi[:, 2]
_, inv = np.unique(key, return_inverse=True)
o = np.argsort(inv, kind='stable'); eb = np.split(o, np.cumsum(np.bincount(inv))[:-1])
# D_E = per-entity block diagonal; symmetric scaling L^-1 M L^-T
n = Mf.shape[0]; rows=[]; cols=[]; vals=[]
for b in eb:
    Li = np.linalg.inv(np.linalg.cholesky(Mf[b][:, b].toarray()))
    rr, cc = np.meshgrid(b, b, indexing='ij'); rows.append(rr.ravel()); cols.append(cc.ravel()); vals.append(Li.ravel())
Linv = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
A = (Linv @ Mf @ Linv.T).tocsc()
w, V = spla.eigsh(A, k=6, sigma=0, which='LM')
names = {(0,1):'curl-edge',(0,2):'curl-face',(0,3):'curl-cell',(1,1):'grad-edge',(1,2):'grad-face',(1,3):'grad-cell'}
print('==', d, 'smallest eigenvalues', np.round(w, 5))
for k in range(2):
    x = V[:, k]   # in entity-orthonormalized coordinates: energy per entity group
    fr = {}
    for (f, dm), nm in names.items():
        m = (fi[:, 0] == f) & (fi[:, 1] == dm)
        fr[nm] = np.sum(x[m]**2)
    # spatial extent: number of entities holding 90% of the norm
    ent = np.array([np.sum(x[b]**2) for b in eb]); s = np.sort(ent)[::-1]
    n90 = np.searchsorted(np.cumsum(s), 0.9) + 1
    print('  mode %d lam=%.4g' % (k, w[k]), ' '.join('%s=%.2f' % (a, b) for a, b in fr.items() if b > 0.005), ' entities for 90%%: %d of %d' % (n90, len(eb)))
