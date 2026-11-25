# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
from analyze import load, restrict_free, cond_sym, Schwarz, schwarz_metrics, blocks_from_color
from basis import assemble
d, p, kind = sys.argv[1], int(sys.argv[2]), sys.argv[3]
alpha = (3, 6)
if len(sys.argv) > 4:
    alpha = tuple(int(a) for a in sys.argv[4].split(','))
Mo, info, colors, cells = load(d)
M = Mo if kind == 'dump' else assemble(d, p, kind, alpha)[0].tocsr()
Mf, finfo, newid = restrict_free(M, info)
keep = newid >= 0
res = {'n': Mf.shape[0]}
k, lmin, lmax = cond_sym(Mf)
res['cond_M'] = k
dg = sp.diags(1.0 / np.sqrt(Mf.diagonal()))
res['cond_diagM'] = cond_sym((dg @ Mf @ dg).tocsr())[0]
key = finfo[:, 1] * 10**7 + finfo[:, 2]
_, inv = np.unique(key, return_inverse=True)
order_ = np.argsort(inv, kind='stable')
eblks = np.split(order_, np.cumsum(np.bincount(inv))[:-1])
res['bj_perentity'] = schwarz_metrics(Mf, Schwarz(Mf, eblks), True)
res['bj_entity'] = schwarz_metrics(Mf, Schwarz(Mf, blocks_from_color(colors['entity'], keep)), True)
for tag in ('ov0', 'ov1p'):
    blks, asgs = cells[tag]
    bl = [newid[b][newid[b] >= 0] for b in blks]
    asg = [newid[a][newid[a] >= 0] for a in asgs]
    res['ras_' + tag] = schwarz_metrics(Mf, Schwarz(Mf, bl, assigned=asg), False)
    res['as_' + tag] = schwarz_metrics(Mf, Schwarz(Mf, bl, additive=True), True)
print('== %s p=%d %s alpha=%s' % (d, p, kind, alpha))
for k_, v in res.items():
    if isinstance(v, dict):
        print('%-12s' % k_, ' '.join('%s=%.4g' % (a, b) for a, b in v.items()))
    else:
        print('%-12s %.4g' % (k_, v))
