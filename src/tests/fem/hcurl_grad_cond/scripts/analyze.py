# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
"""Conditioning / block-Jacobi analysis of dumped H(Curl)+Grad(H^1) mass matrices.

Usage: python analyze.py <dump_dir> [--skip-cond]
"""
import sys
import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla


def load(d, pre=''):
    with open(os.path.join(d, pre + 'mass.txt')) as f:
        nr, nnz = map(int, f.readline().split())
        data = np.loadtxt(f)
    M = sp.csr_matrix((data[:, 2], (data[:, 0].astype(int) - 1, data[:, 1].astype(int) - 1)), shape=(nr, nr))
    M = (M + M.T) * 0.5
    with open(os.path.join(d, pre + 'dofs.txt')) as f:
        nfields, np_, ne, nf, nc = map(int, f.readline().split())
        fields = [list(map(int, f.readline().split())) for _ in range(nfields)]
    # DoF -> (field, entity dim, entity index)
    info = []
    counts = [np_, ne, nf, nc]
    for fi, (nef, g1, g2, g3, g4) in enumerate(fields):
        g = [g1, g2, g3, g4]
        for dim in range(4):
            for e in range(counts[dim]):
                for k in range(g[dim]):
                    info.append((fi, dim, e, k))
    info = np.array(info)
    assert len(info) == nr, (len(info), nr)

    def ints(name):
        v = np.loadtxt(os.path.join(d, pre + name), dtype=int)
        return v[1:]

    colors = {'metis': ints('color_metis.txt'), 'entity': ints('color_entity.txt')}
    cells = {}
    for tag in ('ov0', 'ov1p', 'ov1f'):
        with open(os.path.join(d, pre + 'cell_blocks_%s.txt' % tag)) as f:
            toks = np.array(f.read().split(), dtype=int)
        pos = 0
        nb = toks[pos]; pos += 1
        blks, asgs = [], []
        for _ in range(nb):
            n = toks[pos]; pos += 1
            blks.append(toks[pos:pos + n] - 1); pos += n
            n = toks[pos]; pos += 1
            asgs.append(toks[pos:pos + n] - 1); pos += n
        cells[tag] = (blks, asgs)
    return M, info, colors, cells


def restrict_free(M, info):
    """Remove constrained vertex-gradient DoFs (H^1 field, vertex entities)."""
    free = ~((info[:, 0] == 1) & (info[:, 1] == 0))
    idx = np.nonzero(free)[0]
    newid = -np.ones(len(free), dtype=int)
    newid[idx] = np.arange(len(idx))
    return M[idx][:, idx].tocsr(), info[idx], newid


class Schwarz:
    """Block (overlapping) Schwarz: B^{-1} r = sum_i W_i R_i^T A_i^{-1} R_i r"""
    def __init__(self, M, blocks, assigned=None, additive=False):
        self.n = M.shape[0]
        self.blocks = blocks
        self.facs = []
        self.wmask = []
        for i, b in enumerate(blocks):
            Ai = M[b][:, b].toarray()
            self.facs.append(sla.cho_factor(Ai))
            if additive or assigned is None:
                self.wmask.append(np.ones(len(b), dtype=bool))
            else:
                self.wmask.append(np.isin(b, assigned[i]))

    def apply(self, r):
        z = np.zeros(self.n)
        for b, fac, w in zip(self.blocks, self.facs, self.wmask):
            zi = sla.cho_solve(fac, r[b])
            z[b[w]] += zi[w]
        return z


def cond_sym(M):
    lmax = spla.eigsh(M, k=1, which='LA', return_eigenvectors=False, tol=1e-6)[0]
    lmin = spla.eigsh(M, k=1, sigma=0, which='LM', return_eigenvectors=False, tol=1e-6)[0]
    return lmax / lmin, lmin, lmax


def schwarz_metrics(M, S, symmetric):
    n = M.shape[0]
    BM = spla.LinearOperator((n, n), matvec=lambda x: S.apply(M @ x), dtype=float)
    out = {}
    if symmetric:
        Binv = spla.LinearOperator((n, n), matvec=S.apply, dtype=float)
        lmax = spla.eigsh(M, k=1, M=None, which='LA', return_eigenvectors=False, tol=1e-5,
                          Minv=None) if False else None
        # Generalized: M x = lam B x, with B^{-1} available -> use eigs on B^{-1}M (similar to SPD)
        lmax = np.max(spla.eigs(BM, k=1, which='LR', return_eigenvectors=False, tol=1e-6).real)
        # smallest eigenvalue: largest of I - B^{-1}M / lmax shift
        shift = lmax * 1.05
        Op = spla.LinearOperator((n, n), matvec=lambda x: shift * x - S.apply(M @ x), dtype=float)
        lmin = shift - np.max(spla.eigs(Op, k=1, which='LR', return_eigenvectors=False, tol=1e-6).real)
        out['lmin'] = lmin
        out['lmax'] = lmax
        out['kappa'] = lmax / lmin
        out['rho_undamped'] = max(abs(1 - lmin), abs(lmax - 1))
        out['rho_opt_damped'] = (lmax - lmin) / (lmax + lmin)
    else:
        E = spla.LinearOperator((n, n), matvec=lambda x: x - S.apply(M @ x), dtype=float)
        ev = spla.eigs(E, k=3, which='LM', return_eigenvectors=False, tol=1e-6)
        out['rho_undamped'] = np.max(np.abs(ev))
    return out


def blocks_from_color(color, keep):
    col = color[keep]
    return [np.nonzero(col == c)[0] for c in np.unique(col) if c > 0]


def main():
    d = sys.argv[1]
    skip_cond = '--skip-cond' in sys.argv
    pre = 'lag_' if '--lag' in sys.argv else ''
    M, info, colors, cells = load(d, pre)
    if pre:
        Mf, finfo, newid = M.tocsr(), info, np.arange(M.shape[0])
    else:
        Mf, finfo, newid = restrict_free(M, info)
    keep = newid >= 0
    res = {'n': Mf.shape[0]}
    if not skip_cond:
        k, lmin, lmax = cond_sym(Mf)
        res['cond_M'] = k
        dg = 1.0 / np.sqrt(Mf.diagonal())
        Ds = sp.diags(dg)
        res['cond_diagM'] = cond_sym((Ds @ Mf @ Ds).tocsr())[0]
    for name, color in colors.items():
        blks = blocks_from_color(color, keep)
        S = Schwarz(Mf, blks)
        r = schwarz_metrics(Mf, S, True)
        # check entity alignment: any entity split across blocks?
        col = color[keep]
        key = finfo[:, 1] * 10**7 + finfo[:, 2]
        split = 0
        for kk in np.unique(key):
            if len(np.unique(col[key == kk])) > 1:
                split += 1
        r['split_entities'] = split
        res['bj_' + name] = r
    # Point-block Jacobi with one block per mesh entity (all DoFs on an entity, all fields):
    # measures coupling between entities only (invariant to any within-entity change of basis)
    key = finfo[:, 1] * 10**7 + finfo[:, 2]
    _, inv = np.unique(key, return_inverse=True)
    order_ = np.argsort(inv, kind='stable')
    eblks = np.split(order_, np.cumsum(np.bincount(inv))[:-1])
    res['bj_perentity'] = schwarz_metrics(Mf, Schwarz(Mf, eblks), True)
    for tag, (blks, asgs) in cells.items():
        bl = [newid[b][newid[b] >= 0] for b in blks]
        asg = [newid[a][newid[a] >= 0] for a in asgs]
        if tag == 'ov0' or True:
            Sa = Schwarz(Mf, bl, additive=True)
            res['as_' + tag] = schwarz_metrics(Mf, Sa, True)
        Sr = Schwarz(Mf, bl, assigned=asg)
        res['ras_' + tag] = schwarz_metrics(Mf, Sr, False)
    for k_, v in res.items():
        if isinstance(v, dict):
            print('%-12s' % k_, ' '.join('%s=%.4g' % (a, b) for a, b in v.items()))
        else:
            print('%-12s %.4g' % (k_, v))


if __name__ == '__main__':
    main()
