# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
"""Hierarchy tools for the H(Curl)_p + Grad(H^1_{p+1}) B space (phase 3).

- level(d): mass M, curl-curl K (padded), DoF info, colours, cell blocks, mesh, free-DoF restriction
- injection(fine, coarse): DoF-injection P_{p-1 -> p} on free DoFs (lower-order functions first per entity)
- grad_basis(L): sparse basis of G_p = grad(H^1_{p+1}) on free DoFs (Whitney-carried vertex gradients + Bg)
- gamma(M, U, W): strengthened Cauchy-Schwarz constant between span(U) and span(W) in the M inner product
"""
import os
import sys
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import load, restrict_free  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def read_mat(path, n=None):
    with open(path) as f:
        nr, nnz = map(int, f.readline().split())
        data = np.loadtxt(f)
    n = nr if n is None else n
    A = sp.csr_matrix((data[:, 2], (data[:, 0].astype(int) - 1, data[:, 1].astype(int) - 1)), shape=(n, n))
    return ((A + A.T) * 0.5).tocsr()


def read_mesh(d):
    with open(os.path.join(d, 'mesh.txt')) as f:
        np_, nc = map(int, f.readline().split())
        r = np.array([list(map(float, f.readline().split())) for _ in range(np_)])
        cells = np.array([list(map(int, f.readline().split())) for _ in range(nc)])
    lc = cells[:, 0:4] - 1
    lce = cells[:, 4:10] - 1
    oe = cells[:, 10:22].reshape(nc, 6, 2) - 1
    lcf = cells[:, 22:26] - 1
    return dict(r=r, lc=lc, lce=lce, oe=oe, lcf=lcf, np=np_, nc=nc)


class Level:
    def __init__(self, d):
        self.d = d
        M, info, colors, cells = load(d)
        n = M.shape[0]
        K = read_mat(os.path.join(d, 'curlcurl.txt'), n)
        self.Mfull = M
        self.M, self.info, self.newid = restrict_free(M, info)
        keep = self.newid >= 0
        self.idx = np.nonzero(keep)[0]
        self.K = K[self.idx][:, self.idx].tocsr()
        self.n = self.M.shape[0]
        self.color = colors['entity'][keep]
        self.cells = {}
        for tag, (blks, asgs) in cells.items():
            bl = [self.newid[b][self.newid[b] >= 0] for b in blks]
            asg = [self.newid[a][self.newid[a] >= 0] for a in asgs]
            self.cells[tag] = (bl, asg)
        self.mesh = read_mesh(d)
        with open(os.path.join(d, 'dofs.txt')) as f:
            f.readline()
            self.gs = [list(map(int, f.readline().split()))[1:] for _ in range(2)]
        # entity key per free DoF (dim, entity) and position k within (field, entity)
        self.ekey = self.info[:, 1] * 10**7 + self.info[:, 2]

    def entity_blocks(self):
        _, inv = np.unique(self.ekey, return_inverse=True)
        o = np.argsort(inv, kind='stable')
        return np.split(o, np.cumsum(np.bincount(inv))[:-1])

    def color_blocks(self):
        return [np.nonzero(self.color == c)[0] for c in np.unique(self.color) if c > 0]


def injection(F, C):
    """P (F.n x C.n): coarse free DoF (field, dim, entity, k) -> fine free DoF with the same tuple."""
    def keyof(info):
        return ((info[:, 0] * 4 + info[:, 1]) * 10**7 + info[:, 2]) * 100 + info[:, 3]
    kf = keyof(F.info)
    kc = keyof(C.info)
    pos = {k: i for i, k in enumerate(kf)}
    rows = np.array([pos[k] for k in kc])
    return sp.csr_matrix((np.ones(C.n), (rows, np.arange(C.n))), shape=(F.n, C.n))


def whitney_grad(L):
    """D (free DoFs x np): grad of vertex hats in the Whitney edge basis (sign from oriented local edges)."""
    m = L.mesh
    ne = L.gs[0][1] and m['lce'].max() + 1
    D = {}
    for c in range(m['nc']):
        for k in range(6):
            e = m['lce'][c, k]
            a, b = m['lc'][c, m['oe'][c, k]]
            D[(e, a)] = -1.0
            D[(e, b)] = 1.0
    # Whitney DoF of edge e is free DoF with field 0, dim 1, entity e, k 0
    sel = np.nonzero((L.info[:, 0] == 0) & (L.info[:, 1] == 1))[0]
    eid = np.empty(ne, dtype=int)
    eid[L.info[sel, 2]] = sel
    rows = [eid[e] for (e, v) in D]
    cols = [v for (e, v) in D]
    return sp.csr_matrix((list(D.values()), (rows, cols)), shape=(L.n, m['np']))


def grad_basis(L):
    D = whitney_grad(L)
    D = D[:, 1:]  # drop one vertex (constant has zero gradient)
    gsel = np.nonzero(L.info[:, 0] == 1)[0]
    E = sp.csr_matrix((np.ones(len(gsel)), (gsel, np.arange(len(gsel)))), shape=(L.n, len(gsel)))
    return sp.hstack([D, E]).tocsr()


def gamma(M, U, W):
    """cos of the smallest principal angle between span(U), span(W) in the M inner product (dense)."""
    U = U.toarray() if sp.issparse(U) else U
    W = W.toarray() if sp.issparse(W) else W
    MU = M @ U
    MW = M @ W
    Lu = np.linalg.cholesky(U.T @ MU)
    Lw = np.linalg.cholesky(W.T @ MW)
    C = U.T @ MW
    X = sla.solve_triangular(Lu, C, lower=True)
    X = sla.solve_triangular(Lw, X.T, lower=True).T
    return np.linalg.norm(X, 2)


def unit_cols(n, sel):
    return sp.csr_matrix((np.ones(len(sel)), (sel, np.arange(len(sel)))), shape=(n, len(sel)))


if __name__ == '__main__':
    # validation: P^T M_p P = M_{p-1}, K G = 0, dim G
    for mesh, ps in (('tiny', (1, 2, 3, 4)), ('ref', (1, 2, 3))):
        prev = None
        for p in ps:
            L = Level(os.path.join(HERE, 'dump3', '%s_p%d' % (mesh, p)))
            G = grad_basis(L)
            KG = abs(L.K @ G).max() / abs(L.K).max()
            nh1 = L.gs[1][0] * L.mesh['np'] + int(np.sum(L.info[:, 0] == 1))
            msg = '%s p=%d n=%d  |KG|/|K|=%.1e  dimG=%d (H1_{p+1} - 1 = %d)' % (mesh, p, L.n, KG, G.shape[1], nh1 - 1)
            if prev is not None:
                P = injection(L, prev)
                Mc = (P.T @ L.M @ P).tocsr()
                Kc = (P.T @ L.K @ P).tocsr()
                msg += '  |PtMP-Mc|/|Mc|=%.1e |PtKP-Kc|/|Kc|=%.1e' % (abs(Mc - prev.M).max() / abs(prev.M).max(),
                                                                      abs(Kc - prev.K).max() / abs(prev.K).max())
            print(msg, flush=True)
            prev = L
