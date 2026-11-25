# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
"""Vertex-star (Arnold-Falk-Winther) patch metrics on the B mass block, for comparison with
entity/colour blocks (suggested by the smoothers study). Star(v) = free DoFs on all entities whose
closure contains vertex v (v, edges, faces and cells containing v).

Usage: python star.py <mesh> <pmax> [cfl ...]
"""
import os
import sys
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hier  # noqa: E402
import rebase  # noqa: E402
from analyze import load  # noqa: E402

TET_FC = [[1, 2, 3], [2, 0, 3], [0, 1, 3], [0, 1, 2]]  # face k opposite local vertex k (0-based)


def vertex_entities(mesh):
    ents = {}
    for c in range(mesh['nc']):
        lc = mesh['lc'][c]
        for v in lc:
            s = ents.setdefault(v, set())
            s.add((0, v))
            s.add((3, c))
        for k in range(6):
            for v in lc[mesh['oe'][c, k]]:
                ents[v].add((1, mesh['lce'][c, k]))
        for k in range(4):
            for v in lc[TET_FC[k]]:
                ents[v].add((2, mesh['lcf'][c, k]))
    return ents


def star_blocks(info, mesh):
    ents = vertex_entities(mesh)
    ent_dofs = {}
    for j, row in enumerate(info):
        ent_dofs.setdefault((row[1], row[2]), []).append(j)
    blocks = []
    for v in sorted(ents):
        b = sorted(j for e in ents[v] for j in ent_dofs.get(e, []))
        if b:
            blocks.append(np.array(b))
    return blocks


def star_metrics(A, blocks, coarse_idx=None):
    S = rebase.Smoother(A, blocks)
    k, lmin, lmax = rebase.spd_kappa(A, S.apply)
    out = {'k_star': k}
    n = A.shape[0]
    w = 1.0 / lmax
    sm = lambda y: y - w * S.apply(A @ y)  # noqa: E731
    out['rho_star_w'] = rebase.spectral_radius(n, sm)
    if coarse_idx is not None:
        Cc = rebase.Coarse(A, coarse_idx)
        out['rho_2g_star_w'] = rebase.spectral_radius(n, lambda x: sm(x - Cc.apply(A @ sm(x))))
    return out


def main():
    mesh = sys.argv[1]
    pmax = int(sys.argv[2])
    cfls = [float(c) for c in sys.argv[3:]] or [0.0]
    prev = None
    Ts = {}
    for p in range(1, pmax + 1):
        d = os.path.join(hier.HERE, 'dump3', '%s_p%d' % (mesh, p))
        L = hier.Level(d)
        glob_ok = L.n <= rebase.DENSE_MAX
        newT = {}
        for kind in ('iis', 'iisl'):
            glob = kind == 'iis'
            if glob and not glob_ok:
                continue
            if prev is None:
                newT[kind] = np.eye(L.n) if glob else sp.identity(L.n, format='csr')
            elif kind in Ts:
                newT[kind] = rebase.build_T(L, prev, Ts[kind], 'iis', glob)
        Ts = newT
        C = prev
        prev = L
        coarse_idx = None
        if C is not None:
            coarse_idx = np.asarray(hier.injection(L, C).argmax(axis=0)).ravel()
        blocks = star_blocks(L.info, L.mesh)
        h = np.mean([np.linalg.norm(L.mesh['r'][a] - L.mesh['r'][b]) for c in range(L.mesh['nc'])
                     for a, b in L.mesh['lc'][c][[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]]])
        for cfl in cfls:
            A = (L.M + (cfl * h)**2 * L.K).tocsr()
            for kind in ('i', 'iis', 'iisl'):
                if kind == 'i':
                    AT = A
                elif kind in Ts and not (p == 1):
                    AT = rebase.transformed(A, Ts[kind])
                else:
                    continue
                m = star_metrics(AT, blocks, coarse_idx)
                print('%s p=%d cfl=%g %-5s nstar=%d ' % (mesh, p, cfl, kind, len(blocks)) +
                      ' '.join('%s=%.4g' % kv for kv in m.items()), flush=True)
        # vector Lagrange baseline (same order, mass only), star blocks on its DoFs
        Ml, infol, _, _ = load(d, 'lag_')
        bl = star_blocks(infol, L.mesh)
        m = star_metrics(Ml.tocsr(), bl)
        print('%s p=%d cfl=0 lag   nstar=%d ' % (mesh, p, len(bl)) + ' '.join('%s=%.4g' % kv for kv in m.items()),
              flush=True)


if __name__ == '__main__':
    main()
