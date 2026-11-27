#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only: conditioning of OFT's hexahedral H^1 and H(Curl)
# bases (canonical orientation, unit cube), and comparisons:
#   - H^1: OFT (x J^(1,1) "hpoly_bary", i.e. integrated Legendre / Lobatto) vs the
#     older J^(0,0) definition kept as a comment in lagrange_poly.F90, vs a nodal
#     Gauss-Lobatto-Legendre (GLL) Q_k basis (non-hierarchical);
#   - H(Curl): OFT vs a per-entity, per-level L2-orthogonalized variant of the SAME
#     spaces (each hierarchical level of a face/cell projected off the lower levels of
#     that entity and, for H(Curl), off the same entity's H^1 gradients). The cube's
#     symmetries are isometries and the subtracted sets are invariant, so this keeps the
#     signed-permutation orientation structure, the hierarchy and the gradient-free
#     direct sum. It is a headroom indicator for basis choice within OFT's spaces.
# Transcribed from oft_h1_eval / oft_hcurl_eval (hex branches) and
# src/grid/lagrange_poly.F90, src/grid/hexmesh_type.F90.
# Usage: python hex_conditioning.py   (writes results_hex.txt)
#------------------------------------------------------------------------------
import os
import numpy as np
import sympy as sp

X = sp.symbols('x y z')
OUT = []


def log(msg=''):
    print(msg)
    OUT.append(msg)


# hexmesh_type.F90
HEX_BARY_MAP = (-3, -2, 1, 2, -1, 3)
PFC = [(1, 2, 5), (1, 2, 3), (1, 3, 4), (1, 4, 5), (2, 5, 6), (2, 3, 6), (3, 4, 6), (4, 5, 6)]
EFC = [(1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (2, 3), (3, 4), (4, 5), (2, 6), (3, 6), (4, 6), (5, 6)]
EC = [(5, 3), (2, 4), (3, 5), (4, 2), (1, 6), (1, 6), (1, 6), (1, 6), (5, 3), (2, 4), (3, 5), (4, 2)]
FC = [(5, 2, 3, 4), (1, 5, 6, 3), (1, 2, 6, 4), (1, 3, 6, 5), (2, 1, 4, 6), (2, 5, 4, 3)]

FH = [1 - X[abs(m)-1] if m < 0 else X[abs(m)-1] for m in HEX_BARY_MAP]   # fhex(1..6) -> FH[0..5]


def fh(i):
    return FH[i-1]


def gb(i):
    return sp.Matrix([sp.diff(fh(i), v) for v in X])


def hpoly(ii, x, variant='j11'):
    if variant == 'j11':
        return [1, x, x*(4*x - 2), x*(15*x**2 - 15*x + 3), x*(56*x**3 - 84*x**2 + 36*x - 4),
                x*(210*x**4 - 420*x**3 + 280*x**2 - 70*x + 5)][ii]
    return [1, x, x*(2*x - 1), x*(6*x**2 - 6*x + 1), x*(20*x**3 - 30*x**2 + 12*x - 1),
            x*(70*x**4 - 140*x**3 + 90*x**2 - 20*x + 1)][ii]


def grid2(order):
    out = []
    for k in range(1, order+1):
        for i in range(1, k+1):
            for j in range(1, k+1):
                if i < k and j < k:
                    continue
                out.append((i, j, k))
    return out


def grid3(order):
    out = []
    for k in range(1, order+1):
        for i in range(1, k+1):
            for j in range(1, k+1):
                for l in range(1, k+1):
                    if i < k and j < k and l < k:
                        continue
                    out.append((i, j, l, k))
    return out


def hc_grid2(order):
    out = []
    for k in range(1, order+1):
        for i in range(1, k+1):
            for j in range(1, k+1):
                if i < k and j < k:
                    continue
                if j == 1:
                    out.append((i, i, 1, k))
                    out.append((i, i, 2, k))
                out.append((i, j, 0, k))
    return out


def hc_grid3(order):
    out = []
    for k in range(1, order+1):
        for i in range(1, k+1):
            for j in range(1, k+1):
                for l in range(1, k+1):
                    if i < k and j < k and l < k:
                        continue
                    if l == 1:
                        out += [(0, i, j, 2, k), (i, 0, j, 3, k), (i, j, 0, 4, k)]
                    out.append((i, j, l, 0, k))
                    out.append((i, j, l, 1, k))
    return out


def h1_hex(order, variant='j11'):
    """List of (entity, level, scalar expr). Entities 'v', 'e<i>', 'f<i>', 'c'."""
    out = []
    for i in range(8):
        out.append(('v', 0, sp.Mul(*[fh(c) for c in PFC[i]])))
    for i in range(12):
        e1, e2 = EC[i]
        for j in range(1, order):
            out.append(('e%d' % i, j, hpoly(j, fh(e1), variant)*fh(e2)*fh(EFC[i][0])*fh(EFC[i][1])))
    for i in range(6):
        c = FC[i]
        for (a, b, k) in grid2(order-1):
            out.append(('f%d' % i, k, hpoly(a, fh(c[0]), variant)*hpoly(b, fh(c[1]), variant)*fh(i+1)*fh(c[2])*fh(c[3])))
    for (a, b, d, k) in grid3(order-1):
        out.append(('c', k, hpoly(a, fh(1), variant)*hpoly(b, fh(2), variant)*hpoly(d, fh(3), variant)*fh(4)*fh(5)*fh(6)))
    return out


_T = sp.Symbol('t')


def dh(ii, x):
    return sp.diff(hpoly(ii, _T), _T).subs(_T, x)


def hc_hex(order):
    """List of (entity, level, vector expr) following oft_hcurl_eval (hex, canonical orientation)."""
    out = []
    for i in range(12):
        e1, e2 = EC[i]
        out.append(('e%d' % i, 0, gb(e2)*fh(EFC[i][0])*fh(EFC[i][1])))
    for i in range(6):
        c = FC[i]
        fi = fh(i+1)
        c1, c2, c3, c4 = [fh(v) for v in c]
        for (a, b, form, k) in hc_grid2(order-1):
            if form == 0:
                v1 = hpoly(a, c1)*hpoly(b, c2)
                v2 = dh(a, c1)*hpoly(b, c2)
                v3 = hpoly(a, c1)*dh(b, c2)
                val = gb(c[0])*(v2*c3 - v1)*c4*fi - gb(c[1])*(v3*c4 - v1)*c3*fi
            elif form == 1:
                val = gb(c[0])*hpoly(b, c2)*c4*fi
            else:
                val = gb(c[1])*hpoly(a, c1)*c3*fi
            out.append(('f%d' % i, k, val))
    f1, f2, f3, f4, f5, f6 = FH
    for (a, b, d, form, k) in hc_grid3(order-1):
        if form < 2:
            v1 = hpoly(a, f1)*hpoly(b, f2)*hpoly(d, f3)
            v2 = dh(a, f1)*hpoly(b, f2)*hpoly(d, f3)
            v3 = hpoly(a, f1)*dh(b, f2)*hpoly(d, f3)
            v4 = hpoly(a, f1)*hpoly(b, f2)*dh(d, f3)
            g1 = gb(1)*(v2*f6 - v1)*f4*f5
            g2 = gb(2)*(v3*f4 - v1)*f5*f6
            g3 = gb(3)*(v4*f5 - v1)*f4*f6
            val = g1 - g2 + g3 if form == 0 else g1 - g2 - g3
        elif form == 2:
            val = gb(1)*hpoly(b, f2)*hpoly(d, f3)*f4*f5
        elif form == 3:
            val = gb(2)*hpoly(a, f1)*hpoly(d, f3)*f5*f6
        else:
            val = gb(3)*hpoly(a, f1)*hpoly(b, f2)*f4*f6
        out.append(('c', k, val))
    return out


def gll_lagrange(order):
    """Nodal Q_k basis on GLL points (tensor product), for comparison only."""
    if order == 1:
        nodes = np.array([0., 1.])
    else:
        from numpy.polynomial import legendre as L
        c = np.zeros(order+1)
        c[-1] = 1.
        inner = np.sort(np.real(L.legroots(L.legder(c))))
        nodes = np.concatenate([[-1.], inner, [1.]])*0.5 + 0.5
    nodes = [sp.Float(v, 17) for v in nodes]
    def lag1(i, x):
        e = sp.Integer(1)
        for j, xj in enumerate(nodes):
            if j != i:
                e *= (x - xj)/(nodes[i] - xj)
        return e
    out = []
    for i in range(order+1):
        for j in range(order+1):
            for k in range(order+1):
                out.append(('n', 0, lag1(i, X[0])*lag1(j, X[1])*lag1(k, X[2])))
    return out


#------------------------------------------------------------------------------
# Numerics on the unit cube (tensor Gauss-Legendre, exact for these degrees)
#------------------------------------------------------------------------------
def cube_quad(n=9):
    x, w = np.polynomial.legendre.leggauss(n)
    x = 0.5*(x + 1.)
    w = 0.5*w
    P = np.array(np.meshgrid(x, x, x, indexing='ij')).reshape(3, -1)
    W = np.einsum('i,j,k->ijk', w, w, w).ravel()
    return P, W


def ev(expr, P):
    f = sp.lambdify(X, expr, 'numpy')
    return np.broadcast_to(np.asarray(f(*P), dtype=float), (P.shape[1],)).copy()


def scal_samples(funcs, P):
    V = np.array([ev(e, P) for _, _, e in funcs])                                     # nb, nq
    G = np.array([[ev(sp.diff(e, v), P) for v in X] for _, _, e in funcs])          # nb, 3, nq
    return V, G


def vec_samples(funcs, P):
    V = np.array([[ev(e[c], P) for c in range(3)] for _, _, e in funcs])
    C = []
    for _, _, e in funcs:
        cx = sp.diff(e[2], X[1]) - sp.diff(e[1], X[2])
        cy = sp.diff(e[0], X[2]) - sp.diff(e[2], X[0])
        cz = sp.diff(e[1], X[0]) - sp.diff(e[0], X[1])
        C.append([ev(cx, P), ev(cy, P), ev(cz, P)])
    return V, np.array(C)


def gram(A, B, W):
    if A.ndim == 2:
        return np.einsum('iq,jq,q->ij', A, B, W)
    return np.einsum('icq,jcq,q->ij', A, B, W)


def kappa(M, nzero=0):
    d = np.sqrt(np.abs(np.diag(M)))
    d[d == 0.] = 1.
    Ms = M/np.outer(d, d)
    e = np.sort(np.linalg.eigvalsh(0.5*(Ms + Ms.T)))[nzero:]
    return e[-1]/e[0]


def top_kappa(M, funcs, nzero=0):
    """Jacobi-scaled eigenvalue ratio of the block of the newest hierarchical level
    (what a Jacobi smoother sees in identity-injection p-multigrid)."""
    top = max(f[1] for f in funcs)
    idx = [i for i, f in enumerate(funcs) if f[1] == top and f[0] != 'v']
    if top == 0 or not idx:
        return float('nan')
    return kappa(M[np.ix_(idx, idx)], nzero)


def orthogonalize(funcs, V, D, W, extra=None):
    """Per entity and level: project each level off the lower levels of the same entity
    (and off extra[entity] = list of sample arrays, eg. same-entity H^1 gradients)."""
    V = V.copy()
    D = D.copy()
    ents = [f[0] for f in funcs]
    lev = [f[1] for f in funcs]
    for ent in sorted(set(ents)):
        idx = [i for i in range(len(funcs)) if ents[i] == ent]
        for L in sorted(set(lev[i] for i in idx)):
            cur = [i for i in idx if lev[i] == L]
            low = [i for i in idx if lev[i] < L]
            Bv = [V[i] for i in low]
            Bd = [D[i] for i in low]
            if extra is not None and ent in extra:
                Bv += list(extra[ent][0])
                Bd += list(extra[ent][1])
            if not Bv:
                continue
            Bv = np.array(Bv)
            Bd = np.array(Bd)
            G = gram(Bv, Bv, W)
            for i in cur:
                c = np.linalg.solve(G, gram(Bv, V[i][None], W)[:, 0])
                V[i] = V[i] - np.tensordot(c, Bv, axes=1)
                D[i] = D[i] - np.tensordot(c, Bd, axes=1)
    return V, D


if __name__ == '__main__':
    P, W = cube_quad(9)
    log('## Hex H^1 (orders 2-5): Jacobi-scaled eigenvalue ratio of mass / stiffness (1 zero mode removed)')
    for k in range(2, 6):
        row = []
        for lab, fn in (('OFT J(1,1)', h1_hex(k)), ('old J(0,0)', h1_hex(k, 'j00')), ('GLL nodal', gll_lagrange(k))):
            V, G = scal_samples(fn, P)
            Mm, Kk = gram(V, V, W), gram(G, G, W)
            if lab == 'GLL nodal':
                row.append('%s M %.3e K %.3e' % (lab, kappa(Mm), kappa(Kk, 1)))
            else:
                row.append('%s M %.3e K %.3e (top level M %.3e K %.3e)' % (
                    lab, kappa(Mm), kappa(Kk, 1), top_kappa(Mm, fn), top_kappa(Kk, fn)))
        Vh, Gh = scal_samples(h1_hex(k), P)
        Vo, Go = orthogonalize(h1_hex(k), Vh, Gh, W)
        row.append('OFT level-orth M %.3e K %.3e' % (kappa(gram(Vo, Vo, W)), kappa(gram(Go, Go, W), 1)))
        log('  k=%d (n=%d): ' % (k, len(fn)) + ' | '.join(row))
    log()
    log('## Hex H(Curl) (p=1-4) with grad H^1 (order p+1): Jacobi-scaled eigenvalue ratios')
    log('   (HC mass / HC curl-curl (Whitney-gradient zero modes removed) / combined HC + grad H^1 non-vertex mass)')
    for p in range(1, 5):
        hc = hc_hex(p)
        h1 = [f for f in h1_hex(p+1) if f[0] != 'v']
        V, C = vec_samples(hc, P)
        _, G = scal_samples(h1, P)
        nzero = 7    # grad of the 8 trilinear vertex functions inside the Whitney space
        def stats(Vh, Ch):
            Mc = gram(np.concatenate([Vh, G]), np.concatenate([Vh, G]), W)
            r = np.linalg.matrix_rank(np.concatenate([Vh, G]).reshape(len(Vh)+len(G), -1), tol=1e-9)
            Mh, Ch2 = gram(Vh, Vh, W), gram(Ch, Ch, W)
            return kappa(Mh), kappa(Ch2, nzero), kappa(Mc), r, top_kappa(Mh, hc), top_kappa(Ch2, hc)
        s_old = stats(V, C)
        # same-entity H^1 gradients for the projection
        extra = {}
        for ent in set(f[0] for f in hc):
            gi = [i for i, f in enumerate(h1) if f[0] == ent]
            if gi:
                extra[ent] = (G[gi], np.zeros_like(G[gi]))
        Vo, Co = orthogonalize(hc, V, C, W)
        s_lev = stats(Vo, Co)
        Vg, Cg = orthogonalize(hc, V, C, W, extra=extra)
        s_grd = stats(Vg, Cg)
        log('  p=%d (nHC=%d, nH1 non-vertex=%d): OFT %.3e / %.3e / %.3e | level-orth %.3e / %.3e / %.3e | '
            'level+grad-orth %.3e / %.3e / %.3e | combined rank %d/%d/%d of %d' % (
                p, len(hc), len(h1), *s_old[:3], *s_lev[:3], *s_grd[:3], s_old[3], s_lev[3], s_grd[3],
                len(hc)+len(h1)))
        log('        top-level block (HC mass / curl-curl): OFT %.3e / %.3e | level-orth %.3e / %.3e | '
            'level+grad-orth %.3e / %.3e' % (s_old[4], s_old[5], s_lev[4], s_lev[5], s_grd[4], s_grd[5]))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_hex.txt'), 'w') as fid:
        fid.write('\n'.join(OUT)+'\n')
