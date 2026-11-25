# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
"""Phase 3 steps 3-4: smoother / two-grid metrics for re-based hierarchical bases.

Bases (columns of T in the current free-DoF coefficients, A_T = T^T A T, DoF j keeps its entity):
  i    : current OFT basis (T = I)
  iih  : new level-p functions globally A0-orthogonalized against V_{p-1}           (recursive from p=2)
  iis  : structure preserving: new grads _|_ G_{p-1}; new rotationals _|_ V_{p-1} + G_p   (recursive)
  iisl : patch-local version of iis (projection restricted to the cell-closure patch of each entity)
  iii  : orthonormalization within each (entity, field) block of the current basis
A0 = M (orthogonalization always in the mass inner product); A = M + tau K for the metrics.

Usage: python rebase.py <mesh> <pmax> [cfl ...]
"""
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hier  # noqa: E402

DENSE_MAX = 9000


def dense(A):
    return A.toarray() if sp.issparse(A) else np.asarray(A)


def sub(A, b):
    return dense(A[np.ix_(b, b)]) if not sp.issparse(A) else A[b][:, b].toarray()


# ----------------------------------------------------------------------------- transforms
def orth_global(M, S, Ecols):
    """Columns E - S (S^T M S)^{-1} S^T M E (dense)."""
    S = dense(S)
    E = dense(Ecols)
    MS = M @ S
    G = S.T @ MS
    X = sla.cho_solve(sla.cho_factor(G), MS.T @ E)
    return E - S @ X


def entity_cells(L):
    m = L.mesh
    nc = m['nc']
    cell_ent = []
    for c in range(nc):
        ents = [(0, v) for v in m['lc'][c]] + [(1, e) for e in m['lce'][c]] + \
               [(2, f) for f in m['lcf'][c]] + [(3, c)]
        cell_ent.append(ents)
    ent2cells = {}
    for c, ents in enumerate(cell_ent):
        for de in ents:
            ent2cells.setdefault(de, []).append(c)
    # free DoFs per entity
    ent_dofs = {}
    for j, (fi, dm, e, k) in enumerate(L.info):
        ent_dofs.setdefault((dm, e), []).append(j)
    return cell_ent, ent2cells, ent_dofs


def orth_patch(L, M, S, new):
    """Patch-local version: for the new DoFs on each entity, project out the span of the S columns that
    intersect the closure Q of the entity's cell patch, truncated to Q, using M restricted to Q.
    The correction is supported in Q (bounded fill); S columns live on rows other than `new`."""
    cell_ent, ent2cells, ent_dofs = entity_cells(L)
    S = sp.csc_matrix(S)
    n = L.n
    Sb = S.copy()
    Sb.data[:] = 1.0
    colnnz = np.asarray(Sb.sum(axis=0)).ravel()
    SbT = Sb.T.tocsr()
    newset = np.zeros(n, dtype=bool)
    newset[new] = True
    rows, cols, vals = [], [], []
    for (dm, e), dofs in ent_dofs.items():
        dofs = np.array(dofs)
        jn = dofs[newset[dofs]]
        if len(jn) == 0:
            continue
        Q = set()
        for c in ent2cells[(dm, e)]:
            for de in cell_ent[c]:
                Q.update(ent_dofs.get(de, []))
        Q = np.array(sorted(Q))
        q = np.zeros(n)
        q[Q] = 1.0
        inside = np.nonzero(SbT @ q > 0.5)[0]
        if len(inside) == 0:
            X = np.zeros((len(Q), len(jn)))
        else:
            Sq = S[:, inside][Q].toarray()
            MQ = M[Q][:, Q].toarray()
            MS = MQ @ Sq
            Gm = Sq.T @ MS
            Ej = np.zeros((len(Q), len(jn)))
            pos = np.searchsorted(Q, jn)
            Ej[pos, np.arange(len(jn))] = 1.0
            # rank-revealing solve (patch S may be rank deficient)
            try:
                Y = np.linalg.lstsq(Gm, MS.T @ Ej, rcond=1e-12)[0]
            except np.linalg.LinAlgError:
                # gelsd SVD failed to converge: fall back to complete orthogonal factorization
                Y = sla.lstsq(Gm, MS.T @ Ej, cond=1e-12, lapack_driver='gelsy')[0]
            X = Ej - Sq @ Y
        X[np.abs(X) < 1e-15] = 0.0
        r, cc = np.nonzero(X)
        rows.append(Q[r]); cols.append(jn[cc]); vals.append(X[r, cc])
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))


def build_T(L, C, Tc, kind, glob):
    """Transform at level L given coarse level C with its transform Tc (None for p=1)."""
    n = L.n
    if C is None or kind == 'i':
        return sp.identity(n, format='csr')
    P = hier.injection(L, C)
    coarse_mask = np.asarray(P.sum(axis=1)).ravel() > 0
    new = np.nonzero(~coarse_mask)[0]
    Ccols = P @ Tc  # coarse functions (already re-based) in fine coefficients
    if glob:
        Ccols = dense(Ccols)
    else:
        Ccols = sp.csr_matrix(Ccols)
    M = L.M
    T = sp.lil_matrix((n, n)) if not glob else np.zeros((n, n))
    if glob:
        T[:, coarse_mask] = Ccols @ dense(P.T)[:, coarse_mask] if False else 0.0
        # place coarse columns at their fine DoF index
        cidx = np.asarray(P.argmax(axis=0)).ravel()
        T[:, cidx] = Ccols
    else:
        cidx = np.asarray(P.argmax(axis=0)).ravel()
    if kind == 'iih':
        if glob:
            T[:, new] = orth_global(M, Ccols, hier.unit_cols(n, new))
        else:
            X = orth_patch(L, M, Ccols, new)
    elif kind == 'iis':
        ng = new[L.info[new, 0] == 1]
        nr = new[L.info[new, 0] == 0]
        Gc = P @ (hier.grad_basis(C))
        if glob:
            Xg = orth_global(M, Gc, hier.unit_cols(n, ng))
            T[:, ng] = Xg
            if len(nr):
                T[:, nr] = orth_global(M, np.hstack([Ccols, Xg]), hier.unit_cols(n, nr))
        else:
            Xg = orth_patch(L, M, Gc, ng)
            Xr = orth_patch(L, M, sp.hstack([Ccols, Xg[:, ng]]).tocsc(), nr) if len(nr) else None
            X = Xg if Xr is None else Xg + Xr
    else:
        raise ValueError(kind)
    if glob:
        return T
    Cpart = sp.csr_matrix((n, n))
    Cpart = sp.csr_matrix(Ccols)  # n x nc
    perm = sp.csr_matrix((np.ones(len(cidx)), (np.arange(len(cidx)), cidx)), shape=(len(cidx), n))
    return (Cpart @ perm + X).tocsr()


def block_orthonormal(L, A0, coarse_mask):
    """(iii): T block diagonal over (entity, field, coarse/new) with T_b^T A0_bb T_b = I
    (keeps the coarse span, so two-grid metrics stay comparable)."""
    key = (L.ekey * 2 + L.info[:, 0]) * 2 + coarse_mask.astype(int)
    _, inv = np.unique(key, return_inverse=True)
    rows, cols, vals = [], [], []
    for b in range(inv.max() + 1):
        idx = np.nonzero(inv == b)[0]
        Lb = np.linalg.cholesky(A0[idx][:, idx].toarray())
        Tb = sla.solve_triangular(Lb, np.eye(len(idx)), lower=True).T
        r, c = np.nonzero(np.abs(Tb) > 0)
        rows.append(idx[r]); cols.append(idx[c]); vals.append(Tb[r, c])
    n = L.n
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))


# ----------------------------------------------------------------------------- metrics
class Smoother:
    def __init__(self, A, blocks, assigned=None):
        self.n = A.shape[0]
        self.blocks = blocks
        self.facs = [sla.cho_factor(sub(A, b)) for b in blocks]
        self.w = [np.ones(len(b), bool) if assigned is None else np.isin(b, assigned[i]) for i, b in enumerate(blocks)]

    def apply(self, r):
        z = np.zeros(self.n)
        for b, f, w in zip(self.blocks, self.facs, self.w):
            z[b[w]] += sla.cho_solve(f, r[b])[w]
        return z


class Coarse:
    def __init__(self, A, idx):
        self.idx = idx
        self.n = A.shape[0]
        Ac = sub(A, idx)
        self.f = sla.cho_factor(Ac)

    def apply(self, r):
        z = np.zeros(self.n)
        z[self.idx] = sla.cho_solve(self.f, r[self.idx])
        return z


class SubspaceSolve:
    """Exact A-solve on span(Z) (Z given in the same coordinates as A)."""
    def __init__(self, A, Z):
        self.Z = Z
        self.f = sla.cho_factor(dense(Z.T @ (A @ Z)))

    def apply(self, r):
        return self.Z @ sla.cho_solve(self.f, self.Z.T @ r)


def spd_kappa(A, Binv_apply):
    n = A.shape[0]
    BA = spla.LinearOperator((n, n), matvec=lambda x: Binv_apply(A @ x), dtype=float)
    lmax = np.max(spla.eigs(BA, k=1, which='LR', return_eigenvectors=False, tol=1e-6, maxiter=50000).real)
    sh = 1.05 * lmax
    Op = spla.LinearOperator((n, n), matvec=lambda x: sh * x - Binv_apply(A @ x), dtype=float)
    lmin = sh - np.max(spla.eigs(Op, k=1, which='LR', return_eigenvectors=False, tol=1e-6, maxiter=50000).real)
    return lmax / lmin, lmin, lmax


def spectral_radius(n, E):
    op = spla.LinearOperator((n, n), matvec=E, dtype=float)
    ev = spla.eigs(op, k=2, which='LM', return_eigenvectors=False, tol=1e-6, maxiter=50000)
    return np.max(np.abs(ev))


def metrics(L, A, coarse_idx, Gz=None):
    n = A.shape[0]
    out = {}
    d = 1.0 / np.asarray(A.diagonal() if sp.issparse(A) else np.diag(A)).ravel()
    out['k_diag'] = spd_kappa(A, lambda r: d * r)[0]
    Sbj = Smoother(A, L.color_blocks())
    k, lmin, lmax = spd_kappa(A, Sbj.apply)
    out['k_bj'] = k
    out['rho_bj'] = max(abs(1 - lmin), abs(lmax - 1))
    Sent = Smoother(A, L.entity_blocks())
    out['k_ent'] = spd_kappa(A, Sent.apply)[0]
    if coarse_idx is not None:
        Cc = Coarse(A, coarse_idx)
        # two-level additive (coarse + entity-aligned BJ)
        out['k_2lev'] = spd_kappa(A, lambda r: Cc.apply(r) + Sbj.apply(r))[0]
        w = 1.0 / lmax  # damped BJ (convergent), symmetric two-grid V(1,1)
        sm = lambda y: y - w * Sbj.apply(A @ y)  # noqa: E731
        out['rho_2g_bjw'] = spectral_radius(n, lambda x: sm(x - Cc.apply(A @ sm(x))))
        new_idx = np.setdiff1d(np.arange(n), coarse_idx)
        Acn = A[np.ix_(coarse_idx, new_idx)] if not sp.issparse(A) else A[coarse_idx][:, new_idx]
        out['xc'] = abs(Acn).max() / abs(A).max()
        bl, asg = L.cells['ov1p']
        Sras = Smoother(A, bl, asg)
        out['rho_2g_ras'] = spectral_radius(n, lambda x: (lambda y: y - Sras.apply(A @ y))(x - Cc.apply(A @ x)))
        if Gz is not None:
            Gs = SubspaceSolve(A, Gz)
            out['rho_2g_ras_G'] = spectral_radius(
                n, lambda x: (lambda y: y - Sras.apply(A @ y))((lambda w: w - Gs.apply(A @ w))(x - Cc.apply(A @ x))))
    return out


def transformed(A, T):
    if sp.issparse(T):
        return (T.T @ A @ T).tocsr()
    return T.T @ (A @ T)


def main():
    mesh = sys.argv[1]
    pmax = int(sys.argv[2])
    cfls = [float(c) for c in sys.argv[3:]] or [0.0]
    kinds = sys.argv[3:] and None
    levels = []
    Ts = {'i': None, 'iih': None, 'iis': None, 'iisl': None}
    for p in range(1, pmax + 1):
        L = hier.Level(os.path.join(hier.HERE, 'dump3', '%s_p%d' % (mesh, p)))
        C = levels[-1] if levels else None
        glob_ok = L.n <= DENSE_MAX
        newT = {}
        for kind in ('i', 'iih', 'iis', 'iisl'):
            base = 'iis' if kind == 'iisl' else kind
            glob = kind in ('iih', 'iis')
            if glob and not glob_ok:
                continue
            Tc = Ts.get(kind)
            if C is not None and Tc is None and kind != 'i':
                continue
            t0 = time.time()
            if C is None:
                newT[kind] = sp.identity(L.n, format='csr') if not glob else np.eye(L.n)
            else:
                newT[kind] = build_T(L, C, Tc, base, glob)
            print('  built %s p=%d in %.1fs' % (kind, p, time.time() - t0), flush=True)
        Ts = newT
        levels.append(L)
        if (p == 1 and pmax > 1) or p < int(os.environ.get('PMIN', '1')):
            continue
        coarse_idx = None
        if C is not None:
            P = hier.injection(L, C)
            coarse_idx = np.asarray(P.argmax(axis=0)).ravel()
        h = np.mean([np.linalg.norm(L.mesh['r'][a] - L.mesh['r'][b])
                     for c in range(L.mesh['nc']) for a, b in L.mesh['lc'][c][[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]]])
        G = hier.grad_basis(L)
        for cfl in cfls:
            tau = (cfl * h)**2
            A = (L.M + tau * L.K).tocsr()
            for kind in ('i', 'iii', 'iih', 'iis', 'iisl'):
                if kind == 'iii':
                    cm = np.zeros(L.n, bool)
                    if coarse_idx is not None:
                        cm[coarse_idx] = True
                    T = block_orthonormal(L, L.M, cm)
                elif kind in Ts:
                    T = Ts[kind]
                else:
                    continue
                t0 = time.time()
                AT = transformed(A, T)
                if sp.issparse(AT) and AT.shape[0] <= DENSE_MAX and kind in ('iih', 'iis'):
                    AT = AT.toarray()
                Gz = None
                if cfl > 0:
                    if sp.issparse(T):
                        Gz = spla.splu(T.tocsc()).solve(G.toarray())
                    else:
                        Gz = np.linalg.solve(T, G.toarray())
                m = metrics(L, AT, coarse_idx, Gz)
                print('%s p=%d cfl=%g %-5s ' % (mesh, p, cfl, kind) + ' '.join('%s=%.4g' % kv for kv in m.items()) +
                      '  (%.0fs)' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
