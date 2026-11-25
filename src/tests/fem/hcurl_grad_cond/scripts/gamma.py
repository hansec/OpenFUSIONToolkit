# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
"""Phase 3 step 2: strengthened Cauchy-Schwarz constants, global (assembled M) and local (element M).

Splittings of V_p = H(Curl)_p + Grad(H^1_{p+1}) (free DoFs):
  a : V_{p-1} (injected)            vs  new level-p functions
  b1: Bc (H(Curl) field DoFs)       vs  Bg (H^1 field DoFs)
  b2: R (Bc with the Whitney-carried vertex gradients removed; Euclidean complement)  vs  G_p = grad(H^1_{p+1})
  a_r/a_g: within-family hierarchy, Bc_{p-1} vs new Bc, G_{p-1} vs new grads (Bg)
"""
import os
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hier  # noqa: E402


def gamma_sparse(M, U, W):
    A11 = (U.T @ M @ U).tocsc()
    A22 = (W.T @ M @ W).tocsc()
    A12 = (U.T @ M @ W).tocsr()
    l1 = spla.splu(A11)
    l2 = spla.splu(A22)
    n = A11.shape[0]
    op = spla.LinearOperator((n, n), matvec=lambda x: l1.solve(A12 @ l2.solve(A12.T @ x)), dtype=float)
    ev = spla.eigs(op, k=1, which='LR', return_eigenvectors=False, tol=1e-8, maxiter=20000)
    return np.sqrt(max(ev.real.max(), 0.0))


def whitney_complement(L):
    """R basis: Whitney edge coefficient space minus span(D) (Euclidean complement), plus other Bc DoFs."""
    D = hier.whitney_grad(L)
    wsel = np.nonzero((L.info[:, 0] == 0) & (L.info[:, 1] == 1))[0]
    Dw = D[wsel].toarray()
    U_, s, _ = np.linalg.svd(Dw, full_matrices=True)
    rank = int(np.sum(s > 1e-10 * s[0]))
    Cw = U_[:, rank:]
    Cw[np.abs(Cw) < 1e-14] = 0.0
    Rw = sp.csr_matrix(Cw)
    Rw = sp.csr_matrix((Rw.data, (wsel[Rw.nonzero()[0]], Rw.nonzero()[1])), shape=(L.n, Cw.shape[1]))
    osel = np.nonzero((L.info[:, 0] == 0) & (L.info[:, 1] > 1))[0]
    return sp.hstack([Rw, hier.unit_cols(L.n, osel)]).tocsr()


def global_gammas(L, C):
    n = L.n
    out = {}
    f0 = np.nonzero(L.info[:, 0] == 0)[0]
    f1 = np.nonzero(L.info[:, 0] == 1)[0]
    out['b1'] = gamma_sparse(L.M, hier.unit_cols(n, f0), hier.unit_cols(n, f1))
    out['b2'] = gamma_sparse(L.M, whitney_complement(L), hier.grad_basis(L))
    if C is not None:
        P = hier.injection(L, C)
        coarse = np.asarray(P.sum(axis=1)).ravel() > 0
        new = np.nonzero(~coarse)[0]
        out['a'] = gamma_sparse(L.M, P, hier.unit_cols(n, new))
        cr = np.nonzero(coarse & (L.info[:, 0] == 0))[0]
        nr = np.nonzero(~coarse & (L.info[:, 0] == 0))[0]
        if len(nr):
            out['a_rot'] = gamma_sparse(L.M, hier.unit_cols(n, cr), hier.unit_cols(n, nr))
        # gradients: G_{p-1} (injected) vs new Bg
        Gc = P @ hier.grad_basis(C)
        ng = np.nonzero(~coarse & (L.info[:, 0] == 1))[0]
        out['a_grad'] = gamma_sparse(L.M, Gc, hier.unit_cols(n, ng))
    return out


def local_gammas(d, p):
    """Max over cells of element-level constants (prototype element matrices, OFT basis)."""
    import basis
    from analyze import load
    _, info_full, _, _ = load(d)
    m = hier.read_mesh(d)
    with open(os.path.join(d, 'dofs.txt')) as f:
        f.readline()
        gs = [list(map(int, f.readline().split()))[1:] for _ in range(2)]
    gsc = [[0, 1, max(0, gs[0][2] - (p + 1) + 1 - 0), 0], None]  # placeholder (unused)
    # coarse counts per (field, dim): gstruct of order p-1
    cnt = {1: [[0, 1, 0, 0], [1, 1, 0, 0]], 2: [[0, 1, 2, 0], [1, 2, 1, 0]],
           3: [[0, 1, 5, 3], [1, 3, 3, 1]], 4: [[0, 1, 9, 11], [1, 4, 6, 4]]}
    basis.ELEM_HOOK = []
    basis.assemble(d, p, 'oft')
    res = {'a': [], 'b1': [], 'b2': [], 'a_rot': [], 'a_grad': []}
    for c, (gi, Me) in enumerate(basis.ELEM_HOOK):
        inf = info_full[gi]
        keep = ~((inf[:, 0] == 1) & (inf[:, 1] == 0))
        Me = Me[keep][:, keep]
        inf = inf[keep]
        nl = len(inf)
        f0 = inf[:, 0] == 0
        res['b1'].append(gamma_dense(Me, np.eye(nl)[:, f0], np.eye(nl)[:, ~f0]))
        # local Whitney gradients: grad(l_a) = sum_k s_ka w_k  (s = +1 if a is head, -1 if tail)
        wl = np.nonzero(f0 & (inf[:, 1] == 1))[0]   # 6 local Whitney, in local edge order k
        Dl = np.zeros((nl, 4))
        for k in range(6):
            a, b = m['oe'][c, k]
            Dl[wl[k], a] -= 1.0
            Dl[wl[k], b] += 1.0
        Gl = np.hstack([Dl[:, 1:], np.eye(nl)[:, ~f0]])
        Uw, s, _ = np.linalg.svd(Dl[wl], full_matrices=True)
        Cw = np.zeros((nl, 3))
        Cw[wl] = Uw[:, 3:]
        Rl = np.hstack([Cw, np.eye(nl)[:, f0 & (inf[:, 1] > 1)]])
        res['b2'].append(gamma_dense(Me, Rl, Gl))
        if p > 1:
            cc = cnt[p - 1]
            coarse = np.array([k < cc[fi][dm] for fi, dm, _, k in inf])
            res['a'].append(gamma_dense(Me, np.eye(nl)[:, coarse], np.eye(nl)[:, ~coarse]))
            if np.any(~coarse & f0):
                res['a_rot'].append(gamma_dense(Me, np.eye(nl)[:, coarse & f0], np.eye(nl)[:, ~coarse & f0]))
            # local coarse gradients: vertex grads (via Whitney) + coarse Bg
            Gcl = np.hstack([Dl[:, 1:], np.eye(nl)[:, coarse & ~f0]])
            res['a_grad'].append(gamma_dense(Me, Gcl, np.eye(nl)[:, ~coarse & ~f0]))
    basis.ELEM_HOOK = None
    return {k: max(v) for k, v in res.items() if v}


def gamma_dense(M, U, W):
    return hier.gamma(M, U, W)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'global'
    if which == 'global':
        for mesh, ps in (('tiny', (1, 2, 3, 4)), ('ref', (1, 2, 3))):
            prev = None
            for p in ps:
                L = hier.Level(os.path.join(hier.HERE, 'dump3', '%s_p%d' % (mesh, p)))
                g = global_gammas(L, prev)
                print('global %s p=%d ' % (mesh, p) + ' '.join('%s=%.6f' % kv for kv in g.items()), flush=True)
                prev = L
    else:
        for p in (1, 2, 3):
            g = local_gammas(os.path.join(hier.HERE, 'dump3', 'tiny_p%d' % p), p)
            print('local tiny p=%d ' % p + ' '.join('%s=%.6f' % kv for kv in g.items()), flush=True)
