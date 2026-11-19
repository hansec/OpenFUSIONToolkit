#------------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only: which canonical functions the CURRENT OFT tet
# transforms actually need (nonzero pattern), the tet sort+parity orientation
# classes, and S4-invariance of the interior blocks (for a "virtual sort").
#------------------------------------------------------------------------------
import itertools
import numpy as np
import sympy as sp
from oft_tet_basis import (h1_basis, hcurl_basis, random_tet_points, sample_matrix, TET_ED, TET_FC,
                           F, hc_cell, h1_cell, TetBasis)

pts = random_tet_points(40)
FACE = 13


def groups(basis, kind):
    """Column groups: (label, sample matrix) per entity of a basis."""
    out = []
    for ent in range(15):
        fn = basis.entity_funcs(ent, 1)
        if not fn or (kind == 'grad' and ent < 4):
            continue
        lab = ('v%d' % (ent+1) if ent < 4 else 'e%d' % (ent-3) if ent < 10 else 'f%d' % (ent-9) if ent < 14 else 'int')
        out.append((lab, sample_matrix(fn, kind, pts)))
    return out


def pattern(target, grps, tag):
    A = np.hstack([g[1] for g in grps])
    X, _, _, _ = np.linalg.lstsq(A, target, rcond=None)
    res = np.max(np.abs(A @ X - target))/np.max(np.abs(A))
    used = []
    i0 = 0
    for lab, g in grps:
        n = g.shape[1]
        if np.max(np.abs(X[i0:i0+n])) > 1.e-10:
            used.append('%s%s(%d)' % (tag, lab, n))
        i0 += n
    return used, res


print('## Nonzero regressor pattern of current transforms (face 4 = (1,2,3), all states)')
for k in (4, 5):
    b = h1_basis(k)
    used = set()
    for s in range(2, 7):
        u, r = pattern(sample_matrix(b.entity_funcs(FACE, s), 'h1', pts), groups(b, 'h1'), '')
        used |= set(u)
    print('H1 k=%d face: %s' % (k, sorted(used)))
for p in (2, 3, 4):
    hc = hcurl_basis(p)
    used = set()
    for s in range(2, 7):
        g = [(l, m) for l, m in groups(hc, 'hc')] + [('g' + l, m) for l, m in groups(h1_basis(p+1), 'grad')]
        u, r = pattern(sample_matrix(hc.entity_funcs(FACE, s), 'hc', pts), g, '')
        used |= set(u)
    print('HCurl p=%d face (regressors HC + grad H1_%d): %s, resid %.1e' % (p, p+1, sorted(used), r))

print()
print('## Sort-by-global-index + parity fix: distinct (edge state, face state) patterns')
pats = {}
for gl in itertools.permutations(range(4)):
    order = sorted(range(4), key=lambda i: gl[i])    # local vertices sorted by global index
    lc = [gl[i] for i in order]
    perm_sign = np.linalg.det(np.eye(4)[order])
    if perm_sign < 0:
        lc[2], lc[3] = lc[3], lc[2]   # tetmesh_invert_cell: swap local 3,4
    est = tuple(1 if lc[a-1] < lc[b-1] else 2 for a, b in TET_ED)
    fst = []
    for fc in TET_FC:
        g = [lc[v-1] for v in fc]
        # state = which orient_face maps (1,2,3) positions to sorted order
        srt = tuple(sorted(range(3), key=lambda i: g[i]))
        from oft_tet_basis import orient_face
        st = [s for s in range(1, 7) if tuple(v-1 for v in orient_face((1, 2, 3), s)) == srt][0]
        fst.append(st)
    pats.setdefault((est, tuple(fst)), 0)
    pats[(est, tuple(fst))] += 1
for k, v in pats.items():
    print('  edge states %s face states %s : %d of 24 global orderings' % (k[0], k[1], v))

print()
print('## S4 invariance of the interior blocks (needed by a per-cell "virtual sort")')
for kind, maker, cf, ks in (('h1', h1_basis, h1_cell, (4, 5)), ('hc', hcurl_basis, hc_cell, (3, 4))):
    for k in ks:
        b = maker(k)
        A = np.hstack([g[1] for g in groups(b, kind)])
        n_int = b.ent_ndof[14]
        own0 = b.ent_dof0[14]
        worst = 0.
        worst_own = 0.
        for perm in itertools.permutations(range(4)):
            fp = [F[i] for i in perm]
            if kind == 'h1':
                fn = [cf(j, *fp) for j in range(1, n_int+1)]
            else:
                fn = [dict(zip([i+1 for i in perm], cf(j, *fp))) for j in range(1, n_int+1)]
            B = sample_matrix(fn, kind, pts)
            Xo, _, _, _ = np.linalg.lstsq(A[:, own0:own0+n_int], B, rcond=None)
            worst_own = max(worst_own, np.max(np.abs(A[:, own0:own0+n_int] @ Xo - B))/np.max(np.abs(B)))
            X, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
            worst = max(worst, np.max(np.abs(A @ X - B))/np.max(np.abs(B)))
        print('  %s order %d interior (%d DOFs): max residual over 24 relabelings: own block %.1e, whole cell %.1e' % (
            kind, k, n_int, worst_own, worst))

print()
print('## H(Curl) interior relabelings fitted with H(Curl) + grad H^1_{p+1} (composite space)')
for k in (3, 4):
    b = hcurl_basis(k)
    A = np.hstack([g[1] for g in groups(b, 'hc')] + [g[1] for g in groups(h1_basis(k+1), 'grad')])
    worst = 0.
    for perm in itertools.permutations(range(4)):
        fp = [F[i] for i in perm]
        fn = [dict(zip([i+1 for i in perm], hc_cell(j, *fp))) for j in range(1, b.ent_ndof[14]+1)]
        B = sample_matrix(fn, 'hc', pts)
        X, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
        worst = max(worst, np.max(np.abs(A @ X - B))/np.max(np.abs(B)))
    print('  HCurl p=%d interior: max residual over 24 relabelings %.1e' % (k, worst))
