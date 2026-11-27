#------------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only: S3 characters of the tet face blocks of OFT's
# H^1 and H(Curl) bases (and of each hierarchical increment), computed from the
# action of the 6 face states on the face block modulo all other canonical
# functions (H(Curl): modulo H(Curl) + grad H^1).
#------------------------------------------------------------------------------
import sys
import numpy as np
from oft_tet_basis import h1_basis, hcurl_basis, random_tet_points, sample_matrix
from xform_classify import s3_decompose, monomial_ok

FACE = 13  # entity index of face (1,2,3); all faces are equivalent


def face_action(basis, pts, extra=None):
    kind = 'h1' if basis.kind == 'h1' else 'hc'
    cols = []
    for ent in range(15):
        fn = basis.entity_funcs(ent, 1)
        if fn:
            cols.append(sample_matrix(fn, kind, pts))
    A = np.hstack(cols)
    if extra is not None:
        ex = []
        for ent in range(4, 15):
            fn = extra.entity_funcs(ent, 1)
            if fn:
                ex.append(sample_matrix(fn, 'grad', pts))
        A = np.hstack([A] + ex)
    d0 = basis.ent_dof0[FACE]
    nd = basis.ent_ndof[FACE]
    mats = {}
    for s in range(1, 7):
        B = sample_matrix(basis.entity_funcs(FACE, s), kind, pts)
        X, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
        res = np.max(np.abs(A @ X - B))/np.max(np.abs(A))
        assert res < 1.e-10, res
        mats[s] = X[d0:d0+nd, :]
    return mats


def characters(mats, idx):
    return {s: np.trace(m[np.ix_(idx, idx)]) for s, m in mats.items()}


def report(name, mats, levels):
    """levels: list of (label, dof indices of the block up to that order, new indices)."""
    for label, upto, new in levels:
        mfull = s3_decompose(characters(mats, upto))
        line = '%-10s %-6s full block n=%2d  (T,A,S)=%s %s' % (name, label, len(upto), mfull.astype(int),
                                                             'monomial' if monomial_ok(mfull) else 'NOT monomial')
        if new is not None and len(new) != len(upto):
            # increment: quotient of block(<=order) by block(<order); block(<order) is invariant
            prev = [i for i in upto if i not in new]
            mprev = s3_decompose(characters(mats, prev))
            minc = mfull - mprev
            line += ' | increment n=%d (T,A,S)=%s %s' % (len(new), minc.astype(int),
                                                         'monomial' if monomial_ok(minc) else 'NOT monomial')
        print(line)


if __name__ == '__main__':
    pts = random_tet_points(40)
    # H^1 order 5 contains all lower orders hierarchically (face dofs 1 | 2-3 | 4-6)
    m = face_action(h1_basis(5), pts)
    report('H1 face', m, [('k=3', [0], [0]), ('k=4', [0, 1, 2], [1, 2]), ('k=5', list(range(6)), [3, 4, 5])])
    # H(Curl) order 4 with grad H^1 order 5 (face dofs 1-2 | 3-5 | 6-9)
    m = face_action(hcurl_basis(4), pts, extra=h1_basis(5))
    report('HCurl face', m, [('p=2', [0, 1], [0, 1]), ('p=3', list(range(5)), [2, 3, 4]),
                             ('p=4', list(range(9)), [5, 6, 7, 8])])
    # Hierarchy check: are lower-order blocks invariant on their own?
    for s in range(2, 7):
        mm = m[s]
        leak = max(np.max(np.abs(mm[2:, :2])), np.max(np.abs(mm[5:, 2:5])))
        if leak > 1.e-12:
            print('HCurl face state %d: lower-order block mixes into higher-order face DOFs (max %.2e)' % (s, leak))
