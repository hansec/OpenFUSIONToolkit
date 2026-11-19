#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only (basis study for the REFERENCE-mode quadrature cache).
#
# Python port of the Phase 2 orientation-transform classifier
# (qcache_build_xform in src/fem/fem_quad_cache.F90) for tet bases described in
# oft_tet_basis.py, plus representation-theory helpers for the face symmetry
# group S3.
#------------------------------------------------------------------------------
import numpy as np
from oft_tet_basis import TET_ED, TET_FC, sample_matrix

KIND_NAMES = ('identity', 'signed perm', 'dense', 'none')


def _sub_entities(ent):
    """Sub-entities used as extra regressors (Fortran order: vertices, edges, interior)."""
    if 4 <= ent < 10:
        return [v-1 for v in TET_ED[ent-4]]
    if 10 <= ent < 14:
        fc = TET_FC[ent-10]
        subs = [v-1 for v in fc]
        for e, ed in enumerate(TET_ED):
            if ed[0] in fc and ed[1] in fc:
                subs.append(4+e)
        subs.append(14)
        return subs
    return []


def _fit(A, B):
    X, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    X[np.abs(X) < 1.e-13] = 0.
    m = np.abs(np.abs(X) - 1.) < 1.e-13
    X[m] = np.sign(X[m])
    res = np.max(np.abs(B - A @ X))/max(np.max(np.abs(A)), 1.e-300)
    return X, res


class Classifier:
    """Fit oriented entity blocks against canonical (state 1) functions, as in qcache_build_xform.

    extra: optional list of (TetBasis, kind) pairs whose canonical functions are appended
    as additional regressors in the last-resort fit (eg. grad H^1 for H(Curl)).
    """

    def __init__(self, basis, pts, tol=1.e-10, extra=None):
        self.basis = basis
        self.tol = tol
        kind = 'h1' if basis.kind == 'h1' else 'hc'
        self.kind = kind
        self.pts = pts
        self.canon = {}
        cols = []
        for ent in range(15):
            fn = basis.entity_funcs(ent, 1)
            cols.append(sample_matrix(fn, kind, pts) if fn else None)
        self.A = np.hstack([c for c in cols if c is not None])
        self.extraA = None
        if extra is not None:
            ex = []
            for eb, ekind in extra:
                for ent in range(15):
                    fn = eb.entity_funcs(ent, 1)
                    if fn and not (ekind == 'grad' and ent < 4):
                        ex.append(sample_matrix(fn, ekind, pts))
            self.extraA = np.hstack(ex)

    def block(self, ent, state):
        kind = self.kind
        return sample_matrix(self.basis.entity_funcs(ent, state), kind, self.pts)

    def classify(self, verbose=False):
        b = self.basis
        counts = np.zeros(4, dtype=int)
        details = []
        maxres = 0.
        for ent in range(4, 14):
            nd = b.ent_ndof[ent]
            if nd == 0:
                continue
            d0 = b.ent_dof0[ent]
            own = list(range(d0, d0+nd))
            subs = []
            for se in _sub_entities(ent):
                subs += list(range(b.ent_dof0[se], b.ent_dof0[se]+b.ent_ndof[se]))
            for s in range(2, b.ent_nst[ent]+1):
                B = self.block(ent, s)
                regs = own
                X, res = _fit(self.A[:, regs], B)
                if res > self.tol and subs:
                    regs = own + subs
                    X, res = _fit(self.A[:, regs], B)
                if res > self.tol:
                    regs = list(range(b.nce))
                    X, res = _fit(self.A[:, regs], B)
                kind = None
                xres = None
                if res > self.tol:
                    kind = 3
                    if self.extraA is not None:
                        Aext = np.hstack([self.A, self.extraA])
                        _, xres = _fit(Aext, B)
                else:
                    maxres = max(maxres, res)
                    is_perm = True
                    perm = []
                    for j in range(nd):
                        nz = np.nonzero(X[:, j])[0]
                        if len(nz) != 1 or nz[0] >= nd or abs(abs(X[nz[0], j])-1.) > 1.e-10:
                            is_perm = False
                            break
                        perm.append(nz[0])
                    if is_perm and sorted(perm) != list(range(nd)):
                        is_perm = False
                    if is_perm:
                        kind = 0 if (perm == list(range(nd)) and np.all(np.diag(X[:nd, :nd]) > 0)) else 1
                    else:
                        kind = 2
                counts[kind] += 1
                details.append(dict(ent=ent, state=s, kind=kind, nreg=len(regs), res=res, xres=xres,
                                    X=X if kind != 3 else None))
        self.counts = counts
        self.details = details
        self.maxres = maxres
        return counts


#------------------------------------------------------------------------------
# Representation theory of S3 (face symmetry group)
#------------------------------------------------------------------------------
# Face states 1..6 realize the following permutations of the face's positions;
# classes: identity (state 1), 3-cycles (states 2,3), transpositions (states 4,5,6)
S3_CLASS = {1: 'e', 2: 'c', 3: 'c', 4: 't', 5: 't', 6: 't'}


def s3_decompose(chars):
    """Multiplicities (m_T, m_A, m_S) of trivial, sign and standard irreps from a character
    given as dict state -> trace."""
    ce = chars[1]
    cc = np.mean([chars[2], chars[3]])
    ct = np.mean([chars[4], chars[5], chars[6]])
    mT = (ce + 2*cc + 3*ct)/6.
    mA = (ce + 2*cc - 3*ct)/6.
    mS = (2*ce - 2*cc)/6.
    return np.round([mT, mA, mS], 8)


def monomial_ok(m):
    """A real S3 representation has a basis transforming by signed permutations iff every
    standard irrep can be paired with a 1-d irrep (Ind from Z2 of a +-1 character gives
    T+S or A+S; Ind from Z3 of the trivial character gives T+A)."""
    mT, mA, mS = m
    return mS <= mT + mA + 1.e-9
