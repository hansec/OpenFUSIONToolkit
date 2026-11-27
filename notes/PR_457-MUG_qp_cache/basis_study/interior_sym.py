#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only: feasibility of an S4-invariant (vertex-relabeling
# invariant) gradient-free H(Curl) cell-interior block for tets.
#
# OFT's interior functions (hc_cell) are projected, level by level (degree 3:
# DOFs 1-3, degree 4: DOFs 4-11), L2-orthogonally on a REGULAR tet (S4 acts by
# isometries) off
#   - the gradients of the H^1 interior bubbles up to order p+1, and
#   - the previous (already projected) interior level.
# The result differs from OFT's functions by interior gradients and lower interior
# functions only, so H(Curl)+grad H^1 is unchanged, the direct sum (gradient-free
# property) is kept, and the levels stay hierarchical. Since the projection
# commutes with S4, each level's span becomes S4-invariant. Checks:
#   1. ranks (span / direct sum) with OFT faces and with the sym-proj faces;
#   2. S4 invariance of the interior block (24 relabelings), old vs new;
#   3. Jacobi-scaled conditioning of local matrices (unit and regular tet).
# Usage: python interior_sym.py   (writes results_interior.txt)
#------------------------------------------------------------------------------
import os
import itertools
from functools import lru_cache
import numpy as np
import sympy as sp
import oft_tet_basis as otb
from oft_tet_basis import (F, TetBasis, HC_GSTRUCT, whitney, hc_face, hc_cell, h1_cell, h1_basis,
                           hcurl_basis, random_tet_points, sample_matrix)
from sym_orth import hc_face_orth_formulas, _subs
from verify_sym import canon, rank, cond_spd, local_mats, combined_mass, tet_quadrature, set_geometry

L = sp.symbols('l1:5')
GMET4 = sp.Matrix(4, 4, lambda i, j: 3 if i == j else -1)   # regular tet: grad l_i . grad l_j ~ 4 delta_ij - 1
OUT = []


def log(msg=''):
    print(msg)
    OUT.append(msg)


def _int4(expr):
    p = sp.Poly(sp.expand(expr), *L)
    out = sp.Integer(0)
    for alpha, coef in p.terms():
        num = sp.Integer(1)
        for a in alpha:
            num *= sp.factorial(a)
        out += coef*num*6/sp.factorial(sum(alpha)+3)
    return out


def ip4(y, z):
    out = sp.Integer(0)
    for i in range(4):
        for j in range(4):
            if y[i] != 0 and z[j] != 0:
                out += GMET4[i, j]*_int4(y[i]*z[j])
    return out


def grad4(h):
    return tuple(sp.expand(sp.diff(h, L[i])) for i in range(4))


def project_off(gens, basis_list):
    n = len(basis_list)
    G = sp.Matrix(n, n, lambda i, j: ip4(basis_list[i], basis_list[j]))
    out = []
    for g in gens:
        c = G.LUsolve(sp.Matrix(n, 1, lambda i, j: ip4(basis_list[i], g)))
        out.append(tuple(sp.expand(g[m] - sum(c[k]*basis_list[k][m] for k in range(n))) for m in range(4)))
    return out


@lru_cache(None)
def hc_interior_sym_formulas():
    levels = [list(range(1, 4)), list(range(4, 12))]
    h1lv = [[1], [1, 2, 3, 4]]          # H^1 interior DOFs of order 4 / 5
    done = []
    for il, lev in enumerate(levels):
        gens = [tuple(sp.sympify(v) for v in hc_cell(d, *L)) for d in lev]
        grads = [grad4(h1_cell(d, *L)) for d in h1lv[il]]
        done += project_off(gens, done + grads)
    return done


def hcurl_basis_isym(order, faces='oft'):
    cf = hc_interior_sym_formulas()
    sub = lambda d, a, b, c, e: tuple(sp.sympify(v).subs(dict(zip(L, (a, b, c, e))), simultaneous=True)
                                      for v in cf[d-1])
    if faces == 'oft':
        ffun = hc_face
    else:
        fforms = hc_face_orth_formulas()
        ffun = lambda d, a, b, c: tuple(_subs(v, a, b, c) for v in fforms[d-1])
    return TetBasis('hc', HC_GSTRUCT[order], efun=whitney, ffun=ffun, cfun=sub,
                    name='HCisym(%s faces)_%d' % (faces, order))


def hcurl_basis_projfaces(order):
    fforms = hc_face_orth_formulas()
    return TetBasis('hc', HC_GSTRUCT[order], efun=whitney,
                    ffun=lambda d, a, b, c: tuple(_subs(v, a, b, c) for v in fforms[d-1]), cfun=hc_cell,
                    name='HCorthg_%d' % order)


def interior_s4_residual(basis, pts, extra=None):
    """Max relative residual of the 24 relabeled interior blocks fitted to the interior block
    itself (own) and to all cell functions (+ grad H^1 if extra)."""
    A = canon(basis, pts)
    if extra is not None:
        A_all = np.hstack([A, canon(extra, pts, kind='grad', skip_vertices=True)])
    else:
        A_all = A
    d0 = basis.ent_dof0[14]
    nd = basis.ent_ndof[14]
    Aown = A[:, d0:d0+nd]
    worst_own = worst_all = 0.
    for perm in itertools.permutations(range(4)):
        fp = [F[i] for i in perm]
        fns = []
        for j in range(1, nd+1):
            y = basis.cfun(j, *fp)
            fns.append(dict(zip([i+1 for i in perm], y)))
        B = sample_matrix(fns, 'hc', pts)
        for M, key in ((Aown, 'own'), (A_all, 'all')):
            X, _, _, _ = np.linalg.lstsq(M, B, rcond=None)
            r = np.max(np.abs(M @ X - B))/np.max(np.abs(B))
            if key == 'own':
                worst_own = max(worst_own, r)
            else:
                worst_all = max(worst_all, r)
    return worst_own, worst_all


if __name__ == '__main__':
    pts = random_tet_points(40)
    log('## Projected interior formulas (degree-3 level, exact)')
    for i, y in enumerate(hc_interior_sym_formulas()[:3]):
        log('  interior %d: %s' % (i+1, [sp.factor(v) for v in y]))
    log()
    log('## 1. Spaces')
    for p in (3, 4):
        g = canon(h1_basis(p+1), pts, kind='grad', skip_vertices=True)
        old = canon(hcurl_basis(p), pts)
        for nb in (hcurl_basis_isym(p, 'oft'), hcurl_basis_isym(p, 'proj')):
            new = canon(nb, pts)
            log('%-24s p=%d: rank HC old/new %d/%d, rank(old|new) %d; rank(HC+gradH1) old %d new %d union %d; '
                'NED2 %d; direct-sum count %d' % (
                    nb.name, p, rank(old), rank(new), rank(np.hstack([old, new])), rank(np.hstack([old, g])),
                    rank(np.hstack([new, g])), rank(np.hstack([old, new, g])), (p+1)*(p+2)*(p+3)//2,
                    new.shape[1]+g.shape[1]))
        # hierarchy: p=3 interior functions are the first 3 of p=4 (same formulas by construction)
        b3 = hcurl_basis_isym(3, 'proj').entity_funcs(14)
        b4 = hcurl_basis_isym(4, 'proj').entity_funcs(14)[:3]
        same = all(sp.simplify(sp.sympify(b3[j][v]) - sp.sympify(b4[j][v])) == 0 for j in range(3) for v in b3[j])
        if p == 4:
            log('hierarchy: p=3 interior functions identical to the first 3 of p=4: %s' % same)
    log()
    log('## 2. S4 invariance of the interior block (24 vertex relabelings; max relative fit residual)')
    for p in (3, 4):
        for b in (hcurl_basis(p), hcurl_basis_isym(p, 'oft')):
            own, allr = interior_s4_residual(b, pts, extra=h1_basis(p+1))
            log('%-24s p=%d: own interior block %.1e, HC + grad H1 %.1e' % (b.name, p, own, allr))
    log()
    log('## 3. Conditioning (Jacobi-scaled eigenvalue ratio; zero modes removed)')
    qp, qw = tet_quadrature(8)
    geoms = {'unit': np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
             'regular': np.array([[1., 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])}
    for gname, verts in geoms.items():
        set_geometry(verts)
        log('Geometry: %s tet' % gname)
        for p in (3, 4):
            row = []
            for lab, b in (('old', hcurl_basis(p)), ('old faces + sym int', hcurl_basis_isym(p, 'oft')),
                           ('sym-proj faces', hcurl_basis_projfaces(p)),
                           ('sym-proj faces + sym int', hcurl_basis_isym(p, 'proj'))):
                M, C = local_mats(b, qp, qw, 'all')
                Mc = combined_mass(b, h1_basis(p+1), qp, qw)
                row.append('%s: M %.3e cc %.3e comb %.3e' % (lab, cond_spd(M), cond_spd(C, 3), cond_spd(Mc)))
            log('  p=%d  ' % p + ' | '.join(row))
    set_geometry(geoms['unit'])
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_interior.txt'), 'w') as fid:
        fid.write('\n'.join(OUT)+'\n')
