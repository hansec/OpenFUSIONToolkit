#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only: S3-equivariant orthogonalization of the prototype
# symmetric face families (sym_basis.py).
#
# Each hierarchical level's generators are replaced by their L2-orthogonal
# projection (exact rational arithmetic) off
#   - the previous levels of the same face block, and
#   - for H(Curl): the gradients of the H^1 face functions up to the same degree,
# using the L2 inner product of a REGULAR tetrahedron. The symmetry group S3 of a
# face acts on the regular tet by isometries fixing the opposite vertex, so these
# projections commute with S3: each level stays exactly invariant and keeps its
# signed-permutation / 2x2 structure. The spans required by OFT are unchanged:
# H^1 levels change by lower H^1 face functions; H(Curl) levels change by lower
# rotational face functions plus gradients of H^1 face functions.
#------------------------------------------------------------------------------
import itertools
from functools import lru_cache
import sympy as sp
from oft_tet_basis import TetBasis, H1_GSTRUCT, HC_GSTRUCT, h1_edge, h1_cell, whitney, hc_cell
from sym_basis import hc_face_sym, h1_face_sym

X = sp.symbols('x1:4')  # face barycentrics (formula variables)
# Regular tet: grad lam_i . grad lam_j proportional to 4*delta_ij - 1 (i,j over all 4 vertices);
# formulas only use grad lam_1..3 of the face
GMET = sp.Matrix(3, 3, lambda i, j: 4 if i == j else -1)


def _int_mono(a, b, c):
    """Integral over the tet of x1^a x2^b x3^c (volume normalized to 1/6 scale)."""
    return sp.factorial(a)*sp.factorial(b)*sp.factorial(c)*6/sp.factorial(a+b+c+3)


def _integrate(expr):
    p = sp.Poly(sp.expand(expr), *X)
    out = sp.Integer(0)
    for (a, b, c), coef in p.terms():
        out += coef*_int_mono(a, b, c)
    return out


def ip_scalar(u, v):
    return _integrate(u*v)


def ip_form(y, z):
    out = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            if y[i] != 0 and z[j] != 0:
                out += GMET[i, j]*_integrate(y[i]*z[j])
    return out


def grad_form(h):
    return tuple(sp.expand(sp.diff(h, X[i])) for i in range(3))


def project_off(gens, basis_list, ip):
    """Return gens minus their ip-orthogonal projection onto span(basis_list)."""
    if not basis_list:
        return gens
    n = len(basis_list)
    G = sp.Matrix(n, n, lambda i, j: ip(basis_list[i], basis_list[j]))
    out = []
    for g in gens:
        rhs = sp.Matrix(n, 1, lambda i, j: ip(basis_list[i], g))
        c = G.LUsolve(rhs)
        if isinstance(g, tuple):
            new = list(g)
            for k in range(n):
                for m in range(3):
                    new[m] = new[m] - c[k]*basis_list[k][m]
            out.append(tuple(sp.expand(v) for v in new))
        else:
            out.append(sp.expand(g - sum(c[k]*basis_list[k] for k in range(n))))
    return out


@lru_cache(None)
def h1_face_orth_formulas():
    lv = [[1], [2, 3], [4, 5, 6]]
    done = []
    out = []
    for level in lv:
        gens = [h1_face_sym(d, *X) for d in level]
        gens = project_off(gens, done, ip_scalar)
        done += gens
        out += gens
    return out


@lru_cache(None)
def hc_face_orth_formulas(with_grads=True):
    lv = [[1, 2], [3, 4, 5], [6, 7, 8, 9]]
    h1f = h1_face_orth_formulas()
    h1lv = [[0], [1, 2], [3, 4, 5]]   # H1 face levels k=3,4,5 match H(Curl) degrees 2,3,4
    done = []
    out = []
    grads = []
    for il, level in enumerate(lv):
        gens = [hc_face_sym(d, *X, variant='orth') for d in level]
        if with_grads:
            grads += [grad_form(h1f[i]) for i in h1lv[il]]
        gens = project_off(gens, done + grads, ip_form)
        done += gens
        out += gens
    return out


def _subs(expr, a, b, c):
    return sp.sympify(expr).subs({X[0]: a, X[1]: b, X[2]: c}, simultaneous=True)


def h1_basis_orth(order):
    forms = h1_face_orth_formulas()
    return TetBasis('h1', H1_GSTRUCT[order], efun=h1_edge,
                    ffun=lambda d, a, b, c: _subs(forms[d-1], a, b, c), cfun=h1_cell,
                    name='H1orth_%d' % order)


def hcurl_basis_orth(order, with_grads=True):
    forms = hc_face_orth_formulas(with_grads)
    return TetBasis('hc', HC_GSTRUCT[order], efun=whitney,
                    ffun=lambda d, a, b, c: tuple(_subs(v, a, b, c) for v in forms[d-1]), cfun=hc_cell,
                    name=('HCorthg_%d' if with_grads else 'HCorth_%d') % order)


if __name__ == '__main__':
    for i, y in enumerate(hc_face_orth_formulas()):
        print('HC face', i+1, [sp.factor(v) for v in y])
    for i, h in enumerate(h1_face_orth_formulas()):
        print('H1 face', i+1, sp.factor(h))
