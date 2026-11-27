#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only: prototype "symmetric" (S3-invariant, hierarchical)
# tet face families for OFT's H^1 and gradient-free H(Curl) spaces. Edge, vertex
# and interior functions are OFT's, unchanged. Not used by OFT.
#
# Notation (face with oriented vertices a,b,c -> f1,f2,f3): b = f1*f2*f3,
# w_ij = f_i grad f_j - f_j grad f_i (Whitney form). H(Curl) face functions are
# returned as coefficients (y1,y2,y3) of (grad f1, grad f2, grad f3).
#------------------------------------------------------------------------------
import sympy as sp
from oft_tet_basis import TetBasis, H1_GSTRUCT, HC_GSTRUCT, h1_edge, h1_cell, whitney, hc_cell


def _w(i, j, f):
    """Whitney form w_ij as a coefficient triple."""
    y = [0, 0, 0]
    y[j] += f[i]
    y[i] -= f[j]
    return y


def _add(*terms):
    out = [0, 0, 0]
    for c, y in terms:
        for k in range(3):
            out[k] += c*y[k]
    return tuple(sp.expand(v) for v in out)


#------------------------------------------------------------------------------
# H(Curl) gradient-free face block, hierarchical and S3-invariant per level
#   p=2 (2 DOFs, standard irrep S, one 2x2 block):
#     1: lam_3 w_12 - lam_1 w_23 ... written as 2C - A - B, 2: A - B (OFT DOF 2)
#   p=3 (+3, A+S, signed permutation triple):  lam_k^2 w_ij  (k opposite edge ij)
#   p=4 (+4, T+A+S): triple lam_k^2 (lam_i - lam_j) w_ij (permuted, no signs)
#                    + singleton b (w_12 + w_23 + w_31) (sign character)
#------------------------------------------------------------------------------
def hc_face_sym(dof, f1, f2, f3, variant='int'):
    f = (f1, f2, f3)
    A = (f2*f3, 0, 0)
    B = (0, f1*f3, 0)
    C = (0, 0, f1*f2)
    b = f1*f2*f3
    if dof == 1:
        if variant == 'int':
            # B - C = lam_1 w_23 (with A - B = -lam_3 w_12): integer 2x2 transforms
            return _add((1, B), (-1, C))
        # 'orth': (2C - A - B)/sqrt(3) together with A - B gives orthogonal 2x2 transforms
        return _add((2/sp.sqrt(3), C), (-1/sp.sqrt(3), A), (-1/sp.sqrt(3), B))
    if dof == 2:
        return _add((1, A), (-1, B))
    # cyclic triples: (k; i, j) with k the vertex opposite edge (i,j)
    cyc = [(2, 0, 1), (0, 1, 2), (1, 2, 0)]
    if 3 <= dof <= 5:
        k, i, j = cyc[dof-3]
        return _add((f[k]**2, _w(i, j, f)))
    if 6 <= dof <= 8:
        k, i, j = cyc[dof-6]
        return _add((f[k]**2*(f[i] - f[j]), _w(i, j, f)))
    if dof == 9:
        return _add((b, _w(0, 1, f)), (b, _w(1, 2, f)), (b, _w(2, 0, f)))
    raise ValueError(dof)


#------------------------------------------------------------------------------
# H^1 face block, hierarchical and S3-invariant per level
#   k=3: -2 b (OFT DOF 1, trivial irrep)
#   k=4: 2 b (x2 - x1) (OFT DOF 3) and 2 b (x1 + x2 - 2 x3) (standard irrep, one 2x2 block)
#   k=5: 2 b x_j x_k for the 3 vertex pairs (permutation triple)
#------------------------------------------------------------------------------
def h1_face_sym(dof, x1, x2, x3, variant='int'):
    b = x1*x2*x3
    if dof == 1:
        return -2*b
    if dof == 2:
        if variant == 'int':
            return 2*b*(x3 - x2)
        return 2*b*(x1 + x2 - 2*x3)/sp.sqrt(3)
    if dof == 3:
        return 2*b*(-x1 + x2)
    if dof == 4:
        return 2*b*x2*x3
    if dof == 5:
        return 2*b*x1*x3
    if dof == 6:
        return 2*b*x1*x2
    raise ValueError(dof)


def h1_face_perm(dof, x1, x2, x3):
    """Non-hierarchical alternative (order 4 block = b*x_i, pure permutations at all orders)."""
    b = x1*x2*x3
    return [-2*b*x1, -2*b*x2, -2*b*x3, 2*b*x2*x3, 2*b*x1*x3, 2*b*x1*x2][dof-1]


def h1_basis_sym(order, variant='int'):
    return TetBasis('h1', H1_GSTRUCT[order], efun=h1_edge,
                    ffun=lambda d, a, b, c: h1_face_sym(d, a, b, c, variant), cfun=h1_cell,
                    name='H1sym_%d' % order)


def h1_basis_perm(order):
    if order == 3:
        return TetBasis('h1', H1_GSTRUCT[order], efun=h1_edge, ffun=lambda d, a, b, c: -2*a*b*c,
                        cfun=h1_cell, name='H1perm_3')
    return TetBasis('h1', H1_GSTRUCT[order], efun=h1_edge, ffun=h1_face_perm, cfun=h1_cell,
                    name='H1perm_%d' % order)


def hcurl_basis_sym(order, variant='int'):
    return TetBasis('hc', HC_GSTRUCT[order], efun=whitney,
                    ffun=lambda d, a, b, c: hc_face_sym(d, a, b, c, variant), cfun=hc_cell,
                    name='HCsym_%d' % order)
