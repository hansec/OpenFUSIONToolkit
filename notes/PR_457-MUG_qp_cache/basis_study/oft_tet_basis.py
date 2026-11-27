#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only (basis study for the REFERENCE-mode quadrature cache).
#
# Symbolic transcription of OFT's tetrahedral H^1 (src/fem/h1_basis.F90, oft_h1_eval*)
# and H(Curl) (src/fem/hcurl_basis.F90, oft_hcurl_eval*) bases, plus prototype
# modified ("symmetric") face families. Nothing here is used by OFT itself.
#
# Conventions follow OFT exactly:
#  - logical coordinates f1..f4 (barycentric), unit logical frame of the qcache
#    module: grad f2,f3,f4 = e1,e2,e3 and grad f1 = -(e1+e2+e3);
#  - tet_ed / tet_fc from src/grid/tetmesh_type.F90;
#  - an edge in state 1 keeps cell_ed order, state 2 reverses it; a face in
#    state s uses code = s (s<=3) or 3-s (s>3) and orient_listn (reverse if
#    code<0, then CSHIFT by |code|-1), as in fem_quad_cache.F90.
#------------------------------------------------------------------------------
import itertools
import numpy as np
import sympy as sp

F = sp.symbols('f1:5')
TET_ED = [(1, 4), (2, 4), (3, 4), (2, 3), (3, 1), (1, 2)]
TET_FC = [(2, 3, 4), (3, 1, 4), (1, 2, 4), (1, 2, 3)]
GRADS = np.array([[-1., -1., -1.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])


def orient_edge(ed, state):
    return tuple(ed) if state == 1 else (ed[1], ed[0])


def face_code(state, np_face=3):
    return state if state <= np_face else np_face - state


def orient_face(fc, state):
    code = face_code(state)
    lst = list(fc)
    if code < 0:
        lst = lst[::-1]
    k = abs(code) - 1
    return tuple(lst[k:] + lst[:k])


#------------------------------------------------------------------------------
# H^1 (oft_h1_evale / oft_h1_evalf / oft_h1_evalc)
#------------------------------------------------------------------------------
def h1_edge(dof, x1, x2):
    return [-2*x1*x2,
            2*x1*x2*(-x1 + x2),
            2*x1*x2*(-x1**2 + 3*x1*x2 - x2**2),
            2*x1*x2*(-x1**3 + 6*x1**2*x2 - 6*x1*x2**2 + x2**3)][dof-1]


def h1_face(dof, x1, x2, x3):
    b = x1*x2*x3
    return [-2*b,
            2*b*(x1 + x2 - x3),
            2*b*(-x1 + x2),
            2*b*(-x1**2 - 2*x1*x2 + 4*x1*x3 - x2**2 + 4*x2*x3 - x3**2),
            2*b*(x1**2 - x1*x3 - x2**2 + x2*x3),
            2*b*(-x1**2 + 3*x1*x2 - x2**2)][dof-1]


def h1_cell(dof, x1, x2, x3, x4):
    b = x1*x2*x3*x4
    return [-2*b, 2*b*(x1 + x2 + x3 - x4), 2*b*(x1 + x2 - x3), 2*b*(-x1 + x2)][dof-1]


H1_GSTRUCT = {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1, 2, 1, 0), 4: (1, 3, 3, 1), 5: (1, 4, 6, 4)}

#------------------------------------------------------------------------------
# H(Curl): functions are dicts {vertex: y_v} meaning sum_v y_v grad(f_v)
#------------------------------------------------------------------------------
def hc_face(dof, f1, f2, f3):
    s = f1 + f2 + f3
    y = [(-2*f2*f3, -2*f1*f3, 2*f1*f2),
         (f2*f3, -f1*f3, 0),
         (2*f2*f3*(f2 - f3), 2*f1*f3*(f1 - f3), 2*f1*f2*(-f1 - f2 + 2*f3)),
         (f2*f3*(-f1 - f2 + f3), -f1*f3*(-f1 - f2 + f3), 0),
         (2*f2*f3*(-2*f1 + f2), 2*f1*f3*(-f1 + 2*f2), 2*f1*f2*(f1 - f2)),
         (f2*f3*(2*f1**2 - 2*f2**2 + 8*f2*f3 - 2*f3**2),
          f1*f3*(-2*f1**2 + 8*f1*f3 + 2*f2**2 - 2*f3**2),
          f1*f2*(2*f1**2 + 4*f1*f2 - 16*f1*f3 + 2*f2**2 - 16*f2*f3 + 6*f3**2)),
         (f2*f3*(sp.Rational(3, 2)*(f1 + f2 - f3)**2 - sp.Rational(1, 2)*s**2),
          -f1*f3*(sp.Rational(3, 2)*(f1 + f2 - f3)**2 - sp.Rational(1, 2)*s**2), 0),
         (f2*f3*(2*f1**2 + 4*f1*f2 - 4*f1*f3 - 2*f2**2 + 2*f2*f3),
          f1*f3*(2*f1**2 - 4*f1*f2 - 2*f1*f3 - 2*f2**2 + 4*f2*f3),
          f1*f2*(-2*f1**2 + 4*f1*f3 + 2*f2**2 - 4*f2*f3)),
         (f2*f3*(-6*f1**2 + 12*f1*f2 - 2*f2**2),
          f1*f3*(-2*f1**2 + 12*f1*f2 - 6*f2**2),
          f1*f2*(2*f1**2 - 6*f1*f2 + 2*f2**2))][dof-1]
    return y


def hc_cell(dof, f1, f2, f3, f4):
    y = [(-2*f2*f3*f4, -2*f1*f3*f4, 2*f1*f2*f4, -2*f1*f2*f3),
         (-2*f2*f3*f4, -2*f1*f3*f4, 2*f1*f2*f4, 2*f1*f2*f3),
         (f2*f3*f4, -f1*f3*f4, 0, 0),
         (f2*f3*f4*(4*f1 + 2*f2 + 2*f3 - 2*f4), f1*f3*f4*(2*f1 + 4*f2 + 2*f3 - 2*f4),
          f1*f2*f4*(-2*f1 - 2*f2 + 2*f4), f1*f2*f3*(2*f1 + 2*f2 + 2*f3 - 4*f4)),
         (f2*f3*f4*(2*f2 + 2*f3 - 2*f4), f1*f3*f4*(2*f1 + 2*f3 - 2*f4),
          f1*f2*f4*(-2*f1 - 2*f2 - 4*f3 + 2*f4), f1*f2*f3*(-2*f1 - 2*f2 - 2*f3 + 4*f4)),
         (f2*f3*f4*(-f1 - f2 - f3 + f4), f1*f3*f4*(f1 + f2 + f3 - f4), 0, 0),
         (f2*f3*f4*(2*f2 - 2*f3), f1*f3*f4*(2*f1 - 2*f3),
          f1*f2*f4*(-2*f1 - 2*f2 + 4*f3), f1*f2*f3*(2*f1 + 2*f2 - 2*f3)),
         (f2*f3*f4*(2*f2 - 2*f3), f1*f3*f4*(2*f1 - 2*f3),
          f1*f2*f4*(-2*f1 - 2*f2 + 4*f3), f1*f2*f3*(-2*f1 - 2*f2 + 2*f3)),
         (f2*f3*f4*(-f1 - f2 + f3), f1*f3*f4*(f1 + f2 - f3), 0, 0),
         (f2*f3*f4*(-4*f1 + 2*f2), f1*f3*f4*(-2*f1 + 4*f2),
          f1*f2*f4*(2*f1 - 2*f2), f1*f2*f3*(-2*f1 + 2*f2)),
         (f2*f3*f4*(-4*f1 + 2*f2), f1*f3*f4*(-2*f1 + 4*f2),
          f1*f2*f4*(2*f1 - 2*f2), f1*f2*f3*(2*f1 - 2*f2))][dof-1]
    return y


HC_GSTRUCT = {1: (0, 1, 0, 0), 2: (0, 1, 2, 0), 3: (0, 1, 5, 3), 4: (0, 1, 9, 11)}


#------------------------------------------------------------------------------
# Generic cell basis description
#------------------------------------------------------------------------------
class TetBasis:
    """A tet basis: per-entity DOF formula generators.

    kind 'h1' (scalar expression) or 'hc' (dict vertex -> coefficient of grad f_v).
    vfun(v) / efun(dof, x1, x2) / ffun(dof, x1, x2, x3) / cfun(dof, f1..f4) return
    either a sympy expression ('h1') or a tuple of coefficients on the entity's
    (oriented) vertices ('hc').
    """

    def __init__(self, kind, gstruct, vfun=None, efun=None, ffun=None, cfun=None, name=''):
        self.kind = kind
        self.gstruct = gstruct
        self.vfun, self.efun, self.ffun, self.cfun = vfun, efun, ffun, cfun
        self.name = name
        # entities: 4 vertices, 6 edges, 4 faces, 1 interior
        self.ent_ndof = [gstruct[0]]*4 + [gstruct[1]]*6 + [gstruct[2]]*4 + [gstruct[3]]
        self.ent_dof0 = list(np.cumsum([0] + self.ent_ndof[:-1]))
        self.nce = sum(self.ent_ndof)
        self.ent_nst = [1]*4 + [2]*6 + [6]*4 + [1]

    def entity_funcs(self, ent, state=1):
        """Symbolic functions of entity `ent` (0-based, OFT order) in a given state."""
        out = []
        nd = self.ent_ndof[ent]
        for j in range(1, nd+1):
            if ent < 4:
                v = ent + 1
                expr = F[v-1] if self.kind == 'h1' else None
                out.append(expr)
            elif ent < 10:
                ed = orient_edge(TET_ED[ent-4], state)
                y = self.efun(j, F[ed[0]-1], F[ed[1]-1])
                out.append(y if self.kind == 'h1' else dict(zip(ed, y)))
            elif ent < 14:
                fc = orient_face(TET_FC[ent-10], state)
                y = self.ffun(j, *[F[v-1] for v in fc])
                out.append(y if self.kind == 'h1' else dict(zip(fc, y)))
            else:
                y = self.cfun(j, *F)
                out.append(y if self.kind == 'h1' else dict(zip((1, 2, 3, 4), y)))
        return out


def h1_basis(order):
    return TetBasis('h1', H1_GSTRUCT[order], efun=h1_edge, ffun=h1_face, cfun=h1_cell, name='H1_%d' % order)


def whitney(dof, x1, x2):
    return (-x2, x1)


def hcurl_basis(order, ffun=hc_face, cfun=hc_cell, name=None):
    return TetBasis('hc', HC_GSTRUCT[order], efun=whitney, ffun=ffun, cfun=cfun,
                    name=name or 'HCurl_%d' % order)


#------------------------------------------------------------------------------
# Numerical evaluation
#------------------------------------------------------------------------------
_lamb_cache = {}


def _lamb(expr):
    key = sp.srepr(expr)
    if key not in _lamb_cache:
        _lamb_cache[key] = sp.lambdify(F, expr, 'numpy')
    return _lamb_cache[key]


def _evalexpr(expr, pts):
    v = _lamb(sp.sympify(expr))(*[pts[:, i] for i in range(4)])
    return np.broadcast_to(np.asarray(v, dtype=float), (pts.shape[0],)).copy()


def eval_h1(expr, pts):
    """Value [N] and gradient [N,3] of a scalar formula."""
    val = _evalexpr(expr, pts)
    grad = np.zeros((pts.shape[0], 3))
    for i in range(4):
        d = sp.diff(expr, F[i])
        if d != 0:
            grad += np.outer(_evalexpr(d, pts), GRADS[i])
    return val, grad


def eval_hc(ydict, pts):
    """Value [N,3] and curl [N,3] of sum_v y_v grad f_v."""
    val = np.zeros((pts.shape[0], 3))
    curl = np.zeros((pts.shape[0], 3))
    for v, y in ydict.items():
        y = sp.sympify(y)
        if y == 0:
            continue
        val += np.outer(_evalexpr(y, pts), GRADS[v-1])
        for k in range(4):
            d = sp.diff(y, F[k])
            if d != 0:
                curl += np.outer(_evalexpr(d, pts), np.cross(GRADS[k], GRADS[v-1]))
    return val, curl


def sample_matrix(funcs, kind, pts, grad_only=False):
    """Stack [value; derivative] samples of a list of functions as columns.

    kind 'h1': rows = [val (N), grad (3N)] (or grad only); 'hc': [val (3N), curl (3N)].
    'grad' : H^1 formulas represented by their gradient as a vector field: [grad (3N), 0 (3N)]
    """
    cols = []
    for fn in funcs:
        if kind == 'h1':
            v, g = eval_h1(fn, pts)
            cols.append(g.ravel() if grad_only else np.concatenate([v, g.ravel()]))
        elif kind == 'grad':
            v, g = eval_h1(fn, pts)
            cols.append(np.concatenate([g.ravel(), np.zeros(3*pts.shape[0])]))
        else:
            v, c = eval_hc(fn, pts)
            cols.append(np.concatenate([v.ravel(), c.ravel()]))
    return np.array(cols).T if cols else np.zeros((0, 0))


def random_tet_points(n, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.dirichlet(np.ones(4), size=n)
    return x


def all_face_perms():
    """The 6 face states and the vertex permutation (position map) each realizes."""
    out = []
    for s in range(1, 7):
        out.append((s, orient_face((1, 2, 3), s)))
    return out
