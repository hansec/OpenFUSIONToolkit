# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../README.md).
# Not part of the regression suite.
"""Prototype H(Curl)_p + Grad(H^1_{p+1}) tetrahedral bases on an OFT mesh dump.

- 'oft': transcription of OFT's hierarchical basis (hcurl_basis.F90 / h1_basis.F90),
  used to validate mesh/orientation/numbering against the dumped OFT mass matrix.
- 'sz':  Schoberl-Zaglmayr basis (integrated Legendre / scaled Jacobi, gradient split),
  same DoF counts per entity (gstruct) and the same global layout.

Basis functions are written as coefficient vectors c_a(lambda) of grad(lambda_a), a=1..4
(cell local vertices), so a vector function is sum_a c_a grad(lambda_a).
"""
import sys
import os
import numpy as np
import scipy.sparse as sp
import sympy as sy

L = sy.symbols('l1:5')
ELEM_HOOK = None  # if a list, assemble() appends (global ids, element matrix) per cell


# ----------------------------------------------------------------------------- polynomials
def leg(n, x):
    return sy.legendre(n, x)


def intleg(n, x):
    """Integrated Legendre L_n(x) = int_{-1}^x P_{n-1}, n>=2"""
    return sy.expand((leg(n, x) - leg(n - 2, x)) / (2 * n - 1))


def scaled(expr_fn, n, x, t):
    """t^n f(x/t) as a polynomial"""
    y = sy.Symbol('y_')
    e = sy.expand(expr_fn(n, y))
    return sy.expand(sy.expand(e.subs(y, x / t) * t**n))


def jac(n, a, x):
    return sy.jacobi(n, a, 0, x)


# ----------------------------------------------------------------------------- OFT basis
def oft_hcurl_face(dof, fc):
    f1, f2, f3 = (L[fc[0]], L[fc[1]], L[fc[2]])
    t = {
        1: (-2*f2*f3, -2*f1*f3, 2*f1*f2),
        2: (f2*f3, -f1*f3, 0),
        3: (2*f2*f3*(f2-f3), 2*f1*f3*(f1-f3), 2*f1*f2*(-f1-f2+2*f3)),
        4: (f2*f3*(-f1-f2+f3), -f1*f3*(-f1-f2+f3), 0),
        5: (2*f2*f3*(-2*f1+f2), 2*f1*f3*(-f1+2*f2), 2*f1*f2*(f1-f2)),
        6: (f2*f3*(2*f1**2-2*f2**2+8*f2*f3-2*f3**2), f1*f3*(-2*f1**2+8*f1*f3+2*f2**2-2*f3**2),
            f1*f2*(2*f1**2+4*f1*f2-16*f1*f3+2*f2**2-16*f2*f3+6*f3**2)),
        7: (f2*f3*(sy.Rational(3, 2)*(f1+f2-f3)**2 - sy.Rational(1, 2)*(f1+f2+f3)**2),
            -f1*f3*(sy.Rational(3, 2)*(f1+f2-f3)**2 - sy.Rational(1, 2)*(f1+f2+f3)**2), 0),
        8: (f2*f3*(2*f1**2+4*f1*f2-4*f1*f3-2*f2**2+2*f2*f3), f1*f3*(2*f1**2-4*f1*f2-2*f1*f3-2*f2**2+4*f2*f3),
            f1*f2*(-2*f1**2+4*f1*f3+2*f2**2-4*f2*f3)),
        9: (f2*f3*(-6*f1**2+12*f1*f2-2*f2**2), f1*f3*(-2*f1**2+12*f1*f2-6*f2**2), f1*f2*(2*f1**2-6*f1*f2+2*f2**2)),
    }[dof]
    c = [0, 0, 0, 0]
    for k in range(3):
        c[fc[k]] = t[k]
    return c


def oft_hcurl_cell(dof):
    f1, f2, f3, f4 = L
    return {
        1: [-2*f2*f3*f4, -2*f1*f3*f4, 2*f1*f2*f4, -2*f1*f2*f3],
        2: [-2*f2*f3*f4, -2*f1*f3*f4, 2*f1*f2*f4, 2*f1*f2*f3],
        3: [f2*f3*f4, -f1*f3*f4, 0, 0],
        4: [f2*f3*f4*(4*f1+2*f2+2*f3-2*f4), f1*f3*f4*(2*f1+4*f2+2*f3-2*f4), f1*f2*f4*(-2*f1-2*f2+2*f4),
            f1*f2*f3*(2*f1+2*f2+2*f3-4*f4)],
        5: [f2*f3*f4*(2*f2+2*f3-2*f4), f1*f3*f4*(2*f1+2*f3-2*f4), f1*f2*f4*(-2*f1-2*f2-4*f3+2*f4),
            f1*f2*f3*(-2*f1-2*f2-2*f3+4*f4)],
        6: [f2*f3*f4*(-f1-f2-f3+f4), f1*f3*f4*(f1+f2+f3-f4), 0, 0],
        7: [f2*f3*f4*(2*f2-2*f3), f1*f3*f4*(2*f1-2*f3), f1*f2*f4*(-2*f1-2*f2+4*f3), f1*f2*f3*(2*f1+2*f2-2*f3)],
        8: [f2*f3*f4*(2*f2-2*f3), f1*f3*f4*(2*f1-2*f3), f1*f2*f4*(-2*f1-2*f2+4*f3), f1*f2*f3*(-2*f1-2*f2+2*f3)],
        9: [f2*f3*f4*(-f1-f2+f3), f1*f3*f4*(f1+f2-f3), 0, 0],
        10: [f2*f3*f4*(-4*f1+2*f2), f1*f3*f4*(-2*f1+4*f2), f1*f2*f4*(2*f1-2*f2), f1*f2*f3*(-2*f1+2*f2)],
        11: [f2*f3*f4*(-4*f1+2*f2), f1*f3*f4*(-2*f1+4*f2), f1*f2*f4*(2*f1-2*f2), f1*f2*f3*(2*f1-2*f2)],
    }[dof]


def oft_h1_edge(dof, ed):
    x1, x2 = L[ed[0]], L[ed[1]]
    return {1: -2*x1*x2, 2: 2*x1*x2*(-x1+x2), 3: 2*x1*x2*(-x1**2+3*x1*x2-x2**2),
            4: 2*x1*x2*(-x1**3+6*x1**2*x2-6*x1*x2**2+x2**3)}[dof]


def oft_h1_face(dof, fc):
    x1, x2, x3 = L[fc[0]], L[fc[1]], L[fc[2]]
    b = 2*x1*x2*x3
    return {1: -b, 2: b*(x1+x2-x3), 3: b*(-x1+x2), 4: b*(-x1**2-2*x1*x2+4*x1*x3-x2**2+4*x2*x3-x3**2),
            5: b*(x1**2-x1*x3-x2**2+x2*x3), 6: b*(-x1**2+3*x1*x2-x2**2)}[dof]


def oft_h1_cell(dof):
    x1, x2, x3, x4 = L
    b = 2*x1*x2*x3*x4
    return {1: -b, 2: b*(x1+x2+x3-x4), 3: b*(x1+x2-x3), 4: b*(-x1+x2)}[dof]


# ----------------------------------------------------------------------------- SZ basis
def sz_u(i, a, b):
    """Edge factor u_i = L^S_{i+2}(l_b - l_a, l_a + l_b), degree i+2"""
    return scaled(intleg, i + 2, L[b] - L[a], L[a] + L[b])


def sz_v(i, j, a, b, c, alpha_off=3):
    """Face factor v_ij = l_c P^S_j^{(2i+alpha_off,0)}(2 l_c - t, t), t=l_a+l_b+l_c, degree j+1"""
    t = L[a] + L[b] + L[c]
    al = 2 * i + alpha_off
    return sy.expand(L[c] * scaled(lambda n, y: jac(n, al, y), j, 2 * L[c] - t, t))


def sz_w(i, j, k, d, alpha_off=6):
    """Cell factor w_ijk = l_d P_k^{(2i+2j+alpha_off,0)}(2 l_d - 1), degree k+1 (sum of l = 1)"""
    t = L[0] + L[1] + L[2] + L[3]
    al = 2 * i + 2 * j + alpha_off
    return sy.expand(L[d] * scaled(lambda n, y: jac(n, al, y), k, 2 * L[d] - t, t))


def grad_coefs(phi):
    return [sy.diff(phi, l) for l in L]


def sz_lists(p, jacobi=True):
    """Return per-entity function generators (index lists) for H(curl) order p and H1 order p+1."""
    q = p + 1
    face_curl = []   # (type, i, j)
    for tot in range(0, p - 1):          # hierarchical by total degree
        for i in range(tot + 1):
            face_curl.append(('t2', i, tot - i))
        face_curl.append(('t3', 0, tot))
    cell_curl = []
    for tot in range(0, p - 2):
        for i in range(tot + 1):
            for j in range(tot - i + 1):
                k = tot - i - j
                cell_curl.append(('t2a', i, j, k))
                cell_curl.append(('t2b', i, j, k))
        for j in range(tot + 1):
            cell_curl.append(('t3', 0, j, tot - j))
    face_h1 = [(i, tot - i) for tot in range(0, q - 2) for i in range(tot + 1)]
    cell_h1 = [(i, j, tot - i - j) for tot in range(0, q - 3) for i in range(tot + 1) for j in range(tot - i + 1)]
    return face_curl, cell_curl, face_h1, cell_h1


def build_local(p, kind, alpha=(3, 6)):
    """Symbolic local functions. Returns dict of lists keyed by entity kind; each item is a function
    of orientation returning the 4 coefficient expressions."""
    q = p + 1
    out = {}
    if kind == 'oft':
        g_fc = {1: 0, 2: 2, 3: 5, 4: 9}[p]
        g_cc = {1: 0, 2: 0, 3: 3, 4: 11}[p]
        out['curl_face'] = [lambda fc, d=d: oft_hcurl_face(d, fc) for d in range(1, g_fc + 1)]
        out['curl_cell'] = [lambda d=d: oft_hcurl_cell(d) for d in range(1, g_cc + 1)]
        out['h1_edge'] = [lambda ed, d=d: grad_coefs(oft_h1_edge(d, ed)) for d in range(1, q)]
        g_hf = {1: 0, 2: 0, 3: 1, 4: 3, 5: 6}[q]
        g_hc = {1: 0, 2: 0, 3: 0, 4: 1, 5: 4}[q]
        out['h1_face'] = [lambda fc, d=d: grad_coefs(oft_h1_face(d, fc)) for d in range(1, g_hf + 1)]
        out['h1_cell'] = [lambda d=d: grad_coefs(oft_h1_cell(d)) for d in range(1, g_hc + 1)]
    else:
        fcl, ccl, fh, ch = sz_lists(p)
        af, ac = alpha

        def face_curl(fc, item):
            a, b, c = fc
            typ, i, j = item
            if typ == 't2':
                u = sz_u(i, a, b)
                v = sz_v(i, j, a, b, c, af)
                gu, gv = grad_coefs(u), grad_coefs(v)
                return [gu[m] * v - u * gv[m] for m in range(4)]
            v = sz_v(0, j, a, b, c, af)
            co = [0, 0, 0, 0]
            co[b] = L[a] * v
            co[a] = -L[b] * v
            return co

        def cell_curl(item):
            typ, i, j, k = item
            if typ == 't3':
                v = sz_v(0, j, 0, 1, 2, af)
                w = sz_w(0, j, k, 3, ac)
                co = [0, 0, 0, 0]
                co[1] = L[0] * v * w
                co[0] = -L[1] * v * w
                return co
            u = sz_u(i, 0, 1)
            v = sz_v(i, j, 0, 1, 2, af)
            w = sz_w(i, j, k, 3, ac)
            gu, gv, gw = grad_coefs(u), grad_coefs(v), grad_coefs(w)
            s = (1, -1, 1) if typ == 't2a' else (1, 1, -1)
            return [s[0]*gu[m]*v*w + s[1]*u*gv[m]*w + s[2]*u*v*gw[m] for m in range(4)]

        out['curl_face'] = [lambda fc, it=it: face_curl(fc, it) for it in fcl]
        out['curl_cell'] = [lambda it=it: cell_curl(it) for it in ccl]
        out['h1_edge'] = [lambda ed, i=i: grad_coefs(sz_u(i, ed[0], ed[1])) for i in range(0, q - 1)]
        out['h1_face'] = [lambda fc, it=it: grad_coefs(sy.expand(sz_u(it[0], fc[0], fc[1]) *
                                                                 sz_v(it[0], it[1], fc[0], fc[1], fc[2], af)))
                          for it in fh]
        out['h1_cell'] = [lambda it=it: grad_coefs(sy.expand(sz_u(it[0], 0, 1) * sz_v(it[0], it[1], 0, 1, 2, af) *
                                                             sz_w(it[0], it[1], it[2], 3, ac)))
                          for it in ch]
    return out


# ----------------------------------------------------------------------------- mesh / quadrature
def tet_quad(n):
    """Collapsed Gauss rule on the reference tet, returns barycentric points [nq,4], weights summing to 1."""
    from numpy.polynomial.legendre import leggauss
    x, w = leggauss(n)
    x = 0.5 * (x + 1); w = 0.5 * w
    pts, wts = [], []
    for a, wa in zip(x, w):
        for b, wb in zip(x, w):
            for c, wc in zip(x, w):
                xi = a
                eta = b * (1 - a)
                zeta = c * (1 - a) * (1 - b)
                jac_ = (1 - a)**2 * (1 - b)
                l2, l3, l4 = xi, eta, zeta
                pts.append([1 - l2 - l3 - l4, l2, l3, l4])
                wts.append(wa * wb * wc * jac_ * 6.0)
    return np.array(pts), np.array(wts)


def load_mesh(d):
    with open(os.path.join(d, 'mesh.txt')) as f:
        np_, nc = map(int, f.readline().split())
        r = np.array([list(map(float, f.readline().split())) for _ in range(np_)])
        cells = np.array([list(map(int, f.readline().split())) for _ in range(nc)])
    lc = cells[:, 0:4] - 1
    lce = cells[:, 4:10] - 1
    oe = cells[:, 10:22].reshape(nc, 6, 2) - 1      # Fortran (2,6) column-major -> [edge][2]
    lcf = cells[:, 22:26] - 1
    of = cells[:, 26:38].reshape(nc, 4, 3) - 1      # (3,4) -> [face][3]
    ne = lce.max() + 1
    nf = lcf.max() + 1
    return r, lc, lce, oe, lcf, of, np_, ne, nf, nc


def assemble(d, p, kind, alpha=(3, 6), field=None):
    r, lc, lce, oe, lcf, of, np_, ne, nf, nc = load_mesh(d)
    loc = build_local(p, kind, alpha)
    q = p + 1
    g_c = [0, 1, len(loc['curl_face']), len(loc['curl_cell'])]
    g_h = [1, q - 1, len(loc['h1_face']), len(loc['h1_cell'])]
    n_curl = g_c[1] * ne + g_c[2] * nf + g_c[3] * nc
    n_h1 = g_h[0] * np_ + g_h[1] * ne + g_h[2] * nf + g_h[3] * nc
    pts, wts = tet_quad(p + 3)
    # cache lambdified coefficient functions per orientation tuple
    cache = {}

    def lam(key, fn):
        if key not in cache:
            co = fn()
            cache[key] = sy.lambdify(L, [sy.sympify(cc) for cc in co], 'numpy')
        return cache[key]

    def ev(f):
        vals = f(*pts.T)
        return np.array([np.broadcast_to(np.asarray(v, dtype=float), (len(wts),)) for v in vals])  # [4,nq]

    rows, cols, vals = [], [], []
    bvec = None
    if field is not None:
        n_tot = n_curl + n_h1
        bvec = np.zeros(n_tot)
    for c in range(nc):
        X = r[lc[c]]
        Jm = (X[1:] - X[0]).T
        vol = abs(np.linalg.det(Jm)) / 6.0
        Ginv = np.linalg.inv(Jm)          # rows: grad of l2,l3,l4 wrt x
        gl = np.vstack([-Ginv.sum(axis=0), Ginv])   # [4,3]
        G = gl @ gl.T
        funcs, gidx = [], []
        # H(curl) edges (Whitney)
        for k in range(6):
            e1, e2 = oe[c, k]
            co = np.zeros((4, len(wts)))
            co[e1] = -pts[:, e2]
            co[e2] = pts[:, e1]
            funcs.append(co); gidx.append(lce[c, k])
        for k in range(4):
            fc = tuple(of[c, k])
            for m, fn in enumerate(loc['curl_face']):
                funcs.append(ev(lam(('cf', fc, m), lambda fn=fn, fc=fc: fn(fc))))
                gidx.append(g_c[1] * ne + lcf[c, k] * g_c[2] + m)
        for m, fn in enumerate(loc['curl_cell']):
            funcs.append(ev(lam(('cc', m), fn)))
            gidx.append(g_c[1] * ne + g_c[2] * nf + c * g_c[3] + m)
        # H1 gradients (vertex DoFs included, removed later)
        for k in range(4):
            co = np.zeros((4, len(wts))); co[k] = 1.0
            funcs.append(co); gidx.append(n_curl + lc[c, k])
        for k in range(6):
            ed = tuple(oe[c, k])
            for m, fn in enumerate(loc['h1_edge']):
                funcs.append(ev(lam(('he', ed, m), lambda fn=fn, ed=ed: fn(ed))))
                gidx.append(n_curl + np_ + lce[c, k] * g_h[1] + m)
        for k in range(4):
            fc = tuple(of[c, k])
            for m, fn in enumerate(loc['h1_face']):
                funcs.append(ev(lam(('hf', fc, m), lambda fn=fn, fc=fc: fn(fc))))
                gidx.append(n_curl + np_ + g_h[1] * ne + lcf[c, k] * g_h[2] + m)
        for m, fn in enumerate(loc['h1_cell']):
            funcs.append(ev(lam(('hc', m), fn)))
            gidx.append(n_curl + np_ + g_h[1] * ne + g_h[2] * nf + c * g_h[3] + m)
        F = np.array(funcs)                         # [nf,4,nq]
        V = np.einsum('iaq,ab->ibq', F, G)
        Me = vol * np.einsum('ibq,jbq,q->ij', V, F, wts)
        gi = np.array(gidx)
        if ELEM_HOOK is not None:
            ELEM_HOOK.append((gi.copy(), Me.copy()))
        if field is not None:
            xq = pts @ X                                  # [nq,3]
            Fq = field(xq)                                # [nq,3]
            phys = np.einsum('iaq,ax->iqx', F, gl)        # [nf,nq,3]
            np.add.at(bvec, gi, vol * np.einsum('iqx,qx,q->i', phys, Fq, wts))
        rows.append(np.repeat(gi, len(gi))); cols.append(np.tile(gi, len(gi))); vals.append(Me.ravel())
    n = n_curl + n_h1
    M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    if field is not None:
        return M, (n_curl, n_h1, g_c, g_h, np_, ne, nf, nc), bvec
    return M, (n_curl, n_h1, g_c, g_h, np_, ne, nf, nc)


if __name__ == '__main__':
    d = sys.argv[1]; p = int(sys.argv[2])
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from analyze import load
    Mo, info, colors, cells = load(d)
    Mp, lay = assemble(d, p, 'oft')
    diff = abs(Mp - Mo).max()
    print('n', Mo.shape[0], Mp.shape[0], 'max|M_proto - M_oft| =', diff, 'max|M| =', abs(Mo).max())
