#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# EXPERIMENTAL, analysis only: verify the prototype symmetric face families in
# sym_basis.py against OFT's current tet bases (oft_tet_basis.py):
#   1. spaces: H^1 span unchanged; H(Curl)+grad H^1 span unchanged; new H(Curl)
#      still gradient-free (H(Curl) + grad H^1 is a direct sum apart from the
#      Whitney/vertex overlap); H(Curl)-only span changes by gradients only;
#   2. orientation transforms (Python port of qcache_build_xform) and their block
#      structure per hierarchical level;
#   3. conditioning of local mass/stiffness matrices on the unit and a regular tet.
#
# Usage: python verify_sym.py [--quick]   (writes results_sym.txt next to this file)
#------------------------------------------------------------------------------
import os
import sys
import numpy as np
import oft_tet_basis as otb
from oft_tet_basis import h1_basis, hcurl_basis, random_tet_points, sample_matrix
from xform_classify import Classifier
from sym_basis import h1_basis_sym, h1_basis_perm, hcurl_basis_sym
from sym_orth import h1_basis_orth, hcurl_basis_orth

OUT = []


def log(msg=''):
    print(msg)
    OUT.append(msg)


def canon(basis, pts, kind=None, skip_vertices=False):
    kind = kind or ('h1' if basis.kind == 'h1' else 'hc')
    cols = []
    for ent in range(15):
        if skip_vertices and ent < 4:
            continue
        fn = basis.entity_funcs(ent, 1)
        if fn:
            cols.append(sample_matrix(fn, kind, pts))
    return np.hstack(cols)


def rank(A, tol=1.e-9):
    s = np.linalg.svd(A, compute_uv=False)
    return int(np.sum(s > tol*s[0]))


#------------------------------------------------------------------------------
# 1. Space checks
#------------------------------------------------------------------------------
def space_checks(pts):
    _space_checks(pts)
    log()


def _space_checks(pts):
    log('## 1. Space checks (ranks of sampled function sets, values+derivatives)')
    for k in range(1, 6):
        old = canon(h1_basis(k), pts)
        for nb in (h1_basis_sym(k), h1_basis_orth(k), h1_basis_perm(k) if k >= 3 else None):
            if nb is None:
                continue
            new = canon(nb, pts)
            log('H1 order %d %-9s: nce=%2d rank(old)=%2d rank(new)=%2d rank(old|new)=%2d -> %s' % (
                k, nb.name.split('_')[0], old.shape[1], rank(old), rank(new), rank(np.hstack([old, new])),
                'same span' if rank(old) == rank(new) == rank(np.hstack([old, new])) == old.shape[1] else 'DIFFERENT'))
    for p in range(1, 5):
        h1 = h1_basis(p+1)
        g = canon(h1, pts, kind='grad', skip_vertices=True)
        old = canon(hcurl_basis(p), pts)
        for nbas in (hcurl_basis_sym(p), hcurl_basis_orth(p)):
            hc_space_line(p, old, canon(nbas, pts), g, nbas.name)


def hc_space_line(p, old, new, g, name):
    r_comb_old = rank(np.hstack([old, g]))
    r_comb_new = rank(np.hstack([new, g]))
    r_all = rank(np.hstack([old, new, g]))
    ned2 = (p+1)*(p+2)*(p+3)//2
    log(name + ' p=%d: nce=%2d rank(HC_old)=%2d rank(HC_new)=%2d rank(HC_old|HC_new)=%2d; '
        'rank(HC+gradH1)= old %3d new %3d union %3d (NED2_p dim %3d, direct sum %3d)' % (
            p, old.shape[1], rank(old), rank(new), rank(np.hstack([old, new])), r_comb_old, r_comb_new,
            r_all, ned2, new.shape[1]+g.shape[1]))


#------------------------------------------------------------------------------
# 2. Orientation transforms
#------------------------------------------------------------------------------
def levels_of(basis, ent):
    nd = basis.ent_ndof[ent]
    if ent >= 10 and ent < 14:
        cuts = [0, 1, 3, 6] if basis.kind == 'h1' else [0, 2, 5, 9]
        return [(cuts[i], min(cuts[i+1], nd)) for i in range(3) if cuts[i] < nd]
    return [(0, nd)]


def is_signed_perm(X):
    if X.shape[0] != X.shape[1]:
        return False
    nz = np.abs(X) > 1.e-10
    if not (np.all(nz.sum(axis=0) == 1) and np.all(nz.sum(axis=1) == 1)):
        return False
    return np.allclose(np.abs(X[nz]), 1.)


def xform_report(basis, pts, extra=None):
    c = Classifier(basis, pts, extra=extra)
    counts = c.classify()
    own_only = all(d['nreg'] == basis.ent_ndof[d['ent']] for d in c.details if d['kind'] in (0, 1, 2))
    nblocks = 0
    nperm_levels = 0
    max_dense = 0
    for d in c.details:
        if d['kind'] != 2 or d['nreg'] != basis.ent_ndof[d['ent']]:
            continue
        X = d['X']
        lv = levels_of(basis, d['ent'])
        # block diagonal by hierarchical level?
        mask = np.zeros_like(X, dtype=bool)
        for (a, b) in lv:
            mask[a:b, a:b] = True
        if np.max(np.abs(X[~mask]), initial=0.) > 1.e-10:
            max_dense = max(max_dense, X.shape[0])
            continue
        for (a, b) in lv:
            nblocks += 1
            if is_signed_perm(X[a:b, a:b]):
                nperm_levels += 1
            else:
                max_dense = max(max_dense, b-a)
    line = '%-10s identity/perm/dense/none = %s; own-block only: %s; max regressors %d' % (
        basis.name, list(counts), own_only, max([d['nreg'] for d in c.details if d['kind'] < 3]+[0]))
    if counts[2] > 0:
        line += '; dense pairs: level blocks signed-perm %d/%d, largest non-permutation block %dx%d' % (
            nperm_levels, nblocks, max_dense, max_dense)
    if counts[3] > 0 and extra is not None:
        line += '; kind-3 fit with grad H1 added: max resid %.1e' % max(
            d['xres'] for d in c.details if d['xres'] is not None)
    log(line)
    return c


def example_transforms(basis, pts, ent=13):
    c = Classifier(basis, pts)
    c.classify()
    log('Face transforms of %s (face 4, own block, columns = oriented DOFs):' % basis.name)
    for d in c.details:
        if d['ent'] == ent and d['kind'] < 3:
            X = d['X']
            log('  state %d: %s' % (d['state'], np.array2string(np.round(X, 6), separator=',').replace('\n', '')))


#------------------------------------------------------------------------------
# 3. Conditioning
#------------------------------------------------------------------------------
def tet_quadrature(n):
    x, w = np.polynomial.legendre.leggauss(n)
    x = 0.5*(x+1.)
    w = 0.5*w
    pts, wts = [], []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                u, v, t = x[i], x[j], x[k]
                X = u
                Y = v*(1.-u)
                Z = t*(1.-u)*(1.-v)
                pts.append((1.-X-Y-Z, X, Y, Z))
                wts.append(w[i]*w[j]*w[k]*(1.-u)**2*(1.-v))
    return np.array(pts), np.array(wts)


def set_geometry(verts):
    """Barycentric gradients of a tet with vertices verts[4,3] (OFT logical order)."""
    T = np.array([verts[1]-verts[0], verts[2]-verts[0], verts[3]-verts[0]]).T
    Ti = np.linalg.inv(T)
    g = np.zeros((4, 3))
    g[1:, :] = Ti
    g[0, :] = -Ti.sum(axis=0)
    otb.GRADS[:] = g


def cond_spd(M, nzero=0):
    """Eigenvalue ratio after symmetric Jacobi (diagonal) scaling, zero modes removed."""
    d = np.sqrt(np.abs(np.diag(M)))
    d[d == 0.] = 1.
    M = M/np.outer(d, d)
    e = np.linalg.eigvalsh(0.5*(M+M.T))
    e = np.sort(e)[nzero:]
    return e[-1]/e[0]


def local_mats(basis, qp, qw, which):
    n = qp.shape[0]
    cols = []
    for ent in range(15):
        if which == 'gradnv' and ent < 4:
            continue
        for fn in basis.entity_funcs(ent, 1):
            if basis.kind == 'h1':
                v, g = otb.eval_h1(fn, qp)
                cols.append((v[:, None], g))
            else:
                v, cu = otb.eval_hc(fn, qp)
                cols.append((v, cu))
    V = np.stack([c[0] for c in cols], axis=-1)   # n, ncomp, nb
    D = np.stack([c[1] for c in cols], axis=-1)   # n, 3, nb
    if which == 'gradnv':
        V = D
    M = np.einsum('q,qci,qcj->ij', qw, V, V)
    K = np.einsum('q,qci,qcj->ij', qw, D, D)
    return M, K


def conditioning(qp, qw):
    log('## 3. Conditioning of local matrices (eigenvalue ratio after Jacobi scaling; zero modes removed)')
    geoms = {'unit': np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
             'regular': np.array([[1., 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])}
    for gname, verts in geoms.items():
        set_geometry(verts)
        log('Geometry: %s tet' % gname)
        for k in range(2, 6):
            row = []
            cands = [('old', h1_basis(k)), ('sym-int', h1_basis_sym(k)), ('sym-orth', h1_basis_sym(k, 'orth')),
                     ('sym-proj', h1_basis_orth(k))]
            if k >= 3:
                cands.append(('perm', h1_basis_perm(k)))
            for lab, b in cands:
                M, K = local_mats(b, qp, qw, 'all')
                row.append('%s M %.3e K %.3e' % (lab, cond_spd(M), cond_spd(K, 1)))
            log('  H1 k=%d: ' % k + ' | '.join(row))
        for p in range(1, 5):
            row = []
            row2 = []
            for lab, b in (('old', hcurl_basis(p)), ('sym-int', hcurl_basis_sym(p)), ('sym-orth', hcurl_basis_sym(p, 'orth')),
                           ('sym-proj', hcurl_basis_orth(p))):
                M, C = local_mats(b, qp, qw, 'all')
                row.append('%s M %.3e curlcurl %.3e' % (lab, cond_spd(M), cond_spd(C, 3)))
                row2.append('%s %.3e' % (lab, cond_spd(combined_mass(b, h1_basis(p+1), qp, qw))))
            log('  HCurl p=%d: ' % p + ' | '.join(row))
            log('  HCurl+gradH1 p=%d mass: ' % p + ' | '.join(row2))
    set_geometry(np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    log()


def combined_mass(hc, h1, qp, qw):
    vals = []
    for ent in range(15):
        for fn in hc.entity_funcs(ent, 1):
            v, _ = otb.eval_hc(fn, qp)
            vals.append(v)
    for ent in range(4, 15):
        for fn in h1.entity_funcs(ent, 1):
            _, g = otb.eval_h1(fn, qp)
            vals.append(g)
    V = np.stack(vals, axis=-1)
    return np.einsum('q,qci,qcj->ij', qw, V, V)


if __name__ == '__main__':
    quick = '--quick' in sys.argv
    pts = random_tet_points(40)
    space_checks(pts)
    log('## 2. Orientation transforms (Python port of qcache_build_xform)')
    for k in range(2, 6):
        xform_report(h1_basis(k), pts)
        xform_report(h1_basis_sym(k), pts)
        xform_report(h1_basis_orth(k), pts)
        if k >= 3:
            xform_report(h1_basis_perm(k), pts)
    for p in range(1, 5):
        xform_report(hcurl_basis(p), pts, extra=[(h1_basis(p+1), 'grad')])
        xform_report(hcurl_basis_sym(p), pts, extra=[(h1_basis(p+1), 'grad')])
        xform_report(hcurl_basis_orth(p), pts, extra=[(h1_basis(p+1), 'grad')])
    log()
    example_transforms(hcurl_basis_sym(2), pts)
    example_transforms(hcurl_basis_sym(2, 'orth'), pts)
    example_transforms(h1_basis_sym(4), pts)
    log()
    if not quick:
        qp, qw = tet_quadrature(8)
        conditioning(qp, qw)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_sym.txt'), 'w') as fid:
        fid.write('\n'.join(OUT)+'\n')
