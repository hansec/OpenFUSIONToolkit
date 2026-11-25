# SPDX-License-Identifier: LGPL-3.0-only
# Usage: python3 schur_analysis.py <run_dir containing xmhd_pb_J.h5, xmhd_pb_P.h5, xmhd_pb_fields.txt>
# Offline analysis of PB preconditioner: exact discrete Schur complement vs assembled parabolized S_V
import sys
import h5py
import numpy as np
import scipy.sparse as sp

d = sys.argv[1]


def load(fn, path):
    with h5py.File(fn, "r") as f:
        g = f[path]
        nr = int(np.array(g["nr"]).ravel()[0])
        nc = int(np.array(g["nc"]).ravel()[0])
        kr = np.array(g["kr"]).ravel().astype(np.int64)
        lc = np.array(g["lc"]).ravel().astype(np.int64)
        M = np.array(g["M"]).ravel()
        lg = np.array(g["lg"]).ravel().astype(np.int64)
    A = sp.csr_matrix((M, lc - 1, kr - 1), shape=(nr, nc)).tocoo()
    # fold duplicated (periodic/shared) local rows/cols onto global indices
    ng = lg.max()
    G = sp.coo_matrix((A.data, (lg[A.row] - 1, lg[A.col] - 1)), shape=(ng, ng)).tocsr()
    return G, lg


Jsp, lg = load(d + "/xmhd_pb_J.h5", "J")
Psp, _ = load(d + "/xmhd_pb_P.h5", "P")
J = Jsp.toarray()
P = Psp.toarray()
fl = np.loadtxt(d + "/xmhd_pb_fields.txt", ndmin=2).astype(int)
idx = {1: set(), 2: set(), 3: set()}
for i, off, n, g in fl:
    idx[g] |= set((lg[off : off + n] - 1).tolist())
idx = {k: sorted(v) for k, v in idx.items()}
B, V, S = [np.array(idx[k], dtype=int) for k in (1, 2, 3)]
BS = np.concatenate([B, S])
print("sizes B V S", len(B), len(V), len(S))
print("||P_BB-J_BB||/||J_BB|| =", np.linalg.norm(P[np.ix_(B, B)] - J[np.ix_(B, B)]) / np.linalg.norm(J[np.ix_(B, B)]))
print("||P_SS-J_SS|| =", np.linalg.norm(P[np.ix_(S, S)] - J[np.ix_(S, S)]) if len(S) else 0)
Jbs = J[np.ix_(BS, BS)]
Sex = J[np.ix_(V, V)] - J[np.ix_(V, BS)] @ np.linalg.solve(Jbs, J[np.ix_(BS, V)])
Sv = P[np.ix_(V, V)]
dW_ex = Sex - J[np.ix_(V, V)]
dW_ap = Sv - J[np.ix_(V, V)]
print(
    "||dW_approx||/||dW_exact|| =",
    np.linalg.norm(dW_ap) / np.linalg.norm(dW_ex),
    " rel diff",
    np.linalg.norm(dW_ap - dW_ex) / np.linalg.norm(dW_ex),
)
ev = np.linalg.eigvals(np.linalg.solve(Sv, Sex))
print(
    "eig(Sv^-1 Sex): min|.|=%.3e max|.|=%.3e  real range [%.3e, %.3e]  max|imag| %.2e"
    % (abs(ev).min(), abs(ev).max(), ev.real.min(), ev.real.max(), abs(ev.imag).max())
)
q = np.quantile(abs(ev), [0.01, 0.1, 0.5, 0.9, 0.99])
print("   |eig| quantiles 1/10/50/90/99%:", np.round(q, 3))
# full-system preconditioned spectrum (type 1 PB with exact block solves)
N = J.shape[0]


Pbs = P[np.ix_(BS, BS)]
Mi = np.zeros((N, N))
# M^-1 = block LDU inverse with P blocks: build explicitly
Pbs_inv = np.linalg.inv(Pbs)
Sv_inv = np.linalg.inv(Sv)
Jvb = J[np.ix_(V, BS)]
Jbv = J[np.ix_(BS, V)]
# x_bs* = Pbs^-1 r_bs; x_v = Sv^-1 (r_v - Jvb x_bs*); x_bs = Pbs^-1 (r_bs - Jbv x_v)
Mi_vv = Sv_inv
Mi_vb = -Sv_inv @ Jvb @ Pbs_inv
Mi_bv = -Pbs_inv @ Jbv @ Sv_inv
Mi_bb = Pbs_inv + Pbs_inv @ Jbv @ Sv_inv @ Jvb @ Pbs_inv
Mi[np.ix_(BS, BS)] = Mi_bb
Mi[np.ix_(BS, V)] = Mi_bv
Mi[np.ix_(V, BS)] = Mi_vb
Mi[np.ix_(V, V)] = Mi_vv
evf = np.linalg.eigvals(Mi @ J)
print(
    "eig(M_pb^-1 J): |lambda-1| max=%.3e, |lambda| range [%.3e, %.3e], n(|lambda-1|>0.5)=%d of %d"
    % (abs(evf - 1).max(), abs(evf).min(), abs(evf).max(), (abs(evf - 1) > 0.5).sum(), N)
)


# ---- Alternative Schur approximations and dense GMRES iteration counts
def gmres_its(Mi, A, tol=1e-8, maxit=500):
    # right-preconditioned GMRES on A x = b (random b), count its to rel tol (dense, exact arithmetic-ish)
    n = A.shape[0]
    rng = np.random.default_rng(0)
    b = rng.standard_normal(n)
    AM = A @ Mi
    Q = np.zeros((n, maxit + 1))
    H = np.zeros((maxit + 1, maxit))
    beta = np.linalg.norm(b)
    Q[:, 0] = b / beta
    for k in range(maxit):
        w = AM @ Q[:, k]
        for j in range(k + 1):
            H[j, k] = Q[:, j] @ w
            w -= H[j, k] * Q[:, j]
        H[k + 1, k] = np.linalg.norm(w)
        if H[k + 1, k] > 0:
            Q[:, k + 1] = w / H[k + 1, k]
        e1 = np.zeros(k + 2)
        e1[0] = beta
        y, res, *_ = np.linalg.lstsq(H[: k + 2, : k + 1], e1, rcond=None)
        r = np.linalg.norm(e1 - H[: k + 2, : k + 1] @ y) / beta
        if r < tol:
            return k + 1
    return maxit


def pb_inverse(first, second, A11, A22hat):
    # block LDU inverse eliminating `first` (A11 exact) then Schur approx A22hat on `second`
    Mi = np.zeros((N, N))
    A11i = np.linalg.inv(A11)
    Si = np.linalg.inv(A22hat)
    A21 = J[np.ix_(second, first)]
    A12 = J[np.ix_(first, second)]
    Mi[np.ix_(second, second)] = Si
    Mi[np.ix_(second, first)] = -Si @ A21 @ A11i
    Mi[np.ix_(first, second)] = -A11i @ A12 @ Si
    Mi[np.ix_(first, first)] = A11i + A11i @ A12 @ Si @ A21 @ A11i
    return Mi


print("GMRES its (dense, tol 1e-8):")
Mdiag = np.zeros((N, N))
Mdiag[np.ix_(BS, BS)] = np.linalg.inv(J[np.ix_(BS, BS)])
Mdiag[np.ix_(V, V)] = np.linalg.inv(J[np.ix_(V, V)])
print(" blockdiag(J)=", gmres_its(Mdiag, J))
print("  PB V-Schur continuum S_V (code)  =", gmres_its(pb_inverse(BS, V, Jbs, Sv), J))
print("  PB V-Schur exact                 =", gmres_its(pb_inverse(BS, V, Jbs, Sex), J))
Sdiag = J[np.ix_(V, V)] - J[np.ix_(V, BS)] @ np.diag(1 / np.diag(Jbs)) @ J[np.ix_(BS, V)]
print("  PB V-Schur diag(J_BB)            =", gmres_its(pb_inverse(BS, V, Jbs, Sdiag), J))
Jvv = J[np.ix_(V, V)]
SB_diag = Jbs - J[np.ix_(BS, V)] @ np.diag(1 / np.diag(Jvv)) @ J[np.ix_(V, BS)]
print("  PB B-Schur diag(J_VV)            =", gmres_its(pb_inverse(V, BS, Jvv, SB_diag), J))
# lumped (row-sum abs) variant
lump = np.abs(Jvv).sum(axis=1)
SB_lump = Jbs - J[np.ix_(BS, V)] @ np.diag(1 / lump) @ J[np.ix_(V, BS)]
print("  PB B-Schur |rowsum|(J_VV)        =", gmres_its(pb_inverse(V, BS, Jvv, SB_lump), J))
SBex = Jbs - J[np.ix_(BS, V)] @ np.linalg.solve(Jvv, J[np.ix_(V, BS)])
ev = np.linalg.eigvals(np.linalg.solve(SB_diag, SBex))
print(
    "  eig(SB_diag^-1 SB_ex): |.| range [%.3e, %.3e], real [%.3e, %.3e]"
    % (abs(ev).min(), abs(ev).max(), ev.real.min(), ev.real.max())
)
