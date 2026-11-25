# SPDX-License-Identifier: LGPL-3.0-only
# Offline feasibility study of a NIMROD-style semi-implicit (SI) advance for MUG (xmhd.F90).
#
# Usage: python3 si_analysis.py <dump_dir_dt1> <dump_dir_dt2> <dump_dir_dt3> [--dts 1e-5,2e-5,4e-5]
#   Each dump dir holds xmhd_pb_J.h5 / xmhd_pb_P.h5 / xmhd_pb_fields.txt written by `<pb_pre type="1" dump="T">`
#   (run_alfven.py --kind pb_lu --pbtype 1 --dump --nsteps 2).
#
# Method
#   MUG's Crank-Nicolson Jacobian is J(dt) = M + jdt*L with jdt = dt/2, so two dumps at different dt give the
#   mass matrix M and the linear(ized) operator L (the third dump checks that J is affine in dt). The type-1 P
#   matrix holds S_V = J_VV + jdt^2*W, where W is the continuum parabolization (NIMROD-type SI) operator, so
#   W = (P_VV - J_VV)/jdt^2. With M, L and W in hand, any dt can be analysed without new runs:
#   - CN map:  u+ = (M + jdt L)^-1 (M - jdt L) u
#   - SI map (staggered, leapfrog-like; V at half steps, X = all non-V fields explicit in the V equation):
#       (M_V + jdt L_VV + s jdt^2 W) V+ = (M_V - jdt L_VV + s jdt^2 W) V- - dt L_VX X
#       (M_X + jdt L_XX) X+ = (M_X - jdt L_XX) X - dt L_XV V+
#     s is the SI scale factor (the `alfven` attribute).
#   Reported: affine-fit check, stability (spectral radius) vs dt and s, W vs the exact discrete ideal operator
#   K = L_VX M_X^-1 L_XV mode by mode, and phase/damping of the resolved low-frequency modes for CN vs SI.
import sys
import argparse
import h5py
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

ap = argparse.ArgumentParser()
ap.add_argument("dirs", nargs=3)
ap.add_argument("--dts", default="1e-5,2e-5,4e-5")
ap.add_argument("--scan_dt", default="1e-5,1e-4,1e-3,1e-2")
ap.add_argument("--scan_s", default="0,0.5,1,1.25,1.5,2,4")
ap.add_argument("--nmodes", type=int, default=12)
args = ap.parse_args()
dts = [float(x) for x in args.dts.split(",")]


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
    ng = lg.max()  # fold duplicated (periodic/shared) local rows/cols onto global indices
    G = sp.coo_matrix((A.data, (lg[A.row] - 1, lg[A.col] - 1)), shape=(ng, ng)).toarray()
    return G, lg


Js, Ps = [], []
for d in args.dirs:
    J, lg = load(d + "/xmhd_pb_J.h5", "J")
    P, _ = load(d + "/xmhd_pb_P.h5", "P")
    Js.append(J)
    Ps.append(P)
fl = np.loadtxt(args.dirs[0] + "/xmhd_pb_fields.txt", ndmin=2).astype(int)
grp = np.zeros(Js[0].shape[0], dtype=int)
for i, off, n, g in fl:
    grp[lg[off : off + n] - 1] = g
jd = [dt / 2.0 for dt in dts]

# ---- 1. Extract M, L; check affine dependence on dt
L = (Js[1] - Js[0]) / (jd[1] - jd[0])
M = Js[0] - jd[0] * L
pred = M + jd[2] * L
print("== 1. Affine check: ||J(dt3) - (M + jdt3 L)|| / ||J(dt3)|| = %.3e" % (np.linalg.norm(pred - Js[2]) / np.linalg.norm(Js[2])))
Ws = [(Ps[k] - Js[k]) / jd[k] ** 2 for k in range(3)]
V = np.where(grp == 2)[0]
for k in range(3):
    Ws[k] = Ws[k][np.ix_(V, V)]
print("   W consistency across dumps: ||W2-W1||/||W1|| = %.3e, ||W3-W1||/||W1|| = %.3e"
      % (np.linalg.norm(Ws[1] - Ws[0]) / np.linalg.norm(Ws[0]), np.linalg.norm(Ws[2] - Ws[0]) / np.linalg.norm(Ws[0])))
W_full = Ws[0]

# drop constrained (Dirichlet / identity) rows: L row zero and M row purely diagonal
offd = np.abs(M - np.diag(np.diag(M))).sum(axis=1)
fixed = (np.abs(L).sum(axis=1) == 0) & (offd == 0)
keep = np.where(~fixed)[0]
print("   constrained DOFs dropped: %d of %d" % (fixed.sum(), len(fixed)))
M = M[np.ix_(keep, keep)]
L = L[np.ix_(keep, keep)]
grp = grp[keep]
vkeep = np.isin(V, keep)
W = W_full[np.ix_(vkeep, vkeep)]
V = np.where(grp == 2)[0]
X = np.where(grp != 2)[0]
print("   sizes: V %d, X %d;  ||M_VX||/||M|| = %.2e;  ||W - W^T||/||W|| = %.2e"
      % (len(V), len(X), np.linalg.norm(M[np.ix_(V, X)]) / np.linalg.norm(M), np.linalg.norm(W - W.T) / np.linalg.norm(W)))
# M_VX, M_XV are physically zero; the fitted values are residue of the (weak) state dependence of J between dumps
M[np.ix_(V, X)] = 0.0
M[np.ix_(X, V)] = 0.0
MV, MX = M[np.ix_(V, V)], M[np.ix_(X, X)]
LVV, LVX, LXV, LXX = L[np.ix_(V, V)], L[np.ix_(V, X)], L[np.ix_(X, V)], L[np.ix_(X, X)]

# ---- semi-discrete spectrum: du/dt = -M^-1 L u
mu, E = sla.eig(-L, M)
fin = np.isfinite(mu)
wmax = np.abs(mu[fin].imag).max()
print("== Semi-discrete spectrum: max |Im mu| (fastest wave) = %.3e rad/s; max Re mu = %.2e" % (wmax, mu[fin].real.max()))

# ---- 2. W vs exact discrete ideal operator (V space)
# Eliminating X from the ideal coupling gives M_V V'' = -Kr V with Kr = -L_VX M_X^-1 L_XV (restoring operator).
Kr = -LVX @ np.linalg.solve(MX, LXV)
Ks = 0.5 * (Kr + Kr.T)
Wsym = 0.5 * (W + W.T)
ew = np.linalg.eigvalsh(Wsym)
ek = np.linalg.eigvalsh(Ks)
print("== 2. W (continuum SI operator) vs Kr = -L_VX M_X^-1 L_XV (exact discrete ideal operator)")
print("   ||Kr - Kr^T||/||Kr|| = %.2e;  ||W||/||Kr|| = %.3f" % (np.linalg.norm(Kr - Kr.T) / np.linalg.norm(Kr), np.linalg.norm(W) / np.linalg.norm(Kr)))
print("   eig(W): [%.3e, %.3e];  eig(Kr): [%.3e, %.3e]  (negative parts are round-off relative to max)" % (ew.min(), ew.max(), ek.min(), ek.max()))
# True discrete ideal normal modes: Kr y = w^2 M_V y
w2, Y = sla.eigh(Ks, MV)
wt = np.sqrt(np.maximum(w2, 0))
print("   max discrete ideal frequency = %.3e rad/s (full-system max |Im mu| = %.3e)" % (wt.max(), wmax))
seff = np.einsum("ij,ij->j", Y, Wsym @ Y) / np.maximum(np.einsum("ij,ij->j", Y, Ks @ Y), 1e-300)
print("   Per-mode overestimate s_eff = (y W y)/(y Kr y) for discrete normal modes, binned by w/wmax:")
edges = [1e-3, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.01]
for lo_, hi_ in zip(edges[:-1], edges[1:]):
    sel = (wt / wmax >= lo_) & (wt / wmax < hi_)
    if sel.sum():
        print("     w/wmax in [%.3f,%.2f): n=%4d  s_eff min %.3f  median %.3f  max %.3g" % (lo_, hi_, sel.sum(), seff[sel].min(), np.median(seff[sel]), seff[sel].max()))
nlow = (wt / wmax < 1e-3).sum()
wc_low = np.sqrt(np.maximum(np.einsum("ij,ij->j", Y, Wsym @ Y), 0))[wt / wmax < 1e-3]
print("     w/wmax < 1e-3 (discrete null space, e.g. flows the H(curl) field cannot feel): n=%d; their continuum freq"
      " sqrt(yWy/yMy)/wmax: median %.3f, max %.3f" % (nlow, np.median(wc_low / wmax) if nlow else 0, wc_low.max() / wmax if nlow else 0))


# Ideal-only staggered SI dispersion: (M_V + s jdt^2 W)(V+ - 2V + V-) = -dt^2 Kr V  ->  cos(W dt) = 1 - theta/2,
# theta = eig(dt^2 Kr, M_V + s jdt^2 W); stable iff 0 <= theta <= 4. CN: W dt = 2 atan(w dt / 2).
def si_dispersion(dt, s):
    th, Z = sla.eigh(dt**2 * Ks, MV + s * (dt / 2) ** 2 * Wsym)
    wz = np.sqrt(np.maximum(np.einsum("ij,ij->j", Z, Ks @ Z), 0) / np.einsum("ij,ij->j", Z, MV @ Z))  # true freq of SI mode
    om = np.arccos(np.clip(1 - th / 2, -1, 1)) / dt
    return th, wz, om


# ---- 3. Stability scan
def cn_map(dt):
    jdt = dt / 2
    return np.linalg.solve(M + jdt * L, M - jdt * L)


def si_map(dt, s):
    jdt = dt / 2
    AV = MV + jdt * LVV + s * jdt**2 * W
    BV = MV - jdt * LVV + s * jdt**2 * W
    AX = MX + jdt * LXX
    BX = MX - jdt * LXX
    GVV = np.linalg.solve(AV, BV)
    GVX = -dt * np.linalg.solve(AV, LVX)
    GXV = -dt * np.linalg.solve(AX, LXV @ GVV)
    GXX = np.linalg.solve(AX, BX - dt * LXV @ GVX)
    n = M.shape[0]
    G = np.zeros((n, n))
    G[np.ix_(V, V)] = GVV
    G[np.ix_(V, X)] = GVX
    G[np.ix_(X, V)] = GXV
    G[np.ix_(X, X)] = GXX
    return G


scan_dt = [float(x) for x in args.scan_dt.split(",")]
scan_s = [float(x) for x in args.scan_s.split(",")]
print("== 3. Spectral radius of one-step map (stable if <= 1)")
print("   %-9s %-7s %-10s" % ("dt", "CFL", "CN") + "".join("SI s=%-6g" % s for s in scan_s))
spec = {}
for dt in scan_dt:
    row = "   %-9.0e %-7.1f %-10.6f" % (dt, wmax * dt, np.abs(np.linalg.eigvals(cn_map(dt))).max())
    for s in scan_s:
        ev = np.linalg.eigvals(si_map(dt, s))
        spec[(dt, s)] = ev
        row += "%-11.6f" % np.abs(ev).max()
    print(row)

# ---- 4. Accuracy: ideal-only dispersion of SI vs CN, binned by the true discrete frequency
print("== 4. Ideal-mode dispersion: median numerical/true frequency ratio, binned by w*dt (theta_max: SI stable iff <= 4)")
for dt in scan_dt:
    for s_ in [1.0, 1.25, 2.0]:
        th, wz, om = si_dispersion(dt, s_)
        line = "   dt=%.0e s=%-5g theta_max=%.3f |" % (dt, s_, th.max())
        for lo_, hi_ in [(0.01, 0.3), (0.3, 1), (1, 3), (3, 10), (10, 30), (30, 100), (100, 1e9)]:
            sel = (wz * dt >= lo_) & (wz * dt < hi_)
            if sel.sum():
                x = np.median(wz[sel] * dt)
                line += " wdt[%g,%g) n=%d SI %.3f CN %.3f |" % (lo_, hi_, sel.sum(), np.median(om[sel] / wz[sel]), 2 * np.arctan(x / 2) / x)
        print(line)
sys.stdout.flush()

# ---- 5. Accuracy for well-represented physical modes (s_eff close to 1): Rayleigh estimate of the SI frequency
# For a discrete normal mode y (Kr y = w^2 M_V y): theta_y = dt^2 yKr y / y(M_V + s jdt^2 W)y, cos(W_SI dt) = 1 - theta_y/2.
good = np.where((wt > 1e-3 * wmax) & (seff < 1.5))[0]
good = good[np.argsort(wt[good])]
print("== 5. Physical modes with s_eff < 1.5: n=%d; lowest w = %.4e (pi*v_A = %.4e for the test's v_A=1e4)" % (len(good), wt[good[0]], np.pi * 1e4))
sel = good[:: max(1, len(good) // 8)][:8]
for dt in scan_dt:
    for s_ in [1.0, 1.25, 2.0]:
        line = "   dt=%.0e s=%-5g" % (dt, s_)
        for j in sel:
            y = Y[:, j]
            x = wt[j] * dt
            th = dt**2 * (y @ Ks @ y) / (y @ (MV + s_ * (dt / 2) ** 2 * Wsym) @ y)
            om = np.arccos(np.clip(1 - th / 2, -1, 1)) / dt
            line += " | w=%.2e wdt=%.2f s_eff=%.2f SI/w=%.3f CN/w=%.3f" % (wt[j], x, seff[j], om / wt[j], 2 * np.arctan(x / 2) / x)
        print(line)
sys.stdout.flush()
