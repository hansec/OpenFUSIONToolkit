#!/usr/bin/env python3
"""Generate Benchmark A (cube Alfven wave) case directories and runq list files."""
import os
import subprocess
import sys

ROOT = os.environ.get("OFT_ROOT", os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../..")))
BENCH = os.path.join(ROOT, "notes/mug_large_cfl_preconditioning/scripts/bench_case.py")
PY = os.environ.get("OFT_PYTHON", sys.executable)
OMEGA0 = 1.9e5  # fastest discrete Alfven frequency at R0, p=2
# DOF multipliers relative to p=2 (B, V+N+T)
NDOF0 = {"B": 1856, "VNT": 768 + 512}
PMUL = {2: (1.0, 1.0), 3: (2.6, 3.4), 4: (5.3, 8.0)}


def ndof(R, p):
    return 8**R * (NDOF0["B"] * PMUL[p][0] + NDOF0["VNT"] * PMUL[p][1])


def mem_est(R, p):
    return max(1.0, 0.5 + 2.5e-5 * ndof(R, p) * (p / 2.0))


def gen(out, R, p, dt, s, mode, nsteps, extra=()):
    nlocal = max(1, int(round(ndof(R, p) / 300.0)))
    args = [PY, BENCH, "--geom", "cube", "--out", out, "--nbase", str(R + 1), "--nlevels", str(R + 1),
            "--order", str(p), "--dt", "%.4E" % dt, "--nsteps", str(nsteps), "--si", str(s), "--pre", "pb1",
            "--nlocal", str(nlocal), "--threads", "2"]
    if mode == "lin":
        args.append("--linear")
    elif mode == "cnmf":
        args += ["--mfnk", "T"]
    args += list(extra)
    subprocess.run(args, check=True, stdout=subprocess.DEVNULL)
    return "%s::%.1f" % (out, mem_est(R, p))


def main():
    os.makedirs("A1", exist_ok=True)
    os.makedirs("A2", exist_ok=True)
    a1 = []
    for R in (0, 1, 2):
        for p in (2, 3, 4):
            for cfl in (0.3, 3, 30, 300):
                dt = cfl / (OMEGA0 * (p / 2.0) * 2**R)
                for s in (0, 0.5, 1, 2):
                    for mode in ("lin", "cn", "cnmf"):
                        if R == 2 and (s == 0.5 or mode == "cn"):
                            continue
                        out = "A1/cube_R%d_p%d_cfl%g_s%g_%s" % (R, p, cfl, s, mode)
                        a1.append((R, gen(out, R, p, dt, s, mode, 6)))
    # Run small cases first
    a1.sort(key=lambda x: x[0])
    open("A1.list", "w").write("\n".join(x[1] for x in a1) + "\n")
    a2 = []
    T = 1.0e-4
    for R in (0, 1):
        for p in (2, 3, 4):
            for N in (500, 125, 25, 10, 5):
                for s in (0, 0.5, 1, 2):
                    if N == 500 and s != 0:
                        continue
                    out = "A2/cube_R%d_p%d_N%d_s%g_cn" % (R, p, N, s)
                    a2.append((R, gen(out, R, p, T / N, s, "cn", N)))
    a2.sort(key=lambda x: x[0])
    open("A2.list", "w").write("\n".join(x[1] for x in a2) + "\n")
    print(len(a1), len(a2))


if __name__ == "__main__":
    sys.exit(main())
