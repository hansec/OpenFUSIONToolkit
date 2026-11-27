#!/usr/bin/env python3
"""Timing lane: every case runs alone on the machine with 4 MPI ranks x 2 OpenMP threads (8 P-cores).

D: production-style preconditioner cost, 20-step restart segments from the s=0, dt=2e-7 reference (C2) at three
   phases (one factorization per segment, xmhd_prefreq=20). Block counts are totals over the 4 ranks.
B: spheromak tilt from t=0, 8 steps, R0 p=2/3 and R1 p=2 (weak scaling R0 -> R1 at fixed ranks), plus an R1
   strong-scaling set at fixed 8 cores (2x4, 4x2, 8x1).
"""
import os
import sys
import re
import shutil
import subprocess

ROOT = os.environ.get("OFT_ROOT", os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../..")))
BENCH = os.path.join(ROOT, "notes/mug_large_cfl_preconditioning/scripts/bench_case.py")
PY = os.environ.get("OFT_PYTHON", sys.executable)
G = os.environ.get("OFT_TREE", os.path.join(ROOT, "builds/GCC/MPI_PETSc"))
MPI = G + "/mpich-4_2_3/bin/mpirun -np {ranks}"
REF = "C2/sph_R0_p2_ilu_s0_dt2e-7"
NR, NT = 4, 2


def make(out, args, ranks=NR, threads=NT, nlevels=2, nlocal_total=None, asm=False):
    a = [PY, BENCH, "--geom", "sph", "--out", out, "--nbase", "1", "--nlevels", str(nlevels), "--ranks", str(ranks),
         "--threads", str(threads), "--mpirun", MPI] + args
    if nlocal_total is not None:
        a += ["--nlocal", str(max(1, nlocal_total // ranks))]
    if asm:
        a.append("--overlap")
    subprocess.run(a, check=True, stdout=subprocess.DEVNULL)
    if asm:  # additive Schwarz: overlap inside each rank and across rank boundaries
        p = os.path.join(out, "oft_in.xml")
        x = open(p).read().replace("<loverlap>T</loverlap>", "<loverlap>T</loverlap><boverlap>T</boverlap>")
        open(p, "w").write(x)
    return out


def main():
    lines = []
    # ---- D
    phases = {"lin": 500, "peak": 1100, "late": 2000}
    cfg = [("ilu", 0, "2.e-7", "ilu", None, False)]
    cfg += [("bj%d" % n, 0, "2.e-7", "bj", n, False) for n in (256, 64, 16, 4)]
    cfg += [("asm16", 0, "2.e-7", "bj", 16, True)]
    cfg += [("bj%d" % n, 0, "1.e-6", "bj", n, False) for n in (64, 16, 4)]
    cfg += [("asm16", 0, "1.e-6", "bj", 16, True)]
    cfg += [("pb%d" % n, 1, "1.e-6", "pb1", n, False) for n in (1024, 256, 64, 16)]
    cfg += [("bj16", 1, "1.e-6", "bj", 16, False)]
    cfg += [("pb%d" % n, 0.5, "3.e-6", "pb1", n, False) for n in (256, 64, 16)]
    cfg += [("bj16", 0.5, "3.e-6", "bj", 16, False)]
    cfg += [("pb%d" % n, 1, "3.e-6", "pb1", n, False) for n in (64, 16)]
    os.makedirs("T/D", exist_ok=True)
    for ph, ind in phases.items():
        for tag, s, dt, pre, nl, ov in cfg:
            out = "T/D/%s_%s_s%g_dt%s" % (ph, tag, s, dt.replace(".e", "e"))
            make(out, ["--dt", dt, "--nsteps", "20", "--si", str(s), "--pre", pre, "--rst_freq", "1000"],
                 nlocal_total=nl, asm=ov)
            p = os.path.join(out, "oft.in")
            inp, n = re.subn(r"^ rst_ind=\d+", " rst_ind=%d" % ind, open(p).read(), flags=re.M)
            assert n == 1
            open(p, "w").write(inp)
            shutil.copy(os.path.join(REF, "xMHD_%05d.rst" % ind), out)
            lines.append("%s:8:%.1f" % (out, 4.0 if nl is None else 4.0 + 16.0 / nl ** 0.5))
    # ---- B (from t=0)
    os.makedirs("T/B", exist_ok=True)
    for dt in ("1e-6", "3e-6", "1e-5", "3e-5"):
        for s in (0, 1, 2):
            lines.append(make("T/B/R0_p2_pb1024_s%g_dt%s" % (s, dt), ["--dt", dt, "--nsteps", "8", "--si", str(s), "--pre", "pb1"],
                              nlocal_total=1024) + ":8:5")
        lines.append(make("T/B/R0_p2_pb1024i20_s1_dt%s" % dt, ["--dt", dt, "--nsteps", "8", "--si", "1", "--pre", "pb1",
                                                              "--inner_its", "20"], nlocal_total=1024) + ":8:5")
        lines.append(make("T/B/R0_p2_pb64_s1_dt%s" % dt, ["--dt", dt, "--nsteps", "8", "--si", "1", "--pre", "pb1"],
                          nlocal_total=64) + ":8:8")
        lines.append(make("T/B/R0_p2_ilu_s0_dt%s" % dt, ["--dt", dt, "--nsteps", "8", "--si", "0", "--pre", "ilu"]) + ":8:4")
        for s in (0, 1):
            lines.append(make("T/B/R0_p2_bj64_s%g_dt%s" % (s, dt), ["--dt", dt, "--nsteps", "8", "--si", str(s), "--pre", "bj"],
                              nlocal_total=64) + ":8:8")
        for s in (0, 1):
            lines.append(make("T/B/R0_p3_pb2800_s%g_dt%s" % (s, dt), ["--order", "3", "--dt", dt, "--nsteps", "8", "--si", str(s),
                                                                    "--pre", "pb1"], nlocal_total=2800) + ":8:14")
        lines.append(make("T/B/R0_p3_bj256_s0_dt%s" % dt, ["--order", "3", "--dt", dt, "--nsteps", "8", "--si", "0", "--pre", "bj"],
                          nlocal_total=256) + ":8:18")
    # weak scaling R0 -> R1 at the same layout; R1 has 8x the DOFs and blocks
    for dt in ("1e-6", "1e-5"):
        for s in (0, 1):
            lines.append(make("T/B/R1_p2_pb8192_s%g_dt%s" % (s, dt), ["--dt", dt, "--nsteps", "8", "--si", str(s), "--pre", "pb1"],
                              nlevels=3, nlocal_total=8192) + ":8:30")
    # strong scaling at fixed 8 cores
    for r, t in ((2, 4), (8, 1)):
        lines.append(make("T/B/R1_p2_pb8192_s1_dt1e-5_np%dt%d" % (r, t), ["--dt", "1e-5", "--nsteps", "8", "--si", "1",
                                                                        "--pre", "pb1"], ranks=r, threads=t, nlevels=3,
                          nlocal_total=8192) + ":8:30")
    open("T.list", "w").write("\n".join(lines) + "\n")
    print(len(lines))


if __name__ == "__main__":
    main()
