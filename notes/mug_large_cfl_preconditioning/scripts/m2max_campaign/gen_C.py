#!/usr/bin/env python3
"""Generate the spheromak accuracy cases (C1 linear growth rate, C2 nonlinear tilt to 0.6 ms), R0 p=2, serial, 2 threads.

Equivalent to the bench_case.py calls used for the M2 Max campaign (parameters recovered from each case.json).
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "..", "bench_case.py")
PY = os.environ.get("OFT_PYTHON", sys.executable)


def make(out, args):
    subprocess.run([PY, BENCH, "--geom", "sph", "--out", out, "--threads", "2"] + args, check=True,
                   stdout=subprocess.DEVNULL)
    return out


def main():
    c1, c2 = [], []
    # C1: linear (BDF2) tilt to 300 us; PB type 1 + GMRES(8)/BJ-LU(1024) on P for every s (s=0 included)
    runs = {"1.e-7": (3000, 50), "2e-7": (1500, 50), "1e-6": (300, 10), "3e-6": (100, 5), "1e-5": (30, 2)}
    for dt, (n, rf) in runs.items():
        for s in ((0,) if dt == "1.e-7" else (0, 0.5, 1, 2)):
            c1.append(make("C1/sph_lin_s%g_dt%s" % (s, dt.replace(".e", "e")),
                           ["--linear", "--dt", dt, "--nsteps", str(n), "--rst_freq", str(rf), "--si", str(s),
                            "--pre", "pb1", "--nlocal", "1024"]))
    # C2: nonlinear tilt to 0.6 ms (matrix-free Newton, as in the example)
    c2.append(make("C2/sph_R0_p2_ilu_s0_dt2e-7", ["--dt", "2.e-7", "--nsteps", "3000", "--rst_freq", "50", "--si", "0",
                                                  "--pre", "ilu"]))
    c2.append(make("C2/sph_R0_p2_bj16_s0_dt1e-6", ["--dt", "1.e-6", "--nsteps", "600", "--rst_freq", "50", "--si", "0",
                                                   "--pre", "bj", "--nlocal", "16"]))
    for s, dt, n, extra in ((1, "1.e-6", 600, []), (1, "3.e-6", 200, []), (1, "1.e-5", 60, []), (0.5, "3.e-6", 200, []),
                            (1, "3.e-6", 200, ["--si_full"])):
        tag = "s%g%s" % (s, "full" if extra else "")
        c2.append(make("C2/sph_R0_p2_pb1_%s_dt%s" % (tag, dt.replace(".e", "e")),
                       ["--dt", dt, "--nsteps", str(n), "--rst_freq", "50" if n > 60 else "10", "--si", str(s),
                        "--pre", "pb1", "--nlocal", "1024"] + extra))
    open("C1.list", "w").write("\n".join("%s::3" % d for d in c1) + "\n")
    open("C2.list", "w").write("\n".join("%s::3" % d for d in c2) + "\n")
    print(len(c1), len(c2))


if __name__ == "__main__":
    main()
