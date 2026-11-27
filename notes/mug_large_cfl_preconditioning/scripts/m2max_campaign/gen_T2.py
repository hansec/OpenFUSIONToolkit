#!/usr/bin/env python3
"""Follow-up timing cases: full-Jacobian block-Jacobi LU with the SI term (missing from T), same 4x2 layout."""
import os, re, shutil
from gen_T import make, REF
lines = []
phases = {"lin": 500, "peak": 1100, "late": 2000}
cfg = [("bj64", 1, "1.e-6", 64), ("bj256", 1, "1.e-6", 256), ("bj256", 0, "1.e-6", 256), ("bj64", 0.5, "3.e-6", 64)]
for ph, ind in phases.items():
    for tag, s, dt, nl in cfg:
        out = "T/D/%s_%s_s%g_dt%s" % (ph, tag, s, dt.replace(".e", "e"))
        make(out, ["--dt", dt, "--nsteps", "20", "--si", str(s), "--pre", "bj", "--rst_freq", "1000"], nlocal_total=nl)
        p = os.path.join(out, "oft.in")
        inp, n = re.subn(r"^ rst_ind=\d+", " rst_ind=%d" % ind, open(p).read(), flags=re.M)
        assert n == 1
        open(p, "w").write(inp)
        shutil.copy(os.path.join(REF, "xMHD_%05d.rst" % ind), out)
        lines.append("%s:8:%.1f" % (out, 4.0 + 16.0 / nl ** 0.5))
for dt in ("1e-6", "3e-6", "1e-5", "3e-5"):
    lines.append(make("T/B/R0_p3_bj256_s1_dt%s" % dt, ["--order", "3", "--dt", dt, "--nsteps", "8", "--si", "1", "--pre", "bj"],
                      nlocal_total=256) + ":8:18")
open("T2.list", "w").write("\n".join(lines) + "\n")
print(len(lines))
