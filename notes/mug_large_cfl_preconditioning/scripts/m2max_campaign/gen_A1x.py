import subprocess, os
from gen_A import ndof, mem_est, PY, BENCH, OMEGA0
os.makedirs("A1x", exist_ok=True)
out = []
for R in (0, 1):
    for p in (2, 3, 4):
        for cfl in (3, 30, 300):
            dt = cfl / (OMEGA0 * (p / 2.0) * 2**R)
            n = ndof(R, p)
            for tag, extra in (("lu", ["--pre", "pb1_lu"]),
                               ("bj3k", ["--pre", "pb1", "--nlocal", str(max(1, int(round(n / 3000.0))))]),
                               ("i20", ["--pre", "pb1", "--nlocal", str(max(1, int(round(n / 300.0)))), "--inner_its", "20"])):
                d = "A1x/cube_R%d_p%d_cfl%g_s1_lin_%s" % (R, p, cfl, tag)
                subprocess.run([PY, BENCH, "--geom", "cube", "--out", d, "--nbase", str(R + 1), "--nlevels", str(R + 1),
                                "--order", str(p), "--dt", "%.4E" % dt, "--nsteps", "6", "--si", "1", "--linear",
                                "--threads", "2"] + extra, check=True, stdout=subprocess.DEVNULL)
                out.append("%s::%.1f" % (d, mem_est(R, p) * (3 if tag == "lu" else 1)))
open("A1x.list", "w").write("\n".join(out) + "\n")
print(len(out))
