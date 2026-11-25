#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
# Spheromak-tilt experiment driver (fixed dt, few steps)
import os
import re
import sys
import subprocess
import time
import resource
import argparse
import shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_alfven import make_xml, OFT_BUILD

EX = os.path.join(OFT_BUILD, "examples", "MUG", "Spheromak_tilt")
ap = argparse.ArgumentParser()
ap.add_argument("--tag", required=True)
ap.add_argument("--kind", default="bj")
ap.add_argument("--dt", default="1.e-6")
ap.add_argument("--nsteps", default="3")
ap.add_argument("--nlocal", type=int, default=8)
ap.add_argument("--overlap", action="store_true")
ap.add_argument("--its", type=int, default=8)
ap.add_argument("--inner_its", type=int, default=None)
ap.add_argument("--pbtype", type=int, default=4)
ap.add_argument("--maxits", type=int, default=200)
ap.add_argument("--mfnk", default="T")
ap.add_argument("--inner", default="lu")
ap.add_argument("--timeout", type=int, default=7200)
a = ap.parse_args()
d = os.path.join(os.getcwd(), "runs_sph", a.tag)
os.makedirs(d, exist_ok=True)
shutil.copy(os.path.join(EX, "cyl_tilt.h5"), d)
inp = open(os.path.join(EX, "oft.in")).read()
inp = inp.replace(" nlevels=2", " nlevels=1")
inp = re.sub(r" dt=2.e-7", f" dt={a.dt}", inp)
inp = re.sub(r" nsteps=3000", f" nsteps={a.nsteps}", inp)
inp = re.sub(r" ittarget=40", " ittarget=100000", inp)
inp = re.sub(r" xmhd_mfnk=T", f" xmhd_mfnk={a.mfnk}\n gmres_maxits={a.maxits}", inp)
inp = re.sub(r" rst_freq=10", " rst_freq=1000", inp)
open(os.path.join(d, "oft.in"), "w").write(inp)
open(os.path.join(d, "oft_in.xml"), "w").write(
    make_xml(
        a.kind, nlocal=a.nlocal, overlap=a.overlap, its=a.its, pbtype=a.pbtype, inner_its=a.inner_its, inner=a.inner
    )
)
env = dict(os.environ, OMP_NUM_THREADS="2")
t0 = time.time()
try:
    p = subprocess.run(
        [os.path.join(EX, "MUG_sph_tilt"), "oft.in", "oft_in.xml"],
        cwd=d,
        env=env,
        capture_output=True,
        text=True,
        timeout=a.timeout,
    )
    out = p.stdout + p.stderr
    rc = p.returncode
except subprocess.TimeoutExpired as e:
    out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
    rc = "TIMEOUT"
open(os.path.join(d, "log.txt"), "w").write(out)
steps = re.findall(r"Timestep\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)", out)
rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
print(
    f"{a.tag:34s} rc={rc} steps={len(steps)} lits={[int(s[2]) for s in steps]} nlits={[int(s[3]) for s in steps]} dt={[s[5] for s in steps]} wall={time.time() - t0:.0f}s rss={rss:.0f}MB",
    flush=True,
)
if rc != 0:
    print(out[-1500:])
