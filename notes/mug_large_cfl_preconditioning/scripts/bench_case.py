#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
"""Generate a self-contained MUG benchmark case directory (cube Alfven wave or spheromak tilt).

The script only writes files: oft.in, oft_in.xml, case.json and run.sh (plus the spheromak mesh), so the
case can be submitted with any scheduler. See HANDOFF_BENCHMARKS.md for the benchmark matrix.

Mesh levels follow OFT conventions: `nlevels` counts all grid levels including the MPI transfer level, and
`nbase` levels are built before domain decomposition. Serial runs need nbase == nlevels; MPI runs need
nlevels > nbase. The mesh resolution is therefore R = nlevels-1 refinements (serial) or nlevels-2 (MPI),
each refinement giving 8x the cells.

Preconditioner presets (--pre):
  ilu     GMRES(its) + block-Jacobi(nlocal, default -1) + ILU on J (the spheromak example's default)
  bj      GMRES(its) + block-Jacobi(nlocal) + inner (lu/ilu) on J, optional --overlap
  lu      direct LU of the full J (small cases only)
  mg      polynomial MG on J with a GMRES/block-Jacobi smoother (cube only, needs minlev < FE levels)
  pb1     <pb_pre type=1>: GMRES(inner_its) + block-Jacobi(nlocal) + inner on P
  pb1_lu  <pb_pre type=1> with a direct LU of P
  pb1_mg  <pb_pre type=1> with polynomial MG on P (cube only)

Examples:
  python3 bench_case.py --geom sph --out sph_R0_p2_s1_dt1e-5 --dt 1e-5 --nsteps 10 --si 1 --pre pb1 --nlocal 1024
  python3 bench_case.py --geom sph --out sph_R1_p2 --nbase 1 --nlevels 3 --ranks 16 --threads 2 --dt 3e-6 --si 1 --pre pb1
  python3 bench_case.py --geom cube --out cube_r1_p3 --order 3 --dt 1e-4 --nsteps 10 --si 1 --pre pb1_lu
"""
import argparse
import json
import os
import re
import shutil
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from run_alfven import make_xml, OFT_ROOT, OFT_BUILD  # noqa: E402

sys.path.insert(0, os.path.join(OFT_ROOT, "src", "tests", "physics", "MUG"))
from mug_testing import ALFVEN_TEMPLATE  # noqa: E402

SPH_SRC = os.path.join(OFT_ROOT, "src", "examples", "MUG", "Spheromak_tilt")


def build_xml(a):
    common = dict(nlocal=a.nlocal, overlap=a.overlap, its=a.its, pbtype=1, alfven=a.pb_alfven, sound=a.pb_sound,
                  kcycle=a.kcycle, inner_its=a.inner_its, inner=a.inner)
    if a.pre == "ilu":
        return make_xml("bj", nlocal=a.nlocal if a.nlocal_set else -1, overlap=False, its=a.its, inner="ilu")
    kind = {"bj": "bj", "lu": "lu", "mg": "mg", "pb1": "pb_bj", "pb1_lu": "pb_lu", "pb1_mg": "pb_mg"}[a.pre]
    return make_xml(kind, **common)


def si_lines(a):
    out = ""
    if a.si > 0.0:
        out += " xmhd_si=%g\n" % a.si
    if a.si_sound > 0.0:
        out += " xmhd_si_sound=%g\n" % a.si_sound
    if a.si_full:
        out += " xmhd_si_full=T\n"
    return out


def cube_input(a):
    rep = " bbc='bc'\n nu_visc=1.E-3\n gmres_maxits=%d\n" % a.maxits
    if a.linear:
        rep += " xmhd_lin_euler_start=T\n"
        if a.lin_be:
            rep += " xmhd_lin_bdf2=F\n"
    rep += si_lines(a)
    inp = ALFVEN_TEMPLATE.format(nbase=a.nbase, nlevels=a.nlevels, order=a.order, minlev=a.minlev, nu_mhd="-1",
                                 linear="T" if a.linear else "F", dt=a.dt, nsteps=a.nsteps, lin_tol=a.lin_tol,
                                 mfnk="T" if a.mfnk else "F", petsc="F", mesh_type=2 if a.hex else 1,
                                 rep_options=rep)
    inp = inp.replace("ittarget=100", "ittarget=%d" % a.ittarget)
    return inp, os.path.join(OFT_BUILD, "tests", "physics", "MUG", "test_alfven"), []


def sph_input(a):
    inp = open(os.path.join(SPH_SRC, "oft.in")).read()

    def sub(pattern, repl):
        nonlocal inp
        new, n = re.subn(pattern, repl, inp, count=1, flags=re.M)
        if n != 1:
            raise RuntimeError("Could not update spheromak input: " + pattern)
        inp = new

    sub(r"^ nlevels=\d+", " nlevels=%d" % a.nlevels)
    sub(r"^ nbase=\d+", " nbase=%d" % a.nbase)
    sub(r"^ order=\d+", " order=%d" % a.order)
    sub(r"^ linear=[TF]", " linear=%s" % ("T" if a.linear else "F"))
    sub(r"^ dt=2\.e-7", " dt=%s" % a.dt)
    sub(r"^ nsteps=\d+", " nsteps=%d" % a.nsteps)
    sub(r"^ rst_freq=\d+", " rst_freq=%d" % a.rst_freq)
    sub(r"^ lin_tol=\S+", " lin_tol=%s" % a.lin_tol)
    sub(r"^ ittarget=\d+", " ittarget=%d" % a.ittarget)
    extra = " gmres_maxits=%d\n%s" % (a.maxits, si_lines(a))
    sub(r"^ xmhd_mfnk=[TF].*$", " xmhd_mfnk=%s       ! Use matrix-free Jacobian operator\n%s" % ("T" if a.mfnk else "F",
                                                                                         extra.rstrip("\n")))
    sub(r"^ rst_end=\d+", " rst_end=%d" % a.nsteps)
    inp = re.sub(r"\n\n+", "\n\n", inp)
    return inp, os.path.join(OFT_BUILD, "examples", "MUG", "Spheromak_tilt", "MUG_sph_tilt"), ["cyl_tilt.h5"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geom", choices=("cube", "sph"), required=True)
    ap.add_argument("--out", required=True, help="case directory to create")
    ap.add_argument("--nbase", type=int, default=1)
    ap.add_argument("--nlevels", type=int, default=1)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--minlev", type=int, default=None, help="cube only: lowest FE level for MG (default: top level)")
    ap.add_argument("--hex", action="store_true", help="cube only: hexahedral mesh")
    ap.add_argument("--dt", default="1.e-6")
    ap.add_argument("--nsteps", type=int, default=10)
    ap.add_argument("--rst_freq", type=int, default=1000, help="spheromak restart-file frequency")
    ap.add_argument("--linear", action="store_true", help="linearized run (BDF2 unless --lin_be)")
    ap.add_argument("--lin_be", action="store_true")
    ap.add_argument("--mfnk", choices=("T", "F"), default=None, help="matrix-free Newton (default T sph, F cube)")
    ap.add_argument("--si", type=float, default=0.0, help="xmhd_si (semi-implicit scale s)")
    ap.add_argument("--si_sound", type=float, default=0.0)
    ap.add_argument("--si_full", action="store_true")
    ap.add_argument("--pre", default="pb1", choices=("ilu", "bj", "lu", "mg", "pb1", "pb1_lu", "pb1_mg"))
    ap.add_argument("--nlocal", type=int, default=None, help="block-Jacobi blocks per rank")
    ap.add_argument("--overlap", action="store_true")
    ap.add_argument("--its", type=int, default=8, help="GMRES iterations of the J-level preconditioner/smoother")
    ap.add_argument("--inner_its", type=int, default=8, help="GMRES iterations of the solve on P (pb1)")
    ap.add_argument("--inner", default="lu", choices=("lu", "ilu"))
    ap.add_argument("--kcycle", type=int, default=0)
    ap.add_argument("--pb_alfven", type=float, default=1.0)
    ap.add_argument("--pb_sound", type=float, default=1.0)
    ap.add_argument("--maxits", type=int, default=400, help="outer GMRES iteration limit")
    ap.add_argument("--lin_tol", default="1.E-9")
    ap.add_argument("--ittarget", type=int, default=100000, help="large value keeps dt fixed unless Newton fails")
    ap.add_argument("--ranks", type=int, default=1)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--mpirun", default="mpirun -np {ranks}", help="launcher template for run.sh")
    ap.add_argument("--exe", default=None, help="override the executable path")
    a = ap.parse_args()

    a.nlocal_set = a.nlocal is not None
    if a.nlocal is None:
        a.nlocal = 8
    if a.mfnk is None:
        a.mfnk = "T" if a.geom == "sph" else "F"
    a.mfnk = a.mfnk == "T"
    if a.minlev is None:
        a.minlev = a.nlevels + a.order - 1
    if a.ranks > 1 and a.nlevels <= a.nbase:
        ap.error("MPI runs need nlevels > nbase (the extra level is the distributed transfer level)")
    if a.ranks == 1 and a.nlevels != a.nbase:
        ap.error("serial runs need nbase == nlevels")
    if a.geom == "cube" and a.lin_tol == "1.E-9":
        a.lin_tol = "1.E-12"
    if a.geom == "sph" and a.pre in ("mg", "pb1_mg"):
        ap.error("the spheromak example does not expose minlev, so its preconditioner is single-level")

    os.makedirs(a.out, exist_ok=True)
    inp, exe, extra = (cube_input if a.geom == "cube" else sph_input)(a)
    if a.exe:
        exe = a.exe
    open(os.path.join(a.out, "oft.in"), "w").write(inp)
    open(os.path.join(a.out, "oft_in.xml"), "w").write(build_xml(a) + "\n")
    for f in extra:
        shutil.copy(os.path.join(SPH_SRC, f), a.out)
    refinements = a.nlevels - 1 - (1 if a.ranks > 1 else 0)
    meta = dict(vars(a), exe=exe, refinements=refinements)
    json.dump(meta, open(os.path.join(a.out, "case.json"), "w"), indent=1, sort_keys=True)
    launch = exe + " oft.in oft_in.xml"
    if a.ranks > 1:
        launch = a.mpirun.format(ranks=a.ranks) + " " + launch
    with open(os.path.join(a.out, "run.sh"), "w") as fid:
        fid.write("#!/bin/bash\n# Generated by bench_case.py; run from this directory\n")
        fid.write("cd \"$(dirname \"$0\")\"\nexport OMP_NUM_THREADS=%d\n" % a.threads)
        fid.write("if [ -x /usr/bin/time ]; then TIMECMD=\"/usr/bin/time -v -o time.log\"; else TIMECMD=\"\"; fi\n")
        fid.write("$TIMECMD %s > log.txt 2>&1\n" % launch)
    os.chmod(os.path.join(a.out, "run.sh"), 0o755)
    print("wrote %s (R%d, p=%d, ranks=%d, pre=%s, s=%g)" % (a.out, refinements, a.order, a.ranks, a.pre, a.si))


if __name__ == "__main__":
    main()
