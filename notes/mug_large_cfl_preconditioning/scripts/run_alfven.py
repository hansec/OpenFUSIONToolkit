#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
# Experiment driver: Alfven-wave cube (test_alfven) with configurable dt / order / levels / preconditioner XML
import sys
import os
import re
import subprocess
import time
import argparse
import resource

OFT_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
OFT_BUILD = os.environ.get("OFT_BUILD_DIR", os.path.join(OFT_ROOT, "builds", "build_release"))
sys.path.insert(0, os.path.join(OFT_ROOT, "src", "tests", "physics", "MUG"))
from mug_testing import ALFVEN_TEMPLATE

BUILD = os.path.join(OFT_BUILD, "tests", "physics", "MUG")
XML = {}


def bj(nlocal, overlap, its=2, inner="lu"):
    ov = "<loverlap>T</loverlap>" if overlap else ""
    return f"""<solver type="gmres"><its>{its}</its><nrits>{its}</nrits>
<pre type="block_jacobi"><nlocal>{nlocal}</nlocal>{ov}<solver type="{inner}"></solver></pre></solver>"""


def lu():
    return '<solver type="gmres"><its>1</its><nrits>1</nrits><pre type="block_jacobi"><nlocal>1</nlocal><solver type="lu"></solver></pre></solver>'


def make_xml(
    kind,
    nlocal=8,
    overlap=True,
    its=2,
    pbtype=1,
    alfven=1.0,
    sound=1.0,
    kcycle=0,
    cscale=False,
    mg=False,
    inner_its=None,
    inner="lu",
    dump=False,
):
    attrs = ""
    if kcycle:
        attrs += f' kcycle="{kcycle}"'
    if cscale:
        attrs += ' cscale="T"'
    if kind == "lu":  # direct solve of full J
        body = '<pre type="gmres"><its>1</its><nrits>1</nrits><pre type="block_jacobi"><nlocal>1</nlocal><solver type="lu"></solver></pre></pre>'
        return f"<oft><xmhd>{body}</xmhd></oft>"
    if kind == "bj":  # one-level smoother on full J
        ov = "<loverlap>T</loverlap>" if overlap else ""
        return f'<oft><xmhd><pre type="gmres"><its>{its}</its><nrits>{its}</nrits><pre type="block_jacobi"><nlocal>{nlocal}</nlocal>{ov}<solver type="{inner}"></solver></pre></pre></xmhd></oft>'
    if kind == "mg":  # p-MG on full J
        return f'<oft><xmhd><pre type="mg"{attrs}><smoother direction="up">{bj(nlocal, overlap, its, inner)}</smoother></pre></xmhd></oft>'
    if kind.startswith("pb"):  # physics-based; inner = 'lu' | 'bj' | 'mg'
        pinner = kind.split("_")[1]
        a = f'<pb_pre type="{pbtype}" alfven="{alfven}" sound="{sound}"'
        if dump:  # write J, P and field layout to HDF5 at the first preconditioner update (needs nsteps>=2)
            a += ' dump="T"'
        if pinner == "lu":
            return f"<oft><xmhd>{a}>{lu()}</pb_pre></xmhd></oft>"
        if pinner == "bj":
            return f"<oft><xmhd>{a}>{bj(nlocal, overlap, inner_its or its, inner)}</pb_pre></xmhd></oft>"
        if pinner == "mg":
            return f'<oft><xmhd>{a}{attrs}><smoother direction="up">{bj(nlocal, overlap, its, inner)}</smoother></pb_pre></xmhd></oft>'
    raise ValueError(kind)


def run(args):
    tag = args.tag
    d = os.path.join(os.getcwd(), "runs", tag)
    os.makedirs(d, exist_ok=True)
    rep = " bbc='bc'\n nu_visc=1.E-3\n gmres_maxits=%d\n" % args.maxits
    if args.linear:
        rep += " xmhd_lin_euler_start=T\n"
    inp = ALFVEN_TEMPLATE.format(
        nbase=args.nbase,
        nlevels=args.nlevels,
        order=args.order,
        minlev=args.minlev,
        nu_mhd="-1",
        linear="T" if args.linear else "F",
        dt=args.dt,
        nsteps=args.nsteps,
        lin_tol=args.lin_tol,
        mfnk="F",
        petsc="F",
        mesh_type=1,
        rep_options=rep,
    ).replace("ittarget=100", "ittarget=100000")
    open(os.path.join(d, "oft.in"), "w").write(inp)
    xml = make_xml(
        args.kind,
        nlocal=args.nlocal,
        overlap=args.overlap,
        its=args.its,
        pbtype=args.pbtype,
        alfven=args.alfven,
        sound=args.sound,
        kcycle=args.kcycle,
        cscale=args.cscale,
        inner_its=args.inner_its,
        inner=args.inner,
        dump=args.dump,
    )
    open(os.path.join(d, "oft_in.xml"), "w").write(xml)
    env = dict(os.environ, OMP_NUM_THREADS="2", OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    t0 = time.time()
    cmd = [os.path.join(BUILD, "test_alfven"), "oft.in", "oft_in.xml"]
    if args.np > 1:
        cmd = ["mpirun", "-np", str(args.np), os.path.join(BUILD, "test_alfven"), "oft.in", "oft_in.xml"]
    try:
        p = subprocess.run(cmd, cwd=d, env=env, capture_output=True, text=True, timeout=args.timeout)
        out = p.stdout + p.stderr
        rc = p.returncode
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        rc = "TIMEOUT"
    wall = time.time() - t0
    open(os.path.join(d, "log.txt"), "w").write(out)
    steps = re.findall(r"Timestep\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)", out)
    lits = [int(s[2]) for s in steps]
    nlits = [int(s[3]) for s in steps]
    dts = [float(s[5]) for s in steps]
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
    err = ""
    if os.path.exists(os.path.join(d, "alfven.results")):
        err = " ".join(l.strip() for l in open(os.path.join(d, "alfven.results")))
    fail = "FAIL" if (rc != 0 or not steps) else ""
    print(
        f"{tag:40s} rc={rc} steps={len(steps)} lits={lits} nlits={nlits} dt_end={dts[-1] if dts else None} wall={wall:.1f}s rss={rss:.0f}MB err=[{err}] {fail}",
        flush=True,
    )
    if fail and args.verbose:
        print(out[-3000:])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--kind", default="lu")
    ap.add_argument("--dt", default="8.E-7")
    ap.add_argument("--nsteps", default="5")
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--minlev", type=int, default=2)
    ap.add_argument("--nbase", type=int, default=1)
    ap.add_argument("--nlevels", type=int, default=1)
    ap.add_argument("--nlocal", type=int, default=8)
    ap.add_argument("--overlap", action="store_true")
    ap.add_argument("--its", type=int, default=2)
    ap.add_argument("--inner_its", type=int, default=None)
    ap.add_argument("--pbtype", type=int, default=1)
    ap.add_argument("--alfven", type=float, default=1.0)
    ap.add_argument("--sound", type=float, default=1.0)
    ap.add_argument("--kcycle", type=int, default=0)
    ap.add_argument("--cscale", action="store_true")
    ap.add_argument("--linear", action="store_true")
    ap.add_argument("--lin_tol", default="1.E-12")
    ap.add_argument("--maxits", type=int, default=200)
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--timeout", type=int, default=3000)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--inner", default="lu")
    ap.add_argument("--dump", action="store_true")
    run(ap.parse_args())
