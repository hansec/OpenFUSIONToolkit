#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
"""Summarize MUG benchmark case directories (made by bench_case.py, or any MUG run directory) into a CSV.

For each directory it reads log.txt (time-step lines), xmhd.hist (history file), alfven.results (cube waveform
errors), time.log (GNU time, peak RSS) and case.json (case parameters), and reports iteration, dt, cost and
physics metrics. For runs with a growing kinetic energy (spheromak tilt) it fits the mode growth rate
gamma = 0.5 * d ln(KE)/dt over a window (auto-detected, or --fit t0,t1 in seconds).

Usage:
  python3 analyze_runs.py case_dir [case_dir ...] [--csv summary.csv] [--fit 2e-5,1e-4] [--traces traces.csv]

--traces writes time, kinetic energy, toroidal current and magnetic energy of every case to one long CSV for
plotting physics comparisons (e.g. s=0 small-dt reference vs s=1 large-dt runs).
"""
import argparse
import csv
import json
import os
import re
import struct

import numpy as np


def read_hist(path):
    """Read an OFT history file; uses OpenFUSIONToolkit.io.histfile if importable, else a minimal parser."""
    try:
        from OpenFUSIONToolkit.io import histfile
        h = histfile(path)
        return {k: np.asarray(h[k]) for k in h.keys()}
    except Exception:
        pass
    content = open(path, "rb").read()
    start = content.find(b"--- BEGIN DATA ---")
    if start < 0:
        raise IOError("legacy history format not supported by the fallback reader: " + path)
    header = content[:start].decode(errors="replace")
    fields = re.search(r"^fields:(.*)$", header, re.M).group(1).split()
    types = re.search(r"^field_types:(.*)$", header, re.M).group(1).split()
    fmt = "=" + "".join({"i4": "i", "i8": "q", "r4": "f", "r8": "d"}[t] for t in types)
    rec = struct.calcsize(fmt)
    off = start + 19
    rows = []
    while off + 4 <= len(content):
        if content[off:off + 1] == b"\n":
            off += 1
            continue
        (n,) = struct.unpack_from("i", content, off)
        if n != rec:
            break
        rows.append(struct.unpack_from(fmt, content, off + 4))
        off += 4 + rec + 4
    data = np.array(rows, dtype=float) if rows else np.zeros((0, len(fields)))
    return {f: data[:, i] for i, f in enumerate(fields)}


def parse_log(path):
    out = open(path, errors="replace").read() if os.path.exists(path) else ""
    st = re.findall(r"Timestep\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s*$", out, re.M)
    if st:
        steps = [(int(a), float(b), int(c), int(d), float(f)) for a, b, c, d, e, f in st]
    else:  # linear runs: no Newton count
        st = re.findall(r"Timestep\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s*$", out, re.M)
        steps = [(int(a), float(b), int(c), 0, float(e)) for a, b, c, d, e in st]
    failed = bool(re.search(r"ERROR|oft_abort|failed to converge|MPI_ABORT", out))
    return steps, failed


def fit_growth(t, ke, window):
    """Return (gamma, r2, t0, t1) for gamma = 0.5 d ln(KE)/dt on the window; auto window if None."""
    ok = ke > 0
    t, ke = t[ok], ke[ok]
    if len(t) < 4:
        return (float("nan"),) * 4
    if window is None:
        kmax = ke.max()
        i0 = np.argmax(ke > 10.0 * ke[min(1, len(ke) - 1)]) if np.any(ke > 10.0 * ke[1]) else 1
        above = np.where(ke >= 0.1 * kmax)[0]
        i1 = above[0] if len(above) else len(ke) - 1
        if i1 - i0 < 3:
            i0, i1 = 1, len(ke) - 1
        sel = slice(i0, i1 + 1)
        tt, kk = t[sel], ke[sel]
    else:
        m = (t >= window[0]) & (t <= window[1])
        tt, kk = t[m], ke[m]
    if len(tt) < 3:
        return (float("nan"),) * 4
    y = np.log(kk)
    p = np.polyfit(tt, y, 1)
    res = y - np.polyval(p, tt)
    r2 = 1.0 - np.sum(res**2) / max(np.sum((y - y.mean())**2), 1e-300)
    return 0.5 * p[0], r2, tt[0], tt[-1]


def summarize(d, window):
    row = {"case": os.path.basename(os.path.normpath(d))}
    meta_path = os.path.join(d, "case.json")
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
    for k in ("geom", "refinements", "order", "ranks", "threads", "pre", "si", "linear", "mfnk", "dt", "nsteps",
              "nlocal", "inner", "inner_its"):
        row[k] = meta.get(k, "")
    steps, failed = parse_log(os.path.join(d, "log.txt"))
    row["failed"] = failed
    row["steps_done"] = len(steps)
    if steps:
        lits = np.array([s[2] for s in steps], dtype=float)
        nlits = np.array([s[3] for s in steps], dtype=float)
        dts = np.array([s[4] for s in steps])
        row["lits_mean"] = round(lits.mean(), 1)
        row["lits_max"] = int(lits.max())
        row["nlits_mean"] = round(nlits.mean(), 2)
        row["dt_min"] = dts.min()
        row["dt_last"] = dts[-1]
        try:
            dt_req = float(meta.get("dt", dts[0]))
        except ValueError:
            dt_req = dts[0]
        row["dt_reduced_steps"] = int(np.sum(dts < 0.999 * dt_req))
        row["t_end"] = steps[-1][1] + dts[-1]
    hist = os.path.join(d, "xmhd.hist")
    if os.path.exists(hist):
        h = read_hist(hist)
        t, ke = h.get("time", np.zeros(0)), h.get("ven", np.zeros(0))
        if "stime" in h and len(h["stime"]) > 1:
            st = h["stime"]
            row["wall_per_step_s"] = round(float(np.median(st[1:])), 3)
            row["wall_total_s"] = round(float(np.sum(st)), 1)
            if len(t) > 1 and t[-1] > t[0]:
                row["wall_per_simtime_s_per_us"] = round(float(np.sum(st)) / ((t[-1] - t[0]) * 1e6), 3)
        if len(ke):
            i = int(np.argmax(ke))
            row["ke_peak"] = float(ke[i])
            row["t_ke_peak"] = float(t[i])
            g, r2, t0, t1 = fit_growth(t, ke, window)
            row.update(gamma=g, gamma_r2=r2, fit_t0=t0, fit_t1=t1)
        if "tcurr" in h and len(h["tcurr"]) > 1:
            tc = np.abs(h["tcurr"])
            drop = np.where(tc < 0.5 * tc[0])[0]
            row["t_curr_half"] = float(t[drop[0]]) if len(drop) else ""
        if "men" in h and len(h["men"]):
            row["men_final_over_initial"] = float(h["men"][-1] / h["men"][0]) if h["men"][0] != 0 else ""
    res = os.path.join(d, "alfven.results")
    if os.path.exists(res):
        vals = [float(x) for x in open(res).read().split()]
        if len(vals) >= 2:
            row["berr"], row["verr"] = vals[0], vals[1]
    tlog = os.path.join(d, "time.log")
    if os.path.exists(tlog):
        tl = open(tlog).read()
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", tl)
        if m:
            row["peak_rss_GB"] = round(int(m.group(1)) / 1024**2, 2)
        else:  # BSD/macOS "time -l" reports bytes
            m = re.search(r"(\d+)\s+maximum resident set size", tl)
            if m:
                row["peak_rss_GB"] = round(int(m.group(1)) / 1024**3, 2)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--fit", default=None, help="growth-rate fit window t0,t1 [s]")
    ap.add_argument("--traces", default=None, help="write time traces of all cases to this CSV")
    a = ap.parse_args()
    window = tuple(float(x) for x in a.fit.split(",")) if a.fit else None
    rows = [summarize(d, window) for d in a.dirs]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    if a.csv:
        with open(a.csv, "w", newline="") as fid:
            w = csv.DictWriter(fid, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    for r in rows:
        print(", ".join("%s=%s" % (k, ("%.4g" % v) if isinstance(v, float) else v) for k, v in r.items()))
    if a.traces:
        with open(a.traces, "w", newline="") as fid:
            w = csv.writer(fid)
            w.writerow(["case", "time", "ven", "tcurr", "men", "lits", "nlits"])
            for d in a.dirs:
                hist = os.path.join(d, "xmhd.hist")
                if not os.path.exists(hist):
                    continue
                h = read_hist(hist)
                name = os.path.basename(os.path.normpath(d))
                for i in range(len(h["time"])):
                    w.writerow([name] + [h[k][i] if k in h else "" for k in ("time", "ven", "tcurr", "men", "lits",
                                                                              "nlits")])


if __name__ == "__main__":
    main()
