#!/usr/bin/env python3
"""Run bench_case.py case directories under a core and memory budget.

Usage: runq.py --cores 10 --mem 52 case1[:cores[:memGB]] case2 ...
Cores default to ranks*threads from case.json; memory estimate defaults to 2 GB.
A watchdog kills the most recently started job if the total RSS of all jobs exceeds --mem.
"""
import argparse
import json
import os
import subprocess
import sys
import time


def tree_rss_gb(pid):
    """Total RSS (GB) of pid and all its descendants."""
    out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="], capture_output=True, text=True).stdout
    kids, rss = {}, {}
    for line in out.split("\n"):
        f = line.split()
        if len(f) != 3:
            continue
        p, pp, r = int(f[0]), int(f[1]), int(f[2])
        kids.setdefault(pp, []).append(p)
        rss[p] = r
    tot, stack = 0, [pid]
    while stack:
        p = stack.pop()
        tot += rss.get(p, 0)
        stack += kids.get(p, [])
    return tot / 1024**2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cores", type=int, default=10)
    ap.add_argument("--mem", type=float, default=52.0, help="total RSS cap [GB]")
    ap.add_argument("--timeout", type=float, default=0.0, help="per-job wall limit [s] (0: none)")
    ap.add_argument("--cores_file", default=None, help="file holding the core budget, re-read every loop")
    ap.add_argument("--list", default=None, help="file with one case spec per line")
    ap.add_argument("cases", nargs="*")
    a = ap.parse_args()
    queue = []
    if a.list:
        a.cases += [l.strip() for l in open(a.list) if l.strip() and not l.startswith("#")]
    for c in a.cases:
        f = c.split(":")
        d = f[0]
        meta = json.load(open(os.path.join(d, "case.json")))
        ncore = int(f[1]) if len(f) > 1 and f[1] else meta["ranks"] * meta["threads"]
        mem = float(f[2]) if len(f) > 2 else 2.0
        queue.append((d, ncore, mem))
    running = []  # (proc, d, ncore, mem, t0, peak)
    while queue or running:
        # watchdog / reap
        tot = 0.0
        for j in running:
            j[5] = max(j[5], tree_rss_gb(j[0].pid))
            tot += tree_rss_gb(j[0].pid)
        if tot > a.mem and running:
            j = running[-1]
            print("[runq] RSS %.1f GB > cap: killing %s" % (tot, j[1]), flush=True)
            subprocess.run(["pkill", "-TERM", "-g", str(j[0].pid)])
            open(os.path.join(j[1], "KILLED_MEM"), "w").write("%.1f\n" % tot)
        for j in running:
            if a.timeout > 0 and time.time() - j[4] > a.timeout and j[0].poll() is None:
                print("[runq] timeout: killing %s" % j[1], flush=True)
                subprocess.run(["pkill", "-TERM", "-g", str(j[0].pid)])
                open(os.path.join(j[1], "TIMEOUT"), "w").write("%.0f\n" % a.timeout)
        for j in list(running):
            if j[0].poll() is not None:
                print("[runq] done %s rc=%d %.0fs peak_tree_rss=%.2fGB" % (j[1], j[0].returncode, time.time() - j[4],
                                                                          j[5]), flush=True)
                open(os.path.join(j[1], "runq.txt"), "w").write("rc=%d wall=%.1f peak_tree_rss_GB=%.2f\n" % (
                    j[0].returncode, time.time() - j[4], j[5]))
                running.remove(j)
        # launch
        if a.cores_file and os.path.exists(a.cores_file):
            try:
                a.cores = int(open(a.cores_file).read().split()[0])
            except (ValueError, IndexError):
                pass
        used_c = sum(j[2] for j in running)
        used_m = sum(max(j[3], j[5]) for j in running)
        for q in list(queue):
            d, nc, m = q
            if (used_c + nc <= a.cores and used_m + m <= a.mem) or (not running and a.cores > 0):
                p = subprocess.Popen(["bash", os.path.join(d, "run.sh")], start_new_session=True)
                running.append([p, d, nc, m, time.time(), 0.0])
                used_c += nc
                used_m += m
                queue.remove(q)
                print("[runq] start %s (cores=%d mem~%.1fGB)" % (d, nc, m), flush=True)
        time.sleep(5)


if __name__ == "__main__":
    sys.exit(main())
