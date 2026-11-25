# SPDX-License-Identifier: LGPL-3.0-only
# Summarize run directories: python3 summarize.py runs runs_sph
import os
import re
import sys

for root in sys.argv[1:]:
    for t in sorted(os.listdir(root)):
        f = os.path.join(root, t, "log.txt")
        if not os.path.exists(f):
            continue
        out = open(f, errors="replace").read()
        st = re.findall(r"Timestep\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)", out)
        x = " ".join(open(os.path.join(root, t, "oft_in.xml")).read().split())
        err = ""
        ar = os.path.join(root, t, "alfven.results")
        if os.path.exists(ar):
            err = " ".join(l.strip() for l in open(ar))
        flag = "FAILED" if "failed" in out.lower() or "error" in out.lower() else ""
        print(
            f"{root}/{t:28s} steps={len(st)} lits={[int(s[2]) for s in st]} nlits={[int(s[3]) for s in st]} dt={[s[5] for s in st]} err=[{err}] {flag}"
        )
        print(f"    xml: {x}")
