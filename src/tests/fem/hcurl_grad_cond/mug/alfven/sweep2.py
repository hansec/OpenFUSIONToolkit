# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../../README.md).
# Not part of the regression suite.
import os
import sys
import re
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runalf
from runalf import run
def parse(out):
    return [(int(m.group(1)), int(m.group(2))) for m in re.finditer(r'^ Timestep\s+\d+\s+\S+\s+(\d+)\s+(\d+)', out, re.M)]
rep = sys.argv[1]   # hcurl or lag
results = {}
cases = []
for dt in ('8.E-7', '1.E-5', '1.E-4'):
    for order, mg in ((2, True), (2, False), (3, True)):
        for sm, nparts in (('bj', 1), ('bj', 8), ('bjov', 1), ('cell0', 8), ('cell1', 8)):
            cases.append((order, dt, mg, sm, nparts))
if rep == 'lag':
    runalf.TEMPLATE = runalf.TEMPLATE.replace(" nu_visc=1.E-3\n", " nu_visc=1.E-3\n k_db=1.\n")
exe = '../test_alfven' if rep == 'hcurl' else '../test_alfven_lag'
orig = runalf.subprocess.run
def patched(cmd, **kw):
    cmd = [exe] + cmd[1:]
    return orig(cmd, **kw)
runalf.subprocess.run = patched
for order, dt, mg, sm, nparts in cases:
    tag = '%s_L2_p%d_dt%s_%s_%s_np%d' % (rep, order, dt, 'mg' if mg else '1lev', sm, nparts)
    out = run(tag, order, 2, dt, 6, sm, nparts=nparts, nlocal=8, nlevels=2, mg=mg, timeout=2400)
    its = parse(out)
    lin = [a for a, b in its]
    ok = 'Run Complete' in out
    results[tag] = {'lin': lin, 'ok': ok}
    avg = sum(lin) / len(lin) if lin else float('nan')
    print('%-50s ok=%s steps=%d avg_lin=%.1f max_lin=%s' % (tag, ok, len(lin), avg, max(lin) if lin else '-'), flush=True)
    json.dump(results, open('sweep2_%s.json' % rep, 'w'))
