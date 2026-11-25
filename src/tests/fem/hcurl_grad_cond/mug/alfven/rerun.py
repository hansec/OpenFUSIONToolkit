# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../../README.md).
# Not part of the regression suite.
import os
import sys
import re
import resource
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runalf
cases = [('hcurl', 3, '1.E-4', 'bjov', 1), ('hcurl', 3, '1.E-4', 'cell0', 8), ('hcurl', 3, '1.E-4', 'cell1', 8),
         ('lag', 3, '8.E-7', 'cell1', 8), ('lag', 3, '1.E-5', 'cell1', 8)]
orig = runalf.subprocess.run
base_t = runalf.TEMPLATE
for rep, order, dt, sm, npart in cases:
    exe = '../test_alfven' if rep == 'hcurl' else '../test_alfven_lag'
    runalf.TEMPLATE = base_t if rep == 'hcurl' else base_t.replace(" nu_visc=1.E-3\n", " nu_visc=1.E-3\n k_db=1.\n")
    runalf.subprocess.run = lambda cmd, **kw: orig(['bash', '-c', exe + ' ' + ' '.join(cmd[1:]) + '; echo EXITCODE=$?'], **kw)
    tag = 'rerun_%s_L2_p%d_dt%s_mg_%s_np%d' % (rep, order, dt, sm, npart)
    t0 = time.time()
    out = runalf.run(tag, order, 2, dt, 6, sm, nparts=npart, nlocal=8, nlevels=2, mg=True, timeout=6000)
    lin = [int(m.group(1)) for m in re.finditer(r'^ Timestep\s+\d+\s+\S+\s+(\d+)\s+(\d+)', out, re.M)]
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    sig = re.search(r'EXITCODE=(\d+)', out)
    print('%-50s ok=%s steps=%d avg_lin=%s maxrss_GB=%s exit=%s wall=%.0fs' % (tag, 'Run Complete' in out, len(lin),
          ('%.1f' % (sum(lin)/len(lin))) if lin else '-', '%.2f(max so far)' % (rss/1e6),
          sig.group(0) if sig else ('TIMEOUT' if 'TIMEOUT' in out else '?'), time.time()-t0), flush=True)
