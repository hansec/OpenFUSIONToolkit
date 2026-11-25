# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../../README.md).
# Not part of the regression suite.
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import re
import sys
import json
from runalf import run
def parse(out):
    its = [(int(m.group(1)), int(m.group(2))) for m in re.finditer(r'^ Timestep\s+\d+\s+\S+\s+(\d+)\s+(\d+)', out, re.M)]
    return its
mesh_levels = int(sys.argv[1]) if len(sys.argv) > 1 else 1
results = {}
for order in (2, 3):
    for dt in ('8.E-7', '1.E-5', '1.E-4'):
        for mg in (False, True):
            for sm, nparts in (('bj', 1), ('bj', 8), ('bjov', 1), ('cell0', 8), ('cell1', 8)):
                tag = 'L%d_p%d_dt%s_%s_%s_np%d' % (mesh_levels, order, dt, 'mg' if mg else '1lev', sm, nparts)
                out = run(tag, order, mesh_levels, dt, 10, sm, nparts=nparts, nlocal=8, nlevels=mesh_levels, mg=mg, timeout=1800)
                its = parse(out)
                lin = [a for a, b in its]
                ok = ('Run Complete' in out)
                results[tag] = {'lin': lin, 'nl': [b for a, b in its], 'ok': ok}
                avg = sum(lin) / len(lin) if lin else float('nan')
                print('%-45s ok=%s steps=%d avg_lin=%.1f max_lin=%s' % (tag, ok, len(lin), avg, max(lin) if lin else '-'), flush=True)
                json.dump(results, open('sweep_L%d.json' % mesh_levels, 'w'))
