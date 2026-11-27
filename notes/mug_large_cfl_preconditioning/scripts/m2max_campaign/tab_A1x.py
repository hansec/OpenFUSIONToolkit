import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from analyze_runs import parse_log
print('A1x: linear s=1, outer its/step [wall s, peak GB]; inner = pb1 default (~300-DOF bj LU), lu (exact LU of P), bj3k (~3000-DOF blocks), i20 (20 inner its)')
for R in (0, 1):
    for p in (2, 3, 4):
        for c in ('3', '30', '300'):
            row = 'R%d p%d CFL%-4s' % (R, p, c)
            d0 = 'A1/cube_R%d_p%d_cfl%s_s1_lin' % (R, p, c)
            for d in [d0] + ['A1x/cube_R%d_p%d_cfl%s_s1_lin_%s' % (R, p, c, t) for t in ('lu', 'bj3k', 'i20')]:
                if not os.path.exists(os.path.join(d, 'runq.txt')):
                    row += '%22s' % ('...' if os.path.isdir(d) else '-')
                    continue
                steps, failed = parse_log(os.path.join(d, 'log.txt'))
                rq = dict(kv.split('=') for kv in open(os.path.join(d, 'runq.txt')).read().split())
                its = '%.0f' % np.mean([x[2] for x in steps]) if steps else ''
                if failed or os.path.exists(os.path.join(d, 'TIMEOUT')):
                    its += 'F'
                row += '%22s' % ('%s [%.0f, %.1f]' % (its, float(rq['wall']), float(rq.get('peak_tree_rss_GB', 0))))
            print(row)
