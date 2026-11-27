import sys, os, re, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from analyze_runs import parse_log
R = sys.argv[1] if len(sys.argv) > 1 else '0'
for mode in ('lin', 'cn', 'cnmf'):
    print('== R%s mode=%s  (outer its/step mean [Newton its]; F=failed, T=timeout, dt cuts as *)' % (R, mode))
    print('%-4s %-6s' % ('p', 'CFL') + ''.join('%14s' % ('s=' + s) for s in ('0', '0.5', '1', '2')))
    for p in (2, 3, 4):
        for c in ('0.3', '3', '30', '300'):
            row = '%-4d %-6s' % (p, c)
            for s in ('0', '0.5', '1', '2'):
                d = 'A1/cube_R%s_p%d_cfl%s_s%s_%s' % (R, p, c, s, mode)
                if not os.path.exists(os.path.join(d, 'log.txt')):
                    row += '%14s' % '-'
                    continue
                steps, failed = parse_log(os.path.join(d, 'log.txt'))
                if os.path.exists(os.path.join(d, 'TIMEOUT')):
                    failed = 'T'
                if steps:
                    import numpy as np
                    l = np.mean([x[2] for x in steps]); n = np.mean([x[3] for x in steps])
                    dt0 = steps[0][4]; cut = any(x[4] < 0.999 * max(y[4] for y in steps) for x in steps)
                    txt = '%.0f' % l + ('' if mode == 'lin' else '[%.1f]' % n) + ('*' if cut else '')
                else:
                    txt = ''
                if failed:
                    txt += ('T' if failed == 'T' else 'F')
                row += '%14s' % txt
            print(row)
