import os
import numpy as np
OM0 = 1.9e5
print('A2 cube Alfven wave, T=1e-4: B-error (rel) per case; N=500 s=0 is the reference-dt run')
for R in (0, 1):
    for p in (2, 3, 4):
        print('-- R%d p%d   %-8s %-7s' % (R, p, 'N', 'CFL') + ''.join('%11s' % ('s=' + s) for s in ('0', '0.5', '1', '2')))
        for N in (500, 125, 25, 10, 5):
            cfl = OM0 * (p / 2.0) * 2**R * 1e-4 / N
            row = '           %-8d %-7.3g' % (N, cfl)
            for s in ('0', '0.5', '1', '2'):
                d = 'A2/cube_R%d_p%d_N%d_s%s_cn' % (R, p, N, s)
                r = os.path.join(d, 'alfven.results')
                if not os.path.isdir(d):
                    row += '%11s' % ''
                elif os.path.exists(r):
                    v = [float(x) for x in open(r).read().split()]
                    row += '%11.3e' % v[0]
                elif os.path.exists(os.path.join(d, 'runq.txt')):
                    row += '%11s' % 'F'
                else:
                    row += '%11s' % '...'
            print(row)
