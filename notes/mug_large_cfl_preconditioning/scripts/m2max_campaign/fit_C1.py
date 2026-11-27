import sys, glob, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from analyze_runs import read_hist
T0, T1 = 1.0e-4, 3.0e-4
def gamma(d):
    h = read_hist(os.path.join(d, 'xmhd.hist'))
    t, ke = h['time'], h['ven']
    m = (t >= T0) & (t <= T1) & (ke > 0)
    if m.sum() < 5:
        return None, t[-1], h
    c = np.polyfit(t[m], np.log(ke[m]), 1)
    return 0.5 * c[0], t[-1], h   # KE ~ exp(2 gamma t)
ref, _, _ = gamma('C1/sph_lin_s0_dt1e-7')
print('gamma_ref (s=0, dt=1e-7) = %.5e 1/s' % ref)
for d in sorted(glob.glob('C1/sph_lin_*')):
    if not os.path.exists(os.path.join(d, 'xmhd.hist')):
        continue
    g, tend, h = gamma(d)
    done = os.path.exists(os.path.join(d, 'runq.txt'))
    txt = '-' if g is None else '%.5e  err=%+7.3f%%' % (g, 100 * (g / ref - 1))
    print('%-28s t_end=%.2e %s  its=%.0f  %s' % (d, tend, 'done' if done else 'run ', h['lits'].mean(), txt))
