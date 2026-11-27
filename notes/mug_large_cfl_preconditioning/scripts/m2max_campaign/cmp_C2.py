import sys, glob, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from analyze_runs import read_hist
ref = read_hist('C2/sph_R0_p2_ilu_s0_dt2e-7/xmhd.hist')
tr = ref['time']
print('ref reached t=%.3e' % tr[-1])
for d in sorted(glob.glob('C2/*/')):
    if not os.path.exists(d + 'xmhd.hist'):
        continue
    h = read_hist(d + 'xmhd.hist')
    t = h['time']
    if len(t) < 2:
        continue
    ke, tc, me = h['ven'], np.abs(h['tcurr']), h['men']
    i = np.argmax(ke)
    half = np.where(tc < 0.5 * tc[0])[0]
    # relative KE difference vs reference at common times (log space)
    tm = t[(t > 5e-5) & (t <= min(t[-1], tr[-1]))]
    dke = np.interp(tm, t, np.log(ke)) - np.interp(tm, tr, np.log(np.maximum(ref['ven'], 1e-300)))
    print('%-38s t_end=%.3e ke_peak=%.3e @%.3e  t_curr_half=%s  men_f/i=%.4f  max|dlnKE| vs ref=%.3f  lits_mean=%.0f nlits=%.2f  wall=%.0fs' % (
        d, t[-1], ke[i], t[i], ('%.3e' % t[half[0]]) if len(half) else '-', me[-1] / me[0],
        np.max(np.abs(dke)) if len(tm) else np.nan, h['lits'].mean(), h['nlits'].mean(), h['stime'].sum()))
