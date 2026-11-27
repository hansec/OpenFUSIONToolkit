import os, sys, glob, re
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

def parse_log(path):
    out = open(path, errors="replace").read()
    st = re.findall(r"Timestep\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s*$", out, re.M)
    steps = [(int(a), float(e), int(c), int(d), float(f)) for a, b, c, d, e, f in st]
    return steps, bool(re.search(r"ERROR|oft_abort|MPI_ABORT", out))
CFG = ['ilu_s0_dt2e-7'] + ['bj%d_s0_dt2e-7' % n for n in (256, 64, 16, 4)] + ['asm16_s0_dt2e-7'] + \
      ['bj%d_s0_dt1e-6' % n for n in (64, 16, 4)] + ['asm16_s0_dt1e-6'] + ['pb%d_s1_dt1e-6' % n for n in (1024, 256, 64, 16)] + \
      ['bj16_s1_dt1e-6'] + ['pb%d_s0.5_dt3e-6' % n for n in (256, 64, 16)] + ['bj16_s0.5_dt3e-6'] + ['pb%d_s1_dt3e-6' % n for n in (64, 16)]
print('D: 20-step restart segments, 4 ranks x 2 threads, exclusive. step_s = mean step wall (excl. setup); us/h = simulated us per wall hour (full-dt steps only; restart ramps dt from 2e-7)')
for ph in ('lin', 'peak', 'late'):
    print('== phase %s' % ph)
    print('%-20s %6s %6s %7s %8s %8s %8s %7s' % ('config', 'its', 'newton', 'step_s', 'wall_s', 'us/h', 'xILU', 'GB'))
    base = None
    for c in CFG:
        d = 'T/D/%s_%s' % (ph, c)
        if not os.path.exists(os.path.join(d, 'runq.txt')):
            print('%-20s %s' % (c, '...' if os.path.isdir(d) else '-'))
            continue
        rq = dict(kv.split('=') for kv in open(os.path.join(d, 'runq.txt')).read().split())
        steps, failed = parse_log(os.path.join(d, 'log.txt'))
        dt = float(re.search(r'dt(\S+)$', c).group(1))
        if not steps:
            print('%-20s FAILED rc=%s' % (c, rq['rc'])); continue
        steps = [x for x in steps if x[4] >= 0.999 * dt]
        its = np.mean([x[2] for x in steps]); nl = np.mean([x[3] for x in steps]); st = np.mean([x[1] for x in steps])
        cuts = False; nfull = len(steps)
        simt = sum(x[4] for x in steps)
        ush = simt * 1e6 / (st * len(steps)) * 3600
        if c.startswith('ilu'):
            base = ush
        print('%-20s %6.0f %6.2f %7.2f %8.0f %8.1f %8.2f %7.2f %s%s' % (c, its, nl, st, float(rq['wall']), ush, ush / base if base else np.nan,
              float(rq['peak_tree_rss_GB']), 'n=%d' % nfull, ' F' if failed or rq['rc'] != '0' else ''))
