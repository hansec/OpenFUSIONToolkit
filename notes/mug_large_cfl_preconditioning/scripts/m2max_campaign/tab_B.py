import os, re, glob
import numpy as np
def parse_log(path):
    out = open(path, errors="replace").read()
    st = re.findall(r"Timestep\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s*$", out, re.M)
    return [(int(a), float(e), int(c), int(d), float(f)) for a, b, c, d, e, f in st], bool(re.search(r"ERROR|MPI_ABORT", out))
print('B: spheromak tilt from t=0, 8 steps, exclusive 4x2 unless noted. its/newton/step_s are per-step means; us/h simulated us per wall hour (steps only)')
print('%-34s %5s %6s %8s %7s %8s %6s %s' % ('case', 'its', 'newton', 'step_s', 'wall', 'us/h', 'GB', ''))
for tag in sorted(set(re.sub(r'_dt[^_]+(_np\d+t\d+)?$', '', os.path.basename(d)) for d in glob.glob('T/B/*'))):
    for dt in ('1e-6', '3e-6', '1e-5', '3e-5'):
        for d in sorted(glob.glob('T/B/%s_dt%s*' % (tag, dt))):
            if not os.path.exists(os.path.join(d, 'runq.txt')):
                print('%-34s %s' % (os.path.basename(d), '...')); continue
            rq = dict(kv.split('=') for kv in open(os.path.join(d, 'runq.txt')).read().split())
            steps, failed = parse_log(os.path.join(d, 'log.txt'))
            n = len(steps)
            if n == 0:
                print('%-34s FAILED rc=%s wall=%.0f' % (os.path.basename(d), rq['rc'], float(rq['wall']))); continue
            st = sum(x[1] for x in steps); sim = sum(x[4] for x in steps)
            cut = any(x[4] < 0.999 * float(dt) for x in steps)
            print('%-34s %5.0f %6.2f %8.2f %7.0f %8.1f %6.2f %s' % (os.path.basename(d), np.mean([x[2] for x in steps]), np.mean([x[3] for x in steps]),
                  st / n, float(rq['wall']), sim * 1e6 / st * 3600, float(rq['peak_tree_rss_GB']),
                  ('n=%d ' % n) + ('cut ' if cut else '') + ('F rc=%s' % rq['rc'] if failed or rq['rc'] != '0' else '')))
