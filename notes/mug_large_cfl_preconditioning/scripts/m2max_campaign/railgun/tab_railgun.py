import os, re, subprocess, sys
import numpy as np
CASES = ['coax_base_np3', 'coax_base_dt1.5e-8', 'coax_base_dt6e-8', 'coax_pb1_s0_dt3e-8', 'coax_pb1_s1_dt3e-8',
         'coax_pb1_s1_dt6e-8', 'coax_pb1_s0_dt1.2e-7', 'coax_pb1_s1_dt1.2e-7']
print('Coax railgun (rotational periodicity), 3 MPI ranks x 1 thread, run concurrently with other jobs (wall times indicative only).')
print('base = example preconditioner (GMRES(2) + block-Jacobi LU, 1 block/rank); pb1 = <pb_pre type=1> with GMRES(8)+BJ-LU(440/rank) on P')
print('%-24s %8s %4s %5s %7s %7s %6s %7s  %s' % ('case', 'dt', 's', 'pre', 'its', 'newton', 'GB', 'wall_s', 'sheet speed r=0.11/0.15/0.18 m [km/s]'))
for d in CASES:
    inp = open(os.path.join(d, 'oft_coax.in')).read()
    dt = re.search(r'^ dt=(\S+)', inp, re.M).group(1)
    m = re.search(r'^ xmhd_si=(\S+)', inp, re.M)
    s = m.group(1) if m else '0'
    pre = 'pb1' if 'pb_pre' in open(os.path.join(d, 'oft_coax.xml')).read() else 'base'
    st = re.findall(r'Timestep\s+\d+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+\S+', open(os.path.join(d, 'log.txt')).read())
    its = np.mean([int(a) for a, b in st]); nl = np.mean([int(b) for a, b in st])
    rq = dict(kv.split('=') for kv in open(os.path.join(d, 'runq.txt')).read().split())
    spd = re.findall(r'Coax, r = \S+ m\s+(\S+)', open(os.path.join(d, 'sheet_speed.txt')).read())
    print('%-24s %8s %4s %5s %7.1f %7.2f %6s %7.0f  %s' % (d, dt, s, pre, its, nl, rq.get('peak_tree_rss_GB', '-'), float(rq['wall']), ' / '.join(spd)))
print('\nRelative L2 field differences vs coax_base_np3 (same mesh and decomposition, common plot times)')
sys.stdout.flush()
subprocess.run([sys.executable, 'cmp_fields.py'] + CASES)
