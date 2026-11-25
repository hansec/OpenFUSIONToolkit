# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../../README.md).
# Not part of the regression suite.
"""Evaluate the vertex-star patch smoother (type="patch") on the refined Alfven cube (same setup as sweep2.py).

Usage: python sweep3.py <hcurl|lag> [smoother ...]
Smoothers:
  patch : MG, GMRES(2) smoother with symmetric vertex-star patches (up and down)
  hyb0  : MG, patch (down) + cell Schwarz overlap 0 (up)
  hyb1  : MG, patch (down) + cell Schwarz overlap 1 (up)
  hybov : MG, patch (down) + block Jacobi with loverlap (up)   (configuration of the patch report)
  patch1: one-level GMRES(2) + patch
"""
import os
import sys
import re
import json
import time
import subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runalf  # noqa: E402

PATCH = '<pre type="patch"><mode>symmetric</mode></pre>'
UP = {'hyb0': runalf.SMOOTHERS['cell0'], 'hyb1': runalf.SMOOTHERS['cell1'], 'hybov': runalf.SMOOTHERS['bjov']}


def gm(pre, its=2):
    return '<solver type="gmres"><its>{0}</its><nrits>{0}</nrits>{1}</solver>'.format(its, pre)


def xml_for(sm, nlocal=8):
    if sm == 'patch':
        return '<oft><xmhd><pre type="mg"><smoother direction="up">%s</smoother></pre></xmhd></oft>' % gm(PATCH)
    if sm == 'patch1':
        return '<oft><xmhd><pre type="gmres"><its>2</its><nrits>2</nrits>%s</pre></xmhd></oft>' % PATCH
    up = UP[sm].format(nlocal=nlocal)
    return ('<oft><xmhd><pre type="mg"><smoother direction="down">%s</smoother>'
            '<smoother direction="up">%s</smoother></pre></xmhd></oft>') % (gm(PATCH), gm(up))


def run(rep, tag, order, dt, sm, nsteps=6, nlevels=2, minlev=2, nparts=8, timeout=4000):
    d = 'run_' + tag
    os.makedirs(d, exist_ok=True)
    mg = sm != 'patch1'
    nlev_fe = nlevels + order - 1
    if not mg:
        minlev = nlev_fe
    nu = ','.join(['0'] * (minlev - 1) + ['2'] * (nlev_fe - minlev + 1)) if mg else '-1'
    tmpl = runalf.TEMPLATE
    if rep == 'lag':
        tmpl = tmpl.replace(" nu_visc=1.E-3\n", " nu_visc=1.E-3\n k_db=1.\n")
    with open(os.path.join(d, 'oft.in'), 'w') as f:
        f.write(tmpl.format(nlevels=nlevels, mesh_type=1, order=order, minlev=minlev, dt=dt, nsteps=nsteps, nu=nu,
                            nparts=','.join([str(nparts)] * nlev_fe)))
    with open(os.path.join(d, 'oft_in.xml'), 'w') as f:
        f.write(xml_for(sm))
    exe = '../test_alfven' if rep == 'hcurl' else '../test_alfven_lag'
    env = dict(os.environ, OMP_NUM_THREADS='2')
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, os.path.abspath('memrun.py'), exe, 'oft.in', 'oft_in.xml'], cwd=d,
                           capture_output=True, text=True, timeout=timeout, env=env)
        out = p.stdout + p.stderr
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or '')
        out += '\nTIMEOUT'
    with open(os.path.join(d, 'out.log'), 'w') as f:
        f.write(out)
    return out, time.time() - t0


def main():
    rep = sys.argv[1]
    sms = sys.argv[2:] or ['patch', 'hyb0', 'hyb1', 'hybov']
    results = {}
    fn = 'sweep3_%s.json' % rep
    if os.path.exists(fn):
        results = json.load(open(fn))
    for dt in ('8.E-7', '1.E-5', '1.E-4'):
        for order in (2, 3):
            for sm in sms:
                tag = '%s_L2_p%d_dt%s_%s' % (rep, order, dt, sm)
                if tag in results:
                    continue
                out, wall = run(rep, tag, order, dt, sm)
                lin = [int(m.group(1)) for m in re.finditer(r'^ Timestep\s+\d+\s+\S+\s+(\d+)\s+(\d+)', out, re.M)]
                rss = re.search(r'Maximum resident set size \(kbytes\): (\d+)', out)
                ok = 'Run Complete' in out
                results[tag] = {'lin': lin, 'ok': ok, 'wall': wall, 'rss_GB': int(rss.group(1)) / 1e6 if rss else None}
                avg = sum(lin) / len(lin) if lin else float('nan')
                print('%-40s ok=%s steps=%d avg_lin=%.1f max_lin=%s rss=%sGB wall=%.0fs' % (
                    tag, ok, len(lin), avg, max(lin) if lin else '-',
                    '%.2f' % results[tag]['rss_GB'] if rss else '?', wall), flush=True)
                json.dump(results, open(fn, 'w'))


if __name__ == '__main__':
    main()
