# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../../README.md).
# Not part of the regression suite.
import os
import re
import subprocess
BASE = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'oft.in.base')).read()
SM = {
 'bj': '<pre type="block_jacobi"><nlocal>8</nlocal><solver type="lu"></solver></pre>',
 'bjov': '<pre type="block_jacobi"><nlocal>8</nlocal><loverlap>T</loverlap><solver type="lu"></solver></pre>',
 'cell0': '<pre type="add_schwarz"><overlap_mode>cell</overlap_mode><overlap>0</overlap><solver type="lu"></solver></pre>',
 'cell1': '<pre type="add_schwarz"><overlap_mode>cell</overlap_mode><overlap>1</overlap><solver type="lu"></solver></pre>',
}
def run(tag, sm, nparts, mg, dt='1.e-6', nsteps=3):
    d = 'run_' + tag; os.makedirs(d, exist_ok=True)
    for f in ('MUG_sph_tilt', 'cyl_tilt.h5'):
        if not os.path.exists(os.path.join(d, f)): os.symlink(os.path.abspath(f), os.path.join(d, f))
    nl = 2  # FE levels: mesh1-p1, mesh1-p2
    s = BASE.replace('@MINLEV@', '1' if mg else '2').replace('@DT@', dt).replace('@NSTEPS@', str(nsteps))
    s = s.replace('@NPARTS@', ','.join([str(nparts)] * nl)).replace('@NU@', '0,2' if mg else '-1')
    open(os.path.join(d, 'oft.in'), 'w').write(s)
    inner = SM[sm]
    if mg:
        xml = '<oft><xmhd><pre type="mg"><smoother direction="up"><solver type="gmres"><its>2</its><nrits>2</nrits>%s</solver></smoother></pre></xmhd></oft>' % inner
    else:
        xml = '<oft><xmhd><pre type="gmres"><its>2</its><nrits>2</nrits>%s</pre></xmhd></oft>' % inner
    open(os.path.join(d, 'oft_in.xml'), 'w').write(xml)
    env = dict(os.environ, OMP_NUM_THREADS='2')
    try:
        p = subprocess.run(['./MUG_sph_tilt', 'oft.in', 'oft_in.xml'], cwd=d, capture_output=True, text=True, timeout=3600, env=env)
        out = p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        out = 'TIMEOUT'
    open(os.path.join(d, 'out.log'), 'w').write(out)
    return [int(m.group(1)) for m in re.finditer(r'^ Timestep\s+\d+\s+\S+\s+(-?\d+)', out, re.M)], ('Run Complete' in out)
if __name__ == '__main__':
    for mg in (False, True):
        for sm, npart in (('bj', 1), ('bj', 8), ('bjov', 1), ('cell0', 8), ('cell1', 8)):
            tag = '%s_%s_np%d' % ('mg' if mg else '1lev', sm, npart)
            its, ok = run(tag, sm, npart, mg)
            print('%-22s ok=%s its=%s' % (tag, ok, its), flush=True)
