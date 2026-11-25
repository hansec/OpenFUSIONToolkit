# SPDX-License-Identifier: LGPL-3.0-only
# Offline analysis script for the H(Curl)+Grad B-field conditioning study (see ../../README.md).
# Not part of the regression suite.
import os
import subprocess
TEMPLATE = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc=F
/
&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={nlevels}
 nbase={nlevels}
 grid_order=1
/
&cube_options
 mesh_type={mesh_type}
 ni=1,2,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/
&test_alfven_options
 order={order}
 minlev={minlev}
 delta=1.e-3
 linear=F
/
&xmhd_options
 xmhd_jcb=T
 xmhd_advec=F
 xmhd_hall=F
 xmhd_adv_temp=F
 xmhd_adv_den=F
 bbc='bc'
 nu_visc=1.E-3
 vbc='all'
 dt={dt}
 eta=1.e-3
 nsteps={nsteps}
 rst_freq=1000
 lin_tol=1.E-9
 nl_tol=1.E-5
 nu_xmhd={nu}
 nclean=5000
 rst_ind=0
 maxextrap=2
 ittarget=100
 xmhd_mfnk=F
 xmhd_nparts={nparts}
/
&xmhd_plot_options
 t0=1.
 dt=1.
 rst_start=0
 rst_end={nsteps}
/
"""
SMOOTHERS = {
 'bj': '<pre type="block_jacobi"><nlocal>{nlocal}</nlocal><solver type="lu"></solver></pre>',
 'bjov': '<pre type="block_jacobi"><nlocal>{nlocal}</nlocal><loverlap>T</loverlap><solver type="lu"></solver></pre>',
 'cell0': '<pre type="add_schwarz"><overlap_mode>cell</overlap_mode><overlap>0</overlap><solver type="lu"></solver></pre>',
 'cell1': '<pre type="add_schwarz"><overlap_mode>cell</overlap_mode><overlap>1</overlap><solver type="lu"></solver></pre>',
 'cell1f': '<pre type="add_schwarz"><overlap_mode>cell</overlap_mode><overlap>1</overlap><adjacency>face</adjacency><solver type="lu"></solver></pre>',
}
def run(tag, order, minlev, dt, nsteps, smoother, nparts=1, nlocal=8, mesh_type=1, nlevels=1, nu=None, smits=2, mg=True, timeout=3600):
    d = 'run_' + tag
    os.makedirs(d, exist_ok=True)
    nlev_fe = nlevels + order - 1
    if not mg:
        minlev = nlev_fe
    if nu is None:
        nu = ','.join(['0'] * (minlev - 1) + ['2'] * (nlev_fe - minlev + 1)) if mg else '-1'
    with open(os.path.join(d, 'oft.in'), 'w') as f:
        f.write(TEMPLATE.format(nlevels=nlevels, mesh_type=mesh_type, order=order, minlev=minlev, dt=dt,
                                nsteps=nsteps, nu=nu, nparts=','.join([str(nparts)] * nlev_fe)))
    pre = SMOOTHERS[smoother].format(nlocal=nlocal)
    if mg:
        xml = ('<oft><xmhd><pre type="mg"><smoother direction="up"><solver type="gmres"><its>{0}</its><nrits>{0}</nrits>'
               '{1}</solver></smoother></pre></xmhd></oft>').format(smits, pre)
    else:
        xml = '<oft><xmhd><pre type="gmres"><its>{0}</its><nrits>{0}</nrits>{1}</pre></xmhd></oft>'.format(smits, pre)
    with open(os.path.join(d, 'oft_in.xml'), 'w') as f:
        f.write(xml)
    env = dict(os.environ, OMP_NUM_THREADS='2')
    try:
        p = subprocess.run(['../test_alfven', 'oft.in', 'oft_in.xml'], cwd=d, capture_output=True, text=True,
                           timeout=timeout, env=env)
        out = p.stdout + p.stderr
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or b'').decode() if isinstance(e.stdout, bytes) else (e.stdout or '')
        out += '\nTIMEOUT'
    with open(os.path.join(d, 'out.log'), 'w') as f:
        f.write(out)
    return out
if __name__ == '__main__':
    print(run('probe', 2, 1, '8.E-7', 3, 'bj', nlocal=8)[-3000:])
