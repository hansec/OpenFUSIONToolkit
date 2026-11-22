"""ThinCurr BLAS/LAPACK benchmark: dense and HODLR solves on a torus mesh.

Usage: python bench_thincurr.py <mesh_dir> <ntheta> <nphi> <nthreads> <out_json> [hodlr]
PYTHONPATH must point at the OFT build-tree python directory.
"""
import json
import os
import sys
import time

import numpy as np

from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThinCurr import ThinCurr
from OpenFUSIONToolkit.ThinCurr.meshing import ThinCurr_periodic_toroid

mesh_dir, ntheta, nphi, nthreads, out_json = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
use_hodlr = (len(sys.argv) > 6) and (sys.argv[6] == 'hodlr')

os.makedirs(mesh_dir, exist_ok=True)
mesh_file = os.path.join(mesh_dir, 'torus_{0}x{1}.h5'.format(ntheta, nphi))
xml_file = os.path.join(mesh_dir, 'oft_in.xml')
if not os.path.isfile(mesh_file):
    theta = np.linspace(0.0, 2.0*np.pi, ntheta, endpoint=False)
    phi = np.linspace(0.0, 2.0*np.pi, nphi, endpoint=False)
    r_grid = np.zeros((nphi, ntheta, 3))
    for j in range(nphi):
        R = 1.0 + 0.4*np.cos(theta)
        r_grid[j, :, 0] = R*np.cos(phi[j])
        r_grid[j, :, 1] = R*np.sin(phi[j])
        r_grid[j, :, 2] = 0.4*np.sin(theta)
    ThinCurr_periodic_toroid(r_grid, 1, ntheta, nphi).write_to_file(mesh_file)
if not os.path.isfile(xml_file):
    with open(xml_file, 'w+') as fid:
        fid.write('<oft><thincurr><eta>1.257E-5</eta><icoils><coil_set>'
                  '<coil scale="1.0">1.1, 0.25</coil><coil scale="-1.0">1.1, -0.25</coil>'
                  '</coil_set></icoils></thincurr></oft>\n')

timings = {}


def timed(name, func, *args, **kwargs):
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    timings[name] = time.perf_counter() - t0
    print('BENCH {0:20s} {1:10.3f} s'.format(name, timings[name]), flush=True)
    return result


myOFT = OFT_env(nthreads=nthreads)
tw = ThinCurr(myOFT)
tw.setup_model(mesh_file=mesh_file, xml_filename=xml_file)
tw.setup_io()
timings['nelems'] = int(tw.nelems)
print('nelems = {0}'.format(tw.nelems), flush=True)
results = {}
if use_hodlr:
    tw.set_hodlr_options(target_size=512)
    timed('Lmat_hodlr', tw.compute_Lmat, use_hodlr=True)
    tw.compute_Rmat()
    eig_vals, _ = timed('eigs_iterative', tw.get_eigs, 10, direct=False)
    results['eig_vals'] = eig_vals.tolist()
else:
    timed('Lmat', tw.compute_Lmat)
    tw.compute_Rmat()
    eig_vals, _ = timed('eigs_direct', tw.get_eigs, 10, direct=True)
    results['eig_vals'] = eig_vals.tolist()
    Mc = tw.compute_Mcoil()
    fdriver = np.vstack((Mc[0, :], np.zeros((tw.nelems,))))
    fr = timed('freq_resp_direct', tw.compute_freq_response, fdriver=fdriver, freq=1.E3, direct=True)
    results['freq_resp_norm'] = float(np.linalg.norm(fr))
    coil_currs = np.array([[0.0, 0.0], [1.E-3, 1.E3], [1.0, 1.E3]])  # (time, Icoil)
    timed('td_direct', tw.run_td, 1.E-4, 200, coil_currs=coil_currs, direct=True, status_freq=1000, plot_freq=1000)
with open(out_json, 'w+') as fid:
    json.dump({'timings': timings, 'results': results}, fid, indent=1)
