#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Field line tracing tests (3D tracers on an analytic helical field)
import os
import sys
import h5py
import numpy as np
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', 'fem')))
from test_vlag_periodic import structured_mesh, x_wedge

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={0}
 nbase={1}
 grid_order=1
/

&cube_options
 mesh_type={3}
 ni={10}
 rscale=2.d0,2.d0,2.d0
 shift=-1.d0,-1.d0,-1.d0
 ref_per={7}
/

&test_tracing_options
 tracer_type={2}
 test_poincare={4}
 field_jump={5}
 tol={6}
 field_type={8}
 max_steps={9}
/
"""

# Input file for the native wedge mesh with revolved periodicity
oft_in_wedge_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='wedge'
 cad_type=0
 nlevels={0}
 nbase={1}
 grid_order=1
/

&native_mesh_options
 filename='{3}'
 reflect=T
 ref_periodic=T
/

&test_tracing_options
 tracer_type={2}
 tol=1.E-10
 field_type=3
 max_steps=100000
/
"""

# Step limits for traces on periodic meshes, which do not leave the mesh. LSODE takes much longer
# steps than the RK tracer, whose steps end at cell boundaries
PERIODIC_STEPS = {1: 1000, 4: 10000}


# Common setup function and process handling
def tracing_setup(tracer_type, parallel=False, mesh_type=1, poincare=False, field_jump=0.0, tol=1.E-10,
                  ref_per='F,F,F', field_type=1, max_steps=100000, coarse_cells=1):
    # Both give the same final mesh, as the level where the mesh is distributed is not refined
    if parallel:
        nlevels = 4
        nbase = 3
        nproc = 2
    else:
        nlevels = 3
        nbase = 3
        nproc = 1
    # Similar resolution from a base mesh with 2x2x2 cells refined once less
    if coarse_cells == 2:
        nlevels -= 1
        nbase -= 1
    poincare_flag = 'F'
    if poincare:
        poincare_flag = 'T'
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels, nbase, tracer_type, mesh_type, poincare_flag, field_jump, tol, ref_per,
                                         field_type, max_steps, '{0},{0},{0}'.format(coarse_cells)))
    return run_OFT("./test_tracing", nproc, 120)


def wedge_setup(tracer_type, parallel, cell_type):
    '''Trace a field line rotating about the x-axis through a wedge mesh with revolved periodicity'''
    # Both give the same final mesh, as the level where the mesh is distributed is not refined
    if parallel:
        nlevels = 3
        nbase = 2
        nproc = 2
    else:
        nlevels = 2
        nbase = 2
        nproc = 1
    os.chdir(test_dir)
    # Wedge of half-angle pi/6 about the x-axis (0 <= x <= 1), reflected about z = 0 by OFT. Tet
    # wedges include the axis, hex wedges have an inner radius of 0.5
    r_min = {'tet': 0.0, 'hex': 0.5}[cell_type]
    r, lc, per_nodes = structured_mesh(4, 2, 2, lambda u, v, w: x_wedge(u, v, w, r_min), cell_type)
    mesh_file = 'tracing_wedge_{0}.h5'.format(cell_type)
    with h5py.File(mesh_file, 'w') as h5_file:
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        h5_file.create_dataset('mesh/LC', data=lc+1, dtype='i4')
        h5_file.create_dataset('mesh/REG', data=np.ones((lc.shape[0],)), dtype='i4')
        h5_file.create_dataset('mesh/periodicity/nodes', data=per_nodes+1, dtype='i4')
    with open('oft.in','w+') as fid:
        fid.write(oft_in_wedge_template.format(nlevels, nbase, tracer_type, mesh_file))
    return run_OFT("./test_tracing", nproc, 120)


def jump_trace(parallel, mesh_type, tol):
    '''Trace a field line with jumps in the field strength between cells and return its points'''
    assert tracing_setup(4, parallel, mesh_type, field_jump=0.5, tol=tol)
    with h5py.File('trace_test', 'r') as fid:
        return np.asarray(fid['trace'])

# Tracer type 1 is LSODE (Adams method), type 4 is the cell-limited
# Bogacki-Shampine 3(2) Runge-Kutta tracer
# Mesh type 1 is tetrahedral, 2 is hexahedral (different cell exit tests for type 4)
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_tracing(tracer_type, parallel, mesh_type):
    '''Single field line traced with tracing_line, checked against the analytic helical field line'''
    assert tracing_setup(tracer_type, parallel, mesh_type)

# Poincare sections of several field lines traced concurrently (tracing_poincare),
# which spawns per-line tracers from the template using the tracer's copy method
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_poincare(tracer_type, parallel, mesh_type):
    '''Poincare section of several field lines (tracing_poincare), checked against the analytic crossings'''
    assert tracing_setup(tracer_type, parallel, mesh_type, poincare=True)

# Mesh that is periodic in z: the helical field line is continued across the periodic faces
# (z wraps from 1 to -1) until the step limit, and its end point is compared to the exact
# solution with z mapped into the periodic domain. In parallel, periodic partners of faces can be
# on the same or another domain.
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_tracing_periodic(tracer_type, parallel, mesh_type):
    '''Single field line on a mesh that is periodic in z, continued across the periodic faces'''
    assert tracing_setup(tracer_type, parallel, mesh_type, ref_per='F,F,T', max_steps=PERIODIC_STEPS[tracer_type])


# Mesh that is periodic in x and y (and z): the field line crosses the x and y periodic
# boundaries at the same time (edges of the periodic domain), and for the triply periodic
# mesh also the z boundaries (corners). The base mesh has 2x2x2 cells, as the global grounding
# point of a doubly periodic single-cell hex base mesh would have four periodic copies.
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
@pytest.mark.parametrize("ref_per", ('T,T,F','T,T,T'))
def test_tracing_periodic_xyz(tracer_type, parallel, mesh_type, ref_per):
    '''Single field line on a mesh periodic in several directions, crossing edges and corners of the domain'''
    assert tracing_setup(tracer_type, parallel, mesh_type, ref_per=ref_per, field_type=2,
                         max_steps=PERIODIC_STEPS[tracer_type], coarse_cells=2)


# Wedge about the x-axis with revolved periodicity: the field line rotates about the x-axis,
# repeatedly crossing the periodic planes, until it leaves the mesh at x=1
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("cell_type", ('tet','hex'))
def test_tracing_revolved(tracer_type, parallel, cell_type):
    '''Single field line on a wedge with revolved periodicity, continued across the periodic planes'''
    assert wedge_setup(tracer_type, parallel, cell_type)


# Poincare sections on a mesh that is periodic in z: crossings of x=0 are checked against the
# exact helices with z mapped into the periodic domain, and the jumps across the periodic
# boundaries must not be counted as crossings of the plane z=0
@pytest.mark.coverage
@pytest.mark.parametrize("tracer_type", (1,4))
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_poincare_periodic(tracer_type, parallel, mesh_type):
    '''Poincare sections on a mesh periodic in z, without spurious crossings at the periodic boundaries'''
    assert tracing_setup(tracer_type, parallel, mesh_type, poincare=True, ref_per='F,F,T',
                         max_steps=PERIODIC_STEPS[tracer_type])


# Field whose strength jumps at every cell face, as for finite element fields. The field line
# is still the exact helix (checked by the test program), and the time along it (fourth
# column) is a piecewise linear function of the angle, which is compared to a serial trace with
# a tight tolerance on the same mesh. This covers the retries of steps that overshoot a face
# too far and, in parallel, the steps that end on another domain.
@pytest.mark.coverage
@pytest.mark.parametrize("parallel", (False,True))
@pytest.mark.parametrize("mesh_type", (1,2))
def test_tracing_jump(parallel, mesh_type):
    '''Field line through a field that jumps between cells, checked against a tight-tolerance reference'''
    ref = jump_trace(False, mesh_type, 1.E-11)
    trace = jump_trace(parallel, mesh_type, 1.E-8)
    phase_ref = np.unwrap(np.arctan2(ref[:,1], ref[:,0]))
    phase = np.unwrap(np.arctan2(trace[:,1], trace[:,0]))
    inside = phase <= phase_ref[-1]
    err = np.abs(trace[inside,3] - np.interp(phase[inside], phase_ref, ref[:,3]))
    print('Max time error = {0:.3E}, points = {1}'.format(err.max(), inside.sum()))
    assert err.max() < 2.E-6
