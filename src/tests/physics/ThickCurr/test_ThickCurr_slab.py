#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThickCurr tests: conducting slab between pseudo-vacuum layers inside a perfectly conducting unit cube,
# compared with the analytic decay time of the slowest mode
import pytest
from thickcurr_testing import thickcurr_setup, write_box_mesh, read_results, check_value, slab_tau

ETA = 1.256637E-6  # Slab resistivity [Ohm-m]
SLAB = (0.4, 0.6)
MESH_FILE = 'thickcurr_slab.h5'


def run_slab_eigs(order, nproc=1):
    # Layers in z align with the slab faces; MPI runs add a level that distributes the mesh over the ranks.
    # Note: the pseudo-vacuum makes the LHS ill-conditioned, so a direct solve of the full system is used. In
    # parallel the default block Jacobi + LU preconditioner does not converge for this case.
    write_box_mesh(MESH_FILE, (4, 4, 10), slab=SLAB)
    nlevels = 2 if nproc > 1 else 1
    assert thickcurr_setup('eigs', MESH_FILE, 2, ETA, order=order, nlevels=nlevels, dt=2.E-2, nev=2,
                           nproc=nproc)
    tau = read_results()
    tau_exp = slab_tau(ETA, *SLAB)
    retval = True
    # Slowest mode is 2-fold degenerate (cos(pi*x) and cos(pi*y) variations)
    for i in range(2):
        retval &= check_value("tau_{0}".format(i+1), tau[i], tau_exp, 2.E-2)
    return retval


@pytest.mark.parametrize("order", (2, 3))
def test_slab_eigs(order):
    '''Slowest decay time of a conducting slab surrounded by pseudo-vacuum (thickcurr_eigs), compared with the
    analytic solution from matching potential fields in the vacuum layers to the slab eigenmode.'''
    assert run_slab_eigs(order)

