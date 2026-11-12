#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Marklin tests: lowest Taylor (minimum-energy force-free) eigenstate of a cylinder through the Python
# interface, on straight-sided (linear) and curved (quadratic) mesh geometry [J. B. Taylor,
# Phys. Rev. Lett. 33, 1139 (1974)]; eigenvalues are compared to regression values
import numpy
import pytest
from marklin_testing import mp_run, clear_cache, run_marklin, validate_eigs, run_marklin_interp


def marklin_lin_cell(order):
    eigs = [5.029946,4.958005,4.957746,4.957744]
    if order == 1:
        clear_cache()
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,order,1,None))
    assert validate_eigs(results,eigs[order-1])

@pytest.mark.coverage
def test_marklin_g1_p1():
    '''Linear elements on straight-sided (linear geometry) cells.'''
    marklin_lin_cell(1)

@pytest.mark.coverage
def test_marklin_g1_p2():
    '''Quadratic elements on straight-sided (linear geometry) cells.'''
    marklin_lin_cell(2)

def test_marklin_g1_p3():
    '''Cubic elements on straight-sided (linear geometry) cells.'''
    marklin_lin_cell(3)

def test_marklin_g1_p4():
    '''Quartic elements on straight-sided (linear geometry) cells.'''
    marklin_lin_cell(4)

def marklin_quad_cell(order):
    eigs = [5.027066,4.955207,4.954955,4.954955]
    if order == 1:
        clear_cache()
    results = mp_run(run_marklin,('cyl_Marklin.h5',1,order,2,None))
    assert validate_eigs(results,eigs[order-1])

@pytest.mark.coverage
def test_marklin_g2_p1():
    '''Linear elements on curved (quadratic geometry) cells.'''
    marklin_quad_cell(1)

@pytest.mark.coverage
def test_marklin_g2_p2():
    '''Quadratic elements on curved (quadratic geometry) cells.'''
    marklin_quad_cell(2)

def test_marklin_g2_p3():
    '''Cubic elements on curved (quadratic geometry) cells.'''
    marklin_quad_cell(3)

def test_marklin_g2_p4():
    '''Quartic elements on curved (quadratic geometry) cells.'''
    marklin_quad_cell(4)

@pytest.mark.coverage
def test_marklin_interp():
    '''Vector potential (both gauges) and magnetic field interpolators for the lowest eigenstate on curved
    (quadratic geometry) cells with quadratic elements, checked against regression values of the field
    magnitude at fixed points. Single-point and batched evaluation must agree and points outside the mesh
    must evaluate to NaN.'''
    mag_vals = {
        'A': [0.3605469,0.2121165,0.1862212],
        'A_zn': [0.2172684,0.1539679,0.1328045],
        'B': [1.075002,0.75984,0.6667451],
    }
    results = mp_run(run_marklin_interp,('cyl_Marklin.h5',2,2))
    assert results is not None
    assert results['consistent']
    for key, vals in mag_vals.items():
        assert numpy.all(numpy.isfinite(results[key][:-1]))
        assert numpy.isnan(results[key][-1])
        assert numpy.allclose(results[key][:-1],vals,rtol=1.E-4)
