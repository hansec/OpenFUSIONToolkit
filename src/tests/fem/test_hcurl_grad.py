import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={7}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={1}
 nbase={0}
 grid_order=1
/

&cube_options
 mesh_type={3}
/

&hcurl_grad_op_options
 df_mop={5}
 nu_mop={6}
/

&test_hcurl_grad_options
 order={2}
 mg_test={4}
 bc_test={8}
 poly_test={9}
 hex_enrich={10}
/
"""

# Common setup function and process handling
def hcurl_grad_setup(nbase, nlevels, order, grid_type, mg=False, df='', nu='', petsc=False, bc=False,
                     poly=False, enrich=False):
    petsc_flag=('T' if petsc else 'F')
    mg_flag=('T' if mg else 'F')
    bc_flag=('T' if bc else 'F')
    poly_flag=('T' if poly else 'F')
    enrich_flag=('T' if enrich else 'F')
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, mg_flag, df, nu, petsc_flag, bc_flag,
                                         poly_flag, enrich_flag))
    return run_OFT("./test_hcurl_grad", nproc, 120)

def validate_result(iteration_count,converged_error):
    fid = open('hcurl_grad.results','r')
    its_test = int(fid.readline())
    if iteration_count is not None:
        if abs(iteration_count-its_test) >= max(3,0.05*iteration_count):
            print("FAILED: Iteration count incorrect!")
            print("  Expected = {0}".format(iteration_count))
            print("  Actual =   {0}".format(its_test))
            return False
    error_test = float(fid.readline())
    if abs(converged_error-error_test) > 1.E-6:
        print("FAILED: Residual error incorrect!")
        print("  Expected = {0}".format(converged_error))
        print("  Actual =   {0}".format(error_test))
        return False
    return True

# Single level test function
def single_level(nlevels,order,exp_its,exp_error,grid_type=1,mpi=False,petsc=False,enrich=False):
    if mpi:
        nbase = nlevels-1
        nlev = nlevels+1
    else:
        nbase = nlevels
        nlev = nlevels
    assert hcurl_grad_setup(nbase, nlev, order, grid_type, petsc=petsc, enrich=enrich)
    assert validate_result(exp_its, exp_error)
# Multi-level test function
def multi_level(nlevels,order,exp_its,exp_error,grid_type=1,petsc=False,enrich=False):
    if grid_type == 1:
        df_string='0.0, 0.0, 0.39, 0.36, 0.217, 0.172, 0.142'
        nu_string='0, 0, 80, 16, 8, 4, 2'
    elif enrich:
        df_string='0.000, 0.000, 0.000, 0.3, 0.2, 0.17, 0.15'
        nu_string='0, 0, 0, 40, 16, 8, 4'
    else:
        df_string='0.000, 0.000, 0.000, 0.41, 0.25, 0.24, 0.21'
        nu_string='0, 0, 0, 40, 8, 4, 2'
    assert hcurl_grad_setup(2, nlevels+1, order, grid_type, mg=True, df=df_string, nu=nu_string, petsc=petsc,
                            enrich=enrich)
    assert validate_result(exp_its, exp_error)

# Mass matrix boundary condition test function
def bc_check(order, grid_type=1, mpi=False, enrich=False):
    if mpi:
        nbase, nlev = 2, 3
    else:
        nbase, nlev = 2, 2
    assert hcurl_grad_setup(nbase, nlev, order, grid_type, bc=True, enrich=enrich)
    with open('hcurl_grad_bc.results', 'r') as fid:
        for bc_name in ('none', 'zerob'):
            err = float(fid.readline())
            if err > 1.E-12:
                print('FAILED: Mass matrix BC "{0}" incorrect, max deviation = {1}'.format(bc_name, err))
                assert False

# Polynomial reproduction test function (distorted hexahedral mesh)
def poly_check(order, enrich):
    assert hcurl_grad_setup(2, 2, order, 2, poly=True, enrich=enrich)
    with open('hcurl_grad_poly.results', 'r') as fid:
        err_b = float(fid.readline())
        err_curl = float(fid.readline())
    if enrich:
        # Enriched basis contains the pullback of [P_k]^3 on trilinear cells
        if max(err_b, err_curl) > 1.E-8:
            print('FAILED: [P_k]^3 field not reproduced, errors = {0} {1}'.format(err_b, err_curl))
            assert False
    else:
        # Standard basis cannot reproduce [P_k]^3 on non-affine cells (see doc_fem_hcurl)
        if err_b < 1.E-6:
            print('FAILED: Unexpected reproduction of [P_k]^3 field, error = {0}'.format(err_b))
            assert False

#============================================================================
# Polynomial reproduction on distorted hexahedra (standard and enriched basis)
@pytest.mark.parametrize("order", (1, 2, 3, 4))
@pytest.mark.parametrize("enrich", (False, True))
def test_hex_poly(order, enrich):
    poly_check(order, enrich)

#============================================================================
# Mass matrix boundary conditions (eliminated rows must be identity rows)
@pytest.mark.parametrize("grid_type", (1, 2))
@pytest.mark.parametrize("order", (1, 2, 3, 4))
@pytest.mark.parametrize("mpi", (False, True))
def test_mop_bc(grid_type, order, mpi):
    bc_check(order, grid_type=grid_type, mpi=mpi)
@pytest.mark.parametrize("order", (1, 2, 3, 4))
@pytest.mark.parametrize("mpi", (False, True))
def test_hex_enrich_mop_bc(order, mpi):
    bc_check(order, grid_type=2, mpi=mpi, enrich=True)

#============================================================================
# Test runners for NP=1
def test_r2_p1(petsc_flag=False):
    single_level(2, 1, 57, 1.0, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 61, 1.0, mpi=mpi, petsc=petsc_flag)
def test_r3_p1_mg(petsc_flag=False):
    multi_level(3, 1, 9, 1.0, petsc=petsc_flag)

#============================================================================
# Test runners for NP=2
def test_r2_p2(petsc_flag=False):
    single_level(2, 2, 167, 1.0, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 181, 1.0, mpi=mpi, petsc=petsc_flag)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2_mg(petsc_flag):
    multi_level(3, 2, 29, 1.0, petsc=petsc_flag)

#============================================================================
# Test runners for NP=3
def test_r2_p3(petsc_flag=False):
    single_level(2, 3, 361, 1.0, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 379, 1.0, mpi=mpi, petsc=petsc_flag)
def test_r3_p3_mg(petsc_flag=False):
    multi_level(3, 3, 64, 1.0, petsc=petsc_flag)

#============================================================================
# Test runners for NP=4
def test_r2_p4(petsc_flag=False):
    single_level(2, 4, 712, 1.0, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 743, 1.0, mpi=mpi, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 156, 1.0, petsc=petsc_flag)

#============================================================================
# Test runners for NP=1
def test_hex_r2_p1(petsc_flag=False):
    single_level(2, 1, 7, 1.0, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 7, 1.0, grid_type=2, mpi=mpi, petsc=petsc_flag)
# def test_hex_r3_p1_mg(petsc_flag=False):
#     multi_level(3, 1, 1, 1.0, grid_type=2, petsc=petsc_flag)

#============================================================================
# Test runners for NP=2
def test_hex_r2_p2(petsc_flag=False):
    single_level(2, 2, 192, 1.0, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 278, 1.0, grid_type=2, mpi=mpi, petsc=petsc_flag)
def test_hex_r3_p2_mg(petsc_flag=False):
    multi_level(3, 2, 55, 1.0, grid_type=2, petsc=petsc_flag)

#============================================================================
# Test runners for NP=3
def test_hex_r2_p3(petsc_flag=False):
    single_level(2, 3, 456, 1.0, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 580, 1.0, grid_type=2, mpi=mpi, petsc=petsc_flag)
def test_hex_r3_p3_mg(petsc_flag=False):
    multi_level(3, 3, 122, 1.0, grid_type=2, petsc=petsc_flag)

#============================================================================
# Test runners for NP=4
def test_hex_r2_p4(petsc_flag=False):
    single_level(2, 4, 903, 1.0, grid_type=2, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 1035, 1.0, grid_type=2, mpi=mpi, petsc=petsc_flag)
def test_hex_r3_p4_mg(petsc_flag=False):
    multi_level(3, 4, 248, 1.0, grid_type=2, petsc=petsc_flag)

#============================================================================
# Test runners for the enriched hexahedral basis
@pytest.mark.parametrize("order,exp_its", ((1, 34), (2, 297), (3, 858), (4, 1027)))
def test_hex_enrich_r2(order, exp_its):
    single_level(2, order, exp_its, 1.0, grid_type=2, enrich=True)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("order,exp_its", ((1, 61), (2, 376), (3, 1035), (4, 1110)))
def test_hex_enrich_r3(order, exp_its, mpi):
    single_level(3, order, exp_its, 1.0, grid_type=2, mpi=mpi, enrich=True)
@pytest.mark.parametrize("order,exp_its", ((2, 52), (3, 159)))
def test_hex_enrich_r3_mg(order, exp_its):
    multi_level(3, order, exp_its, 1.0, grid_type=2, enrich=True)
