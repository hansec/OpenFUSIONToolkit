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
 use_petsc={0}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={1}
 nbase=3
 grid_order=1
/

&cube_options
 mesh_type=1
/
"""
#
xml_template = """
<oft>
  {0}
</oft>
"""

# Common setup function and process handling
def xml_setup(xml_file, petsc=False, nprocs=1):
    # Distributed runs decompose the base mesh and refine once more
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    petsc_flag=('T' if petsc else 'F')
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(petsc_flag, nlevels))
    with open('oft_in.xml','w+') as fid:
        fid.write(xml_template.format(xml_file))
    return run_OFT("./test_solver_xml oft.in oft_in.xml", nprocs, 60)

def validate_result(iteration_count,converged_error):
    #
    fid = open('xml.results','r')
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

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_gmres(petsc_flag):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
    </solver>
    """
    iteration_count = 530
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_gmres_jacobi(petsc_flag):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="jacobi"></pre>
    </solver>
    """
    iteration_count = 492
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_cg(petsc_flag):
    xml_file = """
    <solver type="cg">
        <its>-3</its>
    </solver>
    """
    iteration_count = 145
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_cg_jacobi(petsc_flag):
    xml_file = """
    <solver type="cg">
        <its>-3</its>
        <pre type="jacobi"></pre>
    </solver>
    """
    iteration_count = 154
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=petsc_flag)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Native overlapping Schwarz selected from XML: blocks within the process grown by one layer
# (plain additive combination counts overlap rows twice, hence more iterations; partition of
# unity weights avoid this)
@pytest.mark.parametrize("schwarz_type", ("restricted", "additive", "weighted"))
def test_gmres_schwarz(schwarz_type):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <nlocal>4</nlocal>
            <overlap>1</overlap>
            <type>{0}</type>
            <pre type="ilu"></pre>
        </pre>
    </solver>
    """.format(schwarz_type)
    iteration_count = {"restricted": 58, "additive": 102, "weighted": 59}[schwarz_type]
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=False)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Native ILU(0) options from XML (restricted Schwarz sub-solver as above): a small-pivot threshold
# leaves the factorization of this well-conditioned operator unchanged, while a diagonal shift
# factors a perturbed operator and so weakens the preconditioner; the threaded (level-scheduled)
# kernels give identical results, so disabling them changes nothing
@pytest.mark.parametrize("ilu_opts", ("pivot_tol", "diag_shift", "threaded"))
def test_gmres_schwarz_ilu_options(ilu_opts):
    opts = {"pivot_tol": "<pivot_tol>1.E-8</pivot_tol>",
            "diag_shift": "<pivot_tol>0.0</pivot_tol><diag_shift>0.5</diag_shift>",
            "threaded": "<threaded>F</threaded>"}[ilu_opts]
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <nlocal>4</nlocal>
            <overlap>1</overlap>
            <type>restricted</type>
            <pre type="ilu">{0}</pre>
        </pre>
    </solver>
    """.format(opts)
    iteration_count = {"pivot_tol": 58, "diag_shift": 111, "threaded": 58}[ilu_opts]
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=False)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Native ILU(0) on the whole process-local subdomain (Schwarz with one block, applied from serial
# code, where the level-scheduled kernels may use all threads) with and without threading: the
# threaded kernels are bitwise identical to the serial ones, so the iteration count is the same
@pytest.mark.parametrize("threaded", ("T", "F"))
def test_gmres_schwarz_ilu_threaded(threaded):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <nlocal>1</nlocal>
            <overlap>1</overlap>
            <type>restricted</type>
            <pre type="ilu"><threaded>{0}</threaded></pre>
        </pre>
    </solver>
    """.format(threaded)
    iteration_count = 58
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=False)
    assert validate_result(iteration_count,converged_error)

#============================================================================
# Native overlapping Schwarz selected from XML: one subdomain per process grown by one layer
# across processes
@pytest.mark.mpi
@pytest.mark.parametrize("schwarz_type", ("restricted", "additive", "weighted"))
def test_gmres_schwarz_mpi2(schwarz_type):
    xml_file = """
    <solver type="gmres">
        <its>-3</its>
        <nrits>20</nrits>
        <pre type="add_schwarz">
            <nlocal>1</nlocal>
            <overlap>1</overlap>
            <type>{0}</type>
            <pre type="ilu"></pre>
        </pre>
    </solver>
    """.format(schwarz_type)
    iteration_count = {"restricted": 57, "additive": 85, "weighted": 59}[schwarz_type]
    converged_error = 5.1147216970434615
    assert xml_setup(xml_file,petsc=False,nprocs=2)
    assert validate_result(iteration_count,converged_error)
