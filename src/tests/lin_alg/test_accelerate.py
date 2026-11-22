import os
import sys
import platform
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
/

&test_accelerate_options
 nsize={0}
 nrhs={1}
 fact_type={2}
 order={3}
 scaling={4}
 use_xml={5}
/
"""

# Template for XML input file
xml_template = """
<oft>
  <solver type="lu">
    <package>accelerate</package>
    <accelerate>
      <factorization>{0}</factorization>
      <ordering>{1}</ordering>
      <scaling>{2}</scaling>
    </accelerate>
  </solver>
</oft>
"""

# Some options require macOS 26 or later at runtime
def macos_major():
    version = platform.mac_ver()[0]
    if version == '':
        return 0
    return int(version.split('.')[0])
requires_macos26 = pytest.mark.skipif(macos_major() < 26, reason="Requires macOS 26 or later")

# Common setup function and process handling
def accelerate_setup(nsize, nrhs=1, fact_type=-1, order=-1, scaling=-1, xml_options=None):
    os.chdir(test_dir)
    use_xml = 'F'
    command = "./test_accelerate"
    if xml_options is not None:
        use_xml = 'T'
        command = "./test_accelerate oft.in oft_in.xml"
        with open('oft_in.xml','w+') as fid:
            fid.write(xml_template.format(*xml_options))
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nsize, nrhs, fact_type, order, scaling, use_xml))
    return run_OFT(command, 1, 60)

def validate_result(max_error, xml_options=None):
    with open('accelerate.results','r') as fid:
        errors = [float(fid.readline()) for _ in range(4)]
        if xml_options is not None:
            package = fid.readline().strip()
            options = [int(fid.readline()) for _ in range(3)]
    for i, error in enumerate(errors):
        if error > max_error:
            print("FAILED: Solution error too large (solve {0})!".format(i+1))
            print("  Expected < {0}".format(max_error))
            print("  Actual =   {0}".format(error))
            return False
    if xml_options is not None:
        if package != 'accelerate':
            print("FAILED: Incorrect package from XML ({0})".format(package))
            return False
        if options != xml_options:
            print("FAILED: Incorrect Accelerate options from XML!")
            print("  Expected = {0}".format(xml_options))
            print("  Actual =   {0}".format(options))
            return False
    return True

#============================================================================
# Test runner for sparse solves with each factorization type (0 -> LU, 1 -> unpivoted LU,
# 2 -> LU with supernodal pivoting, 3 -> LU with threshold pivoting) and number of RHS
@pytest.mark.coverage
@pytest.mark.parametrize("nrhs", (1, 3))
@pytest.mark.parametrize("fact_type", (0, 1, 2, 3))
def test_sparse(fact_type, nrhs):
    assert accelerate_setup(60, nrhs=nrhs, fact_type=fact_type)
    assert validate_result(1.E-10)

# Test runner for fill-reducing orderings (2 -> AMD, 3 -> METIS) and scaling (2 -> equilibration)
@pytest.mark.parametrize("order", (2, 3))
@pytest.mark.parametrize("scaling", (0, 2))
def test_sparse_options(order, scaling):
    assert accelerate_setup(60, nrhs=2, order=order, scaling=scaling)
    assert validate_result(1.E-10)

# Test runner for options requiring macOS 26 (5 -> MT-METIS ordering, 3 -> Hungarian scaling)
@requires_macos26
@pytest.mark.parametrize("options", ((5, 0), (0, 3)))
def test_sparse_options_macos26(options):
    assert accelerate_setup(60, nrhs=2, order=options[0], scaling=options[1])
    assert validate_result(1.E-10)

# Test runner for options read from XML
@pytest.mark.coverage
@pytest.mark.parametrize("options", (
    (("lu", "default", "default"), [0, 0, 0]),
    (("lu_unpivoted", "amd", "equilibrate"), [1, 2, 2]),
    (("lu_spp", "metis", "default"), [2, 3, 0]),
    (("lu_tpp", "amd", "default"), [3, 2, 0]),
))
def test_sparse_xml(options):
    assert accelerate_setup(60, nrhs=2, xml_options=options[0])
    assert validate_result(1.E-10, xml_options=options[1])
