#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# ThinCurr tests: high-resolution cylinder, dense and with ACA+ (HODLR) matrix compression
import os
import pytest
import numpy as np
from thincurr_testing import ThinCurr_setup, validate_eigs, validate_fr, validate_td, mp_run, test_dir
from oft_testing import approx_list


@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_eig_aca(convert_xml):
    '''L/R eigenvalues of a high-resolution cylinder with the inductance matrix compressed by ACA+ (HODLR), against
    reference values'''
    eigs = approx_list((2.659575E-2, 1.254552E-2, 1.254536E-2, 1.208636E-2), rel=1.E-5)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'eig','F',use_aca=True,jumper_start=2,convert_xml=convert_xml)
    assert validate_eigs(eigs)


@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_aca(convert_xml):
    '''Time-dependent run with a compressed inductance matrix, plus a reduced model built from its slowest modes'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.205529E-4, 6.100309E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([5.453372E3, 5.453372E3], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           jumper_start=2,run_reduced=True)
    assert validate_td(sigs_final,jumpers_final)


@pytest.mark.coverage
@pytest.mark.parametrize("direct_flag", ('F', 'T'))
def test_fr_cyl_hr(direct_flag):
    '''Dense frequency response on the high-resolution cylinder; above 2000 elements the iterative solve uses the
    block-Jacobi preconditioner'''
    fr_real = approx_list((5.884322E-2, 4.882137E-3), rel=1.E-4)
    fr_imag = approx_list((-2.016855E-3, -2.314362E-3), rel=1.E-4)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'fr',direct_flag,freq=5.E3,fr_limit=0,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)


@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_fr_aca(convert_xml):
    '''Frequency response with a compressed inductance matrix, to a looser tolerance than the dense result'''
    fr_real = approx_list((5.888736E-2, 4.881440E-3), rel=1.E-3)
    fr_imag = approx_list((-2.017045E-3, -2.313881E-3), rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'fr','F',use_aca=True,freq=5.E3,fr_limit=0,convert_xml=convert_xml,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)


@pytest.mark.coverage
def test_td_aca_single():
    '''Time-dependent run with the compressed-matrix preconditioner stored in single precision'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.205529E-4, 6.100309E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([5.453372E3, 5.453372E3], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,hodlr_opts={'pre_single': True},
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)


@pytest.mark.coverage
def test_fr_aca_single():
    '''Frequency response with the compressed-matrix preconditioner stored in single precision'''
    fr_real = approx_list((5.888736E-2, 4.881440E-3), rel=1.E-3)
    fr_imag = approx_list((-2.017045E-3, -2.313881E-3), rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'fr','F',use_aca=True,hodlr_opts={'pre_single': True},freq=5.E3,fr_limit=0,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)


@pytest.mark.coverage
@pytest.mark.parametrize("convert_xml", (False, True))
def test_td_volt_aca(convert_xml):
    '''Voltage-driven coil time-dependent run with a compressed inductance matrix'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((1.823946E-4, 1.560471E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1.348327E3, 1348.327, 788.7134, 788.7133], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,convert_xml=convert_xml,
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)


@pytest.mark.parametrize("pre_single,pre_neigs", ((False, 0), (True, 0), (False, 2)))
def test_td_volt_aca_coarse(pre_single,pre_neigs):
    '''Voltage-driven coil run with the two-level (coarse-space) preconditioner, also in single precision and with
    the two slowest L/R eigenmodes added to the coarse space'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((1.823946E-4, 1.560471E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([1.348327E3, 1348.327, 788.7134, 788.7133], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,hodlr_opts={'pre_coarse': True, 'pre_single': pre_single,
                                                                          'pre_neigs': pre_neigs},
                           vcoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           volt_waveform=((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
                           jumper_start=2)
    assert validate_td(sigs_final,jumpers_final)


# HODLR block storage precision options: low-rank factors only, diagonal blocks only, and both (with the
# preconditioner also in single precision)
SINGLE_PREC_OPTS = (
    {'lowrank_single': True},
    {'dense_single': True},
    {'lowrank_single': True, 'dense_single': True, 'pre_single': True},
)
SINGLE_PREC_IDS = ('lowrank', 'dense', 'all')


@pytest.mark.coverage
@pytest.mark.parametrize("hodlr_opts", SINGLE_PREC_OPTS, ids=SINGLE_PREC_IDS)
def test_td_aca_block_single(hodlr_opts):
    '''Time-dependent run (and reduced model with B-field operators) with blocks of the compressed inductance and
    B-field operators stored in single precision'''
    sigs_final = approx_list([4.E-3], rel=1.E-8) +  approx_list((7.205529E-4, 6.100309E-4), rel=1.E-3)
    jumpers_final = approx_list([4.E-3], rel=1.E-8) +  approx_list([5.453372E3, 5.453372E3], rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'td','F',use_aca=True,hodlr_opts=hodlr_opts,
                           icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           curr_waveform=((0.0, 0.0), (1.0, 1.0)),
                           jumper_start=2,run_reduced=True)
    assert validate_td(sigs_final,jumpers_final)


@pytest.mark.coverage
@pytest.mark.parametrize("hodlr_opts", SINGLE_PREC_OPTS, ids=SINGLE_PREC_IDS)
def test_fr_aca_block_single(hodlr_opts):
    '''Frequency response with blocks of the compressed inductance matrix stored in single precision'''
    fr_real = approx_list((5.888736E-2, 4.881440E-3), rel=1.E-3)
    fr_imag = approx_list((-2.017045E-3, -2.313881E-3), rel=1.E-3)
    assert ThinCurr_setup("tw_test-cyl_hr.h5",'fr','F',use_aca=True,hodlr_opts=hodlr_opts,freq=5.E3,
                           fr_limit=0,icoils=((1.1, 0.25), (1.1, -0.25)),
                           floops=((0.9, 0.5), (0.9, 0.0)),
                           jumper_start=2)
    assert validate_fr(fr_real, fr_imag)


def run_hodlr_products(hodlr_opts,cache_tag,mp_q):
    '''Products of the inductance and B-field operators with a fixed potential on the high-resolution cylinder, with
    dense operators (`hodlr_opts=None`) or HODLR operators, and the precision of the diagonal blocks ('DENSE_') and
    low-rank factors ('SPARSE_U_') saved to the HODLR cache files'''
    result = None
    try:
        import h5py
        from OpenFUSIONToolkit import OFT_env
        from OpenFUSIONToolkit.ThinCurr import ThinCurr
        tw_model = ThinCurr(OFT_env(nthreads=-1))
        tw_model.setup_model(mesh_file='tw_test-cyl_hr.h5',jumper_start=2)
        use_hodlr = hodlr_opts is not None
        cache_files = (None, None)
        if use_hodlr:
            tw_model.set_hodlr_options(**hodlr_opts)
            if cache_tag is not None:
                cache_files = ('Lmat_{0}.save'.format(cache_tag), 'Bmat_{0}.save'.format(cache_tag))
        tw_model.compute_Lmat(use_hodlr=use_hodlr,cache_file=cache_files[0])
        tw_model.compute_Bmat(cache_file=cache_files[1])
        pot = np.sin(np.arange(tw_model.nelems)*0.37)
        Lx = np.array(tw_model.apply_Lmat(pot))
        Bx = np.array(tw_model.reconstruct_Bfield(pot))
        dtypes = {'DENSE_': set(), 'SPARSE_U_': set()}
        for cache_file in cache_files:
            if cache_file is None:
                continue
            with h5py.File(cache_file,'r') as fid:
                for key in fid.keys():
                    for prefix in dtypes:
                        if key.startswith(prefix) and not key.startswith(prefix+'hash'):
                            dtypes[prefix].add(str(fid[key].dtype))
        result = (Lx, Bx, dtypes)
    except BaseException as e:
        print(e)
    mp_q.put(result)


def test_hodlr_block_single():
    '''Rounding the low-rank factors and/or diagonal blocks of the HODLR operators to single precision changes
    products by much less than the HODLR approximation itself, cache files with different storage precisions are kept
    apart, and single precision blocks are reused from a cache file'''
    cache_tag = 'block_single'
    for cache_file in ('Lmat_{0}.save'.format(cache_tag), 'Bmat_{0}.save'.format(cache_tag)):
        cache_file = os.path.join(test_dir, cache_file)
        if os.path.exists(cache_file):
            os.remove(cache_file)
    dense = mp_run(run_hodlr_products, (None, None))
    double = mp_run(run_hodlr_products, ({}, cache_tag))
    assert dense is not None and double is not None
    assert double[2] == {'DENSE_': {'float64'}, 'SPARSE_U_': {'float64'}}
    for opts in SINGLE_PREC_OPTS[:2] + ({'lowrank_single': True, 'dense_single': True},):
        # Each setting must not reuse the cache written with the previous one, then must reuse its own
        single = mp_run(run_hodlr_products, (opts, cache_tag))
        single_cached = mp_run(run_hodlr_products, (opts, cache_tag))
        assert single is not None and single_cached is not None
        expected = {'DENSE_': {'float32' if opts.get('dense_single', False) else 'float64'},
                    'SPARSE_U_': {'float32' if opts.get('lowrank_single', False) else 'float64'}}
        assert single[2] == expected
        assert single_cached[2] == expected
        for i in range(2):
            ref = dense[i]
            scale = np.linalg.norm(ref)
            hodlr_err = np.linalg.norm(double[i]-ref)/scale
            single_diff = np.linalg.norm(single[i]-double[i])/scale
            print('{0} operator {1}: HODLR error = {2:.3E}, single vs double = {3:.3E}'.format(opts, i, hodlr_err, single_diff))
            assert hodlr_err < 1.E-3
            assert single_diff < 1.E-6
            assert np.linalg.norm(single[i]-ref)/scale < hodlr_err + 1.E-6
            # Blocks loaded from the cache are identical to those written
            assert np.allclose(single_cached[i], single[i], rtol=0.0, atol=1.E-12*np.abs(ref).max())
    double_after = mp_run(run_hodlr_products, ({}, cache_tag))  # Must not reuse a single precision cache
    assert double_after is not None
    assert double_after[2] == {'DENSE_': {'float64'}, 'SPARSE_U_': {'float64'}}
