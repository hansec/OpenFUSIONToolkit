"""Summarize run_bench.sh output in a results directory as Markdown.

Usage: python summarize.py <results_dir>
"""
import json
import os
import sys

res_dir = sys.argv[1]


def read_micro(name):
    path = os.path.join(res_dir, name)
    if not os.path.isfile(path):
        return None
    out = {}
    with open(path) as fid:
        next(fid)
        for line in fid:
            parts = line.split()
            out[(parts[0], int(parts[1]))] = (float(parts[2]), float(parts[3]), float(parts[4]))
    return out


def read_tc(tag):
    path = os.path.join(res_dir, 'tc_{0}.json'.format(tag))
    if not os.path.isfile(path):
        return None
    with open(path) as fid:
        return json.load(fid)


def fmt(val, spec='{0:.3f}'):
    return '—' if val is None else spec.format(val)


def ratio(base, new):
    if (base is None) or (new is None) or (new == 0.0):
        return '—'
    return '{0:.2f}x'.format(base/new)


lines = []
if os.path.isfile(os.path.join(res_dir, 'sysinfo.txt')):
    lines.append('## System')
    lines.append('```')
    with open(os.path.join(res_dir, 'sysinfo.txt')) as fid:
        lines.append(fid.read().rstrip())
    lines.append('```')
    lines.append('')

# Micro-benchmark
micro_cols = [
    ('OpenBLAS 1t', read_micro('micro_oblas_t1.txt')),
    ('Accel 1t', read_micro('micro_acc_t1.txt')),
    ('OpenBLAS-OMP Nt', read_micro('micro_oblasomp_tN.txt')),
    ('Accel Nt', read_micro('micro_acc_tN.txt')),
]
micro_cols = [col for col in micro_cols if col[1] is not None]
if micro_cols:
    keys = []
    for _, data in micro_cols:
        for key in data:
            if key not in keys:
                keys.append(key)
    lines.append('## BLAS/LAPACK kernels')
    lines.append('')
    lines.append('Time in seconds (dgemv: GFLOP/s, higher is better). '
                 '"Speedup" is OpenBLAS 1t / Accel 1t. Max accuracy check across columns.')
    lines.append('')
    lines.append('| kernel | N | ' + ' | '.join(col[0] for col in micro_cols) + ' | speedup 1t | max check |')
    lines.append('|---|---|' + '---|'*len(micro_cols) + '---|---|')
    base = dict(micro_cols).get('OpenBLAS 1t')
    new = dict(micro_cols).get('Accel 1t')
    for key in keys:
        cells = []
        checks = []
        for _, data in micro_cols:
            if key not in data:
                cells.append('—')
                continue
            if key[0] == 'dgemv':
                cells.append('{0:.1f}'.format(data[key][1]))
            else:
                cells.append('{0:.4f}'.format(data[key][0]))
            checks.append(data[key][2])
        spd = '—'
        if (base is not None) and (new is not None) and (key in base) and (key in new):
            if key[0] == 'dgemv':
                spd = ratio(new[key][1], base[key][1])
            else:
                spd = ratio(base[key][0], new[key][0])
        lines.append('| {0} | {1} | {2} | {3} | {4:.1e} |'.format(key[0], key[1], ' | '.join(cells), spd, max(checks)))
    lines.append('')

# ThinCurr
tc_cases = [
    ('oblas_tN', 'OpenBLAS, OFT N threads'),
    ('acc_tN', 'Accelerate, OFT N threads'),
    ('acc_tN_veclib1', 'Accelerate, OFT N threads, VECLIB_MAXIMUM_THREADS=1'),
    ('oblas_t1', 'OpenBLAS, 1 thread'),
    ('acc_t1', 'Accelerate, 1 thread'),
]
tc_data = [(label, read_tc(tag)) for tag, label in tc_cases]
tc_data = [entry for entry in tc_data if entry[1] is not None]
if tc_data:
    steps = ['Lmat', 'eigs_direct', 'freq_resp_direct', 'td_direct']
    lines.append('## ThinCurr dense (direct) solvers')
    lines.append('')
    lines.append('Torus 48x96 ({0} unknowns), time in seconds.'.format(tc_data[0][1]['timings'].get('nelems')))
    lines.append('')
    lines.append('| case | ' + ' | '.join(steps) + ' |')
    lines.append('|---|' + '---|'*len(steps))
    for label, data in tc_data:
        lines.append('| {0} | {1} |'.format(label, ' | '.join(fmt(data['timings'].get(step)) for step in steps)))
    lines.append('')
    ref = read_tc('oblas_tN') or read_tc('oblas_t1')
    if ref is not None:
        lines.append('Agreement with OpenBLAS (max relative difference):')
        lines.append('')
        for label, data in tc_data:
            eig_diff = max(abs(a - b)/abs(b) for a, b in zip(data['results']['eig_vals'], ref['results']['eig_vals']))
            fr_diff = abs(data['results']['freq_resp_norm'] - ref['results']['freq_resp_norm'])/abs(ref['results']['freq_resp_norm'])
            lines.append('- {0}: eigenvalues {1:.1e}, freq. response norm {2:.1e}'.format(label, eig_diff, fr_diff))
        lines.append('')

hodlr_cases = [
    ('oblas_tN_hodlr', 'acc_tN_hodlr'),
    ('oblas_tN_hodlrL', 'acc_tN_hodlrL'),
]
hodlr_lines = []
for tag_o, tag_a in hodlr_cases:
    data_o = read_tc(tag_o)
    data_a = read_tc(tag_a)
    if (data_o is None) or (data_a is None):
        continue
    eig_diff = max(abs(a - b)/abs(b) for a, b in zip(data_a['results']['eig_vals'], data_o['results']['eig_vals']))
    hodlr_lines.append('| {0} | {1} | {2} | {3} | {4} | {5:.1e} |'.format(
        data_o['timings']['nelems'],
        fmt(data_o['timings']['Lmat_hodlr']), fmt(data_a['timings']['Lmat_hodlr']),
        fmt(data_o['timings']['eigs_iterative']), fmt(data_a['timings']['eigs_iterative']), eig_diff))
if hodlr_lines:
    lines.append('## ThinCurr HODLR (compressed) path, OFT N threads')
    lines.append('')
    lines.append('| unknowns | Lmat_hodlr OpenBLAS | Lmat_hodlr Accel | eigs_iterative OpenBLAS | eigs_iterative Accel | eig rel. diff |')
    lines.append('|---|---|---|---|---|---|')
    lines += hodlr_lines
    lines.append('')

# Sparse direct (LU) solvers
sparse_labels = {
    'acc_accelerate': 'Accelerate sparse LU (default options)',
    'acc_accelerate_metis': 'Accelerate sparse LU, METIS ordering',
    'acc_accelerate_spp_metis': 'Accelerate sparse LU, supernodal pivoting + METIS',
    'acc_strumpack': 'STRUMPACK (Accelerate BLAS)',
    'oblas_strumpack': 'STRUMPACK (OpenBLAS)',
    'acc_umfpack': 'UMFPACK (Accelerate BLAS)',
    'oblas_umfpack': 'UMFPACK (OpenBLAS)',
    'acc_superlu': 'SuperLU (Accelerate BLAS)',
}
sparse_path = os.path.join(res_dir, 'sparse.txt')
if os.path.isfile(sparse_path):
    sparse_rows = []
    with open(sparse_path) as fid:
        for line in fid:
            parts = line.split()
            if len(parts) == 8:
                sparse_rows.append(parts)
    sizes = []
    for row in sparse_rows:
        if (row[1], row[2]) not in sizes:
            sizes.append((row[1], row[2]))
    lines.append('## Sparse direct (LU) solvers, OFT N threads')
    lines.append('')
    lines.append('Order-3 Lagrange Laplacian on a unit cube (`lin_alg/test_native_bjacobi`, GMRES with a '
                 'full LU preconditioner). "solve" is the first solve including factorization, in seconds; '
                 '"wall" also includes mesh and matrix setup. "norm diff" is the relative difference of the '
                 'solution norm from the first successful case.')
    lines.append('')
    for nlevels, ndofs in sizes:
        lines.append('### {0} unknowns (nlevels={1})'.format(ndofs, nlevels))
        lines.append('')
        lines.append('| package | status | solve | wall | GMRES its | norm diff |')
        lines.append('|---|---|---|---|---|---|')
        ref_norm = None
        for row in sparse_rows:
            if row[1] != nlevels:
                continue
            label = sparse_labels.get(row[0], row[0])
            if row[3] != 'ok':
                lines.append('| {0} | {1} | — | {2} | — | — |'.format(label, row[3], row[4]))
                continue
            norm = float(row[7])
            if ref_norm is None:
                ref_norm = norm
            lines.append('| {0} | ok | {1:.2f} | {2} | {3} | {4:.1e} |'.format(
                label, float(row[5]), row[4], row[6], abs(norm - ref_norm)/abs(ref_norm)))
        lines.append('')

print('\n'.join(lines))
