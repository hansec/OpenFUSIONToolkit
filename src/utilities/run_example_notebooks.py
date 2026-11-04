#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
#
# Script for re-executing the example notebooks in "src/examples" against a build
# tree and copying successfully executed notebooks back to the source tree, with
# noisy metadata (timestamps, editor and Python version information) removed and
# all outputs retained.
#
# Notebooks are executed in a fresh copy of the examples tree, so files they
# generate (meshes, results) never reach the source tree; checked-in data files
# that a notebook regenerates with different contents are listed in the summary.
# Within each directory, mesh/geometry generators run first, followed by the order
# used in "test_examples.sh.in" and then the remaining notebooks alphabetically.
# Notebooks that fail are reported and left untouched in the source tree.
#
# Requires "nbformat" and "nbclient", and a Jupyter kernel whose environment can
# import the OFT Python dependencies (numpy, h5py, matplotlib, ...); "--check-kernel"
# tests that. Notebooks needing extra packages can use another kernel through
# "--kernel-for GLOB=KERNEL". A kernel may also be given as the path to a Python
# interpreter (eg. a dedicated virtual environment with "ipykernel" installed), which
# is used through a temporary kernelspec without registering anything globally.
#
# Example:
#   python run_example_notebooks.py --build-dir <build_dir> --work-dir /tmp/nb_run --update
#   python run_example_notebooks.py --build-dir <build_dir> --check-kernel
#   python run_example_notebooks.py --strip-only     # only clean metadata in place
#
#------------------------------------------------------------------------------
import argparse
import atexit
import concurrent.futures
import fnmatch
import glob
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import time

import nbformat

# Notebook-level metadata to keep (everything else is editor/tool noise)
NB_META_KEEP = ('kernelspec', 'language_info')
# Fields of "language_info" that change with the executing environment
LANG_INFO_DROP = ('version',)
# Cell-level metadata to keep (user intent); timestamps, ids, scroll state, etc. are dropped
CELL_META_KEEP = ('tags', 'jupyter', 'editable', 'deletable', 'slideshow')
# Notebook name patterns for generators that other notebooks in a directory depend on
GENERATOR_RE = re.compile(r'(mesh|geom|generator)', re.IGNORECASE)
# Longest file path OFT accepts ("oft_path_slen" in OpenFUSIONToolkit/_core.py), less a
# margin for files the notebooks create with names longer than the checked-in ones
OFT_PATH_LIMIT = 200
PATH_MARGIN = 20


def file_hash(path):
    '''Return the SHA-256 hex digest of a file's contents'''
    with open(path, 'rb') as fid:
        return hashlib.sha256(fid.read()).hexdigest()


def clean_notebook(nb, kernelspec=None):
    '''Remove noisy metadata from a notebook in place, keeping all outputs

    If `kernelspec` is given it replaces the notebook's kernelspec, so that executing
    with a different local kernel does not change the stored one.
    '''
    meta = {}
    for key in NB_META_KEEP:
        if key in nb.metadata:
            meta[key] = nb.metadata[key]
    if kernelspec is not None:
        meta['kernelspec'] = kernelspec
    if 'language_info' in meta:
        meta['language_info'] = {
            key: value for key, value in meta['language_info'].items() if key not in LANG_INFO_DROP
        }
    nb.metadata = nbformat.from_dict(meta)
    for cell in nb.cells:
        cell_meta = {}
        for key in CELL_META_KEEP:
            if cell.metadata.get(key):
                cell_meta[key] = cell.metadata[key]
        cell.metadata = nbformat.from_dict(cell_meta)
    return nb


def notebook_text(nb):
    '''Serialize a notebook exactly as it would be written to disk'''
    text = nbformat.writes(nb)
    if not text.endswith('\n'):
        text += '\n'
    return text


def read_notebook(path):
    '''Read a notebook as nbformat v4'''
    with open(path, 'r', encoding='utf-8') as fid:
        return nbformat.read(fid, as_version=4)


def write_notebook(nb, path):
    '''Write a notebook to disk in the standard Jupyter format'''
    with open(path, 'w', encoding='utf-8') as fid:
        fid.write(notebook_text(nb))


def test_script_order(examples_root):
    '''Return {notebook relative path: index} from the "run_test" lines of test_examples.sh.in'''
    order = {}
    for script in sorted(glob.glob(os.path.join(examples_root, '*', 'test_examples.sh.in'))):
        tool = os.path.basename(os.path.dirname(script))
        with open(script, 'r') as fid:
            for line in fid:
                fields = line.split()
                if len(fields) >= 3 and fields[0] == 'run_test':
                    rel_path = os.path.normpath(os.path.join(tool, fields[1], fields[2]))
                    order.setdefault(rel_path, len(order))
    return order


def plan_notebooks(examples_root, selection):
    '''Group selected notebooks by directory, in execution order'''
    if selection:
        paths = []
        for pattern in selection:
            matches = glob.glob(os.path.join(examples_root, pattern), recursive=True)
            if not matches:
                raise ValueError('No notebooks match "{0}"'.format(pattern))
            paths += matches
    else:
        paths = glob.glob(os.path.join(examples_root, '**', '*.ipynb'), recursive=True)
    rel_paths = sorted(set(
        os.path.relpath(path, examples_root) for path in paths
        if path.endswith('.ipynb') and '.ipynb_checkpoints' not in path
    ))
    order = test_script_order(examples_root)
    groups = {}
    for rel_path in rel_paths:
        groups.setdefault(os.path.dirname(rel_path), []).append(rel_path)
    plan = []
    for dirname in sorted(groups):
        def sort_key(rel_path):
            name = os.path.basename(rel_path)
            return (GENERATOR_RE.search(name) is None, order.get(rel_path, len(order)), name)
        plan.append((dirname, sorted(groups[dirname], key=sort_key)))
    return plan


def kernel_for(rel_path, opts):
    '''Return the kernel for a notebook: the first matching "--kernel-for" pattern, else the default'''
    for pattern, kernel in opts['kernel_overrides']:
        if fnmatch.fnmatch(rel_path, pattern):
            return kernel
    return opts['kernel']


def interpreter_kernels(kernels):
    '''Map kernels given as interpreter paths to temporary kernelspecs

    Kernelspecs are written to a temporary directory added to JUPYTER_PATH (inherited by
    worker processes), and removed at exit. Returns {given value: kernel name}.
    '''
    names = {}
    spec_root = None
    for value in kernels:
        if not (os.path.sep in value and os.path.isfile(value)):
            names[value] = value
            continue
        if spec_root is None:
            spec_root = tempfile.mkdtemp(prefix='oft_nb_kernels_')
            atexit.register(shutil.rmtree, spec_root, True)
            os.environ['JUPYTER_PATH'] = os.pathsep.join([spec_root] + [path for path in [os.environ.get('JUPYTER_PATH')] if path])
        name = 'oft-nb-' + hashlib.sha256(value.encode()).hexdigest()[:12]
        spec_dir = os.path.join(spec_root, 'kernels', name)
        os.makedirs(spec_dir)
        with open(os.path.join(spec_dir, 'kernel.json'), 'w') as fid:
            json.dump({
                'argv': [value, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
                'display_name': value,
                'language': 'python',
            }, fid)
        names[value] = name
    return names


def check_kernel(kernel, label=None):
    '''Start a kernel and report its interpreter and whether the OFT Python dependencies import'''
    from jupyter_client.manager import start_new_kernel
    code = (
        "import sys, importlib\n"
        "missing = []\n"
        "for name in ('numpy', 'scipy', 'h5py', 'matplotlib'):\n"
        "    try:\n"
        "        importlib.import_module(name)\n"
        "    except ImportError:\n"
        "        missing.append(name)\n"
        "print(sys.executable, 'missing:', ','.join(missing) or 'none')\n"
    )
    manager, client = start_new_kernel(kernel_name=kernel)
    output = []
    try:
        client.execute_interactive(
            code, timeout=120,
            output_hook=lambda msg: output.append(msg['content'].get('text', '')) if msg['msg_type'] == 'stream' else None
        )
    finally:
        client.stop_channels()
        manager.shutdown_kernel(now=True)
    text = ''.join(output).strip()
    print('kernel "{0}": {1}'.format(label or kernel, text))
    return 0 if text.endswith('missing: none') else 1


def error_summary(nb):
    '''Return (cell index, "ename: evalue") for the first error output in a notebook'''
    for index, cell in enumerate(nb.cells):
        for output in cell.get('outputs', []):
            if output.get('output_type') == 'error':
                return index, '{0}: {1}'.format(output.get('ename'), output.get('evalue'))
    return None, None


def run_directory(dirname, rel_paths, opts):
    '''Execute the notebooks of one directory, with kernel console output sent to a log file

    Output written directly by compiled code (eg. Fortran solver progress) bypasses the
    notebook and goes to the kernel's stdout/stderr, which the kernel inherits from this
    worker process; redirect it to "<work dir>/<dirname>/kernel_output.log".
    '''
    log_path = os.path.join(opts['work_root'], dirname, 'kernel_output.log')
    log_fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    saved_fds = (os.dup(1), os.dup(2))
    sys.stdout.flush()
    sys.stderr.flush()
    os.dup2(log_fd, 1)
    os.dup2(log_fd, 2)
    try:
        results, modified = execute_directory(dirname, rel_paths, opts)
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved_fds[0], 1)
        os.dup2(saved_fds[1], 2)
        for fd in saved_fds + (log_fd,):
            os.close(fd)
    for result in results:
        if result['status'] == 'failed':
            result['kernel_log'] = log_path
    return results, modified


def execute_directory(dirname, rel_paths, opts):
    '''Execute the notebooks of one directory in the work copy and update the source tree'''
    from nbclient import NotebookClient
    results = []
    work_dir = os.path.join(opts['work_root'], dirname)
    for rel_path in rel_paths:
        src_path = os.path.join(opts['examples_root'], rel_path)
        work_path = os.path.join(opts['work_root'], rel_path)
        result = {'notebook': rel_path, 'status': None, 'seconds': 0.0, 'kernel': kernel_for(rel_path, opts)}
        nb = read_notebook(work_path)
        kernelspec = nb.metadata.get('kernelspec')
        client_args = {
            'timeout': opts['timeout'],
            'kernel_name': opts['kernel_names'][result['kernel']],
            'resources': {'metadata': {'path': work_dir}},
            'record_timing': False,
        }
        if 'coalesce_streams' in NotebookClient.class_trait_names():
            client_args['coalesce_streams'] = True
        start = time.time()
        try:
            NotebookClient(nb, **client_args).execute()
        except Exception as exc:
            result['seconds'] = round(time.time() - start, 1)
            cell, error = error_summary(nb)
            if error is None:
                # Timeouts and kernel deaths leave no error output in the notebook
                lines = str(exc).strip().splitlines()
                error = '{0}: {1}'.format(type(exc).__name__, lines[-1] if lines else '')
            result['status'] = 'failed'
            result['cell'] = cell
            result['error'] = error
            result['failed_copy'] = work_path.replace('.ipynb', '.failed.ipynb')
            write_notebook(nb, result['failed_copy'])
            results.append(result)
            continue
        result['seconds'] = round(time.time() - start, 1)
        clean_notebook(nb, kernelspec=kernelspec)
        write_notebook(nb, work_path)
        old_text = notebook_text(clean_notebook(read_notebook(src_path)))
        if file_hash(src_path) != opts['source_hashes'][rel_path]:
            result['status'] = 'source_changed'
        elif notebook_text(nb) == old_text:
            result['status'] = 'unchanged'
        elif opts['update']:
            write_notebook(nb, src_path)
            result['status'] = 'updated'
        else:
            result['status'] = 'would_update'
        results.append(result)
    # Checked-in data files that the notebooks regenerated with different contents
    modified = []
    src_dir = os.path.join(opts['examples_root'], dirname)
    for name in sorted(os.listdir(src_dir)):
        src_file = os.path.join(src_dir, name)
        work_file = os.path.join(work_dir, name)
        if name.endswith('.ipynb') or not os.path.isfile(src_file) or not os.path.isfile(work_file):
            continue
        if file_hash(src_file) != file_hash(work_file):
            modified.append(os.path.join(dirname, name))
    return results, modified


def main():
    parser = argparse.ArgumentParser(description=__doc__ or 'Re-execute OFT example notebooks')
    parser.add_argument('notebooks', nargs='*', help='Notebook paths or globs relative to src/examples (default: all)')
    parser.add_argument('--build-dir', help='OFT build tree (used as OFT_ROOTPATH)')
    parser.add_argument('--examples-root', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'examples'),
                        help='Examples source tree (default: src/examples next to this script)')
    parser.add_argument('--work-dir', help='Scratch directory for the executed copy of the examples tree '
                        '(must not exist; default: a new temporary directory, kept after the run)')
    parser.add_argument('--kernel', help='Jupyter kernel name or Python interpreter path (default: OFT_PY_KERNEL from the build tree, else python3)')
    parser.add_argument('--kernel-for', action='append', default=[], metavar='GLOB=KERNEL',
                        help='Use KERNEL (a kernel name or Python interpreter path) for notebooks matching GLOB (relative to src/examples; repeatable)')
    parser.add_argument('--check-kernel', action='store_true',
                        help='Check that the kernel(s) start and import the OFT Python dependencies, then exit')
    parser.add_argument('--timeout', type=int, default=1800, help='Per-cell timeout in seconds (default: 1800)')
    parser.add_argument('--jobs', type=int, default=1, help='Directories to run concurrently (default: 1)')
    parser.add_argument('--threads', type=int, default=2, help='OMP_NUM_THREADS for the kernels (default: 2)')
    parser.add_argument('--update', action='store_true', help='Copy successfully executed notebooks back to the source tree')
    parser.add_argument('--summary', help='Write a JSON summary to this path')
    parser.add_argument('--list', action='store_true', help='Print the execution plan and exit')
    parser.add_argument('--strip-only', action='store_true', help='Only remove noisy metadata from the source notebooks, in place')
    args = parser.parse_args()

    examples_root = os.path.abspath(args.examples_root)
    try:
        plan = plan_notebooks(examples_root, args.notebooks)
    except ValueError as exc:
        parser.error(str(exc))
    if args.list:
        for dirname, rel_paths in plan:
            for rel_path in rel_paths:
                print(rel_path)
        return 0
    if args.strip_only:
        for dirname, rel_paths in plan:
            for rel_path in rel_paths:
                path = os.path.join(examples_root, rel_path)
                with open(path, 'r', encoding='utf-8') as fid:
                    old_text = fid.read()
                new_text = notebook_text(clean_notebook(read_notebook(path)))
                if new_text != old_text:
                    with open(path, 'w', encoding='utf-8') as fid:
                        fid.write(new_text)
                    print('stripped   {0}'.format(rel_path))
        return 0

    if args.build_dir is None:
        parser.error('--build-dir is required unless using --list or --strip-only')
    build_dir = os.path.abspath(args.build_dir)
    if not os.path.isdir(os.path.join(build_dir, 'python', 'OpenFUSIONToolkit')):
        parser.error('No OpenFUSIONToolkit Python package in {0}/python (build first)'.format(build_dir))
    kernel = args.kernel
    if kernel is None:
        kernel = 'python3'
        cache_file = os.path.join(build_dir, 'CMakeCache.txt')
        if os.path.isfile(cache_file):
            with open(cache_file, 'r') as fid:
                for line in fid:
                    if line.startswith('OFT_PY_KERNEL:'):
                        kernel = line.split('=', 1)[1].strip()
    kernel_overrides = []
    for entry in args.kernel_for:
        if '=' not in entry:
            parser.error('--kernel-for expects GLOB=KERNEL, got "{0}"'.format(entry))
        kernel_overrides.append(tuple(entry.rsplit('=', 1)))
    # Kernels (including the --check-kernel one) inherit this environment, as in test_examples.sh
    os.environ['OFT_ROOTPATH'] = build_dir
    os.environ['OMP_NUM_THREADS'] = str(args.threads)
    os.environ['PATH'] = os.pathsep.join([os.path.join(build_dir, 'python'), os.path.join(build_dir, 'bin'), os.environ.get('PATH', '')])
    kernel_names = interpreter_kernels([kernel] + [override for _, override in kernel_overrides])
    if args.check_kernel:
        status = 0
        for value in [kernel] + sorted(set(override for _, override in kernel_overrides)):
            status = max(status, check_kernel(kernel_names[value], label=value))
        return status
    if args.work_dir is None:
        work_root = os.path.join(tempfile.mkdtemp(prefix='oftnb_'), 'examples')
    else:
        work_root = os.path.abspath(args.work_dir)
        if os.path.exists(work_root):
            parser.error('Work directory {0} already exists'.format(work_root))
    longest = max(
        (os.path.relpath(os.path.join(root, name), examples_root) for root, _, names in os.walk(examples_root) for name in names),
        key=len
    )
    if len(os.path.join(work_root, longest)) > OFT_PATH_LIMIT - PATH_MARGIN:
        parser.error('Work directory {0} is too long: {1} would exceed OFT\'s {2}-character path limit '
                     '(less a {3}-character margin); use a shorter --work-dir'.format(
                         work_root, os.path.join(work_root, longest), OFT_PATH_LIMIT, PATH_MARGIN))
    shutil.copytree(examples_root, work_root)

    opts = {
        'examples_root': examples_root,
        'work_root': work_root,
        'kernel': kernel,
        'kernel_overrides': kernel_overrides,
        'kernel_names': kernel_names,
        'timeout': args.timeout,
        'update': args.update,
        'source_hashes': {
            rel_path: file_hash(os.path.join(examples_root, rel_path)) for _, rel_paths in plan for rel_path in rel_paths
        },
    }
    nnotebooks = sum(len(rel_paths) for _, rel_paths in plan)
    print('Running {0} notebooks in {1} directories (kernel "{2}", work dir {3})'.format(
        nnotebooks, len(plan), kernel, work_root), flush=True)
    results = []
    modified_data = []
    start = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(run_directory, dirname, rel_paths, opts): dirname for dirname, rel_paths in plan}
        for future in concurrent.futures.as_completed(futures):
            try:
                dir_results, dir_modified = future.result()
            except Exception as exc:
                dir_results = [{'notebook': rel_path, 'status': 'failed', 'seconds': 0.0, 'cell': None,
                                'error': 'Runner error: {0}: {1}'.format(type(exc).__name__, exc)}
                               for dirname, rel_paths in plan if dirname == futures[future] for rel_path in rel_paths]
                dir_modified = []
            for result in dir_results:
                line = '{0:<14} {1} ({2:.0f} s)'.format(result['status'], result['notebook'], result['seconds'])
                if result['status'] == 'failed':
                    line += ' cell {0}: {1}'.format(result.get('cell'), result.get('error'))
                print(line, flush=True)
            results += dir_results
            modified_data += dir_modified

    counts = {}
    for result in results:
        counts[result['status']] = counts.get(result['status'], 0) + 1
    print('SUMMARY {0} in {1:.0f} s'.format(
        ', '.join('{0} {1}'.format(value, key) for key, value in sorted(counts.items())), time.time() - start))
    if modified_data:
        print('REGENERATED DATA (not copied back): {0}'.format(', '.join(sorted(modified_data))))
    if args.summary:
        with open(args.summary, 'w') as fid:
            json.dump({'kernel': kernel, 'kernel_overrides': kernel_overrides, 'build_dir': build_dir, 'work_dir': work_root, 'counts': counts,
                       'results': sorted(results, key=lambda r: r['notebook']),
                       'regenerated_data': sorted(modified_data)}, fid, indent=2)
    return 1 if counts.get('failed', 0) > 0 else 0


if __name__ == '__main__':
    sys.exit(main())
