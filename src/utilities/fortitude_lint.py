#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Lint Fortran sources with fortitude (https://github.com/PlasmaFAIR/fortitude).
#
# Fortitude parses Fortran without running the C preprocessor, so a few OFT
# preprocessor idioms that are valid after preprocessing appear as syntax
# errors, which also cause fortitude to discard most other diagnostics for the
# affected file. To avoid this, each file is copied to a temporary mirror of
# the source tree with the following line-preserving substitutions, so that
# reported line numbers still match the original files:
#
#   1. Lines containing only a debug stack macro (DEBUG_STACK_PUSH or
#      DEBUG_STACK_POP, from include/local.h) are blanked.
#   2. Preprocessor macros used as widths in FORMAT statements
#      (eg. "FORMAT (I XMHD_RST_LEN.XMHD_RST_LEN)") are replaced with a
#      literal width.
#   3. Only for files that still fail to parse: preprocessor directive lines
#      ("#if", "#else", ...) are blanked, so both branches are linted. This is
#      not done by default, as some files rely on the directives to be valid
#      (eg. alternate program bodies in "#if"/"#else" branches).
#
# Fortitude is then run once on the mirror, using the configuration in
# "fortitude.toml" in the source directory, and any temporary paths in its
# output are mapped back to the original files. Any extra arguments are passed
# to "fortitude check" (eg. "--output-format github" or "--statistics"). Note
# that source snippets shown in the output are from the substituted copies.
#
# Usage (from anywhere in the repository):
#   python src/utilities/fortitude_lint.py [paths ...] [-- fortitude args ...]
#
# Paths default to the whole source directory. Automatic fixes ("--fix") are
# not supported, as they would be applied to the temporary copies; run
# "fortitude check --fix <file>" directly instead.
#
#------------------------------------------------------------------------------
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

# Source directory (containing fortitude.toml), relative to this script
SRC_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
CONFIG_FILE = 'fortitude.toml'
# Lines consisting only of a debug stack macro
STACK_MACRO_RE = re.compile(r'^[ \t]*DEBUG_STACK_(PUSH|POP)[ \t]*$', re.MULTILINE)
# FORMAT statements and macros used within them (upper case names with underscores)
FORMAT_RE = re.compile(r'^([ \t]*[0-9]+[ \t]+FORMAT\b)(.*)$', re.MULTILINE | re.IGNORECASE)
FORMAT_MACRO_RE = re.compile(r'(?:\b[A-Z]\s+)?\b[A-Z][A-Z0-9]*_[A-Z0-9_]+\b')
# Preprocessor directive lines
DIRECTIVE_RE = re.compile(r'^[ \t]*#.*$', re.MULTILINE)
# Rules indicating the file could not be parsed
PARSE_ERROR_RULES = ('E000', 'E001')


def fortitude_cmd():
    '''Return the command used to run fortitude'''
    exe = shutil.which('fortitude')
    if exe is None:
        print('Error: fortitude not found, install with "pip install fortitude-lint"')
        sys.exit(2)
    return [exe]


def replace_format_macro(match):
    '''Replace a macro in a FORMAT statement with a literal width'''
    text = match.group(0)
    descriptor = text.split()[0] if (len(text.split()) > 1) else ''
    return descriptor + '5'


def preprocess(text, strip_directives=False):
    '''Apply line-preserving substitutions to make OFT sources parseable'''
    text = STACK_MACRO_RE.sub('', text)
    text = FORMAT_RE.sub(lambda m: m.group(1) + FORMAT_MACRO_RE.sub(replace_format_macro, m.group(2)), text)
    if strip_directives:
        text = DIRECTIVE_RE.sub('', text)
    return text


def write_mirror(mirror_dir, files, strip_directives=False):
    '''Write preprocessed copies of source files (relative to SRC_DIR) to the mirror'''
    for rel_path in files:
        with open(os.path.join(SRC_DIR, rel_path), 'r', encoding='utf-8', errors='surrogateescape') as fid:
            text = fid.read()
        out_path = os.path.join(mirror_dir, rel_path)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, 'w', encoding='utf-8', errors='surrogateescape', newline='') as fid:
            fid.write(preprocess(text, strip_directives))


def list_files(cmd, paths):
    '''Get files to check (relative to SRC_DIR), honoring fortitude.toml file selection'''
    result = subprocess.run(cmd + ['check', '--show-files'] + paths, cwd=SRC_DIR,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if result.returncode != 0:
        print(result.stderr, end='')
        sys.exit(result.returncode)
    files = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if line == '':
            continue
        rel_path = os.path.relpath(os.path.abspath(os.path.join(SRC_DIR, line)), SRC_DIR)
        if rel_path.startswith('..'):
            continue
        files.append(rel_path)
    return sorted(set(files))


def parse_failures(cmd, mirror_dir, files):
    '''Return the set of files in the mirror that fortitude cannot parse'''
    result = subprocess.run(cmd + ['check', '--select', ','.join(PARSE_ERROR_RULES), '--exit-zero',
                                   '--quiet', '--output-format', 'json'] + files,
                            cwd=mirror_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    failed = set()
    for diag in json.loads(result.stdout or '[]'):
        failed.add(os.path.relpath(os.path.abspath(os.path.join(mirror_dir, diag['filename'])), mirror_dir))
    return failed


def main(argv):
    if '--' in argv:
        split = argv.index('--')
        paths, extra_args = argv[:split], argv[split+1:]
    else:
        paths, extra_args = argv, []
    if any(arg.startswith('--fix') or arg in ('--diff', '--unsafe-fixes') for arg in extra_args):
        print('Error: automatic fixes are not supported, run "fortitude check --fix" directly')
        return 2
    # Paths are given relative to the current directory
    paths = [os.path.abspath(path) for path in paths] or [SRC_DIR]
    cmd = fortitude_cmd()
    files = list_files(cmd, paths)
    if len(files) == 0:
        print('No Fortran files to check')
        return 0
    mirror_dir = tempfile.mkdtemp(prefix='oft_fortitude_')
    try:
        shutil.copy(os.path.join(SRC_DIR, CONFIG_FILE), os.path.join(mirror_dir, CONFIG_FILE))
        write_mirror(mirror_dir, files)
        failed = parse_failures(cmd, mirror_dir, files)
        if len(failed) > 0:
            # Retry without preprocessor directives, keeping whichever version parses
            write_mirror(mirror_dir, sorted(failed), strip_directives=True)
            still_failed = parse_failures(cmd, mirror_dir, sorted(failed))
            write_mirror(mirror_dir, sorted(still_failed))
        # Lint mirrored files and map temporary paths back to the source tree
        result = subprocess.run(cmd + ['check'] + extra_args + files, cwd=mirror_dir,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        for stream, out in ((sys.stdout, result.stdout), (sys.stderr, result.stderr)):
            out = out.replace(os.path.realpath(mirror_dir), SRC_DIR).replace(mirror_dir, SRC_DIR)
            stream.write(out)
            stream.flush()
        return result.returncode
    finally:
        shutil.rmtree(mirror_dir, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
