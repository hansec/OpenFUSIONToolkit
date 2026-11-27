# Development notes

This directory keeps the final reports of development sequences: the design notes, benchmark results,
validation records and follow-up lists that were produced while developing a feature, but that do not belong in
the code base on `main`. It lives only on the `ai_dev_notes` branch of the fork and is never merged into `main`.

## Layout

One folder per development sequence, named after the pull request that delivered it, or after the branch if
there is no pull request:

```
notes/
  README.md                         this file
  PR_<number>-<Tool>_<short_topic>/  eg. PR_457-MUG_qp_cache
  <branch_name>/                    for work that never became a PR (replace '/' with '_')
```

Use the PR number of the upstream repository (`OpenFUSIONToolkit/OpenFUSIONToolkit`) and the tool prefix of the
PR title (`MUG`, `TokaMaker`, `ThinCurr`, `Marklin`, `FEM`, `LinAlg`, ...). Each folder must contain a
`README.md` and may contain any other Markdown documents and assets.

## What to put in a folder

- `README.md` (required): what the work is (PR link and title, one-paragraph summary), the branches and commit
  SHAs involved (development branches, the PR head before and after any cleanup, the merge commit once merged),
  an index of the other files, and how to reproduce the results.
- Final reports: design and implementation notes, validation results, benchmark tables with the hardware,
  compiler, libraries and options used. Write them so they can be read without the development history, and
  use the option names and defaults of the merged code.
- Follow-ups: open issues, known limitations, ideas for further work, and problems found in passing that were
  outside the scope of the work.
- Assets: small files needed to reproduce or understand the results, eg. benchmark harnesses or drivers that
  were removed from the PR before merging, input files, plots, or short summary tables. Copy them verbatim from
  a named commit and say which one, so they can be run from a checkout of that commit.

Do not store large raw outputs (per-run JSON/CSV dumps, logs, perf data, binaries). Leave them on the
development branches and reference those branches and commits by name. Keep text files ASCII where practical.

## Adding a new folder

1. Check out `ai_dev_notes` (it is based on `main`, but does not need to track it).
2. Add `notes/<folder>/` with at least a `README.md`, following the layout above.
3. Commit with a message such as `Notes: <Tool> <topic> (PR #<number>)` and push to `ai_dev_notes` on the fork.
