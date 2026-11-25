**Title:** Linear Algebra: Native overlapping (additive) Schwarz preconditioner with graph and cell overlap

**Branch:** `hansec/OpenFUSIONToolkit:claude/native-la-schwarz` (c439741) → `OpenFUSIONToolkit/OpenFUSIONToolkit:main`

> Depends on #474 (`claude/native-la-core`). This branch contains #474; please review it as the diff against `claude/native-la-core`, which holds only the Schwarz work.

## Primary changes

Native overlapping Schwarz preconditioner (`oft_schwarz_precond`, XML `<pre type="add_schwarz">` without PETSc, which previously aborted), with subdomains that overlap across MPI processes.

- **Overlapping subdomains across processes** (`src/lin_alg/native_overlap.F90`, new module `oft_native_overlap_mod`, type `oft_native_overlap`): builds on each process the matrix of a subdomain that may include rows held only by other processes, with persistent scatter plans for numeric updates and for vector gather/write-back. Two entry points: `native_overlap_setup` (owned rows grown by `levels` layers of the global matrix graph) and `native_overlap_setup_dofs` (a given DoF list, one request round). Rows are completed at their owners from the unassembled local matrices (`Mfull`), then sent to the other holders. Periodic meshes are supported (revolved periodicity is untested): the periodic images of a DoF on one process form one subdomain DoF, represented by the owned copy on the owner and the lowest local index elsewhere. Write-back is restricted (owned rows, then the caller's `stitch(0)`), additive (sum over subdomains with `stitch(1)`) or weighted additive (partition of unity, `native_overlap_pu_setup`).
- **Schwarz preconditioner** (`src/lin_alg/native_solvers.F90`, `oft_schwarz_precond` extending `oft_bjprecond`):
  - graph mode, `nlocal=1`: one subdomain per process grown across processes (`overlap=0` is block-Jacobi bit for bit);
  - graph mode, `nlocal>1`/field groups: blocks grown within the owned rows of a process (restricted `overlap=1` is block-Jacobi with `loverlap` bit for bit);
  - cell mode: blocks grown by layers of mesh cells (point/edge/face adjacency) by a cell overlap builder attached to the matrix; the union of a process' blocks is one cross-process subdomain and each block a submatrix of it;
  - restricted, additive and weighted (partition of unity) combinations.
  Block loops follow the block-Jacobi rules of #450: a parallel region only for thread-safe sub-solvers outside active parallel regions (so with #465's STRUMPACK default the blocks are solved serially), and a copy of the initial guess where restricted blocks overlap.
- **Cell-adjacency overlap builder** (`src/fem/fem_cell_overlap.F90`, new module `fem_cell_overlap`, type `oft_fem_cell_overlap`; abstract interface `oft_cell_overlap_builder` and `oft_overlap_dofs` in `lin_alg_base`, attached to a matrix as `oft_matrix%cell_ovl`): partitions the local cells (same METIS partition and DoF coloring as `fem_partition`/`comp_fem_partition`), grows blocks over a halo of cells from neighbouring processes found through the base mesh (`src/grid/mesh_overlap.F90`, new module `oft_mesh_overlap`), including across periodic seams, and returns the DoFs of each block and the DoFs assigned to each block.
- **MUG**: `xmhd` and `xmhd_lag` attach a cell overlap builder to the native Jacobian on every multigrid level; `xmhd_nparts` sets its blocks. For native matrices the block-Jacobi coloring now comes from the builder, which calls the same partitioning routine as before (identical coloring); other back ends keep `comp_fem_partition`.

## Secondary changes

- `oft_scatter`: index sets (`oft_index_set`, `index_set_create/delete/grow`) restored from the core branch's history (they were removed there as unused), plus the collective `oft_exchange_i8` and the debug reporter `oft_exchange_report` (setup time and message volume with `debug>=1`).
- `native_la`: submatrices by arbitrary row and column sets, with the symbolic extraction cached by the parent pattern (`native_extraction`, `native_extraction_get`); block-Jacobi uses the same path.
- `fem_utils`: `fem_partition_cells`/`comp_fem_partition_cells` (the existing partitions, also returning the cell partition); `mesh_local_partition` optionally returns the cell partition.
- Docs: `doc_solvers.md` documents every Schwarz option.
- Tests (new): `lin_alg/test_native_schwarz` (block-Jacobi equivalences, subdomain matrices against the independently assembled global matrix, DoF-list subdomains, updates, iteration trends; plain, x-periodic and fully periodic cubes on 1/2/4 processes), `lin_alg/test_cell_overlap` (builder against brute-force layering), `lin_alg/test_cell_schwarz` (cell-mode exactness, trends, updates; plain and x-periodic), `physics/MUG/test_MUG_schwarz.py` (cell-mode Schwarz as the MG smoother reproducing the block-Jacobi/reference errors, including the triply periodic sound wave and geometric multigrid), new cases in `lin_alg/test_solver_xml.py`.

## API changes

No existing input changes behaviour: block-Jacobi inputs are unchanged (same code path, identical iteration counts in the suite), and `xmhd_nparts` produces the same coloring.

**New XML options** under `<pre type="add_schwarz">` with the native backend (previously an abort without PETSc):

| Option | Values |
|---|---|
| `nlocal` | as for block-Jacobi (graph mode; replaced by the builder's block count in cell mode) |
| `groups` | as for block-Jacobi (graph mode) |
| `boverlap` | as for block-Jacobi (graph mode) |
| `overlap` | layers of overlap, default 1, per level list allowed |
| `type` | `restricted` (default), `additive`, `weighted` |
| `overlap_mode` | `graph` (default), `cell` |
| `adjacency` | `point` (default), `edge`, `face` (cell mode) |
| `base_layers` | base mesh layers bounding the halo in cell mode, default -1 (= `overlap`), per level list allowed |

`loverlap` is rejected for `add_schwarz` (use `overlap`). The nested sub-solver is given as for block-Jacobi.

**New Fortran API:** modules `oft_native_overlap_mod`, `fem_cell_overlap`, `oft_mesh_overlap`; types `oft_schwarz_precond`, `oft_native_overlap`, `oft_fem_cell_overlap`, `oft_cell_overlap_builder` (abstract), `oft_overlap_dofs`, `oft_cell_gid_set`, `native_extraction`, `oft_index_set`; component `oft_matrix%cell_ovl`; constants `OFT_CELL_ADJ_POINT/EDGE/FACE`; procedures `native_overlap_setup`, `native_overlap_setup_dofs`, `native_overlap_update`, `native_overlap_gather`, `native_overlap_scatter`, `native_overlap_pu_setup`, `native_overlap_index`, `native_overlap_delete`, `fem_cell_overlap_create`, `mesh_has_transfer_level`, `mesh_base_gid_inverse`, `mesh_base_candidates`, `mesh_cell_base_index`, `native_extraction_get`, `fem_partition_cells`, `comp_fem_partition_cells`, `index_set_create/delete/grow`, `oft_exchange_i8`, `oft_exchange_report`; new optional arguments `rows`/`cols` of the native submatrix setup and `cpart_out` of `mesh_local_partition`.

## Test results

All runs on Ubuntu 24.04, GCC 14, OpenMP + Open MPI, OpenBLAS/UMFPACK/SuperLU/ARPACK/HDF5/METIS from system packages, 4 cores, `OMP_NUM_THREADS=2`. Builds at 80043fa; the final commit c439741 changes only a test comment and one sentence of `doc_solvers.md`.

| Run | Result |
|---|---|
| #474 baseline (`claude/native-la-core` 1edb7ba), Release, `make test` | 3 failed, 1173 passed, 394 skipped, 229 deselected |
| This branch, Release, `make test` | 3 failed, 1218 passed, 394 skipped, 232 deselected |
| This branch, checked Release (`-DCMAKE_Fortran_FLAGS=-DOFT_DEBUG`, `native_check_shared` active), `make test` | 3 failed, 1218 passed, 394 skipped, 232 deselected; no shared-value check failures |
| This branch, checked build, slow Schwarz cases (`test_native_schwarz` periodic mpi4, `test_cell_schwarz` mpi4_refined and periodic mpi2_refined) | 3 passed |
| This branch, checked build, `physics/MUG/test_alfven.py -m slow` | 23 passed, 12 skipped |
| This branch, checked build, `physics/MUG/test_alfven_lag.py -m slow` | 23 passed, 12 skipped |
| This branch, checked build, `physics/MUG/test_sound.py -m slow` | 32 passed, 24 skipped |
| This branch, checked build, `physics/MUG/test_sound_lag.py -m slow` | 32 passed, 25 skipped |

The MUG skips are PETSc variants (no PETSc in this build). All `test_MUG_schwarz.py` cases (triply periodic sound wave, walls, geometric multigrid) run in the default suite and pass in both builds. The three failures are the same on the baseline and on both builds of this branch and are not related to this PR: `lin_alg/test_native_bjacobi.py::test_ILU_mpi4_nlocal2[0]` and `test_ILU_mpi4_nlocal2_part[0]` (iteration count) and `physics/test_tracing.py::test_poincare[1-True-4]` (non-zero exit). The 45 additional passes are the new tests, the 3 additional deselected cases the new slow ones. Lints: `generate_stack.py -l`, `check_ascii.py` (no new warnings) and `ruff check` are clean.

### Iteration counts

GMRES iterations in `lin_alg/test_native_schwarz` (graph mode). Columns: one subdomain per process with overlap 0 (= block-Jacobi), restricted overlap 1 and 2, additive overlap 1 and 2, weighted overlap 1 and 2; four blocks per process with no overlap, plain additive overlap 1, restricted overlap 1, weighted overlap 1. Overlap 0 reproduces block-Jacobi bit for bit (checked in the test).

| Cube | Procs | BJ | R1 | R2 | A1 | A2 | W1 | W2 | 4 blk BJ | 4 blk A1 | 4 blk R1 | 4 blk W1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plain | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 55 | 28 | 15 | 16 |
| plain | 2 | 49 | 12 | 6 | 12 | 6 | 12 | 6 | 62 | 71 | 50 | 50 |
| plain | 4 | 54 | 16 | 9 | 28 | 20 | 17 | 13 | 68 | 88 | 58 | 59 |
| x-periodic | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 92 | 38 | 24 | 27 |
| x-periodic | 2 | 71 | 21 | 11 | 20 | 11 | 21 | 11 | 96 | 99 | 72 | 72 |
| x-periodic | 4 | 100 | 24 | 14 | 37 | 26 | 27 | 18 | 126 | 140 | 99 | 103 |

(On one process the single subdomain is the whole matrix, solved exactly.) All columns are identical to the follow-up tip (12f539d + ea04bdb) on every case, including the 4-process cases; the first eight columns are identical to the cell-adjacency tip (9eff82f, no weighted mode) on every case, including the 4-process x-periodic case rerun on that tip. The fully periodic cube only checks subdomain matrices (singular Laplacian).

Cell mode, `lin_alg/test_cell_schwarz` (one block per process / three blocks per process). Columns: layers to cover the mesh; point, edge and face adjacency with 0/1/2 cell layers; graph overlap 0/1/2; weighted additive with one point layer (weighted with covering layers converges in 1 iteration in every case).

| Case (procs, levels, order) | Blocks | Cover | Point 0/1/2 | Edge 0/1/2 | Face 0/1/2 | Graph 0/1/2 | Weighted 1 |
|---|---|---|---|---|---|---|---|
| plain (1,3,3) | 1 | 0 | 1/1/1 | 1/1/1 | 1/1/1 | 1/1/1 | 1 |
| plain (1,3,3) | 3 | 3 | 14/7/4 | 14/8/5 | 14/10/8 | 40/13/7 | 11 |
| plain (2,4,3) | 1 | 3 | 21/7/4 | 21/9/6 | 21/12/9 | 39/10/5 | 7 |
| plain (2,4,3) | 3 | 4 | 22/9/6 | 22/10/8 | 22/13/11 | 47/39/39 | 13 |
| plain (4,4,2) | 1 | 4 | 19/8/5 | 19/10/7 | 19/13/10 | 34/11/7 | 11 |
| plain (4,4,2) | 3 | 5 | 19/9/6 | 19/11/8 | 19/13/11 | 39/35/34 | 12 |
| plain (2,5,2) | 1 | 6 | 25/10/7 | 25/12/8 | 25/14/12 | 46/13/8 | 11 |
| plain (2,5,2) | 3 | 8 | 25/12/9 | 25/14/11 | 25/16/14 | 56/47/48 | 15 |
| plain (4,5,2) slow | 1 | 8 | 31/14/9 | 31/15/11 | 31/18/15 | 52/17/10 | 15 |
| plain (4,5,2) slow | 3 | 10 | 33/14/10 | 33/15/12 | 33/19/16 | 62/56/56 | 18 |
| x-periodic (1,3,3) | 3 | 4 | 35/14/11 | 35/16/11 | 35/19/16 | 64/21/12 | 18 |
| x-periodic (2,4,3) | 1 | 3 | 41/15/9 | 41/18/11 | 41/22/18 | 58/17/9 | 15 |
| x-periodic (2,4,3) | 3 | 4 | 43/16/11 | 43/19/14 | 43/25/19 | 76/61/55 | 22 |
| x-periodic (4,4,2) | 1 | 4 | 34/15/9 | 34/17/12 | 34/20/17 | 57/17/11 | 17 |
| x-periodic (4,4,2) | 3 | 5 | 38/16/11 | 38/18/14 | 38/23/18 | 66/59/57 | 20 |
| x-periodic (2,5,2) slow | 1 | 5 | 46/19/13 | 46/23/17 | 46/27/23 | 60/22/14 | 20 |
| x-periodic (2,5,2) slow | 3 | 9 | 52/22/14 | 52/26/18 | 52/31/26 | 82/66/63 | 24 |

All non-weighted columns (and the update/re-setup columns not shown, which equal the fresh setup) are identical to the cell-adjacency tip 9eff82f in every case. The weighted columns are new (cell mode had no weighted combination before). The x-periodic serial single-block row is all ones as for the plain cube.

`lin_alg/test_solver_xml` (`add_schwarz` through XML): restricted/additive/weighted 58/102/59 iterations serial, 30/42/31 on 2 processes, 58/100/60 in cell mode. MUG: `test_MUG_schwarz.py` (cell-mode Schwarz as the multigrid smoother) reproduces the reference errors of the block-Jacobi runs, including the triply periodic sound wave of `test_sound.py`/`test_sound_lag.py` and geometric multigrid on 2 processes.

### Known limitations

- Revolved periodic domains have not been tested with the Schwarz preconditioner (noted in `doc_solvers.md`).
- `test_cell_schwarz` uses a cube periodic in x only: periodic in x and y, its coarse cube aborts during mesh setup (`mesh_global_igrnd`, "Ground point repeated on same domain"), before any Schwarz code runs. Full periodicity is covered by `test_native_schwarz` (fully periodic cube, subdomain matrices) and `test_MUG_schwarz` (triply periodic sound wave, cell mode).

### Timing

Cell-mode setup against the cell-adjacency tip (2 processes, 4 levels, order 3, `debug=1`, medians of 2x3 alternating runs): sum of setup phases 0.528 s vs 0.498 s with one block per process, 0.571 s vs 0.561 s with three; graph-mode setup 0.0360 s vs 0.0359 s. Phases doing identical work vary by a similar amount between runs (for example 0.45 vs 0.415 s and 1.27 vs 1.175 s), so there is no measurable regression.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01NKQvsASycSEqwvD4f9F4eE
