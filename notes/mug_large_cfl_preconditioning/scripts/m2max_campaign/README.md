# M2 Max benchmark campaign scripts

These scripts produced [`../../RESULTS_M2Max.md`](../../RESULTS_M2Max.md). Run them from an empty working directory
(case directories are created relative to the current directory). Environment:

- `OFT_BUILD_DIR`: build tree with `examples/MUG/Spheromak_tilt/MUG_sph_tilt` and `tests/physics/MUG/test_alfven` (read by `../bench_case.py`)
- `OFT_TREE`: build-tree parent holding `mpich-*/bin/mpirun` and `build_railgun/` (used by `gen_T.py` and `railgun/post.sh`)
- `OFT_PYTHON`: Python with numpy/h5py (default: the interpreter running the script)

| script | purpose |
|---|---|
| `gen_A.py` | A1 (iteration matrix: R, p, CFL, s, lin/cn/cnmf) and A2 (accuracy at T=1e-4) cube cases and lists |
| `gen_A1x.py` | A1x inner-solver comparison on P (exact LU, ~3000-DOF blocks, 20 inner its) |
| `gen_C.py` | C1 linear tilt growth-rate and C2 nonlinear tilt cases |
| `gen_T.py` | timing lane: D restart segments (needs `C2/sph_R0_p2_ilu_s0_dt2e-7` restart files) and B scaling, all 4 ranks x 2 threads |
| `gen_T2.py` | follow-up timing cases: full-Jacobian block-Jacobi LU with the SI term |
| `runq.py` | job queue under a core and memory budget (`--cores`, `--mem`, `--timeout`, `--list`, `--cores_file`). A watchdog kills the newest job if total RSS exceeds `--mem` |
| `governor.sh`, `timing_lane.sh` | core budget for the shared accuracy lane; exclusive timing lane that waits until nothing else runs |
| `tab_A1.py [R]`, `tab_A1x.py`, `tab_A2.py`, `fit_C1.py`, `cmp_C2.py`, `tab_D.py`, `tab_B.py` | result tables (as in `../../results/m2max/`) |
| `railgun/mkcase.sh` | railgun case from a base directory `coax_base_np3` holding the example's `oft_coax.in`/`.xml`, mesh, `plot_railgun.py` and a `run.sh` |
| `railgun/post.sh` | plotting run plus `plot_railgun.py` sheet-speed analysis |
| `railgun/cmp_fields.py`, `railgun/tab_railgun.py` | relative L2 field differences vs a base run; correctness table |

Protocol: use the accuracy lane for iteration counts and accuracy only. Every timing case runs alone with the same
ranks x threads layout, because the machine has mixed performance/efficiency cores.
