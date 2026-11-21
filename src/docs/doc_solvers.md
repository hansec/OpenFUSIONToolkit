Available Solvers     {#doc_solvers}
=================

[TOC]

This page outlines the solvers that are currently available through Open FUSION Toolkit solver objects. Additional
solvers may be used through backend implementations, but only internally controlled and supported solvers
are outlined here.

\section doc_solvers_petsc Mixing Native and PETSc Solvers

Solver interfaces have been designed to allow mixing of native and PETSc methods. Most native solvers are
written with generic API calls for linear algebra and as a result can be used when either backend is used
for vectors and matrices. There are many cases however where the PETSc implementation of a particular
solver/preconditioner is more optimal for its matrices and vectors. Additionally, there are some
methods, ex. direct methods, which are only available through PETSc. The
\ref oft_petsc_solvers::oft_petsc_pre_solver "oft_petsc_pre_solver" class is provided to allow using PETSc
solvers as preconditioners for native solvers. This class acts as a wrapper around the PETSc solver, allowing
it to be called using the standard solver calling sequence.  

\section doc_solvers_basic Basic Solvers and Preconditioners

\subsection doc_solvers_basic_jacobi Jacobi Iteration

Jacobi iterations solves a system of linear equations by iterating on the method

\f[ x^{n+1} = x^n + D^{-1}(A*x^n - y) \f]

where \f$ D \f$ is a matrix of only the diagonal elements of \f$ A \f$. Functionality is provided through
either a general native interface \ref oft_solvers::oft_diag_scale "oft_diag_scale" and an optimized PETSc
interface \ref oft_petsc_solvers::oft_petsc_diagprecond "oft_petsc_diagprecond".

**XML Spec:** `<pre type="jacobi">`

\subsection doc_solvers_basic_lu LU Factorization

This method is supported through both the native and PETSc backends through
\ref oft_lu::oft_lusolver "oft_lusolver" and
\ref oft_petsc_solvers::oft_petsc_direct_solver "oft_petsc_direct_solver" for a full solver or
\ref oft_petsc_solvers::oft_petsc_luprecond "oft_petsc_luprecond" as a PETSc preconditioner.

\warning The native solver currently only supports local solves (ie as a sub-solver to Block Jacobi)

**XML Spec:** `<solver type="lu">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `package`  | Solver package to use {super, superd, mumps} | str | No |

PETSc only options

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `type`     | Inverse type {lu, ilu} | str | No |

\subsection doc_solvers_basic_ilu Incomplete LU Factorization

This method is supported through the native backend through \ref oft_lu::oft_ilusolver "oft_ilusolver"
(with the PETSc backend the same XML selects PETSc's ILU instead, and the options below do not apply).
Three native packages are available:
- `native`: ILU(0) (no fill outside the non-zero pattern of the matrix) with native kernels, with an
  OpenMP threaded (level-scheduled) solve; the default unless OFT is built with MKL.
- `mkl`: ILU(0) from MKL; the default when OFT is built with MKL.
- `super`: threshold-based ILU with partial pivoting (ILUTP) from SuperLU (`dgsitrf`, as used by its
  `dgsisx` driver, without the driver's equilibration and MC64 row permutation). Requires OFT to be
  built with SuperLU. Fill is controlled by a drop tolerance and a fill bound, so the factors are
  more accurate (and more expensive) than ILU(0).

The symbolic setup is redone automatically when the non-zero pattern or size of the matrix changes,
and the factorization after the solver is updated. For `super`, a refactorization with an unchanged
pattern reuses the column permutation and elimination tree, which gives the same factors as a new
factorization.

The `super` defaults keep a moderate amount of fill. On the Lagrange Laplacian (quadratic and cubic
elements, 1-30k rows, one process) `drop_tol=1.E-3` gives factors with 3.5-5 times the non-zeros
of the matrix and 4-7 times fewer GMRES iterations than ILU(0), while the factorization is 10-50
times slower than ILU(0) and each application 3-5 times slower. Smaller drop tolerances need a
larger `fill_factor`: once the fill bound is reached entries are dropped by size alone, which can
weaken the preconditioner considerably (`drop_tol=1.E-4` with the default bound needed more
iterations than ILU(0) on the larger cases). The MILU variants needed more iterations than standard
ILU on these operators. ILUTP therefore pays off when a factorization is reused for many solves, or
when ILU(0) converges slowly or not at all. The `native` options (`pivot_tol`, `diag_shift`,
`threaded`) do not apply to `mkl` or `super`, and the `super` options do not apply to the ILU(0)
packages.

\warning The native solver currently only supports local solves (ie as a sub-solver to Block Jacobi or
Additive Schwarz)

**XML Spec:** `<pre type="ilu">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `package`  | Solver package to use {native, mkl, super} | str | No |

`native` only options

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `pivot_tol`  | Relative small-pivot threshold: pivots smaller in magnitude than `pivot_tol` times the largest magnitude in their row are replaced (default: 1.E-12) | real | No |
| `diag_shift` | Relative diagonal shift \f$ \alpha \f$: factor \f$ A + \alpha \mathrm{diag}(A) \f$ (default: 0) | real | No |
| `threaded`   | Use the level-scheduled threaded solve (default: true) | bool | No |

`super` only options

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `drop_tol`     | Drop tolerance (SuperLU `ILU_DropTol`): factor entries smaller than this in relative magnitude are dropped (default: 1.E-3; SuperLU default: 1.E-4) | real | No |
| `fill_factor`  | Bound on the fill of the factors relative to the matrix (SuperLU `ILU_FillFactor`, default: 10) | real | No |
| `milu`         | Modified ILU variant (SuperLU `ILU_MILU`): 0 standard ILU, 1-3 SMILU_1-3 (dropped entries added to the diagonal; default: 0) | int | No |
| `col_perm`     | Column permutation: -1 default (COLAMD), 0 natural, 1 minimum degree on \f$ A^T A \f$, 2 minimum degree on \f$ A^T + A \f$, 3 COLAMD (default: -1) | int | No |
| `fill_tol`     | Zero pivot replacement threshold (SuperLU `ILU_FillTol`, default: 1.E-2) | real | No |
| `pivot_thresh` | Diagonal pivoting threshold (SuperLU `DiagPivotThresh`): 0 always keeps the diagonal pivot when it is non-zero, 1 is partial pivoting (default: 0.1) | real | No |

\subsection doc_solvers_basic_bjacobi Block Jacobi

This method is supported through both the native and PETSc backends through
\ref oft_solvers::oft_bjprecond "oft_bjprecond" and
\ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond".

**XML Spec:** `<pre type="block_jacobi">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `nlocal`   | Number of local blocks (-1 for field based) | int | Yes |
| `groups`   | Group index for each field in multi-field vector [n] | int | Yes |
| `boverlap` | Overlap boundary rows (native only) | bool | false |

\subsection doc_solvers_basic_schwarz Additive Schwarz

This method is provided as a preconditioner through both the native and PETSc backends through
\ref oft_native_solvers::oft_schwarz_precond "oft_schwarz_precond" and
\ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond". Subdomains are grown from a
non-overlapping partition by `overlap` layers along the matrix graph and solved with the
sub-solver given by the nested `pre` element. In the native backend, with `nlocal=1` each process
has one subdomain made of its owned rows, grown across processes; with several local blocks
(`nlocal>1`, or negative for field based blocks) the blocks are grown within the owned rows of
each process. `overlap=0` is block-Jacobi. Periodic meshes are supported: the periodic images of a
DOF held by one process form a single subdomain DOF.

Subdomain solutions are combined (`type`):
- `restricted` (default): each row is taken from the subdomain it was assigned to;
- `additive`: sum of the subdomain solutions (rows in several subdomains are counted once per
  subdomain);
- `weighted`: additive with partition of unity weights, where each subdomain's solution is scaled
  on every row by one over the number of subdomains containing that row (native only). Unlike
  `additive`, rows in the overlap are not over-weighted.

**XML Spec:** `<pre type="add_schwarz">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `nlocal`   | Number of local blocks (-1 for field based) | int | Yes |
| `groups`   | Group index for each field in multi-field vector [n] (native only) | int | Yes |
| `overlap`  | Overlap size (native default: 1) | int | Yes |
| `type`     | Combination of subdomain solutions: `restricted` (default), `additive` or `weighted` (native only) | str | No |

\section doc_solvers_krylov Krylov Solvers

\subsection doc_solvers_krylov_cg Conjugate Gradient

This method is provided through both the native and PETSc backends through \ref oft_cg::oft_cg_solver
"oft_cg_solver" and \ref oft_petsc_solvers::oft_petsc_cg_solver "oft_petsc_cg_solver".

**XML Spec:** `<solver type="cg">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `atol`     | Absolute convergence tolerance | real | Yes |
| `rtol`     | Relative convergence tolerance | real | Yes |
| `its`      | Iteration limit | int | Yes |

\subsection doc_solvers_krylov_gmres Generalized Minimum Residual Method

This method is provided through both the native and PETSc backends through \ref oft_gmres::oft_gmres_solver
"oft_gmres_solver" and \ref oft_petsc_solvers::oft_petsc_gmres_solver "oft_petsc_gmres_solver".

**XML Spec:** `<solver type="gmres">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `atol`     | Absolute convergence tolerance | real | Yes |
| `rtol`     | Relative convergence tolerance | real | Yes |
| `its`      | Iteration limit | int | Yes |
| `nrit`     | Restart size | int | Yes |
