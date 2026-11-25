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

By default the native solver performs local solves (ie as a sub-solver to Block Jacobi). With the
`superd`, `mumps`, and `strumpack` packages, the full distributed system can instead be factored across
all MPI tasks by setting `distributed` to true and using the solver directly on the global system (not
as a sub-solver to Block Jacobi). This option has no effect in serial runs or for local systems.
Distributed factorization with `strumpack` requires STRUMPACK built with MPI support (`--strumpack_mpi=1`
with `build_libs.py`), in which case the STRUMPACK options below also apply.

**XML Spec:** `<solver type="lu">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `package`  | Solver package to use: native {super, superd, umfpack, mumps, mkl, lapack, strumpack}, PETSc {petsc, super, superd, umfpack, mumps}. Native default: first available of mkl, strumpack, umfpack, mumps, super, superd | str | No |

Native only options

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `distributed` | Factor full distributed system across all MPI tasks (`superd`, `mumps`, and `strumpack` only) | bool | No |

Native STRUMPACK options, in a `<strumpack>` group (values of `-1` or omitted use the STRUMPACK defaults)

| Group          | Description | Type      | List |
|----------------|-------------|-----------|------|
| `comp_type`    | Compression of frontal matrices {none, hss, blr} (default: none) | str | No |
| `comp_rel_tol` | Relative compression tolerance | float | No |
| `comp_abs_tol` | Absolute compression tolerance | float | No |
| `leaf_size`    | Size of leaf/tile blocks for compression | int | No |
| `max_rank`     | Maximum rank of compressed blocks | int | No |
| `min_sep_size` | Minimum separator size for compression | int | No |
| `verbose`      | Print STRUMPACK information | bool | No |

\note With compression the STRUMPACK factorization is approximate, which is intended for use as a
preconditioner (eg. Block Jacobi sub-solver). STRUMPACK's iterative refinement is not supported. A
refactorization with an unchanged non-zero pattern reuses the reordering and MC64 matching/scaling from the
last full factorization (local and distributed solves). STRUMPACK only pivots within each front, so after large
changes in matrix values the solution can be inaccurate unless a full factorization is requested
(`update(.TRUE.)`, as for a change in the non-zero pattern), even though the pattern has not changed.

PETSc only options

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `type`     | Inverse type {lu, ilu} | str | No |

\subsection doc_solvers_basic_bjacobi Block Jacobi

This method is supported through both the native and PETSc backends through
\ref oft_solvers::oft_bjprecond "oft_bjprecond" and
\ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond".

**XML Spec:** `<pre type="block_jacobi">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `nlocal`   | Number of local blocks (-1 for per-field, <-1 for custom field `groups`) | int | Yes |
| `groups`   | Group index (1..ngroups) for each field in multi-field vector [n] | int | Yes |
| `boverlap` | Overlap boundary rows (native only) | bool | false |
| `loverlap` | Overlap internal block boundaries by one graph layer (native only) | bool | false |

\subsection doc_solvers_basic_patch Overlapping patch (vertex-star) smoother

This method is supported with the native backend through
\ref oft_patch_solvers::oft_patch_precond "oft_patch_precond". Each patch is a set of matrix rows that is
solved exactly (dense LU) and the patch corrections are combined additively or multiplicatively. The patches
are supplied with the matrix (\ref oft_la_base::oft_patch_info "oft_patch_info"). MUG (`xmhd` and `xmhd_lag`)
attaches vertex-star patches on every multigrid level (including the sub-block matrices \f$ J_{BB} \f$ and
\f$ S_V \f$ when using `block_schur`), built by
\ref fem_utils::comp_fem_vertex_patches "comp_fem_vertex_patches": all degrees of freedom of all fields on a
vertex and on the edges, faces and cells that contain it. Vertex stars contain the gradients of the vertex
hat functions, so they resolve the near-kernel of the curl operator for H(curl) fields.

On distributed meshes, patches of vertices on process boundaries are assembled and solved in full by one
process, so the additive smoother is independent of the number of processes. In the multiplicative modes the
interior patches are swept first, followed by an additive correction from the process-boundary patches.

When `its` is positive the object performs `its` damped Richardson iterations and can be used directly as a
multigrid smoother. Otherwise it acts as a preconditioner (eg. for a GMRES smoother).

**XML Spec:** `<pre type="patch">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `mode`     | Patch combination {additive, multiplicative, symmetric} (default: additive) | str | No |
| `its`      | Richardson iterations (<=0 to act as preconditioner, default: -1) | int | Yes |
| `df`       | Damping factor (default: 1) | real | Yes |
| `single`   | Store patch factors in single precision (default: false) | bool | No |

\note Memory use is the sum of the squared patch sizes. On coupled B-V systems vertex-star smoothing alone
does not converge at Alfven CFL numbers of order one or larger, where B-V wave coupling spans multiple cells;
in that regime use it inside `<b_pre>` of \ref doc_solvers_basic_schur "block_schur" (which also reduces
patch memory significantly) or combine it with an overlapped `block_jacobi` smoother.

\subsection doc_solvers_basic_schur Physics-Based Block-Schur (xMHD)

For 3D extended MHD (\ref xmhd and \ref xmhd_lag), a physics-based Block-Schur preconditioner is provided
through \ref xmhd::oft_xmhd_schur_pre "oft_xmhd_schur_pre" and \ref xmhd_lag::oft_xmhd_schur_pre
"xmhd_lag::oft_xmhd_schur_pre". It partitions the system into magnetic (`Bc+Bg` in \ref xmhd or
`Bx,By,Bz` in \ref xmhd_lag), velocity (`Vx,Vy,Vz`), and optional scalar (`n,Ti,Te`) blocks and
assembles an SPD parabolized Alfven + acoustic velocity Schur complement \f$ S_V \f$ across
all MG levels so that p-MG / h-MG can be applied cleanly to the sub-blocks.

**XML Spec:** `<xmhd><pre type="block_schur">`

| Group        | Description | Type | Default |
|--------------|-------------|------|---------|
| `schur_type` | Factorization type (`lower`, `upper`, or `ldu`) | str | `ldu` |
| `sv_scale`   | Scale factor on parabolized Alfven/acoustic term in \f$ S_V \f$ | real | `1.0` |
| `b_pre`      | Solver or `<pre type="mg">` definition for magnetic block \f$ J_{BB} \f$ | xml | `gmres(2)+BJ(1,ilu)` |
| `v_pre`      | Solver or `<pre type="mg">` definition for velocity Schur \f$ S_V \f$ | xml | MG on \f$ S_V \f$ |
| `s_pre`      | Solver definition for scalar block \f$ J_{SS} \f$ (when active) | xml | `gmres(2)+BJ(-1,ilu)` |

\subsection doc_solvers_basic_mg Multi-Grid Preconditioner

Multi-grid V-cycle preconditioning is provided through \ref oft_native_solvers::oft_ml_precond
"oft_ml_precond".

**XML Spec:** `<pre type="mg">`

| Option         | Description | Type | Default |
|----------------|-------------|------|---------|
| `report`       | Print per-stage residual diagnostics (0=off, 1=total, 2=per-block) | int (attr) | `0` |
| `coarse_damp`  | Damping factor \f$ \omega_c \f$ applied to coarse-grid correction | real (attr/elem) | `1.0` |
| `coarse_shift` | (xMHD only) Mass-matrix shift \f$ \sigma \f$ on coarse levels (\f$ J_c \leftarrow J_c + \sigma M_c \f$) | real (attr/elem) | `0.0` |

\subsection doc_solvers_basic_schwarz Additive Schwarz

This method is currently only supported with the PETSc backend through
\ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond" as a preconditioner.

**XML Spec:** `<pre type="add_schwarz">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `nlocal`   | Number of local blocks (-1 for field based) | int | Yes |
| `overlap`  | Overlap size | int | Yes |

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
