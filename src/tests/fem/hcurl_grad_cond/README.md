# MUG B-field conditioning: H(Curl)+Grad(H^1) bases, smoothers and subspace coupling

This directory records a study of the conditioning of MUG's B-field block. The study:
- characterizes the B block;
- tests whether OFT's hierarchical H(Curl)/H^1 bases should be replaced or modified to make the p-multigrid (p-MG) solves converge better;
- identifies which smoother structure the B block actually needs.

It holds this report, the offline analysis scripts, the MUG run drivers, and the raw results. None of it is part of the regression suite. The only compiled piece is the diagnostic driver `../diag_hcurl_grad_cond.F90`, which writes the matrices and block structures that the scripts read.

**No basis or solver sources were changed.** Branch: `claude/pmg-hcurl-basis-conditioning`.

## Contents

- [Summary and recommendation](#summary-and-recommendation)
- [1. Setup](#1-setup)
- [2. Characterization of the B mass block](#2-characterization-of-the-b-mass-block)
- [3. Alternative basis: Schöberl–Zaglmayr prototype](#3-alternative-basis-schöberlzaglmayr-prototype)
- [4. Nature of the bad modes](#4-nature-of-the-bad-modes)
- [5. MUG time-stepping runs](#5-mug-time-stepping-runs)
- [6. Subspace coupling across p-levels and families](#6-subspace-coupling-across-p-levels-and-families)
- [7. Re-based (stabilized) hierarchical bases](#7-re-based-stabilized-hierarchical-bases)
- [8. Vertex-star (Arnold–Falk–Winther) blocks](#8-vertex-star-arnoldfalkwinther-blocks)
- [9. Interpretation](#9-interpretation)
- [10. Caveats](#10-caveats)
- [11. Reproducing](#11-reproducing)

## Summary and recommendation

1. **Do not rewrite or re-base the hierarchical H(Curl)/H^1 bases.**
   - **An alternative basis does nothing.** A Schöberl–Zaglmayr (SZ) basis spans the same space and keeps the exact grad(H^1) split and DOF-injection nesting. It gives the same block-smoother conditioning as OFT's basis (§3). This is expected: a change of basis inside a smoother block cannot change block Jacobi or Schwarz with exact local solves.
   - **Re-basing helps only in a version that cannot be implemented.** Making the new level-p functions orthogonal to the coarse and gradient spaces (a "stabilized" hierarchical basis) cuts entity-block κ by 5–50× only when done globally (§7). That version has dense basis functions. The patch-local, implementable version recovers little at CFL 0, and nothing or worse at CFL ≥ 1.
   - **Re-basing makes the best smoother worse** (§8).
2. **The coupling hypothesis is confirmed, but it is not a basis-quality problem.**
   - The new functions at each p-level are strongly coupled to the p−1 space, and the rotational and gradient families are strongly coupled to each other.
   - Strengthened Cauchy–Schwarz constants γ are 0.98–0.999, and grow with p (§6). The same values hold within a single element, so they are intrinsic to the H(Curl)+grad(H^1) space with this hierarchy.
   - The near-cancelling (gradient vs rotational) combinations are local to **vertex stars**: the patch of all cells around a vertex.
3. **The recommended direction is a vertex-star (Arnold–Falk–Winther) Schwarz smoother for the B block inside p-MG** (§8). On the B operator A = M + τK:

   | quantity | vertex stars | other smoothers |
   |---|---|---|
   | star-preconditioned κ, CFL 0 | 3.7–5.0, at every p from 1 to 4 | vector Lagrange: 4.8–5.8 |
   | damped two-grid ρ (coarse = p−1) | 0.65–0.75, independent of p (2–4), mesh and CFL (0, 1, 10) | 8 entity-aligned blocks, same two-grid: 0.90–0.97; cell-overlap-1 RAS on the refined mesh: 0.93–1.07 |

   - Star blocks are small: mean 80–450 DOFs, total about 3n.
   - The current cell-overlap-1 Schwarz runs out of memory on realistic problems, because its 8 blocks hold about 4–6n DOFs and up to 1.3·10^4 DOFs each.
4. **Interim practical setting.** Until a star smoother exists, keep the phase-2 recommendation for MUG:
   - p-MG with an FGMRES/GMRES smoother wrapping overlapping Schwarz;
   - cell overlap 1 where memory allows, otherwise cell overlap 0 or `loverlap` (§5).

   Plain block Jacobi, METIS-partitioned or entity-aligned, is 3–5× worse.

### Proposed follow-up (not started; needs approval)

1. **Add an `overlap_mode` "vertex" (star) block builder to the native Schwarz preconditioner.**
   - Blocks: one per vertex, holding the DOFs of all entities whose closure contains the vertex (the vertex, plus the edges, faces and cells around it).
   - Local solves: exact, using LU, or Cholesky for this SPD block.
   - Use it as the p-MG smoother for the B block (or the whole system) in MUG.
2. **Measure MUG iteration counts** on the Alfvén and spheromak cases against cell overlap 0/1.
3. **If CFL ≫ 1 is important**, also test a star smoother combined with an exact or auxiliary solve on the gradient subspace (Hiptmair–Xu-type). §7 shows the gradient solve matters little at CFL ≤ 10 once the coarse p−1 solve is present.

## 1. Setup

### The B space

The B space is V_p = H(Curl)_p + grad(H^1_{p+1}), with two fields:
- **Bc**: H(Curl) hierarchical, order p. It starts with the Whitney edge functions.
- **Bg**: H^1 hierarchical, order p+1, which enters through its gradient.

The vertex DOFs of Bg are constrained (removed). Their gradients are already in span(Whitney). So the gradient subspace is

G_p = grad(H^1_{p+1}) = span{Whitney-carried vertex gradients} + span{free Bg DOFs},

with dim G_p = dim H^1_{p+1} − 1. The Whitney-carried vertex gradient is grad λ_v = Σ_e s_ev w_e, where s_ev = −1 at the tail and +1 at the head of each oriented edge.

### Hierarchy

The p-hierarchy is DOF injection, with lower-order functions first within each entity. Validated:
- P^T M_p P = M_{p−1} and P^T K_p P = K_{p−1} to about 1e-16;
- K G_p = 0 to round-off;
- dim G_p matches.

### Operators

- **Mass:** M is the B mass matrix on free DOFs, giving the CFL → 0 limit of the B block.
- **Finite CFL:** A = M + τK, where K is the H(Curl) curl-curl matrix (zero on Bg) and τ = (CFL·h)², h the mean edge length.

### Meshes (cube, tetrahedra)

| name | description | free B DOFs, p = 1 / 2 / 3 / 4 |
|---|---|---|
| tiny | 96 tets, the Alfvén-test cube (`ni=1,2,4`) | 322 / 1143 / 2788 / 5545 |
| ref | 768 tets (`ni=2,2,2`, one uniform refinement) | 2104 / 8052 / 20336 |

Inputs are in `inputs/`.

### Metrics

| name | meaning |
|---|---|
| κ_diag | κ(D⁻¹A), with D the diagonal |
| κ_bj, ρ_bj | block Jacobi with the 8 entity-aligned blocks MUG uses (`xmhd_nparts`). ρ_bj is undamped; it diverges if > 1. |
| κ_ent | block Jacobi with one block per (entity, field) |
| κ_2lev | additive two-level: exact coarse (p−1) solve + κ_bj blocks |
| ρ_2g,bj | symmetric two-grid V(1,1): coarse = exact p−1 solve; smoother = κ_bj blocks damped by 1/λmax |
| ρ_2g,ras | two-grid: exact p−1 solve, then one undamped RAS sweep. RAS uses the 8 cell blocks with one layer of point overlap (`overlap_mode=cell`, `overlap=1`). |
| ρ_2g,ras+G | as ρ_2g,ras, plus an exact solve on G_p |
| κ_star, ρ_2g,star | as κ_bj and ρ_2g,bj, but with vertex-star blocks (§8) |
| xc | max \|A_coarse,new\| / max \|A\|, the residual coarse–new coupling |

## 2. Characterization of the B mass block

Tiny mesh, free DOFs only:

| p | cond(M) | cond(D⁻¹M) | ρ_bj, METIS (8) | ρ_bj, entity-aligned (8) | κ_bj, entity-aligned | κ_bj, vector Lagrange |
|---|---|---|---|---|---|---|
| 1 | 179 | 73 | 1.21 | 1.02 | 25 | 3.2 |
| 2 | 1.43e4 | 706 | 1.85 | 1.80 | 154 | 5.1 |
| 3 | 8.95e5 | 3271 | 1.70 | 1.99 | 526 | 6.6 |

- A Python transcription of OFT's tetrahedral basis (`scripts/basis.py`) reproduces OFT's assembled B mass matrix to relative error ~1e-15 for p = 1–3. This validates the mesh, orientation, numbering and quadrature handling.
- Partition alignment is not the cause. On the refined mesh at p = 2, entity-aligned blocks give ρ_bj 1.89 and κ 82. RAS with one layer of point overlap gives ρ 1.57.

## 3. Alternative basis: Schöberl–Zaglmayr prototype

`scripts/basis.py` also implements an SZ-type basis on the same mesh, orientation and DOF layout. Its span equals OFT's to 15 digits.

| tiny mesh | OFT p=2 / p=3 | SZ p=2 / p=3 |
|---|---|---|
| κ, entity-aligned 8-block BJ | 154 / 526 | 154 / 524–526 |
| κ, per-entity BJ | 643 / 2338 | 643 / 2338–2529 |
| κ(D⁻¹M) | 706 / 3271 | 706 / 3091–4803 |

The SZ ranges cover Jacobi weights (0,0), (1,2), (3,6) and (5,10). A different hierarchical basis makes no difference: the problem is not within-entity conditioning.

## 4. Nature of the bad modes

- The smallest eigenvectors of the per-entity-scaled M are global, smooth and densely clustered:
  - at p ≥ 2, gradient edge functions nearly cancel face rotationals;
  - at p = 1, the gradient functions nearly cancel the Whitney functions.
- Each such mode spreads over about 10³ entities.
- Two-level additive κ, with a p = 1 coarse level plus blocks, H(Curl)+grad vs vector Lagrange:
  - refined mesh, p = 2: 52 vs 5.7;
  - tiny mesh, p = 3: 347 vs 8.6.
- Two-grid ρ with RAS cell overlap 1 on the refined mesh at p = 2: 0.97 (H(Curl)+grad) vs 0.61 (Lagrange).

## 5. MUG time-stepping runs

### Alfvén cube

Refined mesh, 768 tets, p-MG with a GMRES smoother. Each cell gives average linear iterations per step as H(Curl) / vector-Lagrange B.

| dt | p | BJ METIS | BJ entity | `loverlap` | cell ov0 | cell ov1 |
|---|---|---|---|---|---|---|
| 8e-7 | 2 | 1.3 / 1.3 | 1.7 / 1.3 | 1 / 1 | 1.2 / 1.2 | 1 / 1 |
| 1e-5 | 2 | 10.3 / 5.2 | 11.3 / 5.2 | 3.0 / 2.0 | 5.3 / 3.2 | 2.0 / 2.0 |
| 1e-4 | 2 | 79 / 110 | 108 / 122 | 26 / 19 | 122 / 48 | 11.5 / 9.8 |
| 1e-4 | 3 | 124 / 114 | 66 / 129 | 28.2 / 19 | 67 / 54 | UMFPACK factorization failed / 10.3 |

- The H(Curl) p = 3, cell ov1, 8-block run aborts on the first preconditioner build with `Factorization failed (UMFPACK)`. The overlapping blocks are large; this is probably memory or fill.
- At large dt both representations are limited by the Alfvén coupling, not by the B mass block.

### Spheromak tilt

p = 2, dt = 1e-6, 3 steps, serial, base mesh `nbase=1` (1.66e5 B DOFs). Linear iterations per step:

| smoother | 1-level | p-MG (p=1 coarse) |
|---|---|---|
| block Jacobi, METIS (8) | 109, 108, 99 | 74, 76, 74 |
| block Jacobi, entity-aligned (8) | 119, 119, 104 | 89, 91, 84 |
| block Jacobi + `loverlap` | 24, 26, 25 | 24, 27, 26 |
| cell Schwarz, overlap 0 (8) | 22, 24, 24 | 25, 27, 27 |
| cell Schwarz, overlap 1 (32) | killed (memory, ~14 GB) | killed (memory) |

- With `nbase=2` (1.3e6 B DOFs), every run was killed for lack of memory on the 15 GB test machine.
- Any overlap removes about 75% of the iterations.
- The p = 1 coarse level helps only plain block Jacobi.
- Memory is the practical limit on the cell-overlap-1 smoother.

### Regression check of the merge

`physics/MUG/test_alfven.py::test_hex_nl_r1_p4` failed only when it competed with heavy jobs. Run alone, it passes (363 s).

## 6. Subspace coupling across p-levels and families

γ is the strengthened Cauchy–Schwarz constant in the M inner product: the cosine of the smallest principal angle between two subspaces, γ² = λmax(A11⁻¹A12A22⁻¹A21). γ → 1 means the splitting is unstable: a subspace-correction method on that splitting converges like 1 − (1 − γ²).

Splittings of V_p:

| label | splitting |
|---|---|
| a | V_{p−1} (injected) vs the new level-p functions |
| b1 | raw field split: Bc DOFs vs Bg DOFs |
| b2 | true Helmholtz-type split: R (Bc with the Whitney-carried gradients removed) vs G_p |
| a_rot | within the rotational family: Bc_{p−1} vs new Bc |
| a_grad | within the gradient family: G_{p−1} vs new Bg |

### Global γ (assembled M)

| mesh | p | a | b1 | b2 | a_rot | a_grad |
|---|---|---|---|---|---|---|
| tiny | 1 | – | .956 | .904 | – | – |
| tiny | 2 | .990 | .989 | .988 | .968 | .909 |
| tiny | 3 | .996 | .997 | .998 | .969 | .951 |
| tiny | 4 | .9972 | .9989 | .9990 | .988 | .973 |
| ref | 1 | – | .927 | .890 | – | – |
| ref | 2 | .985 | .981 | .985 | .961 | .897 |
| ref | 3 | .9931 | .9947 | .9955 | .967 | .947 |

### Local γ (maximum over cells, element mass matrices, tiny mesh)

| p | a | b1 | b2 | a_rot | a_grad |
|---|---|---|---|---|---|
| 1 | – | .976 | .897 | – | – |
| 2 | .994 | .994 | .994 | .976 | .957 |
| 3 | .998 | .998 | .999 | .981 | .970 |

### Findings

- Both the hierarchical (a) and the family (b) splittings are close to degenerate, and get worse with p.
- Cross-family coupling (a, b) is stronger than within-family coupling (a_rot, a_grad).
- Separating Bc and Bg naively (b1) is no better than the true split (b2).
- The local γ equal the global ones. The coupling is already present inside a single element, so it cannot come from mesh-scale structure. It is a property of the space and hierarchy.
- This confirms the hypothesis that new rotational and gradient functions at each level couple strongly to each other and to the p−1 space. Any smoother that splits these apart (per-entity or per-field blocks, or additive hierarchical methods) will be slow.

## 7. Re-based (stabilized) hierarchical bases

A basis transform T acts as A_T = Tᵀ A T. Each new DOF keeps its entity, so the same block structures apply. Transforms are applied recursively from p = 1 and are orthogonalized in the M inner product.

| basis | definition |
|---|---|
| i | current OFT basis |
| iii | orthonormalization within each (entity, field, coarse/new) block. By the invariance argument it must equal i for block methods; it is a sanity check. |
| iih | all new level-p functions globally M-orthogonalized against V_{p−1} (dense; ideal bound) |
| iis | structure-preserving, global: new gradients ⟂ G_{p−1}; new rotationals ⟂ (V_{p−1} + new gradients). Keeps G_p spanned by gradient columns. |
| iisl | patch-local iis: each new function is corrected only by coarse or gradient functions that touch its cell-closure patch, truncated to the patch, with patch-restricted M. Sparse and implementable. |

Global transforms are limited to n ≤ 9000, so ref p = 3 has only i, iii and iisl.

### Tiny mesh, CFL 0 (mass only)

| p | basis | κ_diag | κ_bj | ρ_bj | κ_ent | κ_2lev | ρ_2g,bj | xc |
|---|---|---|---|---|---|---|---|---|
| 2 | i | 706 | 155 | 1.82 | 643 | 68.8 | .972 | .15 |
| 2 | iih | 74 | 30.7 | 1.51 | 75 | 25.8 | .953 | 2e-16 |
| 2 | iis | 81 | 27.9 | 1.10 | 80 | 16.1 | .909 | .053 |
| 2 | iisl | 266 | 117 | 1.36 | 277 | 31.7 | .935 | .066 |
| 3 | i | 3270 | 537 | 2.03 | 2338 | 61.4 | .966 | .044 |
| 3 | iih | 88 | 30.7 | 1.51 | 75 | 25.6 | .945 | 3e-16 |
| 3 | iis | 88 | 30.8 | 1.21 | 91 | 14.6 | .902 | .004 |
| 3 | iisl | 605 | 291 | 1.63 | 640 | 39.5 | .944 | .008 |
| 4 | i | 1.39e4 | 1571 | 2.18 | 6647 | 79.6 | .972 | .017 |
| 4 | iih | 120 | 30.7 | 1.51 | 75 | 25.2 | .944 | 2e-16 |
| 4 | iis | 92 | 31.4 | 1.24 | 92 | 14.9 | .903 | .001 |
| 4 | iisl | 1457 | 690 | 1.82 | 1436 | 31.3 | .941 | .003 |

- iii reproduces i exactly for κ_bj, κ_ent, κ_2lev and the two-grid ρ (see `results/`). Only κ_diag moves slightly.
- With iih, κ_bj stays fixed at 30.65 for p ≥ 2. After global orthogonalization, the conditioning is set by the p = 1 → 2 step.

### Tiny mesh, CFL 1 and 10

| p | CFL | basis | κ_bj | κ_2lev | ρ_2g,bj | ρ_2g,ras | ρ_2g,ras+G |
|---|---|---|---|---|---|---|---|
| 2 | 1 | i | 478 | 44.3 | .964 | .0115 | .0058 |
| 2 | 1 | iih | 505 | 89.9 | .978 | .0146 | .0060 |
| 2 | 1 | iis | 264 | 30.4 | .952 | .0145 | .0060 |
| 2 | 1 | iisl | 458 | 34.1 | .958 | .0119 | .0058 |
| 2 | 10 | i | 4.52e4 | 45.3 | .965 | .0115 | .0067 |
| 2 | 10 | iih | 7.32e4 | 1.12e4 | .9998 | .021 | .0070 |
| 2 | 10 | iis | 2.54e4 | 31.2 | .953 | .0145 | .0067 |
| 2 | 10 | iisl | 7.02e4 | 2138 | .999 | .0099 | .0068 |
| 3 | 1 | i | 873 | 28.8 | .938 | .0060 | .0019 |
| 3 | 1 | iih | 2085 | 817 | .998 | .020 | .0019 |
| 3 | 1 | iis | 273 | 21.1 | .932 | .0072 | .0036 |
| 3 | 1 | iisl | 906 | 26.9 | .945 | .0060 | .0020 |
| 3 | 10 | i | 8.24e4 | 29.4 | .940 | .0060 | .0023 |
| 3 | 10 | iih | 2.58e5 | 1.06e5 | 1.000 | .047 | .0026 |
| 3 | 10 | iis | 2.59e4 | 21.3 | .933 | .0072 | .0037 |
| 3 | 10 | iisl | 1.20e5 | 67.2 | .975 | .0059 | .0026 |
| 4 | 1 | i | 1238 | 34.6 | .946 | .0035 | .0007 |
| 4 | 1 | iih | 3215 | 1989 | .999 | .028 | .0008 |
| 4 | 1 | iis | 278 | 21.9 | .934 | .0042 | .0019 |
| 4 | 1 | iisl | 1288 | 26.3 | .946 | .0035 | .0007 |
| 4 | 10 | i | 1.16e5 | 35.3 | .947 | .0035 | .0007 |
| 4 | 10 | iih | 4.36e5 | 2.52e5 | 1.000 | .064 | .0008 |
| 4 | 10 | iis | 2.63e4 | 22.1 | .935 | .0042 | .0019 |
| 4 | 10 | iisl | 1.91e5 | 37.6 | .966 | .0034 | .0007 |

On the tiny mesh the 8 overlapping RAS blocks nearly cover the whole mesh, so ρ_2g,ras ≈ 0.01 says little. The refined mesh is the meaningful test for RAS.

### Refined mesh (768 tets)

| p | CFL | basis | κ_diag | κ_bj | κ_ent | κ_2lev | ρ_2g,bj | ρ_2g,ras | ρ_2g,ras+G |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 0 | i | 363 | 81.7 | 332 | 52.4 | .964 | .974 | – |
| 2 | 0 | iih | 41 | 15.3 | 40 | 15.9 | .918 | 1.05 | – |
| 2 | 0 | iis | 41 | 15.1 | 40 | 11.6 | .878 | .927 | – |
| 2 | 0 | iisl | 103 | 46.9 | 108 | 20.4 | .918 | .929 | – |
| 2 | 1 | i | 901 | 356 | 888 | 22.8 | .932 | .978 | .975 |
| 2 | 1 | iih | 1191 | 315 | 1172 | 57.5 | .965 | 1.06 | 1.05 |
| 2 | 1 | iis | 391 | 178 | 382 | 16.8 | .918 | 1.02 | 1.02 |
| 2 | 1 | iisl | 726 | 324 | 740 | 18.0 | .922 | 1.01 | 1.01 |
| 2 | 10 | i | 7.96e4 | 3.31e4 | 7.68e4 | 22.8 | .933 | .991 | .991 |
| 2 | 10 | iih | 1.32e5 | 4.80e4 | 1.42e5 | 6743 | .9997 | 1.08 | 1.08 |
| 2 | 10 | iis | 3.62e4 | 1.71e4 | 3.50e4 | 17.2 | .920 | 1.02 | 1.02 |
| 2 | 10 | iisl | 1.79e5 | 5.19e4 | 1.77e5 | 645 | .997 | 1.06 | 1.06 |

Refined mesh, p = 3 (bases i, iii and iisl only): still running when this version was committed. Results will be added to `results/rebase_ref_p3.txt` and to this table.

### Findings

- **Global structure-preserving re-basing (iis) is the only variant that helps consistently.** It cuts κ_bj and κ_ent by 5–50× at CFL 0 (more at higher p) and by 2–4× at CFL 1–10.
  - Even then, ρ_2g,bj only moves from about 0.97 to 0.90.
  - RAS two-grid on the refined mesh does not improve at all (0.93–1.02).
  - iis has dense basis functions, so it is not implementable.
- **Orthogonalizing only against V_{p−1} (iih) backfires once curl-curl enters.**
  - At CFL 0 it is as good as iis.
  - At CFL ≥ 1 it is worse than the current basis, by up to 10⁴× in κ_2lev.
  - Reason: it destroys the gradient structure, mixing kernel (gradient) and non-kernel components of K inside each new function.
- **Patch-local iis (iisl) recovers at most 2–5× at CFL 0, and is worse than the current basis at CFL 10.**
  - Truncating the projection to a patch leaves most of the coupling.
  - This matches γ_local ≈ γ_global: the coupling is not removable by local corrections of the basis.
- **The exact G_p solve (ρ_2g,ras+G) changes almost nothing on the refined mesh.** At CFL ≤ 10 the gradient near-kernel is not what limits these smoothers.

## 8. Vertex-star (Arnold–Falk–Winther) blocks

Star(v) contains the free DOFs on every entity whose closure contains vertex v: the vertex itself, the edges, faces and cells around it. This is the classical AFW block for H(Curl) and H(div) problems. With these blocks, local subproblems contain whole local gradients, together with the rotationals that nearly cancel them.

| mesh | p | n | stars | mean / max star size | total size / n | cell ov1 blocks: mean size | cell ov1 blocks: total / n |
|---|---|---|---|---|---|---|---|
| tiny | 2 | 1143 | 38 | 78 / 147 | 2.6 | 818 | 5.7 |
| tiny | 4 | 5545 | 38 | 454 / 885 | 3.1 | 3862 | 5.6 |
| ref | 2 | 8052 | 189 | 111 / 246 | 2.6 | 4171 | 4.1 |
| ref | 3 | 20336 | 189 | 317 / 728 | 2.9 | 10327 | 4.1 |

The star smoother is additive, damped by 1/λmax. Its two-grid is symmetric V(1,1) with an exact p−1 solve: the same form as ρ_2g,bj in §7.

### Tiny mesh (96 tets)

| p | basis | κ_star, CFL 0 / 1 / 10 | ρ_2g,star, CFL 0 / 1 / 10 | vector Lagrange κ_star (CFL 0) |
|---|---|---|---|---|
| 1 | i | 4.70 / 16.6 / 40.3 | – | 5.00 |
| 2 | i | 5.02 / 15.9 / 39.6 | .751 / .729 / .731 | 5.83 |
| 2 | iis | 8.41 / 22.0 / 57.0 | .782 / .762 / .763 | |
| 2 | iisl | 11.4 / 25.7 / 389 | .833 / .911 / .997 | |
| 3 | i | 4.14 / 15.8 / 39.8 | .680 / .689 / .690 | 4.92 |
| 3 | iis | 9.49 / 22.2 / 61.6 | .736 / .748 / .748 | |
| 3 | iisl | 20.5 / 27.0 / 126 | .867 / .840 / .915 | |
| 4 | i | 4.16 / 15.7 / 39.5 | .649 / .652 / .652 | 5.02 |
| 4 | iis | 9.82 / 22.6 / 63.6 | .716 / .726 / .726 | |
| 4 | iisl | 38.5 / 27.6 / 122 | .820 / .819 / .895 | |

### Refined mesh (768 tets)

| p | basis | κ_star, CFL 0 / 1 / 10 | ρ_2g,star, CFL 0 / 1 / 10 | vector Lagrange κ_star (CFL 0) |
|---|---|---|---|---|
| 1 | i | 4.29 / 15.0 / 49.4 | – | 4.77 |
| 2 | i | 4.18 / 16.4 / 54.1 | .715 / .681 / .682 | 5.83 |
| 2 | iis | 6.90 / 20.3 / 71.5 | .749 / .729 / .730 | |
| 2 | iisl | 9.16 / 19.4 / 179 | .814 / .810 / .993 | |
| 3 | i | 3.75 / … / … | .659 / … / … | (pending) |

The rest of the refined-mesh p = 3 star sweep was still running when this version was committed. Results will be added to `results/star_ref.txt`.

### Findings

- **Stars make H(Curl)+grad behave like vector Lagrange.**
  - With vertex stars, the B mass block is as well conditioned as the vector-Lagrange mass block: κ ≈ 4–5 vs 5–6.
  - This holds at every p and on both meshes, versus κ_bj = 82–1571 for entity-aligned blocks.
- **The two-grid rate is p-, mesh- and CFL-robust.** Two-grid ρ with stars and the p−1 coarse solve is 0.65–0.75:
  - it improves slightly with p;
  - it is the same on both meshes;
  - it does not change between CFL 0 and 10.

  With the same damped form, 8 entity-aligned blocks give 0.93–0.97. RAS with cell overlap 1 on the refined mesh gives 0.97–0.99.
- **The single-level κ_star grows with CFL** (about 16 at CFL 1, 40–54 at CFL 10). This is the expected growth from smooth curl-curl modes, which no single-level smoother handles. The coarse p−1 solve removes it, as the flat two-grid ρ shows.
- **Re-basing hurts stars.** iis raises κ_star by about 2× and ρ_2g,star by about 0.05. iisl is worse, and close to 1 at CFL 10. The current basis is the best one to use with star smoothers.

## 9. Interpretation

1. The poor smoother and MG performance of the B block at CFL ≈ 0 comes from strongly coupled pairs:
   - gradient and rotational functions that nearly cancel;
   - the new level-p functions and the p−1 space.

   The coupling is present in every element (γ_local ≈ γ_global ≈ 0.99+), and grows with p.
2. Any smoother whose blocks split such a pair is slow. That includes:
   - per-entity blocks;
   - per-field blocks;
   - entity-aligned partitions where the pair straddles a partition boundary.

   The bad combinations span a vertex patch, so they look global to per-entity smoothers (§4). Large overlapping partitions (cell ov1) capture most of them, but at a high memory cost, and leave a slow two-grid rate.
3. Changing the basis cannot fix this in an implementable way:
   - Within-block basis changes do nothing (§3, iii).
   - Local orthogonalization recovers little, because the coupling is intrinsic to each element's space.
   - Global orthogonalization helps entity blocks but produces dense functions, and backfires at finite CFL unless it preserves the gradient structure.
4. Choosing blocks that contain whole coupled pairs fixes it. Vertex stars are the minimal such blocks; they are the standard AFW construction. They give Lagrange-level conditioning and a two-grid rate of about 0.7, independent of p, mesh and CFL. They are also much cheaper per block than cell-overlap-1 partitions.

## 10. Caveats

- **The analysis covers the B block alone.** Sections 6–8 use the B mass and curl-curl matrices: A = M + τK, with no velocity coupling. The MUG Jacobian adds the Alfvén (B–v) coupling, which dominated at large dt in §5. The star results predict B-block behaviour only; MUG runs with a star smoother are needed to confirm the benefit.
- **Different smoother forms.** The two-grid rates for stars and block Jacobi use a damped additive symmetric smoother. RAS uses an undamped one-sided sweep. Compare like with like (ρ_2g,bj vs ρ_2g,star). The RAS numbers are included because RAS is what OFT's `add_schwarz` does.
- **Small meshes.** Meshes are small (at most 768 tets), tetrahedral only, and periodic boundaries are off. Hex meshes were not analysed.
- **Exact local solves are assumed** everywhere. Memory for star blocks was not measured in OFT.
- **Not computed:**
  - global (dense) transforms for ref p = 3;
  - local γ for the refined mesh (it is the same element family as the tiny mesh).

## 11. Reproducing

1. Build OFT with tests. In `<build_dir>/tests/fem`, run `diag_hcurl_grad_cond` with `inputs/oft.in.tiny` or `inputs/oft.in.ref`, setting `order` to p, in a separate directory per (mesh, p). It writes:
   - `mass.txt`, `curlcurl.txt`, `dofs.txt`;
   - colourings and cell blocks;
   - `mesh.txt`;
   - `lag_*` files for vector Lagrange of the same order.
2. Place the dumps as `scripts/dump3/<mesh>_p<p>/`. `hier.HERE` points to the `scripts/` directory itself.
3. Run from `scripts/`:

   | command | result |
   |---|---|
   | `python hier.py` | hierarchy validation |
   | `python gamma.py global` and `python gamma.py local` | §6 |
   | `python rebase.py tiny 4 0 1 10` and `python rebase.py ref 3 0 1 10` | §7; set `PMIN=3` to report only p = 3 |
   | `python star.py tiny 4 0 1 10` and `python star.py ref 3 0 1 10` | §8 |
   | `analyze.py`, `twolevel.py`, `lowmode.py`, `szcheck.py`, `szmetrics.py`, `val.py` | §2–4 |

4. The MUG drivers are in `mug/`. They expect the test executables (`test_alfven`, `test_alfven_lag`) or the spheromak example binary next to the run directories.

Raw outputs are in `results/`.
