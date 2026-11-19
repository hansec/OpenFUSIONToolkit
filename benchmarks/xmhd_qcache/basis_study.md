# Basis study for the REFERENCE-mode, cell-batched quadrature cache

Branch `claude/xmhd-qcache-basis-study`, from `claude/xmhd-qcache-phase2` (`5328aab`). This is an analysis
study. No OFT source file is changed, and no basis changes its default behavior. All prototypes are
test-only Python scripts in `benchmarks/xmhd_qcache/basis_study/` (see [Reproduction](#reproduction)).

**Why this location.** The report sits next to the Phase 0-2 reports it builds on (`phase2_results.md`). The
scripts it cites live in the same folder. `src/docs` is processed by Doxygen and published as user
documentation, and an internal design study of un-adopted bases does not belong there.

## Summary

1. **The tet face problem is structural for H(Curl) at p = 2, 3, and for H^1 at order 4 if the basis stays
   hierarchical. No choice of basis removes it.** Each face carries a representation of its symmetry group S3:
   - S3 has three irreducible representations: trivial (T), sign (A) and the 2-d standard one (S).
   - A face block can transform by **signed permutations** iff `m_S <= m_T + m_A`, i.e. every S can be paired
     with a 1-d irrep.
   - The multiplicities are intrinsic to the spaces, so they hold for any basis:

   | face block | p=2 / k=3 | p=3 / k=4 | p=4 / k=5 |
   |---|---|---|---|
   | gradient-free H(Curl), whole block | **S (no)** | **A+2S (no)** | T+2A+3S (yes) |
   | H(Curl), hierarchical increment | **S (no)** | A+S (yes) | T+A+S (yes) |
   | H^1 face, whole block | T (yes) | T+S (yes) | 2T+2S (yes) |
   | H^1, hierarchical increment | T (yes) | **S (no)** | T+S (yes) |

   - Edges (Z2), hex quad faces (D4, verified here up to p=4) and Lagrange spaces (symmetric node sets) have no
     obstruction.
   - So "every transform is a signed permutation" cannot be reached on tets by changing only the basis, and
     hierarchical H^1 cannot reach it at order 4.
2. **The defect that costs most today is not the obstruction itself.** It is that OFT's face families are not
   S3-invariant:
   - The current H(Curl) faces need the gradients of H^1 and the cell interior ("kind 3": no transform within the
     space).
   - The current H^1 faces need the interior bubbles: dense blocks with 16-25 regressors.

   Both can be fixed: a **symmetric, hierarchical** choice of the face functions makes every transform an
   **own-entity** transform. That is a signed permutation, plus **one 2x2 block per face** (the single S that
   cannot be paired). The same structure appears independently in MFEM `ND_TetDofTransformation`, Basix N1curl
   and libCEED's curl-oriented restriction.
3. **Prototype: the projected symmetric families** (`sym_orth.py`, "sym-proj"), checked numerically for
   p = 1-4 / k = 1-5:
   - **Spaces:** the H^1 space is unchanged; the combined `H(Curl) ⊕ grad H^1` is unchanged, which is all xMHD's
     B field sees.
   - **Gradient-free and hierarchical:** the new H(Curl) is still gradient-free (direct sum with grad H^1 apart
     from the existing Whitney/vertex overlap), and both spaces stay hierarchical (identity p-injection).
   - **Transforms:** all are own-entity. Kind 3 is eliminated, and so are the 16-25-regressor H^1 blocks.
   - **Conditioning** (local matrices, Jacobi-scaled): within 10% of OFT or better at p <= 3 / k <= 4. At p=4 it is
     within 1.2-2.2x.
4. **What changes:** only the **H(Curl)-only span**, which moves by gradients.
   - Unaffected: nonzero Taylor eigenvalues and `hcurl_divout`, since the gradient-free property holds.
   - Changed: the discrete J2 (hyper-resistivity) field, Taylor vector potentials (gauge), basis-dependent test
     reference values (iteration counts) and old restart files.
5. **New observation.** OFT's H(Curl) **cell-interior** block is not invariant under relabeling a cell's
   vertices. The standalone H(Curl) space therefore already depends on each cell's local vertex order
   (`hcurl_basis.F90:642-715`); the composite does not. This rules out cheap "virtual sort" schemes for H(Curl)
   and is a correctness caveat for any mesh-reordering option. **It is fixable without mesh changes**
   (section 3.7): projecting the interior levels off the H^1 interior gradients makes the interior span exactly
   relabeling-invariant, keeps the gradient-free direct sum and the hierarchy, and needs no cache transforms.
6. **Mesh-side alternative.** Sorting tet vertices by global index, then swapping local 3/4 when the cell is
   inverted (instead of allowing det<0), leaves exactly **2 orientation classes** of 12/24 orderings each. Each
   class has a fixed pattern of states. No transforms remain for any basis: two canonical tables, with det>0
   preserved. It is the strongest result for GPUs, but it is grid-code work across all MG levels and changes the
   H(Curl)-only interior span.

**Recommendation, ranked** (details in [section 6](#6-ranked-recommendations)):
1. Implement the symmetric hierarchical face families (sym-proj) as an **opt-in basis variant** for H(Curl)
   (p >= 2) and H^1 (order >= 4).
2. As a fallback with no basis change, use cross-space transforms for composite (B) fields.
3. Consider the sort+parity mesh option only if per-face 2x2 blocks prove to matter on GPU.
4. Do not adopt a literature basis wholesale.

---

## 1. Where orientation enters (OFT conventions)

**Entity frames are already "global-vertex-ordered".**
- `mesh_global_orient` (`src/grid/mesh_global_util.F90:231-282`) sets `lce` signs and `lcfo` codes from global
  point indices. Periodic ties are broken by local index.
- `find_orient_listn` (`src/base/oft_base.F90:977`) and `mesh_local_orient` (`src/grid/mesh_local_util.F90:33-52`)
  turn these into `oriented_edges`/`oriented_faces`. For triangles this is the vertex list **sorted by global
  index**.
- Every edge/face formula is evaluated on that sorted triple, so conformity holds by construction (as in
  NGSolve/ESEAS).
- Interior functions use the cell's raw local barycentrics `f(1:4)` (`hcurl_basis.F90:642-715`,
  `h1_basis.F90:480-503`), with no orientation step.

**Reference tables.**
- A reference table must hold each entity's functions for each *position* of the sorted order inside the
  reference cell. That gives 2 edge states, 6 triangle face states (S3) and 8 quad face states (D4)
  (`fem_quad_cache.F90:375-399`).
- Phase 2 tabulates every state and fits `phi^s = M phi^1` (`qcache_build_xform`, `fem_quad_cache.F90:904-1110`).
  The cell-batched kernel then maps coefficients to the canonical frame (`qcache_to_canon`, 1696-1740).
- Kind-3 entities are excluded from the canonical contraction and recomputed from their own state block, using
  a short, non-SIMD loop (`qcache_interp_cell_ref`, 1544-1561).

**What a signed permutation needs.** Two conditions, for the block of one entity:
- (a) Its span must be invariant, exactly as cell functions, under the group of the entity. A face block that
  is invariant only "modulo" edges, the interior or gradients gives kind 2 or kind 3.
- (b) The resulting representation must be **monomial**, i.e. realizable by signed permutation matrices.

On (b):
- For Z2 (edges) and D4 (quad faces) every real irrep is monomial; D4's 2-d irrep is the symmetry group of the
  square, i.e. 2x2 signed permutations. So only (a) can fail there.
- For S3, a signed-permutation representation is a sum of representations induced from ±1 characters of
  subgroups: `T`, `A`, `T+S` (Z2, trivial character), `A+S` (Z2, sign character), `T+A` (Z3) and the regular
  representation. Hence the criterion `m_S <= m_T + m_A`.
- The multiplicities of the irreducible representations are invariants of the space. Here "the space" is the
  quotient of the face-trace space by the edge traces, and for H(Curl) also by the gradient face traces, so
  every basis of the same split sees them.

Characters are computed in `rep_analysis.py` from the transform matrices of the current basis. The literature
agent reproduced the same multiplicities independently.

## 2. Survey of the existing bases

### 2.1 Where each space is used

- **xMHD B field.** Composite `H(Curl)_p ⊕ grad H^1_{p+1}` (`oft_hcurl_grad_setup`, `hcurl_basis.F90:203-262`;
  fields 1+2 in `xmhd.F90:4316-4318`).
- **xMHD J2 (hyper-resistivity auxiliary) is standalone H(Curl)** (`xmhd.F90:4339`, `xmhd_lag.F90:1500-1563`). It
  is an L2-type projection, so it depends on the complement.
- **V, N, T** use Lagrange.
- **Standalone H(Curl) users:**
  - Taylor `taylor_hmodes`, `taylor_vac_curr` and the injector routines (`src/physics/taylor.F90:294-527, 672-844,
    849+`). These solve `WOP u = λ KOP u` with KOP the helicity operator (`hcurl_operators.F90:551`). Both
    operators vanish on gradients, and the only cleaning is the linear-Lagrange `hcurl_divout`
    (`hcurl_operators.F90:103-105, 1025-1074`). So they need a **gradient-free** space, but the nonzero λ depend
    only on `span(H(Curl)) + gradients` (a quotient argument).
  - `hcurl_grad_getmop`, which removes the vertex H^1 rows as duplicates of the Whitney gradients
    (`hcurl_grad_operators.F90:961-967, 1052-1054, 1084-1097`).
  - `hcurl_grad_grad` / `gradtp` (683-789).
  - The Marklin wrappers (`src/python/wrappers/marklin_f.F90:118-538`).
  - Restart readers (`src/bin/oft_trace.F90:147-180`, `oft_poincare.F90:204-220`, `marklin_eigs.F90:92-113`).
- **p-multigrid** is a pure 0/1 identity injection for H(Curl) (`hcurl_pinterpmatrix`,
  `hcurl_operators.F90:1359-1507`) and H^1 (`build_pinterpmatrix`, `h1_operators.F90:711-871`). The composite
  prolongation is block-diagonal (`fem_composite.F90:688-739`). All of these **require hierarchical bases**:
  face/cell DOF `j` at order p is the same function at order p+1. Lagrange uses true nodal interpolation
  (`lagrange_operators.F90:1700-2023`).
- **First-DOF assumptions.** Edge DOF 1 is the Whitney function in `hcurl_grad`, `hcurl_grad_grad`, `_mc`/`_bmc`
  and `xmhd.F90:4993-5005`. Nothing proposed here changes edge functions.

### 2.2 Tet bases: results of the Phase 2 classifier (reproduced)

`xform_classify.py` is a line-by-line Python port of `qcache_build_xform` (own block → + sub-entities and
interior → all; snapping to signed permutations), applied to a symbolic transcription of the tet evaluators
(`oft_tet_basis.py`).
- It reproduces the Fortran counts exactly, including the regressor counts.
- The Fortran counts were re-run here for p = 1-4 on a 2x2x2 cube.
- The current transforms use the following functions (`results_structure.txt`, all six face states):

| space | Fortran identity/perm/dense/none | regressors actually used by a face transform |
|---|---|---|
| Lagrange p=1-4 | identity/perm only | own block |
| H^1 k=2, 3 | 6/0/0/0, 20/6/0/0 | own block |
| H^1 k=4 | 0/10/16/0 (16 regressors) | own face (3) + **interior bubble** (1) |
| H^1 k=5 | 0/10/16/0 (25 regressors) | own face (6) + **interior** (4) |
| H(Curl) p=1 | 0/6/0/0 | own edge |
| H(Curl) p=2 | 0/10/0/16 | own face (2) + **grad of the H^1 face function** (1) |
| H(Curl) p=3 | 0/10/0/16 | own face (5) + grad H^1 face (3) + **H(Curl) interior (3) + grad H^1 interior (1)** |
| H(Curl) p=4 | 0/10/0/16 | own face (9) + grad H^1 face (6) + H(Curl) interior (11) + grad H^1 interior (4) |

### 2.3 Why: the construction choices that break equivariance

Notation, for a face with sorted vertices 1,2,3 (barycentrics `x_i`, `b = x1 x2 x3`):
- `A = x2 x3 ∇x1`, `B = x1 x3 ∇x2`, `C = x1 x2 ∇x3`, so that `∇b = A + B + C`;
- `w_ij = x_i ∇x_j − x_j ∇x_i` (Whitney form);
- `s = x1 + x2 + x3`. Inside the cell `s = 1 − x4`, so **`s·φ ≠ φ`**: it differs by `x4 φ`, an interior bubble.

**H^1 faces, order 4** (`h1_basis.F90:450-478`):
- The increment uses `2b(x1+x2−x3)` (line 463) and `2b(x2−x1)`.
- The first linear form has coefficient sum 1, not 0. Under a 3-cycle it becomes `2b(x2+x3−x1) = 2b·s − 4b x1`,
  which is face DOF 1 **plus the interior bubble** `b x4`. That is exactly the 16-regressor fit.
- A sum-zero form such as `2b(x1+x2−2x3) = (3·old + φ_f1 − φ_c1)/2` spans the same space with the interior, and
  makes the increment exactly invariant.
- But sum-zero linear forms carry the irrep S: one 2x2 block is unavoidable for a hierarchical basis.

**H^1 faces, order 5:** the quadratic forms of DOFs 4-6 are likewise not an invariant complement. The fit needs
the 4 interior functions.

**H(Curl) faces, p=2** (`hcurl_basis.F90:597-604`):
- DOF 1 is `−2(A + B − C) = −2∇b + 4C` and DOF 2 is `A − B`.
- The rotational quotient is `span{A,B,C}/span{∇b}`, which is S.
- DOF 1 contains a gradient, `−2∇b`, and is not sum-zero. A 3-cycle sends it to `−2∇b + 4A`, which lies in
  `span{DOF1, DOF2}` only with the help of `∇b`, a gradient that lives in the other space.
- Adding `(2/3)∇b` (i.e. subtracting `(1/3)∇φ^{H1}_{f,1}`) gives `(4/3)(2C − A − B)`. The block then becomes
  exactly invariant, with no other change to the composite space.

**H(Curl) faces, p=3, 4** (DOFs 3-9, lines 606-637):
- The "gradient-free" complement was chosen per DOF, using the sorted order `(1,2,3)` asymmetrically: factors
  like `(−f1−f2+f3)` and `(2f1−f2)`.
- These sets are not invariant modulo lower levels and gradients. Their rotations need grad H^1 face functions,
  and the face function's extension into the cell also differs by H(Curl) and grad H^1 interior functions.

**Interiors:**
- H^1 interiors are S4-invariant (`b4·P_{k−4}`).
- **H(Curl) interiors are not.** Relabeling the 4 vertices takes them outside the H(Curl) span (residual 7.7e-2),
  but stays inside `H(Curl) + grad H^1` (≤1.4e-14; `results_structure.txt`).
- So OFT's standalone H(Curl) space depends on each cell's local vertex order. The composite space does not.
- This is harmless today, because interiors are never oriented, but it matters for the mesh-side options in
  section 4.

### 2.4 Hex and 2-D/boundary bases

**Hex:**
- D4 has no obstruction, and OFT's tensor-product families have definite parity. `hpoly_bary(k,x) =
  x J^{(1,1)}_{k−1}(2x−1)` (`src/grid/lagrange_poly.F90:446-472`), so `x → 1−x` gives a sign and `u ↔ v` permutes
  `hpoly_2d` indices.
- Phase 2 verified p = 1-3. **Re-run here for hex p=4** (H^1 order 5, H(Curl) p=4, Lagrange p=4, and their p−1
  levels): every (entity, state) pair is a signed permutation (54/54 in each space; fit residual ≤ 2e-15). No
  change is needed on hexes.

**2-D/boundary:**
- `oft_bh1_eval` and `oft_bhcurl_eval` (`h1_basis.F90:385-407`, `hcurl_basis.F90:532-564`) reuse the same
  `evalf` formulas with the triangle's **cell** orientation `lco`.
- On a surface mesh the face block *is* the cell interior, so every boundary triangle has one of 6 states, and
  the S3 analysis above applies unchanged. The same symmetric families would fix it. `oft_blag_*` is nodal,
  hence permutations.
- This only matters if boundary integrals are batched. Currently `oft_bhcurl` is used only by `hcurl_bcurl`
  quadrature and `test_hcurl_sop`.

**Lagrange (tet):** nodal on S3-invariant face lattices, with edge reversal by index reversal
(`lagrange_basis.F90:734`). It is always a permutation.

## 3. Simple modifications of the existing bases

### 3.1 What can be reached

From the table in the summary:
- **Hierarchical** (needed for OFT's identity p-injections): **signed permutations everywhere, except one 2x2
  block per triangle face** in H(Curl) (p ≥ 2, the p=2 level) and H^1 (order ≥ 4, the order-4 level). Every
  other increment is monomial and can be given a signed-permutation basis.
- **Non-hierarchical:**
  - H^1 can be pure permutations at every order (Bernstein-like face blocks `b·x^α`).
  - H(Curl) can be signed permutations only at p=4. At p=2 and p=3 it cannot, for any basis.
  - The price is a non-identity p-prolongation.

The 2x2 block costs about 4 flops per face per field in `qcache_to_canon`. It is a fixed matrix per face state
(24 numbers), so on a GPU it is as cheap as a sign flip. The literature uses this choice (MFEM's six 2x2 integer
matrices; Basix's `[[-1,-1],[1,0]]`).

### 3.2 Prototype families (`sym_basis.py`, `sym_orth.py`)

Edge, vertex and interior functions are OFT's, unchanged. Face functions are hierarchical by DOF index, and each
level is exactly S3-invariant as cell functions. Cyclic triples use `(k; i, j)` with `k` the vertex opposite edge
`ij`.

| level | H(Curl) generators | rep | H^1 generators | rep |
|---|---|---|---|---|
| p=2 / k=3 | `(2C−A−B)/√3`, `A−B` (= OFT DOF 2) | S: one 2x2 (orthogonal; integer variant `B−C`, `A−B`) | `−2b` (= OFT DOF 1) | T |
| p=3 / k=4 | `x_k² w_ij` (3) | A+S: signed perm | `2b(x3−x2)`, `2b(x2−x1)` (sum-zero) | S: one 2x2 |
| p=4 / k=5 | `x_k²(x_i−x_j) w_ij` (3) + `b(w12+w23+w31)` | T+S perm + A sign | `2b x_j x_k` (3) | T+S perm |

- **"sym"** uses these generators directly.
- **"sym-proj"** projects each level, in L2 on a **regular** tet (S3 then acts by isometries), off:
  - the previous levels of the same face, and
  - for H(Curl), the gradients of the H^1 face functions up to the same degree.
- The projection commutes with S3, so the structure is kept exactly. The formulas are rational (exact sympy
  arithmetic) and are printed by `python sym_orth.py`.
- In Fortran these would be "generator + a few fixed coefficients times lower face functions / H^1 face
  gradients", evaluated inside `evalf`.
- The non-hierarchical **H1perm** variant (order-4 block `{−2b x_i}`) was also tested.

### 3.3 Space preservation: argument and evidence

**H^1.** Each new face function is a face bubble (`b` factor) of the right degree. The full set has the right
count, and the rank test shows it spans exactly the old space (`rank(old|new) = nce`, k = 1-5, for sym, sym-proj
and perm).

**Composite.**
- New H(Curl) face functions are OFT's functions plus lower face functions plus gradients of H^1 face functions.
  Each is a face bubble.
- So `span(H(Curl)_new) + grad H^1_{p+1} = span(H(Curl)_old) + grad H^1_{p+1} = P_p^3` (NED2_p).
- Independence modulo gradients is preserved, so H(Curl)_new remains a **gradient-free complement**.

Numerical evidence (ranks of sampled values+curls, `results_sym.txt`):

| p | rank HC old / new | rank(HC_old ∪ HC_new) | rank(HC+gradH1) old / new / union | NED2_p dim | direct-sum count |
|---|---|---|---|---|---|
| 1 | 6 / 6 | 6 | 12 / 12 / 12 | 12 | 12 |
| 2 | 14 / 14 | 18 | 30 / 30 / 30 | 30 | 30 |
| 3 | 29 / 29 | 38 (sym), 42 (proj) | 60 / 60 / 60 | 60 | 60 |
| 4 | 53 / 53 | 76 (sym), 80 (proj) | 105 / 105 / 105 | 105 | 105 |

("Direct-sum count" = nce(HC) + nce(H^1) − 4 vertex functions; equality means no gradients beyond Whitney.)

**Preserved:**
- H^1 space;
- composite space (xMHD B, `hcurl_grad` divout, Marklin `get_aint`);
- the gradient-free property (`hcurl_divout`, Taylor's singular solves, `hcurl_grad_getmop` vertex
  elimination);
- nonzero Taylor eigenvalues and helicity/flux quantities (quotient argument);
- hierarchy, so `hcurl_pinterpmatrix` / `build_pinterpmatrix` are unchanged;
- all edge and interior functions, hence all first-DOF / Whitney assumptions;
- hex bases.

**Changed:**
- **the H(Curl)-only span**, which moves by gradients (rank of the union grows by 4/9-13/23-27);
- therefore the discrete J2 field in xMHD (L2 projection in standalone H(Curl)) and the Taylor vector potential
  `hffa` (a gauge);
- the meaning of face DOFs, so old restart/DOF files are invalid;
- basis-dependent test references: `test_hcurl.py` iteration counts and `uu`, `test_hcurl_grad.py` iteration
  counts, the qcache `hcurl_span_check` (which tests the *current* defect), and MG Jacobi damping defaults
  (`hcurl_mloptions`), which should be re-tuned;
- for H1perm only, the 3→4 p-injection, which becomes `b = Σ b x_i + (cell bubbles of both neighbours)`.

### 3.4 Orientation transforms after the change

| space | current | sym / sym-proj |
|---|---|---|
| H^1 k=4 | 0/10/16/0, 16 regressors | 0/10/16/0, **own block only** (3), level blocks: 1x1 identity + one 2x2 |
| H^1 k=5 | 0/10/16/0, 25 regressors | 0/10/16/0, **own block only** (6), level blocks: 32/48 signed perm, rest 2x2 |
| H^1 k=4,5 (perm) | – | **0/26/0/0** (pure permutations; non-hierarchical) |
| H(Curl) p=2 | 0/10/0/**16** | 0/10/16/0, own block (2): one 2x2 |
| H(Curl) p=3 | 0/10/0/**16** | 0/10/16/0, own block (5): 2x2 + 3x3 signed perm |
| H(Curl) p=4 | 0/10/0/**16** | 0/10/16/0, own block (9): 2x2 + two signed-perm levels |

"Dense" here means only own-entity blocks, block-diagonal by level, with the single non-permutation block 2x2.
The 2x2 blocks for the orthogonal p=2 variant are rotations by ±120° and reflections, for example state 2:
`[[-1/2,-√3/2],[√3/2,-1/2]]`.

Effect on the cell-batched kernel:
- About 80% of tet faces are in a state ≠ 1 on real meshes (Phase 2 counts: 1221 of 6144 faces in state 1).
- So today ≈3.2 faces per cell take the kind-3 path, i.e. **≈46/55/54% of the H(Curl) contraction work at
  p = 2/3/4** goes through a short, unvectorized, state-dependent loop that reads extra state blocks.
- With the symmetric bases all of it moves into the single canonical contraction. `fem_qcache_canonical_ok`
  becomes true for every space, which also enables the canonical per-qp path and `fem_qcache_mat_from_ref` for
  Jacobian GEMMs with H(Curl).
- Table memory could also drop to the canonical rows, plus 6 small matrices per face type.
- These are flop-share estimates, not timings. Timing needs a Fortran implementation of the variant.

### 3.5 Conditioning

Local matrices, eigenvalue ratio after symmetric Jacobi scaling (zero modes removed), from `results_sym.txt`.
The regular tet is shown; the unit tet gives the same ranking.

| space | old | sym (int) | sym (orth) | **sym-proj** | perm |
|---|---|---|---|---|---|
| H^1 k=4 mass / stiffness | 4.21e3 / 270 | 3.64e3 / 228 | 3.64e3 / 224 | **3.64e3 / 228** | 4.37e3 / 241 |
| H^1 k=5 mass / stiffness | 1.76e4 / 819 | 4.55e4 / 2667 | 4.55e4 / 2639 | **1.81e4 / 1436** | 2.95e4 / 1672 |
| H(Curl) p=2 mass / curl-curl | 38.9 / 10.6 | 54.1 / 24.3 | 33.0 / 10.6 | **33.0 / 10.6** | – |
| H(Curl) p=3 mass / curl-curl | 192 / 56.8 | 1104 / 372 | 647 / 239 | **209 / 56.7** | – |
| H(Curl) p=4 mass / curl-curl | 616 / 220 | 3002 / 1742 | 2149 / 1327 | **1070 / 487** | – |
| HC+gradH1 mass p=2 / 3 / 4 | 145 / 772 / 4393 | 153 / 2894 / 18910 | 138 / 1919 / 15130 | **138 / 644 / 5290** | – |

The naive symmetric generators are 3-9x worse at p ≥ 3. The S3-equivariant projection recovers OFT's
conditioning (p ≤ 3 within 10% or better; p=4 is 1.2-2.2x, worst for the regular-tet curl-curl matrix). The p=4
generators could likely be tuned further: any S3-equivariant change of the generators is allowed. Adding the
symmetric interior block of section 3.7 also recovers most of the p=4 loss.

*Correction (follow-up):* the first version of `sym_orth.py` used a 3x3 face metric with 4 on the diagonal
instead of the regular-tet value 3 (`∇λ_i·∇λ_j ∝ 4δ_ij − 1`). That matrix is still S3-invariant and positive
definite, so every structural result was unaffected. The sym-proj numbers above are from the corrected metric
and differ from the first version by at most 5% (p=4 combined mass).

### 3.6 Fallback without changing any basis: cross-space transforms for composite fields

The current kind-3 fits are exact in `H(Curl) + grad H^1` (residual ≤ 7e-14), using the face block, the H^1 face
gradients and the interiors of both spaces (table in 2.2). For a field that lives in the composite (B):
- the canonical-frame map can be a **per-cell joint transform** of the H(Curl) and H^1 coefficients;
- curls only need the H(Curl)→H(Curl) part.

This needs no basis change. It does need the two caches to be coupled, and dense blocks that involve interior
DOFs (up to 9x30 at p=4). Standalone H(Curl) fields (J2) and H(Curl)-only test functions still need kind 3, or a
grad-H^1 contraction of the same flux. It is a reasonable fallback if the basis variant is rejected. It is not
better for GPUs.


### 3.7 Symmetric H(Curl) cell interior (follow-up feasibility check)

**Problem.** OFT's H(Curl) interior functions (`hcurl_basis.F90:642-715`, and the inline cell parts of
`eval_all3/4` / `ceval_all3/4`) are evaluated on the cell's raw local barycentrics. Their span is not invariant
under the 24 relabelings of a cell's vertices (section 2.3). So the standalone H(Curl) space, which is what J2
and the Taylor vector potential see, depends on the local vertex order of every cell. The composite and the H^1
interiors are invariant.

**Construction** (`interior_sym.py`, exact sympy arithmetic). This is the face construction applied to the
interior, with S4 instead of S3:
- Take OFT's interior functions level by level: degree 3 (DOFs 1-3) and degree 4 (DOFs 4-11).
- Project each level, in L2 on a regular tet, off the gradients of the H^1 interior bubbles up to order p+1 and
  off the previous projected level.
- S4 acts on the regular tet by isometries and both subtracted spaces are S4-invariant, so the projection
  commutes with relabeling.
- The new functions differ from OFT's only by interior gradients and lower interior functions. They stay
  interior bubbles, so faces, edges and conformity are untouched.

**Result** (`results_interior.txt`):

| check | p=3 | p=4 |
|---|---|---|
| H(Curl)+grad H^1 rank, old / new / union (NED2_p) | 60 / 60 / 60 (60) | 105 / 105 / 105 (105) |
| direct sum with grad H^1 (gradient-free) | yes (60 = 29 + 35 − 4) | yes (105 = 53 + 56 − 4) |
| rank(HC_old ∪ HC_new), OFT faces | 30: one interior function changes by a gradient | 57 |
| hierarchy: degree-3 functions identical in p=3 and p=4 | yes | yes |
| relabeling residual of the interior block, old → new | 7.7e-2 → **1.1e-15** | 7.7e-2 → **2.1e-15** |

At degree 3 the change is small. DOFs 2 and 3 are unchanged, and DOF 1 becomes
`−l2l3l4 ∇l1 − l1l3l4 ∇l2 + 3 l1l2l4 ∇l3 − l1l2l3 ∇l4`, which is OFT DOF 1 corrected by an interior gradient.

Conditioning (Jacobi-scaled eigenvalue ratio, H(Curl) mass / curl-curl / combined HC+gradH1 mass):

| tet | p | old | OFT faces + sym interior | sym-proj faces | sym-proj faces + sym interior |
|---|---|---|---|---|---|
| regular | 3 | 192 / 56.8 / 772 | 194 / 56.8 / 765 | 209 / 56.7 / 644 | 276 / 56.7 / 615 |
| regular | 4 | 616 / 220 / 4393 | **484 / 122 / 2797** | 1070 / 487 / 5290 | 819 / 213 / 2535 |
| unit | 3 | 385 / 97.2 / 1897 | 348 / 97.2 / 1884 | 386 / 106 / 1520 | 449 / 106 / 1403 |
| unit | 4 | 1211 / 376 / 10920 | **878 / 232 / 7340** | 2132 / 681 / 12850 | 1457 / 303 / 6119 |

- With OFT's faces, the symmetric interior is as good as or better than the current one everywhere. At p=4 it
  is better by 1.3-1.8x.
- Combined with the sym-proj faces, the p=4 loss of section 3.5 mostly disappears. The curl-curl and combined
  mass are within 10% of OFT's or better (better at p=4). The H(Curl) mass is 1.2-1.4x higher at p=3,4.

**Consequences:**
- Interiors are never oriented, so `fem_quad_cache` is unaffected: one interior table, no transforms.
  Invariance only makes the space well defined.
- The S4 action on the interior block is a dense own-block representation, not a signed permutation. That does
  not matter here. It would matter only for a "virtual sort" scheme, which with this change would become exact
  for H(Curl) interiors, using 24 small own-block matrices.
- Preserved and changed spaces are the same as for the faces (3.3): the composite, the gradient-free property
  and the hierarchy are preserved. The H(Curl)-only span changes by interior gradients, which affects J2, the
  Taylor A gauge, test references and restart files.
- It touches the same evaluators as the face change, so it is cheapest to do both together as one opt-in
  variant.
- Not done: tuning of the generators, and any Fortran implementation (left for a separate session).

## 4. Mesh-side alternatives

| option | transforms left | det sign | what breaks or changes | effort / risk |
|---|---|---|---|---|
| **UFC sort** (local order = global order; FEniCS ≤2019, MFEM `ReorientTetMesh`) | none | ~50% of cells det<0 | Signed `det` at ~167 `%jacobian(` call sites (`tetmesh_type.F90:461-562, 762-781`; e.g. `hcurl_operators.F90:715`, `xmhd.F90:1851`); `mesh_volumes` re-inverts (`mesh_local.F90:83`); `multigrid_check_ho` rejects vol ≤ 0 | high / high |
| **Sort + parity fix** (sort by `(global lp, local lc)`, then `invert_cell` swap of local 3↔4 if det<0) | none, **two cell classes** with fixed state patterns (12/24 orderings each: edges `(1,1,1,1,2,1)` / `(1,1,2,1,2,1)`, faces `(1,5,1,1)` / `(6,2,1,1)`; `structure_checks.py`) | det>0 kept | h-refinement children come out unsorted (`tetmesh_tessellation.F90:77-135`); fine-level `global%lp` is assigned only after `mesh_local_init` (`multigrid.F90:1244-1249`, `multigrid_build.F90:622-736`), so a re-sort + rebuild of `lce/lcf/...` per level is needed; periodic ties; **interior functions change per cell** (H(Curl)-only interior span changes, see 2.3; H^1 / Lagrange interior DOF meaning changes) → restart files invalid | medium-high / medium |
| **"Virtual sort"** (keep mesh, per-cell vertex permutation: DOF permutation + qp permutation, needs S4-symmetric quadrature; OFT's tet rules are documented as symmetric, `quadrature.F90:292`, not verified here) | none for edges/faces | unchanged | H(Curl) interiors are **not** S4-invariant (2.3), so they must be contracted separately in the cell's own frame (two GEMMs per cell) or transformed per permutation (24 dense blocks); with a per-cell qp gather | medium / low |
| tets only, any option | – | – | hexes don't need it (D4) | – |

The sort+parity option is the only route to zero transforms for every basis and order without negative
determinants. It also makes the H(Curl) basis question moot for the cache. Its cost is mostly in the grid code
(MG levels, MPI global numbering, periodic meshes), and it changes the standalone H(Curl) space through the
interiors. I do not recommend it as the next step. It is worth reconsidering if a GPU implementation finds the
per-face 2x2 blocks or per-cell state lookups costly, which they should not be.

## 5. Literature and reference survey

Sources were read directly where possible: NGSolve `hcurlhofe_impl.hpp`, MFEM `fem/doftrans.cpp`, Basix
`finite-element.h` and its DOF-transformation demo. Paper contents come from abstracts and snippets, so details
marked (u) are unverified.

| approach | spans P_p / NED1 / NED2; grad-rot split | orientation behavior | batched / GPU fit | conditioning p ≤ 4 | p-MG | OFT effort |
|---|---|---|---|---|---|---|
| **Schöberl–Zaglmayr / NGSolve** (COMPEL 24:374, 2005; Zaglmayr thesis 2006) | full NED1/NED2; explicit split (gradients of H^1 + `wuDv−wvDu`-type rotational families); `nograds` flag = OFT-like rotational space | sorted global vertices per entity (`GetEdgeSort` / `GetFaceSort`), no DOF transforms; the rotational face span is **not** S3-invariant at p ≥ 3 (agent check; same defect as OFT) | needs per-cell sorted tables, or a pre-sorted mesh | good (Legendre / Jacobi / Dubiner) | hierarchical | OFT's tet bases are already SZ type (see 5.1); no gain over sym-proj for the cache |
| **Fuentes–Keith–Demkowicz–Nagaraj 2015** (CAMWA 70:353; ESEAS Fortran) | full exact sequence, all shapes; split implicit (u) | "orientation embedded": functions re-evaluated with permuted affine coordinates, i.e. 6/8 tables per face; no transform analysis (u) | same as OFT's current tables | good | hierarchical | moderate (code exists); does not address the problem |
| **Ainsworth–Coyle 2001/2003; Webb 1999** | NED1/NED2 with an explicit split ("complete" vs "mixed" order) | global vertex ordering | as NGSolve | A–C reported worse; Webb (u) | hierarchical | moderate; no advantage |
| **Graglia–Wilton–Peterson 1997** (IEEE TAP 45:329) | NED1 interpolatory; no split | node permutation + 2x2 per point pair (u) | OK | OK | non-hierarchical | poor fit: loses the split OFT needs |
| **Basix / DOLFINx** (Scroggs et al., ACM TOMS 48(2) 2022) | any element | base transformations per edge reversal / face rotation / reflection; permutations when possible, else dense (N1curl degree 2 tet face: `[[-1,-1],[1,0]]`) | `T_apply` per cell; this is OFT's Phase 2 model, formalized | element-dependent | element-dependent | low conceptually; confirms 2x2 is the floor |
| **MFEM `DofTransformation`** (`ND_TetDofTransformation`) | NED1 | signed DOF map for edges + six 2x2 integer matrices per face DOF pair (= the S irrep); replaced `ReorientTetMesh` | yes; libCEED `CreateCurlOriented` restriction (tridiagonal, CUDA kernels) does the same in gather/scatter | OK | non-hierarchical | the runtime model recommended here |
| **libCEED / MFEM PA, deal.II matrix-free (hex)** | tensor-product H^1 / ND | signed permutations (D4); `CeedElemRestrictionCreateOriented` | sum factorization, the GPU standard | good | – | relevant to future hex kernels; OFT hex is already permutation-only |
| **Bernstein–Bézier** (Ainsworth–Andriamaro–Davydov SISC 2011; Kirby 2011; Ainsworth–Fu CMAME 2018) | full spaces; explicit gradient / non-gradient split (A–F) | H^1 pure permutations (multi-indices); H(Curl) non-gradient faces still need 2x2 at p=2,3 (section 1) | excellent (Duffy sum factorization) | grows with p, fine at p ≤ 4 | non-hierarchical (degree elevation, sparse) | high; only worth it with sum factorization |
| **Symmetric FEEC bases** (Licht, FoCM 2023; Berchenko-Kogan, FoCM 2024; AFW CMAME 2009) | P_r Λ^k / P_r^- Λ^k | classify when vertex-permutation-invariant bases exist; consistent with the obstruction found here | – | – | – | theory only |

**Takeaways:**
- No published basis gives signed permutations for the rotational tet face block at p=2,3; it is impossible
  (section 1).
- The production codes that batch on GPUs (MFEM/libCEED, Basix) accept per-face 2x2 integer/orthogonal blocks in
  the gather/scatter. That is exactly what sym-proj gives OFT, while keeping OFT's gradient/rotational split and
  hierarchy, which none of the standard N1curl elements provide.

### 5.1 Context for future exploration: OFT's tet bases are Schöberl–Zaglmayr type

> **Scope note.** This subsection is background for possible *future* work on conditioning. It is **not** a
> recommendation for the current qcache / xMHD changes. It should not change their scope, priorities or
> implementation. Nothing here has been verified against the paper text or measured.

- **Provenance.** Per the OFT developers, the tet H(Curl)/H^1 bases come from Schöberl & Zaglmayr, "High order
  Nédélec elements with local complete sequence properties", COMPEL 24(2):374–384, 2005. The structure matches:
  - OFT's H(Curl) is SZ with the per-entity gradients moved into the separate grad H^1 space (NGSolve's
    `nograds`).
  - The coefficients are scaled Legendre/Jacobi polynomials. For example, p=4 face DOF 7 is `s²·P₂(t/s)` with
    `t = f1+f2−f3` and `s = f1+f2+f3`.
  - Face functions are built on the globally sorted vertex order, used asymmetrically. This is also the origin
    of the S3 defect in section 2.3.
- **Implications.**
  - Adopting "the SZ / NGSolve basis" would not improve conditioning, because OFT already uses that family. This
    fits the observation in 3.5 that OFT's functions are already well tuned.
  - The symmetric families of 3.2 and 3.7 are a symmetric modification of an SZ basis. They keep SZ's local
    complete-sequence property: face functions are corrected only by that face's gradients and lower face
    functions, and the interior only by interior gradients.
- **Possible future direction (conditioning only).** Later refinements of SZ were designed as improvements on
  the 2005 set and could be compared in the Python harness against OFT and sym-proj:
  - Zaglmayr's 2006 thesis (Jacobi-based face and cell families);
  - the Beuchler–Pillwein–Schöberl–Zaglmayr sparsity-optimized tet H(Curl) bases.

  Expected gains are constant factors, most likely at p=4 only. They are unverified; the references are from
  memory and were not checked in the survey.

## 6. Ranked recommendations

**1. Opt-in symmetric hierarchical face families** (sym-proj for H(Curl) p ≥ 2 and H^1 order ≥ 4).
- *How:*
  - Add a per-FE-object variant flag (default off) in `oft_hcurl_fem` / `oft_h1_fem`.
  - Implement the face functions in the scalar evaluators (`oft_hcurl_evalf`/`cevalf`, `oft_h1_evalf`/`gevalf`/
    `d2evalf`) and in the unrolled `eval_all*`/`ceval_all*`/`geval_all*` routines.
  - Do the same for the boundary `bh1`/`bhcurl` evaluators, which share `evalf`.
  - Coefficients come from `sym_orth.py`.
  - No change in `fem_quad_cache`: its classifier will report own-block kind 1/2 only, so `canonical_ok` holds
    everywhere.
  - Optional follow-up: specialize kind 2 for 2x2 blocks and drop the state-block tables.
  - Optional, same evaluators: replace the H(Curl) interior with the relabeling-invariant block of 3.7. This
    makes the standalone H(Curl) space independent of each cell's local vertex order, at no cache cost.
- *Benefit:*
  - kind 3 goes away (≈half of the tet H(Curl) contraction work currently takes the slow correction path);
  - the H^1 k ≥ 4 transforms shrink from 16-25-regressor blocks to own-entity blocks with one 2x2;
  - one canonical table per space, plus 6 fixed 2x2 / signed-perm maps per face;
  - canonical GEMM assembly works for H(Curl);
  - cleanest GPU form.
- *Risks:*
  - Solver accuracy for B: none (composite space unchanged).
  - Changes the discretization of standalone H(Curl) fields (J2, Taylor A) at the truncation-error level.
  - Test reference values, MG damping defaults and old restart files must be updated or versioned.
  - p=4 local conditioning is 1.2-2.2x higher with the face change alone; adding the symmetric interior (3.7)
    brings curl-curl and combined mass back to OFT's level or better.
- *Effort:* about 1 week, including evaluator variants, unit tests (span / curl / gradient-free checks like
  `verify_sym.py`, run in Fortran through the existing qcache classifier), and re-baselining `test_hcurl*` /
  MUG / Marklin.
- *Next step:* implement H(Curl) p=2,3 first (most of the kind-3 cost is there; the p=2 change is a single
  `−(1/3)∇φ^{H1}_{f,1}` term). Measure `refcel` with `fem_qcache_bench` against Phase 2.

**2. Fallback, no basis change: joint (H(Curl), grad H^1) transforms for composite fields.**
- Exact per-face blocks coupling face and interior DOFs of both spaces (3.6).
- Fixes B / curl B, but not standalone H(Curl) fields.
- Effort: medium (coupled caches in `fem_quad_cache`). Risk to accuracy: none.

**3. Later, GPU-motivated: sort + parity tet ordering (two cell classes, zero transforms).**
- Only if measurements show the per-face small maps or state lookups matter.
- High grid-code effort.
- Changes H(Curl)-only interior spans.

**4. Not recommended:**
- adopting NGSolve / ESEAS / Ainsworth–Coyle / Webb / GWP bases (same or worse orientation behavior, loss of
  the split or the hierarchy, high porting cost);
- pure UFC sort with det<0 (about 167 call sites assume det>0);
- Bernstein–Bézier, unless a sum-factorized simplex kernel is planned.

**Independent of the choice:**
- Hex needs nothing (signed permutations verified through p=4).
- The standalone H(Curl) interior's dependence on the local vertex order (2.3) is worth a note in
  `hcurl_basis.F90`, and in any future mesh-reordering work.

## Reproduction

```bash
# Python env with numpy, scipy, sympy (the CI venv plus `pip install sympy`)
cd benchmarks/xmhd_qcache/basis_study
python verify_sym.py        # spaces, transforms, conditioning -> results_sym.txt (~20 s)
python rep_analysis.py      # S3 characters of current face blocks -> results_rep.txt
python structure_checks.py  # regressor patterns, sort+parity classes, interior S4 checks -> results_structure.txt
python sym_orth.py          # prints the exact sym-proj face formulas
python interior_sym.py      # symmetric H(Curl) interior feasibility -> results_interior.txt
```

| file | content |
|---|---|
| `oft_tet_basis.py` | symbolic transcription of OFT's tet H^1 (orders 1-5) and H(Curl) (p 1-4) evaluators, with OFT's orientation conventions |
| `xform_classify.py` | port of `qcache_build_xform`, plus S3 character tools; reproduces the Fortran identity/perm/dense/none counts and regressor counts exactly |
| `sym_basis.py`, `sym_orth.py` | prototype symmetric face families (sym, sym-int / orth, H1perm, sym-proj) |
| `verify_sym.py`, `rep_analysis.py`, `structure_checks.py` | the checks quoted above |
| `interior_sym.py` | symmetric H(Curl) interior feasibility (3.7) -> `results_interior.txt` (~30 s) |

The Fortran classifier was re-run with `builds/build_release/tests/fem/test_fem_qcache` on a 2x2x2 cube. Tets
p = 1-4 match the Python port; hex p=4 is new (all signed permutations). Configuration: Ubuntu 24.04 system
libraries, GCC 13, OpenMPI, oft-build skill section 5. No OFT source or test file was modified.
