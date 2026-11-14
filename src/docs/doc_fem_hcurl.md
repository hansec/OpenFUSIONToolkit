H(Curl) Finite Element Representations      {#doc_fem_hcurl}
========================================

[TOC]

This page describes the H(Curl) finite element space implemented in \ref oft_hcurl_basis and the composite
H(Curl) + Grad(\f$ H^1 \f$) space built from it, including the exact basis definitions on tetrahedra and hexahedra, the
polynomial spaces they span, and known limitations. For general information on how finite element representations
are defined in OFT see \ref doc_oft_fem, and for composite representations see \ref doc_comp_fem.

\section doc_fem_hcurl_overview Overview

The order \f$ k \f$ H(Curl) space (\f$ 1 \le k \le 4 \f$, tetrahedral and hexahedral cells) is a hierarchical set of
tangentially continuous vector functions consisting of
 - the lowest order (Whitney) edge functions, and
 - only the rotational (non-gradient) higher order face and cell functions.

Higher order gradients are deliberately omitted. As a result, within each cell the curl operator is only zero for
the gradients of the lowest order (vertex) functions, which are spanned by the edge functions, and the curl of the
order \f$ k \f$ space spans the same space as the curl of the corresponding full Nedelec space. This space is used
on its own where only the curl of a field matters (eg. force-free eigenmodes in \ref doc_taylor_main "Marklin").

When the full vector field is needed (eg. the magnetic field in \ref doc_mhd_main "MUG") the gradient components are
added by pairing H(Curl) of order \f$ k \f$ with \f$ \nabla H^1 \f$ of order \f$ k+1 \f$, using
\ref oft_hcurl_basis::oft_hcurl_grad_setup "oft_hcurl_grad_setup". Fields in this composite space are represented as
\f$ \mathbf{B} = \mathbf{B}_c + \nabla \phi \f$, with \f$ \mathbf{B}_c \f$ in H(Curl) (field tag `c`) and \f$ \phi \f$
in \f$ H^1 \f$ (field tag `g`).

\section doc_fem_hcurl_dofs Degrees of freedom

Number of basis functions associated with each geometric entity (the `gstruct` array) and per cell (`nce`).

| Order \f$ k \f$ | Tet edge | Tet face | Tet cell | Tet `nce` | Hex edge | Hex face | Hex cell | Hex `nce` |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 1 | 0 | 0  | 6  | 1 | 0  | 0  | 12  |
| 2 | 1 | 2 | 0  | 14 | 1 | 3  | 5  | 35  |
| 3 | 1 | 5 | 3  | 29 | 1 | 8  | 28 | 88  |
| 4 | 1 | 9 | 11 | 53 | 1 | 15 | 81 | 183 |

On hexahedra the face and cell counts are \f$ (k-1)^2 + 2(k-1) \f$ and \f$ 2(k-1)^3 + 3(k-1)^2 \f$ respectively.

With the optional enriched hexahedral basis (see \ref doc_fem_hcurl_hex_enrich) the counts become

| Order \f$ k \f$ | Enriched hex edge | Enriched hex face | Enriched hex cell | Enriched hex `nce` |
|:---:|:---:|:---:|:---:|:---:|
| 1 | 1 | 2  | 3   | 27  |
| 2 | 1 | 7  | 26  | 80  |
| 3 | 1 | 14 | 79  | 175 |
| 4 | 1 | 23 | 174 | 324 |

with face and cell counts \f$ k^2 + 2k - 1 \f$ and \f$ 2k^3 + 3k^2 - 2 \f$ respectively.

\section doc_fem_hcurl_orient Orientation and hierarchy

All higher order functions are evaluated using the global (oriented) vertex ordering of their edge or face, so neighboring
cells evaluate identical functions on shared entities and tangential continuity is preserved. For edges the orientation is
from the lower to the higher global vertex index. For triangular faces the vertices are sorted by global index. For quadrilateral
faces the ordering starts at the lowest global index vertex and proceeds toward its lower-indexed neighbor.

The order \f$ k \f$ functions are a subset (prefix) of the order \f$ k+1 \f$ functions, so interpolation between polynomial
levels in p-multigrid is simple injection.

\section doc_fem_hcurl_tet Tetrahedral basis

Using barycentric coordinates \f$ \lambda_i \f$ and their gradients \f$ \nabla \lambda_i \f$:
 - **Edge** (1 per edge): \f$ \lambda_a \nabla \lambda_b - \lambda_b \nabla \lambda_a \f$, where \f$ a \rightarrow b \f$ is the
   oriented edge.
 - **Face** (2, 5, 9 per face for \f$ k = \f$ 2, 3, 4): polynomial combinations
   \f$ \sum_{i \in \{a,b,c\}} p_i(\lambda_a,\lambda_b,\lambda_c) \nabla \lambda_i \f$ over the oriented face vertices
   \f$ (a,b,c) \f$, which vanish tangentially on all other faces. The functions are listed explicitly in
   \ref oft_hcurl_basis::oft_hcurl_evalf "oft_hcurl_evalf", with their curls in \ref oft_hcurl_basis::oft_hcurl_cevalf "oft_hcurl_cevalf".
   For example, the two quadratic face functions are
   \f$ -2 (\lambda_b \lambda_c \nabla \lambda_a + \lambda_a \lambda_c \nabla \lambda_b - \lambda_a \lambda_b \nabla \lambda_c) \f$ and
   \f$ \lambda_b \lambda_c \nabla \lambda_a - \lambda_a \lambda_c \nabla \lambda_b \f$.
 - **Cell** (0, 3, 11 per cell for \f$ k = \f$ 2, 3, 4): interior bubbles listed in \ref oft_hcurl_basis::oft_hcurl_evalc "oft_hcurl_evalc",
   with their curls in \ref oft_hcurl_basis::oft_hcurl_cevalc "oft_hcurl_cevalc".

\section doc_fem_hcurl_hex Hexahedral basis

The hexahedral basis is written in terms of six face coordinates built from the logical coordinates \f$ x_i \in [0,1] \f$
\f[ \mu = (1-x_3,\ 1-x_2,\ x_1,\ x_2,\ 1-x_1,\ x_3), \f]
where \f$ \mu_i \f$ is one on face \f$ i \f$ and zero on the opposing face (opposing pairs 1-6, 2-4 and 3-5). One dimensional
bubbles are built from \f$ L_n(t) = t\, J^{(1,1)}_{n-1}(2t-1) \f$ (\ref oft_lag_poly::hpoly_bary "hpoly_bary") as
\f$ L_n(\mu_i)\, \mu_{i'} \f$, where \f$ \mu_{i'} = 1 - \mu_i \f$ is the coordinate of the opposing face.

**Edge** (1 per edge): \f$ \mu_a \mu_b \nabla \mu_t \f$, where \f$ \mu_a = \mu_b = 1 \f$ on the edge (`hex_bary_efcoords`)
and \f$ \mu_t \f$ increases along the oriented edge (`hex_bary_ecoords`).

**Face**: let \f$ (c_1,c_2,c_3,c_4) \f$ be the oriented face coordinates, with \f$ c_3 = 1-c_1 \f$, \f$ c_4 = 1-c_2 \f$ and
\f$ c_1 = c_2 = 1 \f$ at the lowest global index vertex, and define \f$ u_i = L_i(c_1)\, c_3 \f$ and \f$ v_j = L_j(c_2)\, c_4 \f$.
Each entry \f$ (i,j,\mathrm{form}) \f$ of the index list (`indsf`, built by \ref oft_hcurl_basis::hcurl_2d_grid "hcurl_2d_grid"),
with \f$ i,j \ge 1 \f$, defines one function. Each is multiplied by the coordinate \f$ \mu_f \f$ of face \f$ f \f$ to extend it into the cell.

| form | function |
|:---:|:---|
| 0 | \f$ \mu_f \left( v_j \nabla u_i - u_i \nabla v_j \right) \f$ |
| 1 | \f$ \mu_f\, v_j \nabla c_1 \f$ |
| 2 | \f$ \mu_f\, u_i \nabla c_2 \f$ |

**Cell**: define \f$ U_i = L_i(\mu_1)\, \mu_6 \f$, \f$ V_j = L_j(\mu_2)\, \mu_4 \f$ and \f$ W_l = L_l(\mu_3)\, \mu_5 \f$. Each entry
\f$ (i,j,l,\mathrm{form}) \f$ of the index list (`indsc`, built by \ref oft_hcurl_basis::hcurl_3d_grid "hcurl_3d_grid") defines one function.

| form | function |
|:---:|:---|
| 0 | \f$ \nabla U_i\, V_j W_l - U_i \nabla V_j\, W_l + U_i V_j \nabla W_l \f$ |
| 1 | \f$ \nabla U_i\, V_j W_l - U_i \nabla V_j\, W_l - U_i V_j \nabla W_l \f$ |
| 2 | \f$ V_j W_l \nabla \mu_1 \f$ |
| 3 | \f$ U_i W_l \nabla \mu_2 \f$ |
| 4 | \f$ U_i V_j \nabla \mu_3 \f$ |

Forms 0 and 1 are the two non-gradient sign combinations of \f$ \nabla (U_i V_j W_l) \f$; the gradient itself is supplied by \f$ H^1 \f$.

\subsection doc_fem_hcurl_hex_enrich Enriched hexahedral basis

An enriched hexahedral basis can be selected with `hex_enrich=.TRUE.` in \ref oft_hcurl_basis::oft_hcurl_setup "oft_hcurl_setup"
(default `.FALSE.`). It restores optimal convergence on hexahedral meshes whose cells remain non-affine at all resolutions
(see \ref doc_fem_hcurl_limits_hex). It has no effect on tetrahedral meshes.

The order \f$ k \f$ enriched space uses the face and cell functions of the standard order \f$ k+1 \f$ space (same functions,
index lists and orientation as above), except for the last entry of the face list, \f$ (i,j,\mathrm{form}) = (k,k,0) \f$, and the last
two entries of the cell list, \f$ (i,j,l,\mathrm{form}) = (k,k,k,0) \f$ and \f$ (k,k,k,1) \f$. Equivalently, the standard order
\f$ k \f$ space is enlarged by the order \f$ k+1 \f$ face and cell functions with bubble index \f$ \max(i,j) = k \f$ or
\f$ \max(i,j,l) = k \f$, minus these three. The composite space still pairs it with \f$ H^1 \f$ of order \f$ k+1 \f$, so it is a
subspace of the first kind Nedelec space \f$ N_{k+1} = Q_{k,k+1,k+1} \times Q_{k+1,k,k+1} \times Q_{k+1,k+1,k} \f$ (the composite
order \f$ k+1 \f$ H(Curl) plus \f$ \nabla Q_{k+1} \f$ space). As the enriched functions are a subset of an existing tangentially
continuous basis evaluated in the global orientation, tangential continuity is unchanged. The order \f$ k \f$ standard functions are a
prefix of the order \f$ k \f$ enriched functions, which are in turn a prefix of the order \f$ k+1 \f$ enriched functions, so p-multigrid
interpolation remains injection. The vertex functions of \f$ H^1 \f$ are the only redundant functions in the composite space, as for
the standard basis. As the enriched functions have one higher polynomial degree, the quadrature order used by the H(Curl) space (and
the composite operators) is increased from \f$ 2k+1 \f$ to \f$ 2k+3 \f$.

**Why this enrichment.** Let \f$ F \f$ be the trilinear map of a cell. H(Curl) functions map covariantly, \f$ \hat{\mathbf{u}} = DF^T (\mathbf{u} \circ F) \f$,
and curls map with the contravariant Piola transform, \f$ \hat{\nabla} \times \hat{\mathbf{u}} = \mathrm{adj}(DF) \left( (\nabla \times \mathbf{u}) \circ F \right) \f$.
By the same argument as Arnold, Boffi and Falk (quadrilaterals) and Falk, Gatto and Monk (hexahedra, "Hexahedral H(div) and H(curl) finite elements",
ESAIM: M2AN 45, 2011), approximation of order \f$ k+1 \f$ for the field requires the reference space to contain
\f$ S_k = \{ DF^T (\mathbf{p} \circ F) : \mathbf{p} \in [P_k]^3 \} \f$ for every trilinear \f$ F \f$, and approximation of order \f$ k \f$ for the
curl requires the reference curls to contain \f$ \hat{\nabla} \times S_k \f$. Since \f$ DF \f$ has entries in \f$ Q_1 \f$, \f$ S_k \subset N_{k+1} \f$, and the
gradient part of \f$ S_k \f$ is in \f$ \nabla Q_{k+1} \f$. For a subspace \f$ V \f$ of \f$ N_{k+1} \f$ paired with \f$ \nabla Q_{k+1} \f$ both conditions
therefore reduce to \f$ \hat{\nabla} \times V \supseteq \hat{\nabla} \times S_k \f$, so the same enrichment serves the H(Curl) space on its own
(curl) and the composite space (field and curl). Computing the span of \f$ S_k \f$ over random trilinear maps gives, per cell,

| Order \f$ k \f$ | Standard composite | Span of \f$ S_k \f$ | \f$ N_{k+1} \f$ | Enriched composite |
|:---:|:---:|:---:|:---:|:---:|
| 1 | 31  | 46  | 54  | 46  |
| 2 | 91  | 133 | 144 | 136 |
| 3 | 205 | 286 | 300 | 292 |
| 4 | 391 | 521 | 540 | 532 |

where all counts are dimensions of the spanned spaces (redundant functions counted once). The span of \f$ S_k \f$ contains the standard composite space, so the
missing content lies entirely in \f$ N_{k+1} \f$. On a face the normal component of \f$ \hat{\nabla} \times S_k \f$ is
\f$ (\partial_1 F \times \partial_2 F) \cdot (\mathbf{q} \circ F) \f$ with \f$ \mathbf{q} \in [P_{k-1}]^3 \f$. The cross product has no
\f$ x_1 x_2 \f$ term, so this lies in \f$ Q_{k,k} \f$ without the \f$ x_1^k x_2^k \f$ term. The only face function that is not needed is
therefore the one whose curl has degree \f$ (k,k) \f$, entry \f$ (k,k,0) \f$. Similarly, the cofactors of \f$ DF \f$ have no
\f$ x_1^2 x_2 x_3 \f$-type terms, which removes the two cell functions of top degree. For \f$ k=1 \f$ the enriched composite space equals
the span of \f$ S_1 \f$ exactly. For \f$ k \ge 2 \f$ a few more cell functions (3, 6 and 11 for \f$ k = \f$ 2, 3, 4) are not needed, but
they cannot be removed by dropping basis functions, so they are kept. This analysis is for trilinear cells. Higher order (eg. quadratic) cell
maps whose deviation from a trilinear map vanishes under refinement, such as curved boundary cells of a smooth geometry, are expected to
behave the same way, but this has not been tested; cells with \f$ \mathcal{O}(1) \f$ curvature at every resolution need even more content.

**Cost.** Per cell, counting shared entities once on a large mesh (about 3 edges, 3 faces and 1 vertex per cell), the number of
unknowns is

| Order \f$ k \f$ | Standard composite | Enriched composite | Standard composite of order \f$ k+1 \f$ |
|:---:|:---:|:---:|:---:|
| 1 | 10  | 19  | 43  |
| 2 | 43  | 76  | 118 |
| 3 | 118 | 187 | 253 |
| 4 | 253 | 370 | -   |

The enriched order \f$ k \f$ space gives on distorted meshes the rates the standard order \f$ k \f$ space gives on regular meshes, at
about 1.5 to 1.9 times the unknowns and with \f$ (k+2)^3/(k+1)^3 \f$ times the quadrature points. Raising the standard order to
\f$ k+1 \f$ to recover the same rates costs more and still loses one order on distorted meshes. On regular meshes the enriched space
contains the standard space, so errors are never larger, and its curl converges one order faster (as \f$ \mathcal{O}(h^{k+1}) \f$).

**Compatibility.** The DOF layout (`gstruct`, `nce`) of the enriched space differs from the standard one, so vectors and restart files
cannot be exchanged between the two; loading a restart file written with the other option fails with a size mismatch. Any code that
reconstructs a saved field (eg. post-processing tools) must construct the space with the same option. Solver iteration counts also
change, and multigrid smoother settings tuned for the standard space may need to be adjusted.

\section doc_fem_hcurl_spaces Spanned polynomial spaces

Here \f$ P_k \f$ denotes polynomials of total degree \f$ \le k \f$ and \f$ Q_k \f$ polynomials of degree \f$ \le k \f$ in
each variable separately. All statements are per cell, on the reference element.

**Tetrahedra**
 - H(Curl) of order \f$ k \f$ plus \f$ \nabla P_{k+1} \f$ (the composite space) is exactly \f$ [P_k]^3 \f$, the second kind Nedelec space.
 - H(Curl) of order \f$ k \f$ plus \f$ \nabla P_k \f$ has the dimension of the first kind Nedelec space, but for \f$ k \ge 2 \f$ it is
   not equal to it, as some face and cell functions are not first kind Nedelec functions.
 - The span of the face and cell functions alone depends on the vertex ordering of the cell; only their span modulo gradients is
   independent of it. This is invisible in the composite space and in curl-only operators, but means the H(Curl) space used on
   its own (eg. with a mass matrix) differs slightly from cell to cell.

**Hexahedra**
 - H(Curl) of order \f$ k \f$ plus \f$ \nabla Q_k \f$ is exactly the first kind Nedelec space
   \f$ Q_{k-1,k,k} \times Q_{k,k-1,k} \times Q_{k,k,k-1} \f$, independent of orientation.
 - The composite space (H(Curl) plus \f$ \nabla Q_{k+1} \f$) contains \f$ [P_k]^3 \f$ but not \f$ [Q_k]^3 \f$. For example
   \f$ (x y z)^k \hat{x} \f$ cannot be represented. The composite space is not smaller than \f$ [Q_k]^3 \f$ (91 vs 81 functions per
   cell for \f$ k=2 \f$); it has different content.

**Hexahedra, enriched basis**
 - H(Curl) of order \f$ k \f$ plus \f$ \nabla Q_{k+1} \f$ is a subspace of \f$ N_{k+1} \f$ of codimension 8 that contains \f$ [P_k]^3 \f$ and,
   on every trilinear cell, the covariant pullback of \f$ [P_k]^3 \f$ (see \ref doc_fem_hcurl_hex_enrich).

**Composite redundancy**: on both cell types the gradients of the \f$ H^1 \f$ vertex functions are also spanned by the H(Curl) edge
functions, so the composite space contains exactly one redundant function per mesh vertex. These must be removed, see
\ref doc_fem_hcurl_comp_use.

\section doc_fem_hcurl_comp_use Using the composite space

The vertex DOFs of the \f$ H^1 \f$ component are redundant (see above) and must be eliminated:
 - Apply \ref oft_hcurl_grad_operators::oft_hcurl_grad_gzerop "oft_hcurl_grad_gzerop" to right hand side vectors and initial guesses.
 - \ref oft_hcurl_grad_operators::hcurl_grad_getmop "hcurl_grad_getmop" with `bc="none"` eliminates these DOFs by replacing the
   vertex rows of the \f$ H^1 \f$ block with identity rows, but leaves the corresponding columns in place. The resulting matrix is
   therefore **not symmetric**, and solves with symmetric solvers (eg. CG) are only valid if the vertex entries of the right hand side
   and initial guess are zero (ie. after applying `oft_hcurl_grad_gzerop`). With `bc="zerob"` the boundary DOFs of both components are
   also eliminated in the same way, matching \ref oft_hcurl_grad_operators::oft_hcurl_grad_zerob "oft_hcurl_grad_zerob". Any `bc` value
   other than `"none"` or `"zerob"` returns the unmodified, singular matrix.

\section doc_fem_hcurl_limits Limitations

\subsection doc_fem_hcurl_limits_hex Approximation on hexahedral meshes

Because the composite space on hexahedra contains \f$ [P_k]^3 \f$, it has the same convergence rates as on tetrahedra when cells are
affine or asymptotically affine: \f$ \mathcal{O}(h^{k+1}) \f$ for the field and \f$ \mathcal{O}(h^k) \f$ for its curl. Not containing
\f$ [Q_k]^3 \f$ only affects the error constant and means tensor product polynomial fields are not reproduced exactly, so manufactured
solutions for verification on hexahedra should use \f$ P_k \f$ or smooth fields.

On hexahedral meshes whose cells remain non-affine at all resolutions, convergence is reduced by roughly one order for both the field and
its curl, and for \f$ k=1 \f$ the curl does not converge. This is a general property of mapped hexahedral H(Curl) elements: on a trilinear
cell a physical \f$ P_k \f$ field maps to a reference field of degree \f$ k+1 \f$ in some directions, which neither this space nor
\f$ [Q_k]^3 \f$ contains. Meshes produced by uniform refinement of a coarse hexahedral mesh, including all meshes generated with the
`cube` and `sphere` mesh types, become asymptotically affine and are not affected. Fine unstructured hexahedral meshes used without
refinement may be.

The optional enriched hexahedral basis (`hex_enrich`, see \ref doc_fem_hcurl_hex_enrich) restores the rates of affine meshes on
such meshes, at the cost of more unknowns. It is off by default.

The table below shows observed convergence rates (from the two finest meshes, up to 4096 hexahedra or 6144 tetrahedra, and up to 512
hexahedra for \f$ k=4 \f$) for the L2 projection
of the smooth field \f$ \mathbf{B} = \left( \sin(2y+z)\cos x,\ e^{x/2}\sin z,\ \cos(x+y)\, e^{z/3} \right) \f$ onto the composite space
on a unit cube. "Smoothly mapped" applies a smooth global distortion to a regular mesh (curved, non-affine cells that become affine under
refinement). "Distorted" randomly moves each interior vertex by up to \f$ 0.25 h \f$, so cells remain non-affine at every resolution.
Each entry gives the rate for \f$ \mathbf{B} \f$ / \f$ \nabla \times \mathbf{B} \f$.

| Mesh | \f$ k=1 \f$ | \f$ k=2 \f$ | \f$ k=3 \f$ | \f$ k=4 \f$ |
|:---|:---:|:---:|:---:|:---:|
| Tetrahedra                          | 2.07 / 1.06  | 2.99 / 2.05 | 4.15 / 3.14 | -           |
| Regular hexahedra                   | 2.00 / 1.00  | 2.89 / 2.04 | 4.00 / 3.00 | 4.85 / 4.07 |
| Smoothly mapped hexahedra           | 1.99 / 1.00  | 2.88 / 2.04 | 4.00 / 2.99 | 4.77 / 4.04 |
| Distorted hexahedra                 | 0.89 / -0.09 | 1.18 / 0.87 | 2.76 / 1.83 | 3.77 / 2.79 |
| Regular hexahedra, enriched         | 2.01 / 2.02  | 3.02 / 3.00 | 4.02 / 4.07 | 5.38 / 5.00 |
| Smoothly mapped hexahedra, enriched | 2.00 / 1.27  | 3.00 / 2.86 | 3.99 / 3.15 | 4.93 / 4.81 |
| Distorted hexahedra, enriched       | 1.90 / 0.87  | 2.81 / 0.11 (a) | 3.88 / 2.93 | 4.85 / 4.07 |

With the enriched basis the distorted-mesh rates are \f$ k+1 \f$ for the field and \f$ k \f$ for the curl, as on affine meshes. On regular
meshes the enriched curl converges one order faster; on smoothly mapped meshes this extra order is lost gradually as the mapping is
resolved, so the measured curl rates there lie between \f$ k \f$ and \f$ k+1 \f$. Absolute errors with the enriched basis are lower in
all cases (on distorted meshes by 1 to 2 orders of magnitude at equal mesh size).

(a) The curl error for this case stagnates between the two finest meshes (\f$ 1.50 \times 10^{-4} \f$ and \f$ 1.39 \times 10^{-4} \f$)
after converging at rate 2.00 between the previous pair; tightening the solver tolerance does not change it. With a milder distortion
(vertices moved by up to \f$ 0.15 h \f$) the same case converges at 2.86 / 1.94, so this appears to come from a few strongly
distorted cells of the random perturbation on the finest mesh rather than from the space. The L2 projection does not control the
curl directly, so it is sensitive to such cells.

\subsection doc_fem_hcurl_limits_div Divergence of composite fields

Higher order (\f$ k \ge 2 \f$) H(Curl) functions are not divergence free. The divergence interpolator
\ref oft_hcurl_grad_operators::oft_hcurl_grad_dinterp "oft_hcurl_grad_dinterp" (used for the \f$ \nabla \cdot \mathbf{B} \f$ diagnostic
in MUG) only includes the divergence of the \f$ \nabla \phi \f$ component, so for these orders it does not return the full divergence
of the field. The weak divergence operator \ref oft_hcurl_grad_operators::hcurl_grad_div "hcurl_grad_div" includes both components.

\subsection doc_fem_hcurl_limits_other Other limitations

 - Geometric (h-refinement) multigrid interpolation is only implemented for tetrahedral meshes; on hexahedral meshes only polynomial
   (p-refinement) levels may be used.
 - The surface (2D) H(Curl) space can be constructed on triangular and quadrilateral surface meshes, but basis evaluation
   (\ref oft_hcurl_basis::oft_bhcurl_eval "oft_bhcurl_eval", \ref oft_hcurl_basis::oft_bhcurl_ceval "oft_bhcurl_ceval") is only
   implemented for triangles.
 - The maximum supported order is \f$ k=4 \f$.
