/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* Apple Accelerate sparse direct (LU) solver interface for the Open FUSION Toolkit
*
* The Accelerate sparse headers use clang-only language extensions, so this file
* is always compiled with Apple clang (see lin_alg/CMakeLists.txt), independent
* of the compiler used for the rest of OFT.
*----------------------------------------------------------------------------*/
#ifdef HAVE_ACCELERATE_SPARSE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/sysctl.h>
#include <pthread.h>
#include <Accelerate/Accelerate.h>

/* The macOS version needed for newer options (LU >= 15.5, MT-METIS and Hungarian
   scaling >= 26.0) is checked at runtime with oft_accelerate_os_atleast, so a clear
   error is reported if they are unavailable. The Accelerate headers themselves use
   __builtin_available, which requires the clang runtime library; it is pre-linked into
   this object by lin_alg/CMakeLists.txt so other compilers can link OFT. */
#pragma clang diagnostic ignored "-Wunguarded-availability-new"
#pragma clang diagnostic ignored "-Wunguarded-availability"

typedef struct {
  int n;                                    /* Number of rows */
  int nnz;                                  /* Number of non-zeros */
  long *col_starts;                         /* CSR row pointer as CSC column starts (0-based) */
  int *row_indices;                         /* CSR column indices as CSC row indices (0-based) */
  bool have_symbolic;                       /* Symbolic factorization is valid? */
  bool have_numeric;                        /* Numeric factorization is valid? */
  SparseOpaqueSymbolicFactorization symbolic;
  SparseOpaqueFactorization_Double numeric;
} oft_accelerate_factors_t;

/* Serializes symbolic factorizations: concurrent LU symbolic factorizations of
   independent matrices from different threads crash with the default (MT-METIS)
   ordering (observed on macOS 27.2). Numeric factorization and solves are not affected. */
static pthread_mutex_t oft_accelerate_symbolic_lock = PTHREAD_MUTEX_INITIALIZER;

/* Is the running macOS version at least major.minor? */
static bool oft_accelerate_os_atleast(int major, int minor) {
  char version[32];
  size_t len = sizeof(version);
  int vmajor = 0, vminor = 0;
  if (sysctlbyname("kern.osproductversion", version, &len, NULL, 0) != 0) return false;
  if (sscanf(version, "%d.%d", &vmajor, &vminor) < 1) return false;
  return (vmajor > major) || (vmajor == major && vminor >= minor);
}

/* Report parameter errors instead of trapping (default behavior) */
static void oft_accelerate_report_error(const char *message) {
  fprintf(stderr, "Accelerate Sparse: %s\n", message);
}

/* Matrix A (OFT CSR storage) as a transposed view of the CSC matrix A^T */
static SparseMatrix_Double oft_accelerate_matrix(oft_accelerate_factors_t *f, double *values) {
  SparseMatrix_Double A;
  A.structure.rowCount = f->n;
  A.structure.columnCount = f->n;
  A.structure.columnStarts = f->col_starts;
  A.structure.rowIndices = f->row_indices;
  A.structure.attributes = (SparseAttributes_t){ .transpose = true, .kind = SparseOrdinary };
  A.structure.blockSize = 1;
  A.data = values;
  return A;
}

/* Does the stored graph differ from the supplied CSR graph (1-based)? */
static bool oft_accelerate_graph_changed(oft_accelerate_factors_t *f, int n, int nnz, int *colind, int *rowptr) {
  if (f->n != n || f->nnz != nnz) return true;
  for (int i = 0; i <= n; i++) if (f->col_starts[i] != (long) rowptr[i] - 1) return true;
  for (int i = 0; i < nnz; i++) if (f->row_indices[i] != colind[i] - 1) return true;
  return false;
}

static void oft_accelerate_free_factors(oft_accelerate_factors_t *f) {
  if (f->have_numeric) SparseCleanup(f->numeric);
  if (f->have_symbolic) SparseCleanup(f->symbolic);
  f->have_numeric = false;
  f->have_symbolic = false;
}

/* Map OFT option values to Accelerate enums (-1 -> default) */
static int oft_accelerate_options(int fact_type, int order, int scaling,
                                  SparseFactorization_t *type, SparseOrder_t *order_method,
                                  SparseScaling_t *scaling_method) {
  switch (fact_type) {
    case -1: case 0: *type = SparseFactorizationLU; break;
    case 1: *type = SparseFactorizationLUUnpivoted; break;
    case 2: *type = SparseFactorizationLUSPP; break;
    case 3: *type = SparseFactorizationLUTPP; break;
    default:
      fprintf(stderr, "Accelerate Sparse: Unknown factorization type %d\n", fact_type);
      return -10;
  }
  switch (order) {
    case -1: case 0: *order_method = SparseOrderDefault; break;
    case 2: *order_method = SparseOrderAMD; break;
    case 3: *order_method = SparseOrderMetis; break;
    case 5:
      if (oft_accelerate_os_atleast(26, 0)) {
        *order_method = SparseOrderMTMetis;
      } else {
        fprintf(stderr, "Accelerate Sparse: MT-METIS ordering requires macOS 26.0 or later\n");
        return -11;
      }
      break;
    default:
      fprintf(stderr, "Accelerate Sparse: Unknown ordering %d\n", order);
      return -11;
  }
  switch (scaling) {
    case -1: case 0: *scaling_method = SparseScalingDefault; break;
    case 2: *scaling_method = SparseScalingEquilibriationInf; break;
    case 3:
      if (oft_accelerate_os_atleast(26, 0)) {
        *scaling_method = SparseScalingHungarianScalingOnly;
      } else {
        fprintf(stderr, "Accelerate Sparse: Hungarian scaling requires macOS 26.0 or later\n");
        return -12;
      }
      break;
    default:
      fprintf(stderr, "Accelerate Sparse: Unknown scaling %d\n", scaling);
      return -12;
  }
  return 0;
}

/*
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs full LU decomposition (symbolic and numeric)
 *      = 2, performs triangular solve
 *      = 3, performs numerical factorization only (same non-zero pattern), a full
 *           factorization is performed if the non-zero pattern has changed
 *      = 4, free all the storage in the end
 *
 * rowptr, colind (input) int*
 *      Matrix graph in CSR format with 1-based (Fortran) indices
 *
 * f_factors (input/output) oft_accelerate_factors_t**
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 * fact_type, order, scaling (input) int
 *      Factorization type, fill-reducing ordering, and scaling (-1 -> default)
 *
 * info (output) int
 *      0 on success, Accelerate SparseStatus_t or negative bridge error otherwise
 */
void oft_accelerate_dgssv_c(int iopt, int n, int nnz, int nrhs,
                            double *values, int *colind, int *rowptr,
                            double *b, int ldb, oft_accelerate_factors_t **f_factors,
                            int fact_type, int order, int scaling, int *info)
{
  oft_accelerate_factors_t *f = *f_factors;
  *info = 0;

  if ( iopt == 1 || iopt == 3 ) { /* Factorization */

    if ( !oft_accelerate_os_atleast(15, 5) ) {
      fprintf(stderr, "Accelerate Sparse: LU factorization requires macOS 15.5 or later\n");
      *info = -1;
      return;
    }
    if ( iopt == 3 ) {
      /* Full factorization if no previous factorization or the pattern has changed */
      if ( f == NULL || !f->have_numeric || oft_accelerate_graph_changed(f, n, nnz, colind, rowptr) ) iopt = 1;
    }

    if ( iopt == 1 ) {
      SparseFactorization_t type;
      SparseOrder_t order_method;
      SparseScaling_t scaling_method;
      *info = oft_accelerate_options(fact_type, order, scaling, &type, &order_method, &scaling_method);
      if ( *info != 0 ) return;
      /* Create or reset storage */
      if ( f == NULL ) {
        f = (oft_accelerate_factors_t *) calloc(1, sizeof(oft_accelerate_factors_t));
        *f_factors = f;
      } else {
        oft_accelerate_free_factors(f);
        free(f->col_starts);
        free(f->row_indices);
      }
      /* Store 0-based copy of graph (CSR of A == CSC of A^T) */
      f->n = n;
      f->nnz = nnz;
      f->col_starts = (long *) malloc((n+1) * sizeof(long));
      f->row_indices = (int *) malloc(nnz * sizeof(int));
      for (int i = 0; i <= n; i++) f->col_starts[i] = (long) rowptr[i] - 1;
      for (int i = 0; i < nnz; i++) f->row_indices[i] = colind[i] - 1;
      SparseMatrix_Double A = oft_accelerate_matrix(f, values);
      /* Symbolic factorization */
      SparseSymbolicFactorOptions sfoptions = {
        .control = SparseDefaultControl,
        .orderMethod = order_method,
        .order = NULL,
        .ignoreRowsAndColumns = NULL,
        .malloc = malloc,
        .free = free,
        .reportError = oft_accelerate_report_error
      };
      pthread_mutex_lock(&oft_accelerate_symbolic_lock);
      f->symbolic = SparseFactor(type, A.structure, sfoptions);
      pthread_mutex_unlock(&oft_accelerate_symbolic_lock);
      if ( f->symbolic.status != SparseStatusOK ) {
        fprintf(stderr, "Accelerate Sparse: Symbolic factorization failed (%d)\n", (int) f->symbolic.status);
        *info = (int) f->symbolic.status;
        SparseCleanup(f->symbolic);
        return;
      }
      f->have_symbolic = true;
      /* Numeric factorization */
      SparseNumericFactorOptions nfoptions = {
        .control = SparseDefaultControl,
        .scalingMethod = scaling_method,
        .scaling = NULL,
        .pivotTolerance = 0.01, /* Accelerate defaults (SparseNumericFactorOptionsDefault) */
        .zeroTolerance = 1e-4 * __DBL_EPSILON__
      };
      f->numeric = SparseFactor(f->symbolic, A, nfoptions);
      f->have_numeric = true;
    } else {
      /* Refactor with same pattern and options */
      SparseMatrix_Double A = oft_accelerate_matrix(f, values);
      SparseRefactor(A, &f->numeric);
    }
    if ( f->numeric.status != SparseStatusOK ) {
      if ( f->numeric.status == SparseMatrixIsSingular ) {
        fprintf(stderr, "Accelerate Sparse: Singular matrix detected in factorization\n");
      } else {
        fprintf(stderr, "Accelerate Sparse: Numeric factorization failed (%d)\n", (int) f->numeric.status);
      }
      *info = (int) f->numeric.status;
    }

  } else if ( iopt == 2 ) { /* Perform solve */

    if ( f == NULL || !f->have_numeric || f->numeric.status != SparseStatusOK ) {
      fprintf(stderr, "Accelerate Sparse: Solve requested without valid factorization\n");
      *info = -1;
      return;
    }
    /* Solve one RHS at a time: with a transposed matrix (as used here) Accelerate
       returns incorrect results for all but the first column of a multi-RHS solve
       (observed on macOS 27.2) */
    for (int j = 0; j < nrhs; j++) {
      DenseVector_Double xb = { .count = f->n, .data = b + (size_t) j * ldb };
      SparseSolve(f->numeric, xb);
    }

  } else if ( iopt == 4 ) { /* free storage */

    if ( f != NULL ) {
      oft_accelerate_free_factors(f);
      free(f->col_starts);
      free(f->row_indices);
      free(f);
    }
    *f_factors = NULL;

  }
}
#endif
