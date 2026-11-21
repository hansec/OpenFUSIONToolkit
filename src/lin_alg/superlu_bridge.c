/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* SuperLU LU/iLU solver interface for the Open FUSION Toolkit
*----------------------------------------------------------------------------*/
#ifdef HAVE_SUPERLU
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Modified for use with Open FUSION Toolkit by Chris Hansen (December 2014)
 */
#include <stdbool.h>
#include <stdint.h>
#include "slu_ddefs.h"

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
    int *etree;
    superlu_options_t options;
#if SUPERLU_VER_MAJOR >= 5
    GlobalLU_t Glu;
#endif
} factors_t;

void
oft_superlu_dgssv_c(int iopt, int n, int nnz, int nrhs,
                    double *values, int *colind, int *rowptr,
                    double *b, int ldb,
                    factors_t **f_factors,
                    int perm_spec, bool iter_refine, int *info)
/*
* This routine can be called from Fortran.
*
* iopt (input) int
*      Specifies the operation:
*      = 1, performs LU decomposition for the first time
*      = 2, performs triangular solve
*      = 3, refactor matrix with the same non-zero pattern
*      = 4, free all the storage in the end
*
* f_factors (input/output) ptr*
*      If iopt == 1, it is an output and contains the pointer pointing to
*                    the structure of the factored matrices.
*      Otherwise, it it an input.
*
*/
{
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    int panel_size, permc_spec, relax;
    trans_t trans;
    superlu_options_t options;
    #if SUPERLU_VER_MAJOR >= 5
    GlobalLU_t Glu;
    #endif
    SuperLUStat_t stat;
    factors_t *LUfactors;

    trans = TRANS; // Row storage format is passed so we solve the transposed system

    if ( iopt == 1 || iopt == 3 ) { /* LU decomposition */
        if ( !*f_factors ) {
            /* Set the default input options. */
            set_default_options(&options);
            /*
            * Get column permutation vector perm_c[], according to permc_spec:
            *   permc_spec = 0: natural ordering
            *   permc_spec = 1: minimum degree on structure of A'*A
            *   permc_spec = 2: minimum degree on structure of A'+A
            *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
            */
            if ( perm_spec == 0 ) {
		            options.ColPerm = NATURAL;
            } else if (perm_spec == 1) {
                options.ColPerm = MMD_ATA;
            } else if (perm_spec == 2) {
                options.ColPerm = MMD_AT_PLUS_A;
            } else if (perm_spec == 3) {
                options.ColPerm = COLAMD;
            } else if (perm_spec == -1) {
                // Use default
            } else {
                ABORT("Invalid column permutation type.");
            }
            if (!iter_refine) options.IterRefine = NOREFINE;
            /* */
            L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
            if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
            if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
        } else {
            LUfactors = *f_factors;
            Destroy_SuperNode_Matrix(LUfactors->L);
            SUPERLU_FREE(LUfactors->L);
            Destroy_CompCol_Matrix(LUfactors->U);
            SUPERLU_FREE(LUfactors->U);
            L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            perm_c = LUfactors->perm_c;
            perm_r = LUfactors->perm_r;
            etree = LUfactors->etree;
            options = LUfactors->options;
#if SUPERLU_VER_MAJOR >= 5
            Glu = LUfactors->Glu;
#endif
        }
        //
        options.Fact = DOFACT;
        if ( iopt == 3 ) options.Fact = SamePattern;
        /* Initialize the statistics variables. */
        StatInit(&stat);
        //
        dCreate_CompCol_Matrix(&A, n, n, nnz, values, colind, rowptr,
                               SLU_NC, SLU_D, SLU_GE);
        if ( iopt == 3 ) {
            //
            sp_preorder(&options, &A, perm_c, etree, &AC);
            panel_size = sp_ienv(1);
            relax = sp_ienv(2);
            //
#if SUPERLU_VER_MAJOR >= 5
            dgstrf(&options, &AC, relax, panel_size, etree,
                   NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, info);
#else
            dgstrf(&options, &AC, relax, panel_size, etree,
                   NULL, 0, perm_c, perm_r, L, U, &stat, info);
#endif
        } else {
            permc_spec = options.ColPerm;
            get_perm_c(permc_spec, &A, perm_c);
            //
            sp_preorder(&options, &A, perm_c, etree, &AC);
            //
            panel_size = sp_ienv(1);
            relax = sp_ienv(2);
            //
#if SUPERLU_VER_MAJOR >= 5
            dgstrf(&options, &AC, relax, panel_size, etree,
            	   NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, info);
#else
            dgstrf(&options, &AC, relax, panel_size, etree,
            	   NULL, 0, perm_c, perm_r, L, U, &stat, info);
#endif
        }
        //
//        if ( *info == 0 ) {
//            SCformat *Lstore = (SCformat *) L->Store;
//            NCformat *Ustore = (NCformat *) U->Store;
//            //printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
//            //printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
//            //printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
//            //mem_usage_t mem_usage;
//            //dQuerySpace(L, U, &mem_usage);
//            //printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
//            //       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
//        } else {
//            printf("dgstrf() error returns INFO= %d\n", *info);
//            if ( *info <= n ) { /* factorization completes */
//                dQuerySpace(L, U, &mem_usage);
//                printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
//                       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
//            }
//        }

        if ( !*f_factors ) {
            /* Save the LU factors in the factors handle */
        	LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
            LUfactors->perm_c = perm_c;
            LUfactors->perm_r = perm_r;
            LUfactors->etree = etree;
            LUfactors->options = options;
#if SUPERLU_VER_MAJOR >= 5
            LUfactors->Glu = Glu;
#endif
            *f_factors = LUfactors;
        }
        LUfactors->L = L;
        LUfactors->U = U;

        /* Free un-wanted storage */
        Destroy_SuperMatrix_Store(&A);
        Destroy_CompCol_Permuted(&AC);
        StatFree(&stat);

    } else if ( iopt == 2 ) { /* Triangular solve */
        /* Initialize the statistics variables. */
        StatInit(&stat);

        /* Extract the LU factors in the factors handle */
        LUfactors = *f_factors;
        L = LUfactors->L;
        U = LUfactors->U;
        perm_c = LUfactors->perm_c;
        perm_r = LUfactors->perm_r;

        dCreate_Dense_Matrix(&B, n, nrhs, b, ldb, SLU_DN, SLU_D, SLU_GE);

        /* Solve the system A*X=B, overwriting B with X. */
        dgstrs (trans, L, U, perm_c, perm_r, &B, &stat, info);

        Destroy_SuperMatrix_Store(&B);
        StatFree(&stat);

    } else if ( iopt == 4 ) { /* Free storage */
        /* Free the LU factors in the factors handle */
        LUfactors = *f_factors;
        SUPERLU_FREE (LUfactors->etree);
        SUPERLU_FREE (LUfactors->perm_r);
        SUPERLU_FREE (LUfactors->perm_c);
        Destroy_SuperNode_Matrix(LUfactors->L);
        Destroy_CompCol_Matrix(LUfactors->U);
        SUPERLU_FREE (LUfactors->L);
        SUPERLU_FREE (LUfactors->U);
        SUPERLU_FREE (LUfactors);
        *f_factors = NULL;
    } else {
        fprintf(stderr,"Invalid iopt=%d passed to oft_superlu_dgssv_c()\n",iopt);
        *info = -1;
    }
}

/* Free an ILU factor handle and everything it owns (L and U are only destroyed
 * once their storage has been created by dgsitrf) */
static void
oft_superlu_ilu_free(factors_t *LUfactors)
{
    if ( !LUfactors ) return;
    if ( LUfactors->L ) {
        if ( LUfactors->L->Store ) Destroy_SuperNode_Matrix(LUfactors->L);
        SUPERLU_FREE(LUfactors->L);
    }
    if ( LUfactors->U ) {
        if ( LUfactors->U->Store ) Destroy_CompCol_Matrix(LUfactors->U);
        SUPERLU_FREE(LUfactors->U);
    }
    if ( LUfactors->etree ) SUPERLU_FREE(LUfactors->etree);
    if ( LUfactors->perm_r ) SUPERLU_FREE(LUfactors->perm_r);
    if ( LUfactors->perm_c ) SUPERLU_FREE(LUfactors->perm_c);
    SUPERLU_FREE(LUfactors);
}

/* Empty SuperMatrix for dgsitrf to fill (Store stays NULL until it does) */
static SuperMatrix *
oft_superlu_ilu_new_matrix(void)
{
    SuperMatrix *M = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    if ( !M ) ABORT("Malloc fails for SuperMatrix.");
    M->Store = NULL;
    return M;
}

void
oft_superlu_dgsisx_c(int iopt, int n, int nnz, int nrhs,
                     double *values, int *colind, int *rowptr,
                     double *b, int ldb,
                     factors_t **f_factors,
                     int perm_spec, double drop_tol, double fill_factor,
                     double fill_tol, double pivot_thresh, int milu,
                     int64_t *nnz_lu, int *info)
/*
* Incomplete LU factorization (ILUTP, dgsitrf) and solve with SuperLU.
* This routine can be called from Fortran.
*
* The matrix is given in CSR format (0-based), which SuperLU sees as the CSC
* format of its transpose: A^T is factored and the transposed system is solved,
* giving the solution of A x = b. Equilibration and the MC64 row permutation of
* the dgsisx driver are not applied.
*
* iopt (input) int
*      Specifies the operation:
*      = 1, symbolic and numeric factorization (frees an existing factorization first)
*      = 2, triangular solve of nrhs right hand sides (b is ldb x nrhs)
*      = 3, refactor a matrix with the same non-zero pattern (reuses the column
*           permutation and elimination tree, options->Fact = SamePattern)
*      = 4, free all the storage
*
* perm_spec, drop_tol, fill_factor, fill_tol, pivot_thresh, milu (input, iopt 1/3)
*      Column permutation (-1 default (COLAMD), 0 natural, 1 MMD on A'*A,
*      2 MMD on A'+A, 3 COLAMD; used by iopt 1 only), options->ILU_DropTol,
*      options->ILU_FillFactor, options->ILU_FillTol, options->DiagPivotThresh
*      and options->ILU_MILU (0 SILU, 1 SMILU_1, 2 SMILU_2, 3 SMILU_3).
*
* f_factors (input/output) ptr*
*      Factor handle: output of iopt 1, input otherwise (NULL after iopt 4 and
*      after a failed factorization).
*
* nnz_lu (output, iopt 1/3) int64_t*
*      Non-zeros of L+U (diagonal counted once).
*
* info (output) int*
*      = 0: success
*      < 0: invalid argument or option (nothing done)
*      1..n (iopt 1/3): number of zero pivots replaced (factorization complete)
*      > n (iopt 1/3): memory allocation failure (handle freed)
*/
{
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    int panel_size, relax;
    trans_t trans;
    superlu_options_t options;
#if SUPERLU_VER_MAJOR >= 5
    GlobalLU_t Glu;
#endif
    SuperLUStat_t stat;
    factors_t *LUfactors;
#if SUPERLU_VER_MAJOR >= 5
    int_t finfo;
#else
    int finfo;
#endif

    trans = TRANS; // Row storage format is passed so we solve the transposed system
    *info = 0;

    if ( iopt == 1 || iopt == 3 ) { /* iLU decomposition */
        if ( drop_tol < 0.0 || fill_factor <= 0.0 || fill_tol < 0.0
             || pivot_thresh < 0.0 || pivot_thresh > 1.0 || milu < 0 || milu > 3
             || perm_spec < -1 || perm_spec > 3 ) {
            fprintf(stderr,"Invalid ILU option passed to oft_superlu_dgsisx_c()\n");
            *info = -1;
            return;
        }
        if ( iopt == 1 ) {
            /* Fresh factorization: options and permutation arrays are sized and set up anew */
            oft_superlu_ilu_free(*f_factors);
            *f_factors = NULL;
            LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
            if ( !LUfactors ) ABORT("Malloc fails for factors_t.");
            LUfactors->L = NULL;
            LUfactors->U = NULL;
            LUfactors->etree = NULL;
            LUfactors->perm_r = NULL;
            LUfactors->perm_c = NULL;
            if ( !(LUfactors->perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
            if ( !(LUfactors->perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
            if ( !(LUfactors->etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
            ilu_set_default_options(&options);
            /*
            * Column permutation vector perm_c[], according to permc_spec:
            *   permc_spec = 0: natural ordering
            *   permc_spec = 1: minimum degree on structure of A'*A
            *   permc_spec = 2: minimum degree on structure of A'+A
            *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
            */
            if ( perm_spec == 0 ) {
                options.ColPerm = NATURAL;
            } else if (perm_spec == 1) {
                options.ColPerm = MMD_ATA;
            } else if (perm_spec == 2) {
                options.ColPerm = MMD_AT_PLUS_A;
            } else if (perm_spec == 3) {
                options.ColPerm = COLAMD;
            } /* perm_spec == -1: SuperLU default (COLAMD) */
            options.Fact = DOFACT;
        } else {
            if ( !*f_factors ) {
                fprintf(stderr,"Refactor without a factorization passed to oft_superlu_dgsisx_c()\n");
                *info = -1;
                return;
            }
            LUfactors = *f_factors;
            /* Previous factors are replaced (dgsitrf allocates new storage) */
            if ( LUfactors->L->Store ) Destroy_SuperNode_Matrix(LUfactors->L);
            LUfactors->L->Store = NULL;
            if ( LUfactors->U->Store ) Destroy_CompCol_Matrix(LUfactors->U);
            LUfactors->U->Store = NULL;
            options = LUfactors->options;
            options.Fact = SamePattern;
        }
        /* Numerical options apply to every factorization */
        options.ILU_DropTol = drop_tol;
        options.ILU_FillFactor = fill_factor;
        options.ILU_FillTol = fill_tol;
        options.DiagPivotThresh = pivot_thresh;
        if ( milu == 0 ) {
            options.ILU_MILU = SILU;
        } else if ( milu == 1 ) {
            options.ILU_MILU = SMILU_1;
        } else if ( milu == 2 ) {
            options.ILU_MILU = SMILU_2;
        } else {
            options.ILU_MILU = SMILU_3;
        }
        if ( !LUfactors->L ) LUfactors->L = oft_superlu_ilu_new_matrix();
        if ( !LUfactors->U ) LUfactors->U = oft_superlu_ilu_new_matrix();
        L = LUfactors->L;
        U = LUfactors->U;
        perm_c = LUfactors->perm_c;
        perm_r = LUfactors->perm_r;
        etree = LUfactors->etree;
        /* Initialize the statistics variables. */
        StatInit(&stat);
        //
        dCreate_CompCol_Matrix(&A, n, n, nnz, values, colind, rowptr,
                               SLU_NC, SLU_D, SLU_GE);
        if ( iopt == 1 ) get_perm_c(options.ColPerm, &A, perm_c);
        /* Elimination tree and postordered perm_c (reused with SamePattern) */
        sp_preorder(&options, &A, perm_c, etree, &AC);
        panel_size = sp_ienv(1);
        relax = sp_ienv(2);
        finfo = 0;
#if SUPERLU_VER_MAJOR >= 5
        dgsitrf(&options, &AC, relax, panel_size, etree,
                NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, &finfo);
#else
        dgsitrf(&options, &AC, relax, panel_size, etree,
                NULL, 0, perm_c, perm_r, L, U, &stat, &finfo);
#endif
        *info = (int) finfo;
        /* Free un-wanted storage */
        Destroy_SuperMatrix_Store(&A);
        Destroy_CompCol_Permuted(&AC);
        StatFree(&stat);
        if ( *info < 0 || *info > n || !L->Store || !U->Store ) {
            /* Failed factorization: no usable factors, release the handle */
            if ( *info == 0 ) *info = n + 1;
            oft_superlu_ilu_free(LUfactors);
            *f_factors = NULL;
            *nnz_lu = 0;
            return;
        }
        /* Save the factors (and options for refactoring) in the factors handle */
        LUfactors->options = options;
#if SUPERLU_VER_MAJOR >= 5
        LUfactors->Glu = Glu;
#endif
        *f_factors = LUfactors;
        *nnz_lu = (int64_t) ((SCformat *) L->Store)->nnz
                + (int64_t) ((NCformat *) U->Store)->nnz - (int64_t) n;

    } else if ( iopt == 2 ) { /* Triangular solve */
        if ( !*f_factors ) {
            fprintf(stderr,"Solve without a factorization passed to oft_superlu_dgsisx_c()\n");
            *info = -1;
            return;
        }
        /* Initialize the statistics variables. */
        StatInit(&stat);

        /* Extract the LU factors in the factors handle */
        LUfactors = *f_factors;
        L = LUfactors->L;
        U = LUfactors->U;
        perm_c = LUfactors->perm_c;
        perm_r = LUfactors->perm_r;

        dCreate_Dense_Matrix(&B, n, nrhs, b, ldb, SLU_DN, SLU_D, SLU_GE);

        /* Solve the system A*X=B, overwriting B with X. */
        dgstrs(trans, L, U, perm_c, perm_r, &B, &stat, info);

        Destroy_SuperMatrix_Store(&B);
        StatFree(&stat);

    } else if ( iopt == 4 ) { /* Free storage */
        oft_superlu_ilu_free(*f_factors);
        *f_factors = NULL;
    } else {
        fprintf(stderr,"Invalid iopt=%d passed to oft_superlu_dgsisx_c()\n",iopt);
        *info = -1;
    }
}
#endif
