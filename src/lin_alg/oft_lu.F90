!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_lu.F90
!
!> Native module for direct solves.
!! Available interfaces:
!! - LAPACK (dgetrf, dgetrs)
!! - SuperLU
!! - MUMPS
!! - PARDISO (MKL interface)
!!
!! @note Currently supports local solves only.
!!
!! @authors Chris Hansen
!! @date May 2014
!! @ingroup doxy_oft_lin_alg
!------------------------------------------------------------------------------
MODULE oft_lu
USE, INTRINSIC :: iso_c_binding, ONLY: c_bool, c_int, c_double, c_ptr, c_null_ptr, c_f_pointer, c_loc, &
  c_associated
USE oft_local
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_graph, OFT_ACCESS_READ, OFT_ACCESS_WRITE
USE oft_solver_base, ONLY: oft_solver
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast, native_pattern_diag
IMPLICIT NONE
#include "local.h"
#ifdef HAVE_MUMPS
#include "dmumps_struc.h"
#endif
#ifdef HAVE_MKL_PARDISO
#define HAVE_MKL
#endif
#ifdef HAVE_MKL
#define DEF_LU_PACK "mkl"
#define DEF_ILU_PACK "mkl"
#else
#define DEF_ILU_PACK "native"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_UMFPACK)
#define DEF_LU_PACK "umfpack"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_MUMPS)
#define DEF_LU_PACK "mumps"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_SUPERLU)
#define DEF_LU_PACK "super"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_SUPERLU_DIST)
#define DEF_LU_PACK "superd"
#endif
#if !defined(DEF_LU_PACK)
#define DEF_LU_PACK "none"
#endif
PRIVATE
PUBLIC lapack_matinv, lapack_cholesky
!------------------------------------------------------------------------------
!> Native ILU(0) preconditioner information object
!!
!! Storage is allocatable so that copies of a solver (eg. `ALLOCATE(..., source=)` in
!! block-Jacobi) are independent deep copies.
!------------------------------------------------------------------------------
TYPE :: native_ilu_struc
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ju !< Start of the U part of each row in `alu`/`jlu` [nr]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: jlu !< Row pointer (1:nr+1) and column indices (nr+2:) of the factors [nnz+1]
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: alu !< Inverted diagonal (1:nr) and off-diagonal (nr+2:) factor values [nnz+1]
END TYPE native_ilu_struc
!------------------------------------------------------------------------------
!> Matrix structure a factorization was set up for (symbolic setup is redone when it changes)
!!
!! Identified by the matrix's non-zero pattern (@ref oft_native_pattern::id, which is
!! never reused) when the matrix has one, otherwise by the row pointer and column
!! list arrays it uses, together with the number of rows and non-zeros. The pointers
!! are references only and are never deallocated through this object.
!------------------------------------------------------------------------------
TYPE :: lu_symbolic_key
  INTEGER(i8) :: pattern_id = -1 !< Id of the matrix pattern (0 if the matrix has none, -1 if not set up)
  INTEGER(i4) :: nr = -1 !< Number of rows
  INTEGER(i4) :: nnz = -1 !< Number of non-zeros
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kr => NULL() !< Row pointer used (matrices without a pattern)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lc => NULL() !< Column list used (matrices without a pattern)
END TYPE lu_symbolic_key
!------------------------------------------------------------------------------
!> SuperLU solver information
!------------------------------------------------------------------------------
TYPE :: superlu_struc
  INTEGER(i4) :: col_perm = -1 !< Column permutation to use (-1 for default)
#ifdef OFT_MPI_F08
  TYPE(mpi_comm) :: comm = MPI_COMM_NULL !< Private MPI communicator for SuperLU-DIST
#else
  INTEGER(i4) :: comm = MPI_COMM_NULL !< Private MPI communicator for SuperLU-DIST
#endif
  TYPE(C_PTR) :: f_factors = C_NULL_PTR !< SuperLU factors pointer
  TYPE(C_PTR) :: grid_handle = C_NULL_PTR !< SuperLU-DIST processor grid pointer
  INTEGER(i4), POINTER, DIMENSION(:) :: kr => NULL() !< Row pointer (0-based indices)
  INTEGER(i4), POINTER, DIMENSION(:) :: lc => NULL() !< Col pointer (0-based indices)
  INTEGER(i4), POINTER, DIMENSION(:) :: csc_map => NULL() !< Col pointer (0-based indices)
  REAL(r8), POINTER, DIMENSION(:) :: csc_vals => NULL() !< Col pointer (0-based indices)
  INTEGER(i4), POINTER, DIMENSION(:) :: lc_orig => NULL() !< Unpermuted copy of `lc` (SuperLU-DIST only)
END TYPE superlu_struc
!------------------------------------------------------------------------------
!> Pardiso solver information
!------------------------------------------------------------------------------
TYPE :: pardiso_struc
  INTEGER(i4) :: mtype = 11 !< Matrix type (1 -> real/symmetric, 11 -> real/nonsymmetric)
  INTEGER(i4) :: msglvl = 0 !< Controls message output (0 -> no output, 1 -> print some stats)
  INTEGER(i8) :: pt(64) = 0 !< Handle to internal data structure
  INTEGER(i4) :: iparm(64) = 0 !< This array is used to pass various parameters to Intel MKL PARDISO
  INTEGER(i4), POINTER, DIMENSION(:) :: perm => NULL() !< Permutation vector
  TYPE(C_PTR) :: f_factors = C_NULL_PTR !< ILU0 factors pointer
END TYPE pardiso_struc
!------------------------------------------------------------------------------
!> LU solver class
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_lusolver
  LOGICAL :: refactor = .TRUE. !< Refactor solution on next application
  LOGICAL :: update_graph = .TRUE. !< Perform full factorization including symbolic steps
  !> Factor the full local matrix of a distributed matrix (`Mfull`, refreshed with
  !! `update_slice`) instead of its local values. Set by @ref oft_native_solvers::oft_bjprecond
  !! "block-Jacobi" with `boundary_overlap`; distributed matrices are rejected otherwise.
  LOGICAL :: full_local = .FALSE.
  LOGICAL(kind=c_bool) :: iter_refine = .FALSE. !< Perform iterative refinement
  INTEGER(i4) :: nrhs = 1 !< Number of right hand sides
  INTEGER(i4), POINTER, DIMENSION(:) :: ipiv => NULL() !< Row pivot array ("lapack" only)
  REAL(r8), POINTER, DIMENSION(:,:) :: atmp => NULL() !< Local dense matrix ("lapack" only)
  REAL(r8), POINTER, DIMENSION(:,:) :: sec_rhs => NULL() !< Storage for additional RHS data
  CHARACTER(LEN=7) :: package = DEF_LU_PACK !< Factorization package
#if defined( HAVE_SUPERLU ) || defined( HAVE_SUPERLU_DIST ) || defined( HAVE_UMFPACK )
  TYPE(superlu_struc) :: superlu_struct !< SuperLU solver information
#endif
#ifdef HAVE_MKL
  TYPE(pardiso_struc) :: pardiso_struct !< MKL/PARDISO solver information
#endif
#ifdef HAVE_MUMPS
  TYPE(DMUMPS_STRUC) :: mumps_struct !< MUMPS solver information
#endif
CONTAINS
  !> Solve linear system
  PROCEDURE :: apply => lusolver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => lusolver_setup_xml
  !> Update solver with new settings/operators and refactor
  PROCEDURE :: update => lusolver_update
  !> Check thread safety
  PROCEDURE :: check_thread => lusolver_check_thread
  !> Clean-up internal storage
  PROCEDURE :: delete => lusolver_delete
END TYPE oft_lusolver
!------------------------------------------------------------------------------
!> ILU solver class
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_ilusolver
  LOGICAL :: refactor = .TRUE. !< Refactor solution on next application
  LOGICAL :: update_graph = .TRUE. !< Perform full factorization including symbolic steps
  !> Factor the full local matrix of a distributed matrix (`Mfull`, refreshed with
  !! `update_slice`) instead of its local values. Set by @ref oft_native_solvers::oft_bjprecond
  !! "block-Jacobi" with `boundary_overlap`; distributed matrices are rejected otherwise.
  LOGICAL :: full_local = .FALSE.
  INTEGER(i4) :: nrhs = 1 !< Number of right hand sides
  !> Relative small-pivot threshold ("native" only): a pivot smaller in magnitude than
  !! `pivot_tol` times the largest magnitude in its row of the matrix is replaced by
  !! that value, keeping its sign (0 to only catch exact zero pivots)
  REAL(r8) :: pivot_tol = 1.d-12
  !> Diagonal shift ("native" only): factor \f$ A + \alpha \mathrm{diag}(A) \f$ instead of \f$ A \f$
  REAL(r8) :: diag_shift = 0.d0
  INTEGER(i4) :: npivot_replaced = 0 !< Number of pivots replaced in the last factorization ("native" only)
  REAL(r8), POINTER, DIMENSION(:,:) :: sec_rhs => NULL() !< Storage for additional RHS data
  CHARACTER(LEN=7) :: package = DEF_ILU_PACK !< Factorization package
  TYPE(lu_symbolic_key) :: sym_key !< Matrix structure the factor storage was set up for
  TYPE(native_ilu_struc) :: native_struct !< Data for native iLU solver
#if defined( HAVE_SUPERLU )
  TYPE(superlu_struc) :: superlu_struct !< SuperLU solver information
#endif
#ifdef HAVE_MKL
  TYPE(pardiso_struc) :: pardiso_struct !< MKL/PARDISO solver information
#endif
CONTAINS
  !> Solve linear system
  PROCEDURE :: apply => ilusolver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => ilusolver_setup_xml
  !> Update solver with new settings/operators and refactor
  PROCEDURE :: update => ilusolver_update
  !> Check thread safety
  PROCEDURE :: check_thread => ilusolver_check_thread
  !> Clean-up internal storage
  PROCEDURE :: delete => ilusolver_delete
END TYPE oft_ilusolver
INTERFACE
#ifdef HAVE_SUPERLU
!------------------------------------------------------------------------------
!> Interface to dgssv from SuperLU
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dgssv(iopt,n,nnz,nrhs,values,colind, &
    rowptr,b,ldb,f_factors,col_perm,iter_refine,info) BIND(C,NAME="oft_superlu_dgssv_c")
  IMPORT c_bool, c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve (nrhs=1 only)
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to SuperLU internal data storage
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  LOGICAL(c_bool), VALUE, INTENT(in) :: iter_refine !< Perform iterative refinement?
  INTEGER(c_int), INTENT(out) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dgssv
!------------------------------------------------------------------------------
!> Interface to dgsisx (iLUT) from SuperLU
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dgsisx(iopt,n,nnz,nrhs,values,colind, &
    rowptr,b,ldb,f_factors,col_perm,fill_tol,info) BIND(C,NAME="oft_superlu_dgsisx_c")
  IMPORT c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve (nrhs=1 only)
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to SuperLU internal data storage
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  REAL(c_double), VALUE, INTENT(in) :: fill_tol !< Fill tolerance for tresholded iLU
  INTEGER(c_int), INTENT(out) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dgsisx
#endif
#ifdef HAVE_SUPERLU_DIST
!------------------------------------------------------------------------------
!> Interface to dgssv from SuperLU-DIST (local only)
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dist_dgssv(iopt,n,nnz,nrhs,values,colind, &
    rowptr,b,ldb,grid_handle,f_factors,col_perm,iter_refine,info) BIND(C,NAME="oft_superlu_dist_dgssv_c")
  IMPORT c_bool, c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve (nrhs=1 only)
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  TYPE(c_ptr), VALUE, INTENT(in) :: grid_handle !< Pointer to SuperLU-DIST communication object
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to SuperLU-DIST internal data storage
  LOGICAL(c_bool), VALUE, INTENT(in) :: iter_refine !< Perform iterative refinement?
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  INTEGER(c_int), INTENT(inout) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dist_dgssv
!------------------------------------------------------------------------------
!> Interface to dgssv from SuperLU
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dist_slugrid(iopt,slu_comm,nprow,npcol,grid_handle,info) BIND(C,NAME="oft_superlu_dist_slugrid_c")
  IMPORT c_int, c_double, c_ptr
#ifdef OFT_MPI_F08
  IMPORT mpi_comm
#endif
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> Setup grid, 2 -> Clear/Release grid)
#ifdef OFT_MPI_F08
  TYPE(mpi_comm), INTENT(in) :: slu_comm !< MPI communicator for processor grid (must be MPI_COMM_SELF)
#else
  INTEGER(c_int), INTENT(in) :: slu_comm
#endif
  INTEGER(c_int), VALUE, INTENT(in) :: nprow !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: npcol !< Number of non-zeros
  TYPE(c_ptr), INTENT(inout) :: grid_handle !< Pointer to SuperLU-DIST processor grid
  INTEGER(c_int), INTENT(inout) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dist_slugrid
#endif
!------------------------------------------------------------------------------
!> Interface to dgssv from UMFPACK
!------------------------------------------------------------------------------
#ifdef HAVE_UMFPACK
  SUBROUTINE oft_umfpack_dgssv(iopt,n,nnz,nrhs,values,colind, &
    rowptr,b,ldb,f_factors,col_perm,iter_refine,info) BIND(C,NAME="oft_umfpack_dgssv_c")
  IMPORT c_bool, c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !<  Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve (nrhs=1 only)
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to UMFPACK internal data storage
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  LOGICAL(c_bool), VALUE, INTENT(in) :: iter_refine !< Perform iterative refinement?
  INTEGER(c_int), INTENT(inout) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_umfpack_dgssv
#endif
#ifdef HAVE_MKL
!------------------------------------------------------------------------------
!> Interface to MKL iLU0 solver
!------------------------------------------------------------------------------
  SUBROUTINE oft_mkl_ilu(iopt,n,nnz,values,colind, &
    rowptr,b,f_factors,info) BIND(C,NAME="oft_mkl_ilu_c")
  IMPORT c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs iLU0 decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to MKL internal data storage
  INTEGER(c_int), INTENT(inout) :: info !< Pointer to MKL internal data storage
  END SUBROUTINE oft_mkl_ilu
#endif
END INTERFACE
INTERFACE lapack_matinv
  MODULE PROCEDURE lapack_matinv_real
  MODULE PROCEDURE lapack_matinv_complex
  MODULE PROCEDURE lapack_matinv_complex4
END INTERFACE lapack_matinv
INTERFACE lapack_cholesky
  MODULE PROCEDURE lapack_cholesky_real
  MODULE PROCEDURE lapack_cholesky_real4
  ! MODULE PROCEDURE lapack_cholesky_complex
END INTERFACE lapack_cholesky
CONTAINS
!------------------------------------------------------------------------------
!> Convert csr graph to csc storage
!------------------------------------------------------------------------------
SUBROUTINE convert_csr_to_csc(csr_graph,csc_graph,map)
TYPE(oft_graph), INTENT(inout) :: csr_graph !< Original graph in csr format
TYPE(oft_graph), INTENT(inout) :: csc_graph !< Resulting graph in csc format
INTEGER(i4), INTENT(inout) :: map(:) !< Mapping from csr entries to csc entries
INTEGER(i4) :: i,j,k
INTEGER(i4), POINTER :: kr_tmp(:)
DEBUG_STACK_PUSH
!---Convert to CSC format
csc_graph%nr=csr_graph%nr; csc_graph%nc=csr_graph%nc
csc_graph%nrg=-1; csc_graph%ncg=-1
csc_graph%nnz=csr_graph%nnz
!---
ALLOCATE(csc_graph%kr(csc_graph%nr+1))
csc_graph%kr=0
DO i=1,csr_graph%nr
  DO j=csr_graph%kr(i),csr_graph%kr(i+1)-1
    csc_graph%kr(csr_graph%lc(j))=csc_graph%kr(csr_graph%lc(j))+1
  END DO
END DO
!---
csc_graph%kr(csc_graph%nr+1)=csc_graph%nnz+1
DO i=csc_graph%nr,1,-1 ! cumulative point to point count
  csc_graph%kr(i)=csc_graph%kr(i+1)-csc_graph%kr(i)
END DO
IF(csc_graph%kr(1)/=1)CALL oft_abort('Bad element to element count', &
'convert_csr_to_csc',__FILE__)
!---
ALLOCATE(kr_tmp(csc_graph%nr))
ALLOCATE(csc_graph%lc(csc_graph%nnz))
kr_tmp=0
DO i=1,csr_graph%nr
  DO j=csr_graph%kr(i),csr_graph%kr(i+1)-1
    k=csc_graph%kr(csr_graph%lc(j))+kr_tmp(csr_graph%lc(j))
    csc_graph%lc(k)=i
    map(k)=j!csr_map(j)
    kr_tmp(csr_graph%lc(j))=kr_tmp(csr_graph%lc(j))+1
  END DO
END DO
DEALLOCATE(kr_tmp)
DEBUG_STACK_POP
END SUBROUTINE convert_csr_to_csc
!------------------------------------------------------------------------------
!> Solve a linear system using a direct solve
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE lusolver_apply(self,u,g)
CLASS(oft_lusolver), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_vector), INTENT(inout) :: g !< RHS/Residual field
!---
INTEGER(i4) :: mode,nrhs,ldb,ierr,i,j,k,info
INTEGER(i4), POINTER :: csr_map(:),kr_tmp(:)
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: mat_vals
REAL(r8), POINTER, DIMENSION(:) :: csc_vals,vtmp
REAL(r8), POINTER, DIMENSION(:,:) :: vals,b
CLASS(oft_native_matrix), POINTER :: A_native
TYPE(oft_timer) :: mytimer
TYPE(oft_graph) :: csr_graph
TYPE(oft_graph) :: csc_graph
DEBUG_STACK_PUSH
IF(TRIM(self%package)=='pardiso')self%package='mkl'
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting LU solver: ',self%package,self%refactor
IF(.NOT.native_matrix_cast(A_native,self%A))CALL oft_abort('Native matrix required', &
  'lusolver_apply',__FILE__)
!---Check call threading for MUMPS
IF(self%package(1:5)=='mumps')THEN
  IF(omp_get_num_threads()>1)CALL oft_abort('MUMPS used by multiple threads', &
    'lusolver_apply',__FILE__)
END IF
CALL lu_get_values(A_native,self%full_local,self%refactor,mat_vals,'lusolver_apply')
!---Initialize solvers
IF(.NOT.self%initialized)THEN
  !---
  SELECT CASE(TRIM(self%package))
    CASE DEFAULT
      CALL oft_abort('Unknown LU package','lusolver_apply', __FILE__)
    CASE("lapack")
      ! Do nothing
    CASE("super")
#ifdef HAVE_SUPERLU
      ALLOCATE(self%superlu_struct%kr(A_native%nr+1))
      ALLOCATE(self%superlu_struct%lc(A_native%nnz))
      self%superlu_struct%kr=A_native%kr-1
      self%superlu_struct%lc=A_native%lc-1
#else
      CALL oft_abort('OFT not compiled with SUPERLU','lusolver_apply',__FILE__)
#endif
    CASE("superd")
#ifdef HAVE_SUPERLU_DIST
      ! CSC graph is built at the first factorization
      mode=1
      nrhs=1
      CALL MPI_Comm_dup(MPI_COMM_SELF,self%superlu_struct%comm,ierr)
      CALL oft_superlu_dist_slugrid(mode,self%superlu_struct%comm,nrhs,nrhs,self%superlu_struct%grid_handle,ierr)
#else
      CALL oft_abort('OFT not compiled with SUPERLU-DIST','lusolver_apply',__FILE__)
#endif
    CASE("umfpack")
#ifdef HAVE_UMFPACK
      ALLOCATE(self%superlu_struct%kr(A_native%nr+1))
      ALLOCATE(self%superlu_struct%lc(A_native%nnz))
      self%superlu_struct%kr=A_native%kr-1
      self%superlu_struct%lc=A_native%lc-1
#else
      CALL oft_abort('OFT not compiled with UMFPACK','lusolver_apply',__FILE__)
#endif
    CASE("mkl")
#ifdef HAVE_MKL
      self%pardiso_struct%pt=0
      self%pardiso_struct%iparm=0
      call pardisoinit(self%pardiso_struct%pt,self%pardiso_struct%mtype,self%pardiso_struct%iparm)
      ALLOCATE(self%pardiso_struct%perm(A_native%nr))
      self%pardiso_struct%iparm(1)=1
      self%pardiso_struct%iparm(27)=1
#else
      CALL oft_abort('OFT not compiled with MKL-PARDISO','lusolver_apply',__FILE__)
#endif
    CASE("mumps")
#ifdef HAVE_MUMPS
      self%mumps_struct%job = -1
#ifdef OFT_MPI_F08
      self%mumps_struct%comm = MPI_COMM_SELF%MPI_VAL
#else
      self%mumps_struct%comm = MPI_COMM_SELF
#endif
      self%mumps_struct%sym = 0
      self%mumps_struct%par = 1
      CALL dmumps(self%mumps_struct)
      self%mumps_struct%icntl((/5,18/))=0 ! Set centralized assembly
      self%mumps_struct%icntl((/20,21/))=0 ! Set centralized dense RHS
      self%mumps_struct%icntl(3:4)=-1 ! Surpress information messages
      self%mumps_struct%n=A_native%nr
      self%mumps_struct%nz=A_native%nnz
      self%mumps_struct%jcn=>A_native%lc
      self%mumps_struct%A=>mat_vals
      ALLOCATE(self%mumps_struct%irn(A_native%nnz))
      DO i=1,A_native%nr
        self%mumps_struct%irn(A_native%kr(i):A_native%kr(i+1)-1)=i
      END DO
#else
      CALL oft_abort('OFT not compiled with MUMPS','lusolver_apply',__FILE__)
#endif
  END SELECT
  self%initialized=.TRUE.
END IF
!---Call factorization
ALLOCATE(vals(u%n,self%nrhs))
nrhs=1
vtmp=>vals(:,1)
CALL g%get_local(vtmp)
IF(self%nrhs>1)THEN
  nrhs=self%nrhs
  DO i=2,self%nrhs
    vals(:,i) = self%sec_rhs(:,i-1)
  END DO
END IF
ldb=A_native%nr
! NULLIFY(vals)
! CALL g%get_local(vals)
SELECT CASE(TRIM(self%package))
  CASE("super")
#ifdef HAVE_SUPERLU
    IF(self%refactor)THEN
      mode=3
      IF(self%update_graph)mode=1
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
      !$omp critical (superlu_solve)
#endif
      CALL oft_superlu_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
        mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
        self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (SuperLU)','lusolver_apply',__FILE__)
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
      !$omp end critical (superlu_solve)
#endif
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=2
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
    !$omp critical (superlu_solve)
#endif
    CALL oft_superlu_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
      mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
      self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (SuperLU)','lusolver_apply',__FILE__)
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
    !$omp end critical (superlu_solve)
#endif
#else
    CALL oft_abort('OFT not compiled with SuperLU','lusolver_apply',__FILE__)
#endif
CASE("superd")
#ifdef HAVE_SUPERLU_DIST
    IF(self%refactor)THEN
      IF(self%update_graph.OR.(.NOT.ASSOCIATED(self%superlu_struct%lc_orig)))THEN
        !---Build CSC graph and perform full factorization
        IF(ASSOCIATED(self%superlu_struct%lc_orig))THEN
          DEALLOCATE(self%superlu_struct%csc_map,self%superlu_struct%csc_vals)
          DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc,self%superlu_struct%lc_orig)
        END IF
        csr_graph%nr=A_native%nr; csr_graph%nc=A_native%nc; csr_graph%nnz=A_native%nnz
        csr_graph%kr=>A_native%kr; csr_graph%lc=>A_native%lc
        ALLOCATE(self%superlu_struct%csc_map(A_native%nnz))
        ALLOCATE(self%superlu_struct%csc_vals(A_native%nnz))
        CALL convert_csr_to_csc(csr_graph,csc_graph,self%superlu_struct%csc_map)
        csc_graph%kr=csc_graph%kr-1; csc_graph%lc=csc_graph%lc-1
        self%superlu_struct%kr=>csc_graph%kr; self%superlu_struct%lc=>csc_graph%lc
        ALLOCATE(self%superlu_struct%lc_orig(A_native%nnz))
        self%superlu_struct%lc_orig=self%superlu_struct%lc
        mode=1
      ELSE
        !---Restore row indices, which the previous factorization (pdgssvx_ABglobal)
        !---overwrote with their row-permuted values, and refactor with the same pattern
        self%superlu_struct%lc=self%superlu_struct%lc_orig
        mode=3
      END IF
      self%superlu_struct%csc_vals=mat_vals(self%superlu_struct%csc_map)
      !!$omp critical (superlu_solve)
      CALL oft_superlu_dist_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
        self%superlu_struct%csc_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
        self%superlu_struct%grid_handle,self%superlu_struct%f_factors, &
        self%superlu_struct%col_perm,self%iter_refine,ierr)
      !!$omp end critical (superlu_solve)
      IF(ierr/=0)CALL oft_abort('Factorization failed (SuperLU-DIST)','lusolver_apply',__FILE__)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=2
    !!$omp critical (superlu_solve)
    CALL oft_superlu_dist_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
      self%superlu_struct%csc_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
      self%superlu_struct%grid_handle,self%superlu_struct%f_factors, &
      self%superlu_struct%col_perm,self%iter_refine,ierr)
    !!$omp end critical (superlu_solve)
    IF(ierr/=0)CALL oft_abort('Solve failed (SuperLU-DIST)','lusolver_apply',__FILE__)
#else
    CALL oft_abort('OFT not compiled with SuperLU-DIST','lusolver_apply',__FILE__)
#endif
  CASE("umfpack")
#ifdef HAVE_UMFPACK
    IF(self%refactor)THEN
      mode=3
      IF(self%update_graph)mode=1
      !!$omp critical (umfpack_solve)
      CALL oft_umfpack_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
        mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
        self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (UMFPACK)','lusolver_apply',__FILE__)
      !!$omp end critical (umfpack_solve)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=2
    !!$omp critical (umfpack_solve)
    CALL oft_umfpack_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
      mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
      self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (UMFPACK)','lusolver_apply',__FILE__)
    !!$omp end critical (umfpack_solve)
#else
    CALL oft_abort('OFT not compiled with UMFPACK','lusolver_apply',__FILE__)
#endif
  CASE("mkl")
#ifdef HAVE_MKL
    ALLOCATE(b(A_native%nr,self%nrhs))
    b=vals
    IF(self%refactor)THEN
      !---Update matrix
      mode=22
      IF(self%update_graph)mode=12
      CALL pardiso(self%pardiso_struct%pt,1,1,self%pardiso_struct%mtype,mode,A_native%nr, &
        mat_vals,A_native%kr,A_native%lc,self%pardiso_struct%perm, &
        nrhs,self%pardiso_struct%iparm,self%pardiso_struct%msglvl,b,vals,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (MKL-PARDISO)','lusolver_apply',__FILE__)
      self%pardiso_struct%iparm(1)=0
      self%pardiso_struct%iparm(27)=0
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=33
    IF(self%iter_refine)THEN
      self%pardiso_struct%iparm(8)=0
    ELSE
      self%pardiso_struct%iparm(8)=1
    END IF
    CALL pardiso(self%pardiso_struct%pt,1,1,self%pardiso_struct%mtype,mode,A_native%nr, &
      mat_vals,A_native%kr,A_native%lc,self%pardiso_struct%perm, &
      nrhs,self%pardiso_struct%iparm,self%pardiso_struct%msglvl,b,vals,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (MKL-PARDISO)','lusolver_apply',__FILE__)
    DEALLOCATE(b)
#else
    CALL oft_abort('OFT not compiled with MKL-PARDISO','lusolver_apply',__FILE__)
#endif
  CASE("mumps")
#ifdef HAVE_MUMPS
    IF(self%refactor)THEN
      ! (1 for analysis, 2 for factorization, 4 for both)
      self%mumps_struct%job = 2
      IF(self%update_graph)self%mumps_struct%job = 4
      CALL dmumps(self%mumps_struct)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    CALL c_f_pointer(C_LOC(vals),self%mumps_struct%rhs,[nrhs*g%n])
    ! self%mumps_struct%rhs => vals(1:g%n,1:nrhs)
    self%mumps_struct%lrhs = g%n
    self%mumps_struct%nrhs = nrhs
    !---
    self%mumps_struct%job = 3 ! Solve phase
    CALL dmumps(self%mumps_struct)
#else
    CALL oft_abort('OFT not compiled with MUMPS','lusolver_apply',__FILE__)
#endif
  CASE("lapack")
    IF(.NOT.ASSOCIATED(self%atmp))THEN
      ALLOCATE(self%ipiv(A_native%nr))
      ALLOCATE(self%atmp(A_native%nr,A_native%nc))
    END IF
    IF(self%refactor)THEN
      !---Update matrix
      self%atmp=0.d0
      DO i=1,A_native%nr
        DO j=A_native%kr(i),A_native%kr(i+1)-1
          self%atmp(i,A_native%lc(j))=mat_vals(j)
        END DO
      END DO
      CALL dgetrf(A_native%nr,A_native%nr,self%atmp,A_native%nr,self%ipiv,info)
      IF(info/=0)CALL oft_abort('Factorization failed (LAPACK)','lusolver_apply',__FILE__)
      self%refactor=.FALSE.
    END IF
    CALL dgetrs('N',A_native%nr,nrhs,self%atmp,A_native%nr,self%ipiv,vals,ldb,info)
    IF(info/=0)CALL oft_abort('Solve failed (LAPACK)','lusolver_apply',__FILE__)
  CASE DEFAULT
    CALL oft_abort('Unknown factorization package','lusolver_apply',__FILE__)
END SELECT
CALL u%restore_local(vals(:,1))
IF(self%nrhs>1)THEN
  DO i=2,self%nrhs
    self%sec_rhs(:,i-1) = vals(:,i)
  END DO
END IF
DEALLOCATE(vals)
DEBUG_STACK_POP
END SUBROUTINE lusolver_apply
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators
!------------------------------------------------------------------------------
recursive subroutine lusolver_update(self,new_pattern)
class(oft_lusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix pattern (optional)
self%refactor=.TRUE.
IF(PRESENT(new_pattern))self%update_graph=new_pattern
end subroutine lusolver_update
!------------------------------------------------------------------------------
!> Setup solver from XML definition
!------------------------------------------------------------------------------
subroutine lusolver_setup_xml(self,solver_node,level)
CLASS(oft_lusolver), INTENT(inout) :: self
TYPE(xml_node), INTENT(in) :: solver_node !< XML node containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node) :: current_node
!---
CHARACTER(LEN=:), ALLOCATABLE :: factor_package
CHARACTER(LEN=3) :: fac_type
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
!---
CALL xml_get_element(solver_node,"package",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,factor_package,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `package` node","lusolver_setup_xml",__FILE__)
  IF(ALLOCATED(factor_package))THEN
    IF(LEN(factor_package)>7)CALL oft_abort('Factorization package name too long','lusolver_setup_xml',__FILE__)
    self%package=factor_package
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'LU solver setup:'
  WRITE(*,'(2X,2A)')'- Package:  ',self%package
END IF
DEBUG_STACK_POP
end subroutine lusolver_setup_xml
!------------------------------------------------------------------------------
!> Check for thread safety
!------------------------------------------------------------------------------
recursive function lusolver_check_thread(self) result(thread_safe)
class(oft_lusolver), intent(inout) :: self
logical :: thread_safe
thread_safe=.TRUE.
IF(self%package(1:5)=='mumps')thread_safe=.FALSE.
IF(self%package(1:6)=='superd')thread_safe=.FALSE.
end function lusolver_check_thread
!------------------------------------------------------------------------------
!> Destroy direct solver and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine lusolver_delete(self,propogate)
class(oft_lusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
INTEGER(i4) :: mode,nrhs,ierr,ldb,nr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ivals
REAL(r8), ALLOCATABLE, DIMENSION(:) :: rvals
self%refactor=.TRUE.
self%update_graph=.TRUE.
IF(.NOT.self%initialized)THEN
  NULLIFY(self%A)
  RETURN
END IF
nrhs=0
ldb=0
ALLOCATE(ivals(1),rvals(1))
SELECT CASE(TRIM(self%package))
#ifdef HAVE_SUPERLU
  CASE("super")
    mode=4
    CALL oft_superlu_dgssv(mode,nrhs,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
    DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc)
    self%superlu_struct%f_factors=C_NULL_PTR
#endif
#ifdef HAVE_SUPERLU_DIST
  CASE("superd")
    mode=4
    CALL oft_superlu_dist_dgssv(mode,self%A%nr,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%grid_handle,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
    self%superlu_struct%f_factors=C_NULL_PTR
    mode=2
    nrhs=1
    CALL oft_superlu_dist_slugrid(mode,self%superlu_struct%comm,nrhs,nrhs,self%superlu_struct%grid_handle,ierr)
    self%superlu_struct%grid_handle=C_NULL_PTR
    CALL MPI_COMM_FREE(self%superlu_struct%comm, ierr)
    IF(ASSOCIATED(self%superlu_struct%lc_orig))THEN
      DEALLOCATE(self%superlu_struct%csc_vals,self%superlu_struct%csc_map)
      DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc,self%superlu_struct%lc_orig)
    END IF
    self%superlu_struct%comm=MPI_COMM_NULL
#endif
#ifdef HAVE_UMFPACK
  CASE("umfpack")
    mode=4
    CALL oft_umfpack_dgssv(mode,nrhs,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
    DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc)
    self%superlu_struct%f_factors=C_NULL_PTR
#endif
#ifdef HAVE_MUMPS
  CASE("mumps")
    self%mumps_struct%job = -2
    CALL dmumps(self%mumps_struct)
#endif
#ifdef HAVE_MKL
  CASE("mkl")
    mode=-1
    CALL pardiso(self%pardiso_struct%pt,1,1,self%pardiso_struct%mtype,mode,nr, &
      rvals,ivals,ivals,self%pardiso_struct%perm, &
      1,self%pardiso_struct%iparm,self%pardiso_struct%msglvl,rvals,rvals,ierr)
#endif
  CASE("lapack")
    DEALLOCATE(self%ipiv,self%atmp)
END SELECT
DEALLOCATE(ivals,rvals)
IF(ASSOCIATED(self%sec_rhs))DEALLOCATE(self%sec_rhs)
NULLIFY(self%A)
self%initialized=.FALSE.
end subroutine lusolver_delete
!------------------------------------------------------------------------------
!> Get the values to factor for a direct or incomplete LU solver
!!
!! With `full_local` the full local matrix (`Mfull`) is used, refreshed with
!! `update_slice` if out of date. Otherwise the matrix must not be distributed
!! across processes, as its local values do not include contributions from other
!! processes on shared rows.
!------------------------------------------------------------------------------
SUBROUTINE lu_get_values(A_native,full_local,refresh,mat_vals,caller)
CLASS(oft_native_matrix), INTENT(inout) :: A_native !< Matrix to factor
LOGICAL, INTENT(in) :: full_local !< Use the full local matrix?
LOGICAL, INTENT(in) :: refresh !< Refresh the full local matrix if out of date?
REAL(r8), CONTIGUOUS, POINTER, INTENT(out) :: mat_vals(:) !< Values to factor
CHARACTER(LEN=*), INTENT(in) :: caller !< Calling routine (for error messages)
IF(full_local)THEN
  IF(refresh.AND.(.NOT.A_native%full_current))CALL A_native%update_slice
  IF(.NOT.ASSOCIATED(A_native%Mfull))CALL oft_abort('Full local matrix not available', &
    caller,__FILE__)
  mat_vals=>A_native%Mfull
ELSE
  IF((.NOT.A_native%force_local).AND.(INT(A_native%nr,i8)/=A_native%nrg))CALL oft_abort( &
    'Matrix is distributed across processes; use this solver inside a "block_jacobi" '// &
    'or "add_schwarz" preconditioner',caller,__FILE__)
  mat_vals=>A_native%M
END IF
END SUBROUTINE lu_get_values
!------------------------------------------------------------------------------
!> Check if a factorization's symbolic setup matches the structure of a matrix
!!
!! @result Structure matches `key`?
!------------------------------------------------------------------------------
FUNCTION lu_key_matches(key,A_native) RESULT(match)
TYPE(lu_symbolic_key), INTENT(in) :: key !< Structure used for the existing setup
CLASS(oft_native_matrix), INTENT(in) :: A_native !< Matrix to factor
LOGICAL :: match
match=(key%nr==A_native%nr).AND.(key%nnz==A_native%nnz)
IF(.NOT.match)RETURN
IF(ASSOCIATED(A_native%pattern))THEN
  match=(key%pattern_id==A_native%pattern%id)
ELSE
  match=(key%pattern_id==0).AND.ASSOCIATED(key%kr,A_native%kr).AND.ASSOCIATED(key%lc,A_native%lc)
END IF
END FUNCTION lu_key_matches
!------------------------------------------------------------------------------
!> Record the structure of a matrix for a factorization's symbolic setup
!------------------------------------------------------------------------------
SUBROUTINE lu_key_set(key,A_native)
TYPE(lu_symbolic_key), INTENT(inout) :: key !< Structure used for the setup
CLASS(oft_native_matrix), INTENT(in) :: A_native !< Matrix to factor
key%nr=A_native%nr
key%nnz=A_native%nnz
IF(ASSOCIATED(A_native%pattern))THEN
  key%pattern_id=A_native%pattern%id
  NULLIFY(key%kr,key%lc)
ELSE
  key%pattern_id=0
  key%kr=>A_native%kr
  key%lc=>A_native%lc
END IF
END SUBROUTINE lu_key_set
!------------------------------------------------------------------------------
!> Check that a CSR matrix has a diagonal entry in every row and sorted columns
!!
!! Required by the ILU(0) factorizations, which store the factors in modified
!! sparse row format (one slot per non-zero plus one) and eliminate in column order.
!! Aborts naming the first offending row.
!------------------------------------------------------------------------------
SUBROUTINE ilu_check_structure(A_native)
CLASS(oft_native_matrix), INTENT(in) :: A_native !< Matrix to factor
INTEGER(i4) :: i,j,ndiag
INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: dloc
CHARACTER(LEN=OFT_ERROR_SLEN) :: msg
IF(A_native%nr/=A_native%nc)CALL oft_abort('ILU requires a square matrix','ilu_check_structure',__FILE__)
NULLIFY(dloc)
IF(ASSOCIATED(A_native%pattern))dloc=>native_pattern_diag(A_native%pattern)
DO i=1,A_native%nr
  !---Sorted (strictly increasing) columns
  DO j=A_native%kr(i),A_native%kr(i+1)-2
    IF(A_native%lc(j)>=A_native%lc(j+1))THEN
      WRITE(msg,'(A,I0,A)')'ILU requires sorted, unique column indices (row ',i,')'
      CALL oft_abort(msg,'ilu_check_structure',__FILE__)
    END IF
  END DO
  !---Diagonal entry
  IF(ASSOCIATED(dloc))THEN
    ndiag=MIN(dloc(i),1)
  ELSE
    ndiag=COUNT(A_native%lc(A_native%kr(i):A_native%kr(i+1)-1)==i)
  END IF
  IF(ndiag<1)THEN
    WRITE(msg,'(A,I0,A)')'ILU requires a diagonal entry in every row (missing in row ',i,')'
    CALL oft_abort(msg,'ilu_check_structure',__FILE__)
  END IF
END DO
END SUBROUTINE ilu_check_structure
!------------------------------------------------------------------------------
!> Apply ILU(0) preconditioner (approximate solve of the linear system)
!!
!! The factor storage is set up (symbolic setup) on first use, on
!! `update(new_pattern=.TRUE.)` and whenever the non-zero pattern or size of the
!! matrix changes; the factorization is redone after any call to `update`.
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE ilusolver_apply(self,u,g)
CLASS(oft_ilusolver), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_vector), INTENT(inout) :: g !< RHS/Residual field
!---
INTEGER(i4) :: i,ierr,nrhs
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: mat_vals,uvals,gvals
REAL(r8), POINTER, DIMENSION(:) :: vtmp
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: vals
CLASS(oft_native_matrix), POINTER :: A_native
CHARACTER(LEN=OFT_ERROR_SLEN) :: msg
DEBUG_STACK_PUSH
IF(TRIM(self%package)=='pardiso')self%package='mkl'
IF(oft_debug_print(3))WRITE(*,'(2X,3A,L1)')'Starting ILU solver: ',TRIM(self%package),' ',self%refactor
IF(.NOT.native_matrix_cast(A_native,self%A))CALL oft_abort('Native matrix required', &
  'ilusolver_apply',__FILE__)
!---Check package
IF(.NOT.self%initialized)THEN
  SELECT CASE(TRIM(self%package))
    CASE("native")
      ! Nothing to do
    CASE("mkl")
#if !defined(HAVE_MKL)
      CALL oft_abort('OFT not compiled with MKL-PARDISO','ilusolver_apply',__FILE__)
#endif
    CASE DEFAULT
      CALL oft_abort('Unknown ILU package','ilusolver_apply', __FILE__)
  END SELECT
  self%initialized=.TRUE.
END IF
!---Symbolic setup for a new matrix structure
IF(self%update_graph.OR.(.NOT.lu_key_matches(self%sym_key,A_native)))THEN
  CALL ilu_check_structure(A_native)
  SELECT CASE(TRIM(self%package))
    CASE("native")
      IF(ALLOCATED(self%native_struct%alu))DEALLOCATE(self%native_struct%alu, &
        self%native_struct%jlu,self%native_struct%ju)
      ALLOCATE(self%native_struct%ju(A_native%nr))
      ALLOCATE(self%native_struct%jlu(A_native%nnz+1))
      ALLOCATE(self%native_struct%alu(A_native%nnz+1))
    CASE("mkl")
#ifdef HAVE_MKL
      !---Factor storage is sized for the previous matrix
      CALL ilusolver_mkl_free(self)
#endif
  END SELECT
  CALL lu_key_set(self%sym_key,A_native)
  self%update_graph=.FALSE.
  self%refactor=.TRUE.
END IF
!---Numerical factorization
CALL lu_get_values(A_native,self%full_local,self%refactor,mat_vals,'ilusolver_apply')
IF(self%refactor)THEN
  SELECT CASE(TRIM(self%package))
    CASE("native")
      CALL ilu0(A_native%nr,mat_vals,A_native%lc,A_native%kr,self%native_struct%alu, &
        self%native_struct%jlu,self%native_struct%ju,self%pivot_tol,self%diag_shift, &
        self%npivot_replaced,ierr)
      IF(ierr/=0)THEN
        IF(ierr>0)THEN
          WRITE(msg,'(A,I0,A)')'Zero pivot in row ',ierr,' (try increasing "pivot_tol" or "diag_shift")'
        ELSE
          WRITE(msg,'(A,I0)')'Non-finite pivot in row ',-ierr
        END IF
        CALL oft_abort(msg,'ilusolver_apply',__FILE__)
      END IF
      IF(self%npivot_replaced>0)CALL ilu_pivot_warning(self%npivot_replaced,A_native%nr)
    CASE("mkl")
#ifdef HAVE_MKL
      ALLOCATE(vtmp(A_native%nr)) ! RHS argument (not used by the factorization)
      CALL oft_mkl_ilu(1,A_native%nr,A_native%nnz,mat_vals,A_native%lc, &
        A_native%kr,vtmp,self%pardiso_struct%f_factors,ierr)
      DEALLOCATE(vtmp)
      IF(ierr/=0)THEN
        WRITE(msg,'(A,I0,A)')'Factorization failed (MKL error ',ierr,')'
        CALL oft_abort(msg,'ilusolver_apply',__FILE__)
      END IF
#endif
  END SELECT
  self%refactor=.FALSE.
END IF
!---Approximate solve
IF(self%nrhs>1)THEN
  IF(.NOT.ASSOCIATED(self%sec_rhs))CALL oft_abort('Additional RHS storage ("sec_rhs") not set', &
    'ilusolver_apply',__FILE__)
  nrhs=self%nrhs
  ALLOCATE(vals(A_native%nr,nrhs))
  vtmp=>vals(:,1)
  CALL g%get_local(vtmp)
  DO i=2,nrhs
    vals(:,i) = self%sec_rhs(:,i-1)
  END DO
  CALL ilusolver_solve(self,A_native%nr,nrhs,vals)
  CALL u%restore_local(vals(:,1))
  DO i=2,nrhs
    self%sec_rhs(:,i-1) = vals(:,i)
  END DO
  DEALLOCATE(vals)
ELSE
  !---Solve in place in the local values of u (no copies for native vectors)
  NULLIFY(uvals,gvals)
  CALL g%borrow_local(gvals,OFT_ACCESS_READ)
  CALL u%borrow_local(uvals,OFT_ACCESS_WRITE)
  uvals=gvals
  CALL g%release_local(gvals,OFT_ACCESS_READ,sync=.FALSE.)
  vals(1:A_native%nr,1:1)=>uvals
  CALL ilusolver_solve(self,A_native%nr,1,vals)
  CALL u%release_local(uvals,OFT_ACCESS_WRITE,sync=.TRUE.)
END IF
DEBUG_STACK_POP
END SUBROUTINE ilusolver_apply
!------------------------------------------------------------------------------
!> Approximate solve with existing ILU(0) factors (in place)
!------------------------------------------------------------------------------
SUBROUTINE ilusolver_solve(self,n,nrhs,x)
CLASS(oft_ilusolver), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: n !< Number of rows
INTEGER(i4), INTENT(in) :: nrhs !< Number of right hand sides
REAL(r8), INTENT(inout) :: x(n,nrhs) !< Input: RHS, Output: Solution
INTEGER(i4) :: i,ierr
CLASS(oft_native_matrix), POINTER :: A_native
SELECT CASE(TRIM(self%package))
  CASE("native")
    CALL lusol(n,nrhs,x,self%native_struct%alu,self%native_struct%jlu,self%native_struct%ju)
  CASE("mkl")
#ifdef HAVE_MKL
    IF(.NOT.native_matrix_cast(A_native,self%A))CALL oft_abort('Native matrix required', &
      'ilusolver_solve',__FILE__)
    DO i=1,nrhs
      CALL oft_mkl_ilu(2,A_native%nr,A_native%nnz,A_native%M,A_native%lc, &
        A_native%kr,x(:,i),self%pardiso_struct%f_factors,ierr)
      IF(ierr/=0)CALL oft_abort('Solve failed (MKL)','ilusolver_solve',__FILE__)
    END DO
#endif
  CASE DEFAULT
    CALL oft_abort('Unknown factorization package','ilusolver_solve',__FILE__)
END SELECT
END SUBROUTINE ilusolver_solve
#ifdef HAVE_MKL
!------------------------------------------------------------------------------
!> Free MKL ILU(0) factor storage
!------------------------------------------------------------------------------
SUBROUTINE ilusolver_mkl_free(self)
CLASS(oft_ilusolver), INTENT(inout) :: self
INTEGER(i4) :: ierr,ivals(1)
REAL(r8) :: rvals(1)
IF(.NOT.c_associated(self%pardiso_struct%f_factors))RETURN
ivals=0; rvals=0.d0
CALL oft_mkl_ilu(3,1,1,rvals,ivals,ivals,rvals,self%pardiso_struct%f_factors,ierr)
self%pardiso_struct%f_factors=C_NULL_PTR
END SUBROUTINE ilusolver_mkl_free
#endif
!------------------------------------------------------------------------------
!> Report replaced ILU(0) pivots (once per process)
!------------------------------------------------------------------------------
SUBROUTINE ilu_pivot_warning(nreplaced,nr)
INTEGER(i4), INTENT(in) :: nreplaced !< Number of replaced pivots
INTEGER(i4), INTENT(in) :: nr !< Number of rows
LOGICAL, SAVE :: warned = .FALSE.
LOGICAL :: do_warn
CHARACTER(LEN=OFT_ERROR_SLEN) :: msg
do_warn=.FALSE.
!$omp critical (oft_ilu_pivot_warning)
IF(.NOT.warned)THEN
  warned=.TRUE.
  do_warn=.TRUE.
END IF
!$omp end critical (oft_ilu_pivot_warning)
IF(do_warn)THEN
  WRITE(msg,'(A,I0,A,I0,A)')'ILU(0): replaced ',nreplaced,' small pivot(s) of ',nr, &
    ' rows (further replacements are not reported)'
  CALL oft_warn(msg)
END IF
END SUBROUTINE ilu_pivot_warning
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators
!!
!! The factorization is redone on the next application, including the symbolic
!! setup if `new_pattern` is `.TRUE.`. Changes to the matrix structure are also
!! detected automatically.
!------------------------------------------------------------------------------
recursive subroutine ilusolver_update(self,new_pattern)
class(oft_ilusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix pattern (optional)
self%refactor=.TRUE.
IF(PRESENT(new_pattern))THEN
  IF(new_pattern)self%update_graph=.TRUE.
END IF
end subroutine ilusolver_update
!------------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! Options: `package` (factorization package), `pivot_tol` (relative small-pivot
!! threshold) and `diag_shift` (relative diagonal shift); the last two are used by
!! the "native" package only.
!------------------------------------------------------------------------------
subroutine ilusolver_setup_xml(self,solver_node,level)
CLASS(oft_ilusolver), INTENT(inout) :: self
TYPE(xml_node), INTENT(in) :: solver_node !< XML node containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
!---
TYPE(xml_node) :: current_node
CHARACTER(LEN=:), ALLOCATABLE :: factor_package
INTEGER(i4) :: ierr
REAL(r8) :: rtmp
DEBUG_STACK_PUSH
!---
CALL xml_get_element(solver_node,"package",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,factor_package,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `package` node","ilusolver_setup_xml",__FILE__)
  IF(ALLOCATED(factor_package))THEN
    IF(LEN(factor_package)>7)CALL oft_abort('Factorization package name too long','ilusolver_setup_xml',__FILE__)
    self%package=factor_package
  END IF
END IF
CALL xml_get_element(solver_node,"pivot_tol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,rtmp,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `pivot_tol` node","ilusolver_setup_xml",__FILE__)
  IF(rtmp<0.d0)CALL oft_abort('"pivot_tol" must be non-negative','ilusolver_setup_xml',__FILE__)
  self%pivot_tol=rtmp
END IF
CALL xml_get_element(solver_node,"diag_shift",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,rtmp,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `diag_shift` node","ilusolver_setup_xml",__FILE__)
  self%diag_shift=rtmp
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'ILU solver setup:'
  WRITE(*,'(2X,2A)')'- Package:    ',self%package
  WRITE(*,'(2X,A,ES11.3)')'- Pivot tol:  ',self%pivot_tol
  WRITE(*,'(2X,A,ES11.3)')'- Diag shift: ',self%diag_shift
END IF
DEBUG_STACK_POP
end subroutine ilusolver_setup_xml
!------------------------------------------------------------------------------
!> Check for thread safety
!------------------------------------------------------------------------------
recursive function ilusolver_check_thread(self) result(thread_safe)
class(oft_ilusolver), intent(inout) :: self
logical :: thread_safe
thread_safe=.TRUE.
end function ilusolver_check_thread
!------------------------------------------------------------------------------
!> Destroy ILU solver and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine ilusolver_delete(self,propogate)
class(oft_ilusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
TYPE(lu_symbolic_key) :: key_reset
self%refactor=.TRUE.
self%update_graph=.TRUE.
self%sym_key=key_reset
IF(.NOT.self%initialized)THEN
  NULLIFY(self%A)
  RETURN
END IF
SELECT CASE(TRIM(self%package))
#ifdef HAVE_MKL
  CASE("mkl")
    CALL ilusolver_mkl_free(self)
#endif
  CASE("native")
    IF(ALLOCATED(self%native_struct%alu))THEN
      DEALLOCATE(self%native_struct%alu,self%native_struct%jlu,self%native_struct%ju)
    END IF
END SELECT
NULLIFY(self%A)
self%initialized=.FALSE.
end subroutine ilusolver_delete
!---------------------------------------------------------------------------------
!> Invert a general real matrix in place using LAPACK (double precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_matinv_real(nrows,Amat,error)
INTEGER(4), INTENT(in) :: nrows !< Number of rows/columns
REAL(8), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(4) :: N,LWORK,info
REAL(8) :: elapsed_time
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: ipiv
REAL(8), ALLOCATABLE, DIMENSION(:) :: rwork
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting real matrix'
  CALL stimer%tick
END IF
N = nrows
ALLOCATE(ipiv(N),rwork(1))
CALL dgetrf(N,N,Amat,N,ipiv,info)
IF(info/=0)THEN
  WRITE(*,*)'DGETRF',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=-1
CALL dgetri(N,Amat,N,ipiv,rwork,lwork,info)
IF(info/=0)THEN
  WRITE(*,*)'DGETRI workspace query',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=MAX(INT(rwork(1),4),N)
IF(oft_debug_print(1).AND.oft_env%pm)WRITE(*,*)'  Block size = ',lwork/N
DEALLOCATE(rwork)
ALLOCATE(rwork(lwork))
CALL dgetri(N,Amat,N,ipiv,rwork,lwork,info)
error=info
IF(info/=0)WRITE(*,*)'DGETRI',info
DEALLOCATE(ipiv,rwork)
IF(oft_env%pm)THEN
  elapsed_time=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_time
END IF
END SUBROUTINE lapack_matinv_real
!---------------------------------------------------------------------------------
!> Invert a general complex matrix in place using LAPACK (double precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_matinv_complex(nrows,Amat,error)
INTEGER(i4), INTENT(in) :: nrows !< Number of rows/columns
DOUBLE COMPLEX, INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(i4) :: N,LWORK,info
REAL(r8) :: elapsed_timer
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ipiv
DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: rwork
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting complex matrix'
  CALL stimer%tick
END IF
N = nrows
ALLOCATE(ipiv(N),rwork(1))
CALL zgetrf(N,N,Amat,N,ipiv,info)
IF(info/=0)THEN
  WRITE(*,*)'ZGETRF',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=-1
CALL zgetri(N,Amat,N,ipiv,rwork,lwork,info)
IF(info/=0)THEN
  WRITE(*,*)'ZGETRI workspace query',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=MAX(INT(REAL(rwork(1),r8),4),N)
IF(oft_debug_print(1).AND.oft_env%pm)WRITE(*,*)'  Block size = ',lwork/N
DEALLOCATE(rwork)
ALLOCATE(rwork(lwork))
CALL zgetri(N,Amat,N,ipiv,rwork,lwork,info)
error=info
IF(info/=0)WRITE(*,*)'ZGETRI',info
DEALLOCATE(ipiv,rwork)
IF(oft_env%pm)THEN
  elapsed_timer=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_timer
END IF
END SUBROUTINE lapack_matinv_complex
!---------------------------------------------------------------------------------
!> Invert a general complex matrix in place using LAPACK (single precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_matinv_complex4(nrows,Amat,error)
INTEGER(i4), INTENT(in) :: nrows !< Number of rows/columns
COMPLEX(c4), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(i4) :: N,LWORK,info
REAL(r8) :: elapsed_timer
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ipiv
COMPLEX(c4), ALLOCATABLE, DIMENSION(:) :: rwork
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting complex matrix (single precision)'
  CALL stimer%tick
END IF
N = nrows
ALLOCATE(ipiv(N),rwork(1))
CALL cgetrf(N,N,Amat,N,ipiv,info)
IF(info/=0)THEN
  WRITE(*,*)'CGETRF',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=-1
CALL cgetri(N,Amat,N,ipiv,rwork,lwork,info)
IF(info/=0)THEN
  WRITE(*,*)'CGETRI workspace query',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=MAX(INT(REAL(rwork(1),r4),4),N)
IF(oft_debug_print(1).AND.oft_env%pm)WRITE(*,*)'  Block size = ',lwork/N
DEALLOCATE(rwork)
ALLOCATE(rwork(lwork))
CALL cgetri(N,Amat,N,ipiv,rwork,lwork,info)
error=info
IF(info/=0)WRITE(*,*)'CGETRI',info
DEALLOCATE(ipiv,rwork)
IF(oft_env%pm)THEN
  elapsed_timer=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_timer
END IF
END SUBROUTINE lapack_matinv_complex4
!---------------------------------------------------------------------------------
!> Invert a symmetric positive-definite real matrix in place using a Cholesky
!! factorization (double precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_cholesky_real(nrows,Amat,error)
INTEGER(4), INTENT(in) :: nrows !< Number of rows/columns
REAL(8), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(4) :: N,info,i,j
REAL(8) :: elapsed_time
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting real matrix'
  CALL stimer%tick
END IF
N = nrows
CALL dpotrf('U',N,Amat,N,info)
IF(info/=0)THEN
  error=info
  RETURN
END IF
CALL dpotri('U',N,Amat,N,info)
DO i=1,N
  DO j=1,i-1
    Amat(i,j)=Amat(j,i)
  END DO
END DO
error=info
IF(oft_env%pm)THEN
  elapsed_time=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_time
END IF
END SUBROUTINE lapack_cholesky_real
!---------------------------------------------------------------------------------
!> Invert a symmetric positive-definite real matrix in place using a Cholesky
!! factorization (single precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_cholesky_real4(nrows,Amat,error)
INTEGER(4), INTENT(in) :: nrows !< Number of rows/columns
REAL(r4), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(4) :: N,info,i,j
REAL(8) :: elapsed_time
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting real matrix (single precision)'
  CALL stimer%tick
END IF
N = nrows
CALL spotrf('U',N,Amat,N,info)
IF(info/=0)THEN
  error=info
  RETURN
END IF
CALL spotri('U',N,Amat,N,info)
DO i=1,N
  DO j=1,i-1
    Amat(i,j)=Amat(j,i)
  END DO
END DO
error=info
IF(oft_env%pm)THEN
  elapsed_time=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_time
END IF
END SUBROUTINE lapack_cholesky_real4
!---------------------------------------------------------------------------------
!> Native incomplete LU factorization with no fill (ILU(0))
!!
!! Derived from the `ilu0` routine of SPARSKIT (Y. Saad). The factors are stored in
!! modified sparse row format: the inverted diagonal of U in `alu(1:n)`, the row
!! pointer in `jlu(1:n+1)` and the off-diagonal entries of L (unit diagonal) and U,
!! with their columns, from `alu(n+2)`/`jlu(n+2)`; `ju(i)` is the start of the U part
!! of row `i`.
!!
!! Pivot protection: a pivot whose magnitude is below `pivot_tol` times the largest
!! magnitude in its row of the matrix is replaced by that value, keeping its sign
!! (positive for an exact zero). With `diag_shift` the factorization is of
!! \f$ A + \alpha \mathrm{diag}(A) \f$.
!!
!! @note Requires a CSR graph with sorted column indices and a diagonal entry in every
!! row (see @ref ilu_check_structure); `alu` and `jlu` must have `nnz+1` entries.
!---------------------------------------------------------------------------------
subroutine ilu0(n,a,lc,kr,alu,jlu,ju,pivot_tol,diag_shift,nreplaced,ierr)
integer(i4), intent(in) :: n !< Size of matrix
real(r8), contiguous, dimension(:), intent(in) :: a !< Matrix entries
integer(i4), contiguous, dimension(:), intent(in) :: lc !< Column indices
integer(i4), contiguous, dimension(:), intent(in) :: kr !< Row pointer
real(r8), contiguous, dimension(:), intent(inout) :: alu !< Entries for L/U factors
integer(i4), contiguous, dimension(:), intent(inout) :: jlu !< Column indices for alu
integer(i4), contiguous, dimension(:), intent(inout) :: ju !< Diagonal indices for alu
real(r8), intent(in) :: pivot_tol !< Relative small-pivot threshold (0 to only catch exact zeros)
real(r8), intent(in) :: diag_shift !< Relative diagonal shift
integer(i4), intent(out) :: nreplaced !< Number of pivots replaced
integer(i4), intent(out) :: ierr !< Error flag (0: success, i>0: zero pivot in row i, i<0: non-finite pivot in row -i)
integer(i4), allocatable, dimension(:) :: iw
integer(i4) :: i,ii,j,jj,jcol,jf,jm,jrow,js,ju0,jw
real(r8) :: tl,rowmax,pthresh
!
ju0 = n+2
jlu(1) = ju0
! initialize work vector to zero's
ALLOCATE(iw(n))
do i=1, n
    iw(i) = 0
end do
!---
ierr = 0
nreplaced = 0
do ii=1,n
  js=ju0
  rowmax=0.d0
  !---Generating row ii of L and U
  do j=kr(ii),kr(ii+1)-1
    ! copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
    jcol=lc(j)
    rowmax=MAX(rowmax,ABS(a(j)))
    if(jcol==ii)then
      alu(ii)=a(j)*(1.d0+diag_shift)
      iw(jcol)=ii
      ju(ii)=ju0
    else
      alu(ju0)=a(j)
      jlu(ju0)=lc(j)
      iw(jcol)=ju0
      ju0=ju0+1
    endif
  end do
  jlu(ii+1)=ju0
  jf=ju0-1
  jm=ju(ii)-1
  !---Move up to diagonal
  do j=js,jm
    jrow=jlu(j)
    tl=alu(j)*alu(jrow)
    alu(j)=tl
    do jj=ju(jrow),jlu(jrow+1)-1
      jw=iw(jlu(jj))
      if(jw/=0)alu(jw)=alu(jw)-tl*alu(jj)
    end do
  end do
  !---Check pivot
  if(.NOT.(ABS(alu(ii))<=HUGE(alu(ii))))then
    ierr = -ii
    exit
  end if
  pthresh=pivot_tol*rowmax
  if(ABS(alu(ii))<pthresh)then
    alu(ii)=SIGN(pthresh,alu(ii))
    nreplaced=nreplaced+1
  else if(alu(ii)==0.d0)then
    ierr = ii
    exit
  endif
  alu(ii)=1.d0/alu(ii)
  !---Reset iw pointer
  iw(ii)=0
  do i=js,jf
    iw(jlu(i))=0
  end do
end do
DEALLOCATE(iw)
END SUBROUTINE ilu0
!---------------------------------------------------------------------------------
!> Approximate solve using ILU(0) factorization produced by ilu0
!---------------------------------------------------------------------------------
subroutine lusol(n, nrhs, x, alu, jlu, ju)
integer(4), intent(in) :: n !< Size of matrix
integer(4), intent(in) :: nrhs !< Number of RHS to solve
real(8), intent(inout) :: x(n,nrhs) !< Input: RHS, Output: Solution
real(8), contiguous, dimension(:), intent(in) :: alu !< Entries for L/U factors
integer(4), contiguous, dimension(:), intent(in) :: jlu !< Column indices for alu
integer(4), contiguous, dimension(:), intent(in) :: ju !< Diagonal indices for alu
integer(4) :: i,k
! forward solve
do i = 1, n
  do k=jlu(i),ju(i)-1
    x(i,:)=x(i,:)-alu(k)*x(jlu(k),:)
  end do
end do
! backward solve
do i=n,1,-1
  do k=ju(i),jlu(i+1)-1
      x(i,:)=x(i,:)-alu(k)*x(jlu(k),:)
  end do
  x(i,:) = alu(i)*x(i,:)
end do
end subroutine
END MODULE oft_lu
