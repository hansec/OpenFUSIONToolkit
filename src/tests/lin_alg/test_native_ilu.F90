!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_native_ilu.F90
!
!> Regression tests for the ILU solver (@ref oft_lu::oft_ilusolver): native ILU(0) and
!! SuperLU ILUTP.
!!
!! The current test cases are (selected by `case_id` in `test_ilu_options`):
!! - 1: Local matrices built directly: ILU(0) of banded matrices with a full band
!!   (no fill outside the pattern) is an exact LU factorization, refactoring after
!!   value changes and after switching to matrices with a different pattern or
!!   size (with and without `update(new_pattern=.TRUE.)`), matrices without a
!!   pattern, small-pivot replacement (sign preserving), diagonal shift, deep copies
!!   of a factored solver, and the real and complex ILU(0) kernels called directly.
!! - 2: Poisson equation on a unit cube with GMRES: ILU (package `ilu_pkg`) applied
!!   directly to the operator (single process only) and inside block-Jacobi with
!!   `boundary_overlap` (full local matrix of the distributed operator).
!! - 3: ILU (package `ilu_pkg`) applied directly to a distributed operator (aborts on
!!   more than one process, solves on one).
!! - 4: Timing of the ILU(0) factorization and solve for the Lagrange Laplacian of
!!   order `order` on a single process (written to `ilu_timing.results`), and of the
!!   factorization and the serial and level-scheduled (threaded) solves for 1, 2, 4,
!!   ... threads.
!! - 5: Level-scheduled (threaded) ILU(0) solve and solver on a non-symmetric variant
!!   of the Lagrange Laplacian: schedule structure, bitwise identity with the serial
!!   solve (real and complex, with pivot replacement and diagonal shift, for several
!!   thread counts), the serial path inside parallel regions, and the threaded solver
!!   against the serial one.
!! - 6: GMRES on the Laplacian of order `order` preconditioned by Schwarz (one
!!   subdomain per process, overlap 1) with ILU(0), with and without threading:
!!   iteration counts and solutions must be identical; wall times of the first
!!   solve (including setup and factorization) and of a repeated solve are written
!!   to `ilu_timing.results`.
!! - 7: SuperLU ILUTP ("super") on local matrices built directly: complete LU with no
!!   dropping (also for several right hand sides and matrices without a pattern),
!!   refactoring after value and option changes (bitwise equal to a new factorization),
!!   after pattern and size changes and after `update(new_pattern=.TRUE.)`, column
!!   permutation changes, zero pivot replacement, copies and repeated deletion.
!! - 8: Comparison of native ILU(0) and SuperLU ILUTP for several drop tolerances on the
!!   Lagrange Laplacian of order `order` on a single process: factorization and solve
!!   times, non-zeros of the factors and GMRES iterations (written to `ilu_compare.results`).
!!
!! Cases 7 and 8, and case 2 or 3 with `ilu_pkg='super'`, print `SKIP TEST` without SuperLU.
!!
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_native_ilu
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_stitching, ONLY: oft_seam
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_native_la, ONLY: oft_native_matrix, oft_native_vector, oft_native_pattern, &
  native_pattern_create, native_matrix_set_pattern, native_matrix_cast, native_pattern_diag, &
  native_pattern_ilu_sched
USE oft_native_solvers, ONLY: oft_bjprecond, oft_schwarz_precond, oft_native_gmres_solver
USE oft_solver_base, ONLY: oft_solver
USE oft_lu, ONLY: oft_ilusolver
USE oft_native_kernels, ONLY: nk_ilu0_factor, nk_ilu0_solve, nk_ilu0_sched_values, &
  native_ilu_sched, native_level_sched, native_ilu_sched_build, native_ilu_sched_free, &
  native_ilu_sched_stats, native_ilu_sched_active, native_ilu_level_rows
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_zerob, oft_lag_getlop, oft_lag_getmop
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr,nfail
INTEGER(i4) :: case_id = 1
INTEGER(i4) :: order = 2
INTEGER(i4) :: its(2) = -1
CHARACTER(LEN=7) :: ilu_pkg = 'native'
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
TYPE(oft_seam), TARGET :: local_seam
NAMELIST/test_ilu_options/case_id,order,ilu_pkg
IF(use_petsc)THEN
  WRITE(*,*)'SKIP TEST'
  STOP
END IF
!---Initialize enviroment
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_ilu_options,IOSTAT=ierr)
CLOSE(io_unit)
local_seam%skip=.TRUE.
nfail=0
#if !defined( HAVE_SUPERLU )
IF((case_id==7).OR.(case_id==8).OR.(TRIM(ilu_pkg)=='super'))THEN
  WRITE(*,*)'SKIP TEST'
  CALL oft_finalize
END IF
#endif
SELECT CASE(case_id)
  CASE(1)
    CALL test_exact()
    CALL test_refactor()
    CALL test_pivots()
    CALL test_shift()
    CALL test_copy()
    CALL test_no_pattern()
    CALL test_kernels()
  CASE(2,3)
    CALL multigrid_construct(mg_mesh)
    IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
    CALL oft_lag_setup(mg_mesh,2,ML_oft_lagrange,minlev=-1)
    lag_zerob%ML_lag_rep=>ML_oft_lagrange
    CALL test_lap(case_id)
  CASE(4)
    CALL multigrid_construct(mg_mesh)
    CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
    lag_zerob%ML_lag_rep=>ML_oft_lagrange
    CALL test_timing()
  CASE(5)
    CALL multigrid_construct(mg_mesh)
    CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
    lag_zerob%ML_lag_rep=>ML_oft_lagrange
    CALL test_threaded()
  CASE(6)
    CALL multigrid_construct(mg_mesh)
    CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
    lag_zerob%ML_lag_rep=>ML_oft_lagrange
    CALL test_schwarz_timing()
  CASE(7)
    CALL test_super_exact()
    CALL test_super_refactor()
    CALL test_super_options()
    CALL test_super_pivots()
    CALL test_super_copy()
  CASE(8)
    CALL multigrid_construct(mg_mesh)
    CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
    lag_zerob%ML_lag_rep=>ML_oft_lagrange
    CALL test_compare()
  CASE DEFAULT
    CALL oft_abort('Unknown test case','main',__FILE__)
END SELECT
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='ilu.results')
  WRITE(io_unit,*)nfail
  WRITE(io_unit,*)its
  CLOSE(io_unit)
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
! Record a failed check
SUBROUTINE check(name,pass)
CHARACTER(LEN=*), INTENT(in) :: name
LOGICAL, INTENT(in) :: pass
IF(.NOT.pass)THEN
  nfail=nfail+1
  WRITE(*,'(2A)')'FAILED: ',name
END IF
END SUBROUTINE check
! Local (single process) square CSR matrix from its structure and values (optionally
! without a pattern object, as built by code that sets kr/lc directly)
FUNCTION make_matrix(n,kr,lc,vals,no_pattern) RESULT(A)
INTEGER(i4), INTENT(in) :: n,kr(:),lc(:)
REAL(r8), INTENT(in) :: vals(:)
LOGICAL, OPTIONAL, INTENT(in) :: no_pattern
TYPE(oft_native_matrix), POINTER :: A
INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: krp,lcp
TYPE(oft_native_pattern), POINTER :: pattern
ALLOCATE(A)
A%force_local=.TRUE.
A%nr=n; A%nc=n
A%nrg=n; A%ncg=n
A%nrslice=n; A%ncslice=n
A%ni=1; A%nj=1
ALLOCATE(A%i_map(1),A%j_map(1))
A%i_map(1)%n=n
A%j_map(1)%n=n
ALLOCATE(krp(n+1),lcp(SIZE(lc)))
krp=kr; lcp=lc
IF(PRESENT(no_pattern))THEN
  IF(no_pattern)THEN
    A%kr=>krp; A%lc=>lcp
    A%nnz=SIZE(lc)
    ALLOCATE(A%M(A%nnz))
    A%M=vals
    RETURN
  END IF
END IF
CALL native_pattern_create(pattern,n,n,krp,lcp)
CALL native_matrix_set_pattern(A,pattern)
ALLOCATE(A%M(A%nnz))
A%M=vals
END FUNCTION make_matrix
! Non-symmetric, diagonally dominant banded matrix with all entries of the band
! (half-bandwidth bw) present, so ILU(0) has no dropped fill
FUNCTION make_banded(n,bw,scale,no_pattern) RESULT(A)
INTEGER(i4), INTENT(in) :: n,bw
REAL(r8), INTENT(in) :: scale
LOGICAL, OPTIONAL, INTENT(in) :: no_pattern
TYPE(oft_native_matrix), POINTER :: A
INTEGER(i4) :: i,j,k
INTEGER(i4), ALLOCATABLE :: kr(:),lc(:)
REAL(r8), ALLOCATABLE :: vals(:)
ALLOCATE(kr(n+1),lc(n*(2*bw+1)),vals(n*(2*bw+1)))
k=0
kr(1)=1
DO i=1,n
  DO j=MAX(1,i-bw),MIN(n,i+bw)
    k=k+1
    lc(k)=j
    IF(i==j)THEN
      vals(k)=scale*(2.d0*bw+1.5d0+0.1d0*SIN(REAL(i,r8)))
    ELSE
      vals(k)=scale*(-1.d0+0.3d0*COS(REAL(i+2*j,r8)))
    END IF
  END DO
  kr(i+1)=k+1
END DO
A=>make_matrix(n,kr,lc(1:k),vals(1:k),no_pattern)
DEALLOCATE(kr,lc,vals)
END FUNCTION make_banded
! Local vector of size n that skips all stitching
FUNCTION make_vec(n) RESULT(v)
INTEGER(i4), INTENT(in) :: n
TYPE(oft_native_vector), POINTER :: v
ALLOCATE(v)
ALLOCATE(v%v(n))
v%n=n
v%ng=INT(n,i8)
v%nslice=n
v%stitch_info=>local_seam
v%v=0.d0
END FUNCTION make_vec
! Free a local matrix or vector
SUBROUTINE free_mat(A)
TYPE(oft_native_matrix), POINTER, INTENT(inout) :: A
CALL A%delete
IF(ASSOCIATED(A%i_map))DEALLOCATE(A%i_map,A%j_map)
DEALLOCATE(A)
END SUBROUTINE free_mat
SUBROUTINE free_vec(v)
TYPE(oft_native_vector), POINTER, INTENT(inout) :: v
DEALLOCATE(v%v)
DEALLOCATE(v)
END SUBROUTINE free_vec
! Relative residual ||A x - b||_inf / ||b||_inf, computed directly from the CSR data
FUNCTION rel_resid(A,x,b) RESULT(res)
TYPE(oft_native_matrix), INTENT(in) :: A
REAL(r8), INTENT(in) :: x(:),b(:)
REAL(r8) :: res,ax
INTEGER(i4) :: i,j
res=0.d0
DO i=1,A%nr
  ax=0.d0
  DO j=A%kr(i),A%kr(i+1)-1
    ax=ax+A%M(j)*x(A%lc(j))
  END DO
  res=MAX(res,ABS(ax-b(i)))
END DO
res=res/MAXVAL(ABS(b))
END FUNCTION rel_resid
! Right hand side with entries of mixed sign
SUBROUTINE fill_rhs(b)
TYPE(oft_native_vector), INTENT(inout) :: b
INTEGER(i4) :: i
DO i=1,b%n
  b%v(i)=1.d0+SIN(0.7d0*REAL(i,r8))
END DO
END SUBROUTINE fill_rhs
! One application of ILU(0) solves systems with no dropped fill to round-off
SUBROUTINE test_exact()
TYPE(oft_native_matrix), POINTER :: A
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu
INTEGER(i4) :: bw
DO bw=1,2
  A=>make_banded(40,bw,1.d0)
  u=>make_vec(40); g=>make_vec(40)
  CALL fill_rhs(g)
  ilu%A=>A
  CALL ilu%apply(u,g)
  CALL check('ILU(0) exact for banded matrix',rel_resid(A,u%v,g%v)<1.d-13)
  CALL check('No pivots replaced',ilu%npivot_replaced==0)
  !---Second application reuses the factors
  u%v=0.d0
  CALL ilu%apply(u,g)
  CALL check('ILU(0) reapplied',rel_resid(A,u%v,g%v)<1.d-13)
  CALL ilu%delete
  CALL free_mat(A); CALL free_vec(u); CALL free_vec(g)
END DO
END SUBROUTINE test_exact
! Refactor after value changes and after switching to matrices of different
! pattern and size (larger matrices overflow factor storage sized for the first)
SUBROUTINE test_refactor()
TYPE(oft_native_matrix), POINTER :: A1,A2,A3
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu
!---Tridiagonal, 30 rows
A1=>make_banded(30,1,1.d0)
u=>make_vec(30); g=>make_vec(30)
CALL fill_rhs(g)
ilu%A=>A1
CALL ilu%apply(u,g)
CALL check('Refactor: initial solve',rel_resid(A1,u%v,g%v)<1.d-13)
!---New values, same pattern: refactor on update
A1%M=2.d0*A1%M
CALL ilu%update
CALL ilu%apply(u,g)
CALL check('Refactor: new values',rel_resid(A1,u%v,g%v)<1.d-13)
CALL free_vec(u); CALL free_vec(g)
!---Larger pentadiagonal matrix (different pattern, more rows and non-zeros), no update call
A2=>make_banded(60,2,1.d0)
u=>make_vec(60); g=>make_vec(60)
CALL fill_rhs(g)
ilu%A=>A2
CALL ilu%apply(u,g)
CALL check('Refactor: larger matrix, new pattern detected',rel_resid(A2,u%v,g%v)<1.d-13)
!---Same size and number of non-zeros but a different pattern (new pattern object)
A3=>make_banded(60,2,0.5d0)
ilu%A=>A3
CALL ilu%apply(u,g)
CALL check('Refactor: new pattern object detected',rel_resid(A3,u%v,g%v)<1.d-13)
!---Explicit request for a new symbolic setup
CALL ilu%update(new_pattern=.TRUE.)
CALL ilu%apply(u,g)
CALL check('Refactor: update(new_pattern)',rel_resid(A3,u%v,g%v)<1.d-13)
CALL free_vec(u); CALL free_vec(g)
!---Back to the smaller matrix
u=>make_vec(30); g=>make_vec(30)
CALL fill_rhs(g)
ilu%A=>A1
CALL ilu%apply(u,g)
CALL check('Refactor: smaller matrix',rel_resid(A1,u%v,g%v)<1.d-13)
CALL ilu%delete
CALL free_mat(A1); CALL free_mat(A2); CALL free_mat(A3)
CALL free_vec(u); CALL free_vec(g)
END SUBROUTINE test_refactor
! Small pivot replacement: [[1,1],[1,1+eps]] has second pivot eps, replaced by
! sign(eps)*pivot_tol*max_j |a_2j|
SUBROUTINE test_pivots()
TYPE(oft_native_matrix), POINTER :: A
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu
REAL(r8) :: eps,piv,x2
INTEGER(i4) :: k
DO k=1,3
  SELECT CASE(k)
    CASE(1)
      eps=0.d0 ! Exact zero pivot (aborted before replacement was added)
    CASE(2)
      eps=4.d0*EPSILON(1.d0)
    CASE(3)
      eps=-4.d0*EPSILON(1.d0)
  END SELECT
  A=>make_matrix(2,[1,3,5],[1,2,1,2],[1.d0,1.d0,1.d0,1.d0+eps])
  u=>make_vec(2); g=>make_vec(2)
  g%v=[1.d0,0.d0]
  ilu%A=>A
  ilu%pivot_tol=1.d-8
  CALL ilu%apply(u,g)
  CALL check('Pivot replaced',ilu%npivot_replaced==1)
  !---Second pivot (1+eps)-1 replaced by +/- pivot_tol*(1+eps); x2 = -1/pivot
  piv=SIGN(ilu%pivot_tol*MAX(1.d0,ABS(1.d0+eps)),(1.d0+eps)-1.d0)
  x2=-1.d0/piv
  CALL check('Replaced pivot value and sign',ABS(u%v(2)-x2)<=1.d-12*ABS(x2))
  CALL check('Solve with replaced pivot',ABS(u%v(1)-(1.d0-x2))<=1.d-12*ABS(x2))
  !---Pivots above the threshold are kept
  IF(eps/=0.d0)THEN
    ilu%pivot_tol=0.1d0*ABS(eps)
    CALL ilu%update
    CALL ilu%apply(u,g)
    CALL check('Pivot above threshold kept',ilu%npivot_replaced==0)
  END IF
  CALL ilu%delete
  CALL free_mat(A); CALL free_vec(u); CALL free_vec(g)
END DO
END SUBROUTINE test_pivots
! Diagonal shift: ILU(0) of a banded matrix with shift alpha solves (A + alpha diag(A)) exactly
SUBROUTINE test_shift()
TYPE(oft_native_matrix), POINTER :: A,As
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu
INTEGER(i4) :: i,j
REAL(r8), PARAMETER :: alpha=0.25d0
A=>make_banded(40,1,1.d0)
As=>make_banded(40,1,1.d0)
DO i=1,As%nr
  DO j=As%kr(i),As%kr(i+1)-1
    IF(As%lc(j)==i)As%M(j)=(1.d0+alpha)*As%M(j)
  END DO
END DO
u=>make_vec(40); g=>make_vec(40)
CALL fill_rhs(g)
ilu%A=>A
ilu%diag_shift=alpha
CALL ilu%apply(u,g)
CALL check('Shifted ILU(0) solves shifted matrix',rel_resid(As,u%v,g%v)<1.d-13)
CALL check('Shifted ILU(0) does not solve original matrix',rel_resid(A,u%v,g%v)>1.d-3)
CALL ilu%delete
CALL free_mat(A); CALL free_mat(As); CALL free_vec(u); CALL free_vec(g)
END SUBROUTINE test_shift
! Copies of a factored solver (as made by block-Jacobi with ALLOCATE(source=)) are independent
SUBROUTINE test_copy()
TYPE(oft_native_matrix), POINTER :: A
TYPE(oft_native_vector), POINTER :: u,g
CLASS(oft_solver), POINTER :: ilu,ilu_copy
A=>make_banded(40,1,1.d0)
u=>make_vec(40); g=>make_vec(40)
CALL fill_rhs(g)
ALLOCATE(oft_ilusolver::ilu)
ilu%A=>A
CALL ilu%apply(u,g)
ALLOCATE(ilu_copy, source=ilu)
!---Deleting the original must leave the copy's factors intact
CALL ilu%delete
DEALLOCATE(ilu)
u%v=0.d0
CALL ilu_copy%apply(u,g)
CALL check('Copy of factored solver',rel_resid(A,u%v,g%v)<1.d-13)
CALL ilu_copy%delete
DEALLOCATE(ilu_copy)
CALL free_mat(A); CALL free_vec(u); CALL free_vec(g)
END SUBROUTINE test_copy
! Matrices without a pattern use a diagonal index private to the solver
SUBROUTINE test_no_pattern()
TYPE(oft_native_matrix), POINTER :: A
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu
A=>make_banded(40,2,1.d0,no_pattern=.TRUE.)
CALL check('Matrix without pattern',.NOT.ASSOCIATED(A%pattern))
u=>make_vec(40); g=>make_vec(40)
CALL fill_rhs(g)
ilu%A=>A
CALL ilu%apply(u,g)
CALL check('ILU(0) exact for matrix without pattern',rel_resid(A,u%v,g%v)<1.d-13)
A%M=2.d0*A%M
CALL ilu%update
CALL ilu%apply(u,g)
CALL check('ILU(0) refactor for matrix without pattern',rel_resid(A,u%v,g%v)<1.d-13)
CALL ilu%delete
CALL free_mat(A); CALL free_vec(u); CALL free_vec(g)
END SUBROUTINE test_no_pattern
! Square banded CSR structure with all entries of the band (half-bandwidth bw)
SUBROUTINE banded_pattern(n,bw,kr,lc,dloc)
INTEGER(i4), INTENT(in) :: n,bw
INTEGER(i4), ALLOCATABLE, INTENT(out) :: kr(:),lc(:),dloc(:)
INTEGER(i4) :: i,j,k
ALLOCATE(kr(n+1),dloc(n))
kr(1)=1
DO i=1,n
  kr(i+1)=kr(i)+MIN(n,i+bw)-MAX(1,i-bw)+1
END DO
ALLOCATE(lc(kr(n+1)-1))
DO i=1,n
  k=kr(i)
  DO j=MAX(1,i-bw),MIN(n,i+bw)
    lc(k)=j
    IF(j==i)dloc(i)=k
    k=k+1
  END DO
END DO
END SUBROUTINE banded_pattern
! Relative residual of a complex CSR system, ||A x - b||_inf / ||b||_inf
FUNCTION cresid(n,kr,lc,a,x,b) RESULT(res)
INTEGER(i4), INTENT(in) :: n,kr(:),lc(:)
COMPLEX(c8), INTENT(in) :: a(:),x(:),b(:)
REAL(r8) :: res
COMPLEX(c8) :: ax
INTEGER(i4) :: i,j
res=0.d0
DO i=1,n
  ax=(0.d0,0.d0)
  DO j=kr(i),kr(i+1)-1
    ax=ax+a(j)*x(lc(j))
  END DO
  res=MAX(res,ABS(ax-b(i)))
END DO
res=res/MAXVAL(ABS(b))
END FUNCTION cresid
! ILU(0) kernels called directly: exact factorization of complex (and real) banded
! matrices with no dropped fill, with and without a diagonal shift, and replacement of
! small pivots keeping the phase (complex) or sign (real) of the pivot
SUBROUTINE test_kernels()
INTEGER(i4), PARAMETER :: n=30
INTEGER(i4) :: bw,i,j,k,nrep,ierr
INTEGER(i4), ALLOCATABLE :: kr(:),lc(:),dloc(:)
COMPLEX(c8), ALLOCATABLE :: a(:),as(:),fac(:),x(:),b(:)
REAL(r8), ALLOCATABLE :: ar(:),facr(:),xr(:)
COMPLEX(c8) :: d,piv,c2(4),cfac(4),cx(2)
REAL(r8) :: pthresh,r2(4),rfac(4),rx(2)
DO bw=1,2
  CALL banded_pattern(n,bw,kr,lc,dloc)
  ALLOCATE(a(kr(n+1)-1),as(kr(n+1)-1),fac(kr(n+1)-1),x(n),b(n))
  DO i=1,n
    DO j=kr(i),kr(i+1)-1
      IF(lc(j)==i)THEN
        a(j)=CMPLX(2.d0*bw+1.5d0+0.1d0*SIN(REAL(i,r8)),0.5d0*COS(REAL(i,r8)),c8)
      ELSE
        a(j)=CMPLX(-1.d0+0.3d0*COS(REAL(i+2*lc(j),r8)),0.4d0*SIN(REAL(i-lc(j),r8)),c8)
      END IF
    END DO
    b(i)=CMPLX(1.d0+SIN(0.7d0*REAL(i,r8)),COS(0.3d0*REAL(i,r8)),c8)
  END DO
  CALL nk_ilu0_factor(n,kr,lc,dloc,a,1.d-12,0.d0,fac,nrep,ierr)
  CALL check('Complex kernel: factorization',(ierr==0).AND.(nrep==0))
  x=b
  CALL nk_ilu0_solve(n,kr,lc,dloc,fac,x)
  CALL check('Complex kernel: exact for banded matrix',cresid(n,kr,lc,a,x,b)<1.d-13)
  !---Shifted factorization solves (A + alpha diag(A)) exactly
  as=a
  as(dloc)=1.25d0*a(dloc)
  CALL nk_ilu0_factor(n,kr,lc,dloc,a,1.d-12,0.25d0,fac,nrep,ierr)
  x=b
  CALL nk_ilu0_solve(n,kr,lc,dloc,fac,x)
  CALL check('Complex kernel: shifted matrix',cresid(n,kr,lc,as,x,b)<1.d-13)
  !---Real instantiation on the real part
  ALLOCATE(ar(SIZE(a)),facr(SIZE(a)),xr(n))
  ar=REAL(a,r8)
  CALL nk_ilu0_factor(n,kr,lc,dloc,ar,1.d-12,0.d0,facr,nrep,ierr)
  xr=REAL(b,r8)
  CALL nk_ilu0_solve(n,kr,lc,dloc,facr,xr)
  CALL check('Real kernel: exact for banded matrix', &
    cresid(n,kr,lc,CMPLX(ar,0.d0,c8),CMPLX(xr,0.d0,c8),CMPLX(REAL(b,r8),0.d0,c8))<1.d-13)
  DEALLOCATE(kr,lc,dloc,a,as,fac,x,b,ar,facr,xr)
END DO
!---Pivot replacement on [[1,0],[1,d]]: the second pivot is exactly d (no cancellation)
! and the largest magnitude in row 2 is 1, so |d| < pivot_tol replaces it
ALLOCATE(kr(3),lc(4),dloc(2))
kr=[1,3,5]; lc=[1,2,1,2]; dloc=[1,4]
pthresh=1.d-8
DO k=1,3
  SELECT CASE(k)
    CASE(1)
      d=(0.d0,0.d0)
    CASE(2)
      d=CMPLX(3.d-12,4.d-12,c8)
    CASE(3)
      d=CMPLX(-1.2d-11,5.d-12,c8)
  END SELECT
  c2=[(1.d0,0.d0),(0.d0,0.d0),(1.d0,0.d0),d]
  CALL nk_ilu0_factor(2,kr,lc,dloc,c2,pthresh,0.d0,cfac,nrep,ierr)
  CALL check('Complex kernel: pivot replaced',(ierr==0).AND.(nrep==1))
  !---Expected pivot: pivot_tol*d/|d|, or +pivot_tol for d=0
  IF(d==(0.d0,0.d0))THEN
    piv=CMPLX(pthresh,0.d0,c8)
  ELSE
    piv=pthresh*d/ABS(d)
  END IF
  CALL check('Complex kernel: replaced pivot keeps phase',ABS(cfac(4)*piv-(1.d0,0.d0))<1.d-14)
  !---Solve with b = [1,0]: x2 = -1/piv, x1 = 1
  cx=[(1.d0,0.d0),(0.d0,0.d0)]
  CALL nk_ilu0_solve(2,kr,lc,dloc,cfac,cx)
  CALL check('Complex kernel: solve with replaced pivot', &
    (ABS(cx(2)+1.d0/piv)<=1.d-14*ABS(1.d0/piv)).AND.(ABS(cx(1)-1.d0)<1.d-14))
  !---Real instantiation keeps the sign
  r2=[1.d0,0.d0,1.d0,REAL(d,r8)]
  CALL nk_ilu0_factor(2,kr,lc,dloc,r2,pthresh,0.d0,rfac,nrep,ierr)
  IF(REAL(d,r8)==0.d0)THEN
    CALL check('Real kernel: zero pivot replaced by +tol',(nrep==1).AND.(rfac(4)==1.d0/pthresh))
  ELSE
    CALL check('Real kernel: replaced pivot keeps sign',(nrep==1).AND. &
      (rfac(4)==1.d0/SIGN(pthresh,REAL(d,r8))))
  END IF
  rx=[1.d0,0.d0]
  CALL nk_ilu0_solve(2,kr,lc,dloc,rfac,rx)
  CALL check('Real kernel: solve with replaced pivot',ABS(rx(2)+rfac(4))<=1.d-14*ABS(rfac(4)))
END DO
!---Exact zero pivot without a threshold, and a non-finite pivot
c2=[(1.d0,0.d0),(0.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0)]
CALL nk_ilu0_factor(2,kr,lc,dloc,c2,0.d0,0.d0,cfac,nrep,ierr)
CALL check('Complex kernel: zero pivot error',ierr==2)
!---|HUGE+i HUGE| overflows
c2(4)=CMPLX(HUGE(1.d0),HUGE(1.d0),c8)
CALL nk_ilu0_factor(2,kr,lc,dloc,c2,0.d0,0.d0,cfac,nrep,ierr)
CALL check('Complex kernel: non-finite pivot error',ierr==-2)
DEALLOCATE(kr,lc,dloc)
END SUBROUTINE test_kernels
! Poisson equation with GMRES, preconditioned by ILU(0) directly (case 2 on one process,
! case 3) and by block-Jacobi with boundary overlap and ILU(0) (case 2)
SUBROUTINE test_lap(case_id)
INTEGER(i4), INTENT(in) :: case_id
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_bjprecond), TARGET :: bj_pre
TYPE(oft_ilusolver), TARGET :: ilu_pre,ilu_sub
CLASS(oft_vector), POINTER :: u,v,b
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
REAL(r8) :: res
INTEGER(i4) :: k
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
CALL ML_oft_lagrange%vec_create(b)
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
CALL u%set(1.d0)
CALL mop%apply(u,b)
CALL lag_zerob%apply(b)
linv%A=>lop
linv%its=-3
linv%nrits=20
DO k=1,2
  IF(case_id==2.AND.k==1.AND.oft_env%nprocs>1)CYCLE
  IF(case_id==3.AND.k==2)CYCLE
  ilu_pre%package=ilu_pkg
  ilu_sub%package=ilu_pkg
  IF(k==1)THEN
    ilu_pre%A=>lop
    linv%pre=>ilu_pre
  ELSE
    bj_pre%A=>lop
    bj_pre%boundary_overlap=.TRUE.
    bj_pre%nlocal=1
    bj_pre%pre=>ilu_sub
    linv%pre=>bj_pre
  END IF
  CALL u%set(0.d0)
  CALL v%add(0.d0,1.d0,b)
  CALL linv%apply(u,v)
  its(k)=linv%cits
  !---Check true residual
  CALL lop%apply(u,v)
  CALL v%add(1.d0,-1.d0,b)
  res=v%dot(v)
  CALL check('GMRES converged',SQRT(res)<1.d-6*SQRT(b%dot(b)))
  CALL linv%pre%delete
  NULLIFY(linv%pre)
END DO
CALL u%delete; CALL v%delete; CALL b%delete
DEALLOCATE(u,v,b)
CALL lop%delete; CALL mop%delete
DEALLOCATE(lop,mop)
CALL linv%delete
END SUBROUTINE test_lap
! Time ILU(0) factorization and solve for the Laplacian: the minimum over several
! batches of the time per solve, and of the time per refactor plus solve minus the
! time per solve for the factorization
SUBROUTINE test_timing()
INTEGER(i4), PARAMETER :: nbatch=5
TYPE(oft_ilusolver) :: ilu
CLASS(oft_vector), POINTER :: u,g
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_native_matrix), POINTER :: A_native
TYPE(oft_timer) :: timer
INTEGER(i4) :: i,k,nfactor,nsolve
REAL(r8) :: tfactor,tsolve,tsetup
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(g)
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
IF(.NOT.native_matrix_cast(A_native,lop))CALL oft_abort('Native matrix required','test_timing',__FILE__)
CALL g%set(1.d0)
ilu%A=>lop
CALL timer%tick
CALL ilu%apply(u,g)
tsetup=timer%tock()
!---Batches of ~0.2 s
nsolve=MAX(10,INT(0.2d0/MAX(tsetup,1.d-6)))
nfactor=MAX(2,nsolve/4)
tsolve=HUGE(1.d0)
tfactor=HUGE(1.d0)
DO k=1,nbatch
  CALL timer%tick
  DO i=1,nsolve
    CALL ilu%apply(u,g)
  END DO
  tsolve=MIN(tsolve,timer%tock()/REAL(nsolve,r8))
  CALL timer%tick
  DO i=1,nfactor
    CALL ilu%update
    CALL ilu%apply(u,g)
  END DO
  tfactor=MIN(tfactor,timer%tock()/REAL(nfactor,r8))
END DO
tfactor=tfactor-tsolve
CALL check('Timing solution finite',u%dot(u)<HUGE(1.d0))
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,I10,A,I12)')'ILU timing: nr = ',A_native%nr,', nnz = ',A_native%nnz
  WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4)')'ILU timing (s): first apply = ',tsetup, &
    ', factor = ',tfactor,', solve = ',tsolve
  OPEN(NEWUNIT=io_unit,FILE='ilu_timing.results')
  WRITE(io_unit,*)A_native%nr,A_native%nnz,nfactor,nsolve
  WRITE(io_unit,*)tsetup,tfactor,tsolve
  CLOSE(io_unit)
END IF
CALL timing_kernels(A_native,MAX(2,nfactor/2),nsolve)
CALL ilu%delete
CALL u%delete; CALL g%delete
DEALLOCATE(u,g)
CALL lop%delete
DEALLOCATE(lop)
END SUBROUTINE test_timing
! Time the ILU(0) kernels for 1, 2, 4, ... threads (up to the OpenMP maximum): the
! (serial) factorization, and the serial and level-scheduled solves with the profitability
! heuristic overridden (min over batches, per call). Also prints the level structure of
! the matrix.
SUBROUTINE timing_kernels(A,nfactor,nsolve)
CLASS(oft_native_matrix), INTENT(in) :: A
INTEGER(i4), INTENT(in) :: nfactor,nsolve
INTEGER(i4), PARAMETER :: nbatch=5
TYPE(native_ilu_sched) :: sched
TYPE(oft_timer) :: timer
INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: dloc
REAL(r8), ALLOCATABLE, DIMENSION(:) :: fac,sfac,x,b
INTEGER(i4) :: i,k,nt,nthreads_max,nrep,ierr
REAL(r8) :: tf,ts(2)
LOGICAL :: use_solve
dloc=>native_pattern_diag(A%pattern)
sched=native_pattern_ilu_sched(A%pattern)
use_solve=sched%use_solve
sched%use_solve=.TRUE.
ALLOCATE(fac(A%nnz),sfac(A%nnz),x(A%nr),b(A%nr))
b=1.d0
nthreads_max=omp_get_max_threads()
IF(oft_env%head_proc)THEN
  CALL native_ilu_sched_stats(sched,'ILU schedule: ')
  WRITE(*,'(A,L2)')'ILU schedule: heuristic would use threaded solve: ',use_solve
  WRITE(*,'(A)')'ILU kernels:  threads          factor (s)    solve serial (s)    solve sched (s)'
END IF
nt=1
DO WHILE(nt<=nthreads_max)
  CALL omp_set_num_threads(nt)
  tf=HUGE(1.d0); ts=HUGE(1.d0)
  DO k=1,nbatch
    CALL timer%tick
    DO i=1,nfactor
      CALL nk_ilu0_factor(A%nr,A%kr,A%lc,dloc,A%M,1.d-12,0.d0,fac,nrep,ierr)
    END DO
    tf=MIN(tf,timer%tock()/REAL(nfactor,r8))
    CALL timer%tick
    DO i=1,nsolve
      x=b
      CALL nk_ilu0_solve(A%nr,A%kr,A%lc,dloc,fac,x)
    END DO
    ts(1)=MIN(ts(1),timer%tock()/REAL(nsolve,r8))
    CALL nk_ilu0_sched_values(A%nr,A%kr,dloc,fac,sched,sfac)
    CALL timer%tick
    DO i=1,nsolve
      x=b
      CALL nk_ilu0_solve(A%nr,A%kr,A%lc,dloc,fac,x,sched,sfac)
    END DO
    ts(2)=MIN(ts(2),timer%tock()/REAL(nsolve,r8))
  END DO
  IF(oft_env%head_proc)WRITE(*,'(A,I8,3ES20.4)')'ILU kernels: ',nt,tf,ts
  nt=2*nt
END DO
CALL omp_set_num_threads(nthreads_max)
CALL native_ilu_sched_free(sched)
DEALLOCATE(fac,sfac,x,b)
END SUBROUTINE timing_kernels
! Bitwise equality of real and complex arrays
FUNCTION same_r(x,y) RESULT(same)
REAL(r8), INTENT(in) :: x(:),y(:)
LOGICAL :: same
same=(SIZE(x)==SIZE(y))
IF(same)same=ALL(TRANSFER(x,0_i8,SIZE(x))==TRANSFER(y,0_i8,SIZE(y)))
END FUNCTION same_r
FUNCTION same_c(x,y) RESULT(same)
COMPLEX(c8), INTENT(in) :: x(:),y(:)
LOGICAL :: same
same=(SIZE(x)==SIZE(y))
IF(same)same=ALL(TRANSFER(x,0_i8,2*SIZE(x))==TRANSFER(y,0_i8,2*SIZE(y)))
END FUNCTION same_c
! Level of each row in one sweep of a schedule
SUBROUTINE row_levels(nr,ls,level)
INTEGER(i4), INTENT(in) :: nr
TYPE(native_level_sched), INTENT(in) :: ls
INTEGER(i4), INTENT(out) :: level(nr)
INTEGER(i4) :: lev,k
DO lev=1,ls%nlevels
  DO k=ls%lptr(lev),ls%lptr(lev+1)-1
    level(ls%perm(k))=lev
  END DO
END DO
END SUBROUTINE row_levels
! Structure of level schedules: each sweep is a permutation grouped by level with every
! dependency in an earlier level, segments cover the levels in order with single-level
! parallel segments of at least min_rows rows, and the reordered layout matches the rows
FUNCTION sched_valid(A,sched,min_rows) RESULT(valid)
CLASS(oft_native_matrix), INTENT(in) :: A
TYPE(native_ilu_sched), INTENT(in) :: sched
INTEGER(i4), INTENT(in) :: min_rows
LOGICAL :: valid
INTEGER(i4) :: i,j,k,nl,nu
INTEGER(i4), ALLOCATABLE :: lvl(:),lvu(:),cnt(:)
valid=(sched%nr==A%nr).AND.(sched%nnz==A%nnz)
IF(.NOT.valid)RETURN
ALLOCATE(lvl(A%nr),lvu(A%nr),cnt(A%nr))
lvl=0; lvu=0
CALL row_levels(A%nr,sched%lower,lvl)
CALL row_levels(A%nr,sched%upper,lvu)
valid=ALL(lvl>0).AND.ALL(lvu>0)
!---Permutations with inverse
cnt=0
DO k=1,A%nr
  cnt(sched%lower%perm(k))=cnt(sched%lower%perm(k))+1
  valid=valid.AND.(sched%lower%ipos(sched%lower%perm(k))==k).AND.(sched%upper%ipos(sched%upper%perm(k))==k)
END DO
valid=valid.AND.ALL(cnt==1)
!---Dependencies in earlier levels, and level 1 rows have none
nl=0; nu=0
DO i=1,A%nr
  DO j=A%kr(i),A%kr(i+1)-1
    IF(A%lc(j)<i)THEN
      valid=valid.AND.(lvl(A%lc(j))<lvl(i))
      nl=nl+1
    ELSE IF(A%lc(j)>i)THEN
      valid=valid.AND.(lvu(A%lc(j))<lvu(i))
      nu=nu+1
    END IF
  END DO
END DO
valid=valid.AND.(sched%nnzl==nl).AND.(sched%lower%kp(A%nr+1)-1==nl).AND.(sched%upper%kp(A%nr+1)-1==A%nnz)
DO k=1,A%nr
  i=sched%lower%perm(k)
  valid=valid.AND.(sched%lower%kp(k+1)-sched%lower%kp(k)==COUNT(A%lc(A%kr(i):A%kr(i+1)-1)<i))
  i=sched%upper%perm(k)
  valid=valid.AND.(sched%col(sched%upper%kp(k))==sched%lower%ipos(i))
END DO
!---Segments
valid=valid.AND.seg_valid(A%nr,sched%lower,min_rows).AND.seg_valid(A%nr,sched%upper,min_rows)
DEALLOCATE(lvl,lvu,cnt)
END FUNCTION sched_valid
! Segments of one sweep cover its levels in order, with single-level parallel segments
! of at least min_rows rows
FUNCTION seg_valid(nr,ls,min_rows) RESULT(ok)
INTEGER(i4), INTENT(in) :: nr
TYPE(native_level_sched), INTENT(in) :: ls
INTEGER(i4), INTENT(in) :: min_rows
LOGICAL :: ok
INTEGER(i4) :: lev,npar,s
ok=(ls%sptr(1)==1).AND.(ls%sptr(ls%nseg+1)==nr+1).AND.(ls%lptr(ls%nlevels+1)==nr+1)
npar=0
DO s=1,ls%nseg
  ok=ok.AND.(ls%sptr(s+1)>ls%sptr(s)).AND.ANY(ls%lptr==ls%sptr(s))
  IF(ls%spar(s))THEN
    npar=npar+ls%sptr(s+1)-ls%sptr(s)
    ok=ok.AND.(ls%sptr(s+1)-ls%sptr(s)>=min_rows)
    DO lev=1,ls%nlevels
      IF(ls%lptr(lev)==ls%sptr(s))ok=ok.AND.(ls%lptr(lev+1)==ls%sptr(s+1))
    END DO
  END IF
END DO
ok=ok.AND.(npar==ls%npar)
END FUNCTION seg_valid
! Factor with the (serial) real kernel, then solve with the serial and the scheduled
! kernels and compare bitwise; returns the number of replaced pivots and the error flag
SUBROUTINE compare_real(A,dloc,vals,b,sched,ptol,shift,name,nrep,ierr)
CLASS(oft_native_matrix), INTENT(in) :: A
INTEGER(i4), INTENT(in) :: dloc(:)
REAL(r8), INTENT(in) :: vals(:),b(:)
TYPE(native_ilu_sched), INTENT(in) :: sched
REAL(r8), INTENT(in) :: ptol,shift
CHARACTER(LEN=*), INTENT(in) :: name
INTEGER(i4), INTENT(out) :: nrep,ierr
REAL(r8), ALLOCATABLE :: fac(:),sfac(:),x0(:),x1(:)
ALLOCATE(fac(A%nnz),sfac(A%nnz),x0(A%nr),x1(A%nr))
CALL nk_ilu0_factor(A%nr,A%kr,A%lc,dloc,vals,ptol,shift,fac,nrep,ierr)
CALL check(name//': factorization succeeded',ierr==0)
IF(ierr==0)THEN
  x0=b
  CALL nk_ilu0_solve(A%nr,A%kr,A%lc,dloc,fac,x0)
  CALL nk_ilu0_sched_values(A%nr,A%kr,dloc,fac,sched,sfac)
  x1=b
  CALL nk_ilu0_solve(A%nr,A%kr,A%lc,dloc,fac,x1,sched,sfac)
  CALL check(name//': identical solution',same_r(x0,x1))
END IF
DEALLOCATE(fac,sfac,x0,x1)
END SUBROUTINE compare_real
! Complex version of compare_real
SUBROUTINE compare_cplx(A,dloc,vals,b,sched,ptol,shift,name,nrep,ierr)
CLASS(oft_native_matrix), INTENT(in) :: A
INTEGER(i4), INTENT(in) :: dloc(:)
COMPLEX(c8), INTENT(in) :: vals(:),b(:)
TYPE(native_ilu_sched), INTENT(in) :: sched
REAL(r8), INTENT(in) :: ptol,shift
CHARACTER(LEN=*), INTENT(in) :: name
INTEGER(i4), INTENT(out) :: nrep,ierr
COMPLEX(c8), ALLOCATABLE :: fac(:),sfac(:),x0(:),x1(:)
ALLOCATE(fac(A%nnz),sfac(A%nnz),x0(A%nr),x1(A%nr))
CALL nk_ilu0_factor(A%nr,A%kr,A%lc,dloc,vals,ptol,shift,fac,nrep,ierr)
CALL check(name//': factorization succeeded',ierr==0)
IF(ierr==0)THEN
  x0=b
  CALL nk_ilu0_solve(A%nr,A%kr,A%lc,dloc,fac,x0)
  CALL nk_ilu0_sched_values(A%nr,A%kr,dloc,fac,sched,sfac)
  x1=b
  CALL nk_ilu0_solve(A%nr,A%kr,A%lc,dloc,fac,x1,sched,sfac)
  CALL check(name//': identical solution',same_c(x0,x1))
END IF
DEALLOCATE(fac,sfac,x0,x1)
END SUBROUTINE compare_cplx
! Level-scheduled ILU(0) solve and solver against the serial ones (bitwise)
SUBROUTINE test_threaded()
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_native_matrix), POINTER :: A
CLASS(oft_vector), POINTER :: u1,u2,g
TYPE(native_ilu_sched) :: sched
TYPE(native_ilu_sched), POINTER :: psched
TYPE(oft_ilusolver) :: ilu1,ilu2
INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: dloc
REAL(r8), ALLOCATABLE :: av(:),b(:),fac0(:),sfac0(:),xt(:),x0(:)
REAL(r8), POINTER, DIMENSION(:) :: v1,v2
COMPLEX(c8), ALLOCATABLE :: cav(:),cb(:)
INTEGER(i4) :: i,j,k,nt,nthreads,nrep,ierr
INTEGER(i4) :: minr(3)
LOGICAL :: ok,in_par
CHARACTER(LEN=40) :: tag
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
IF(.NOT.native_matrix_cast(A,lop))CALL oft_abort('Native matrix required','test_threaded',__FILE__)
dloc=>native_pattern_diag(A%pattern)
nthreads=omp_get_max_threads()
!---Non-symmetric real and complex values with the pattern of the operator
ALLOCATE(av(A%nnz),cav(A%nnz),b(A%nr),cb(A%nr))
DO i=1,A%nr
  DO j=A%kr(i),A%kr(i+1)-1
    av(j)=A%M(j)*(1.d0+0.2d0*SIN(REAL(i+3*A%lc(j),r8)))
    cav(j)=av(j)*CMPLX(1.d0,0.3d0*COS(REAL(j,r8)),c8)
  END DO
  b(i)=1.d0+SIN(0.7d0*REAL(i,r8))
  cb(i)=CMPLX(b(i),COS(0.3d0*REAL(i,r8)),c8)
END DO
!---Schedules with all levels parallel, mixed parallel and serial segments, and the default
minr=[1,8,native_ilu_level_rows]
DO k=1,3
  WRITE(tag,'(A,I0)')'min_rows=',minr(k)
  CALL native_ilu_sched_build(A%nr,A%kr,A%lc,sched,min_rows=minr(k))
  CALL check('Schedule structure, '//TRIM(tag),sched_valid(A,sched,minr(k)))
  !---Force the scheduled solve (the heuristic would not use it on small matrices)
  sched%use_solve=.TRUE.
  CALL check('Scheduled solve used iff threads>1, '//TRIM(tag), &
    native_ilu_sched_active(sched,A%nr).EQV.(nthreads>1))
  CALL check('Schedule rejected for other sizes',.NOT.native_ilu_sched_active(sched,A%nr+1))
  CALL compare_real(A,dloc,av,b,sched,1.d-12,0.d0,'Real '//TRIM(tag),nrep,ierr)
  CALL check('Real: no error or replaced pivots',(ierr==0).AND.(nrep==0))
  !---Large pivot threshold (many replaced pivots) and diagonal shift
  CALL compare_real(A,dloc,av,b,sched,0.9d0,0.1d0,'Real pivots '//TRIM(tag),nrep,ierr)
  CALL check('Real: pivots replaced',(ierr==0).AND.(nrep>0).AND.(nrep<A%nr))
  IF(oft_env%head_proc)WRITE(*,'(3A,I0,A,I0)')'Replaced pivots (',TRIM(tag),'): ',nrep,' of ',A%nr
  CALL compare_cplx(A,dloc,cav,cb,sched,1.d-12,0.d0,'Complex '//TRIM(tag),nrep,ierr)
  CALL compare_cplx(A,dloc,cav,cb,sched,0.9d0,0.1d0,'Complex pivots '//TRIM(tag),nrep,ierr)
  CALL check('Complex: pivots replaced',(ierr==0).AND.(nrep>0).AND.(nrep<A%nr))
END DO
!---Other thread counts in the same run (oversubscribed if there are fewer cores)
DO nt=2,4
  CALL omp_set_num_threads(nt)
  WRITE(tag,'(A,I0)')'threads=',nt
  CALL compare_real(A,dloc,av,b,sched,0.9d0,0.1d0,'Real '//TRIM(tag),nrep,ierr)
  CALL compare_cplx(A,dloc,cav,cb,sched,0.9d0,0.1d0,'Complex '//TRIM(tag),nrep,ierr)
END DO
CALL omp_set_num_threads(nthreads)
!---Inside a parallel region (eg. block-Jacobi with several blocks) the solve runs serially
ALLOCATE(fac0(A%nnz),sfac0(A%nnz),x0(A%nr))
CALL nk_ilu0_factor(A%nr,A%kr,A%lc,dloc,av,0.9d0,0.1d0,fac0,nrep,ierr)
CALL nk_ilu0_sched_values(A%nr,A%kr,dloc,fac0,sched,sfac0)
x0=b
CALL nk_ilu0_solve(A%nr,A%kr,A%lc,dloc,fac0,x0)
ok=.TRUE.
!$omp parallel num_threads(2) private(xt,in_par) reduction(.AND.:ok)
in_par=omp_in_parallel()
ALLOCATE(xt(A%nr))
ok=.NOT.(in_par.AND.native_ilu_sched_active(sched,A%nr))
xt=b
CALL nk_ilu0_solve(A%nr,A%kr,A%lc,dloc,fac0,xt,sched,sfac0)
ok=ok.AND.same_r(xt,x0)
DEALLOCATE(xt)
!$omp end parallel
CALL check('Serial path inside parallel region',ok)
DEALLOCATE(fac0,sfac0,x0)
!---No threading for patterns without a diagonal entry in every row (3x3, none in row 2)
CALL native_ilu_sched_build(3,[1,3,4,6],[1,3,1,2,3],sched,min_rows=1)
sched%use_solve=.TRUE.
CALL check('Schedule unused without diagonal',.NOT.native_ilu_sched_active(sched,3))
CALL native_ilu_sched_free(sched)
!---Pattern schedules are built once and shared
psched=>native_pattern_ilu_sched(A%pattern)
CALL check('Pattern schedule shared',ASSOCIATED(psched,native_pattern_ilu_sched(A%pattern)))
CALL check('Pattern schedule structure',sched_valid(A,psched,native_ilu_level_rows))
!---Threaded solver against the serial one (scheduled solve forced on the pattern schedule)
psched%use_solve=.TRUE.
CALL ML_oft_lagrange%vec_create(u1)
CALL ML_oft_lagrange%vec_create(u2)
CALL ML_oft_lagrange%vec_create(g)
NULLIFY(v1)
CALL g%get_local(v1)
v1=b
CALL g%restore_local(v1)
ilu1%A=>lop
ilu2%A=>lop
ilu2%threaded=.FALSE.
ilu1%pivot_tol=0.9d0
ilu2%pivot_tol=0.9d0
CALL ilu1%apply(u1,g)
CALL ilu2%apply(u2,g)
CALL check('Threaded solver: same replaced pivots',(ilu1%npivot_replaced==ilu2%npivot_replaced) &
  .AND.(ilu1%npivot_replaced>0))
CALL check('Threaded solver: schedule used iff threads>1', &
  ALLOCATED(ilu1%native_struct%sfac).EQV.(nthreads>1))
CALL check('Serial solver: no schedule',.NOT.ALLOCATED(ilu2%native_struct%sfac))
!---Second application (reused factors) and a refactorization
CALL ilu1%apply(u1,g)
CALL ilu1%update
CALL ilu1%apply(u1,g)
NULLIFY(v1,v2)
CALL u1%get_local(v1)
CALL u2%get_local(v2)
CALL check('Threaded solver: identical solution',same_r(v1,v2))
DEALLOCATE(v1,v2)
psched%use_solve=.FALSE.
CALL ilu1%delete; CALL ilu2%delete
CALL u1%delete; CALL u2%delete; CALL g%delete
DEALLOCATE(u1,u2,g)
DEALLOCATE(av,cav,b,cb)
CALL lop%delete
DEALLOCATE(lop)
END SUBROUTINE test_threaded
! GMRES with Schwarz (one subdomain per process, overlap 1) and ILU(0), with threading
! enabled and disabled: identical iterations and solutions, and wall times of a first
! solve (setup, factorization and iterations) and of a repeated solve
SUBROUTINE test_schwarz_timing()
CLASS(oft_vector), POINTER :: u,v,b
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
REAL(r8), POINTER, DIMENSION(:) :: vals
REAL(r8), ALLOCATABLE :: u1(:)
REAL(r8) :: times(2,2)
INTEGER(i4) :: k
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
CALL ML_oft_lagrange%vec_create(b)
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
CALL u%set(1.d0)
CALL mop%apply(u,b)
CALL lag_zerob%apply(b)
DO k=1,2
  CALL schwarz_solve(lop,u,v,b,k==1,times(:,k),its(k))
  NULLIFY(vals)
  CALL u%get_local(vals)
  IF(k==1)THEN
    ALLOCATE(u1(SIZE(vals)))
    u1=vals
  ELSE
    CALL check('Schwarz: identical solution with and without threading',same_r(u1,vals))
  END IF
  DEALLOCATE(vals)
END DO
CALL check('Schwarz: same iterations with and without threading',its(1)==its(2))
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,I0,A,I0,A,I0)')'Schwarz timing: procs = ',oft_env%nprocs,', threads = ', &
    omp_get_max_threads(),', iterations = ',its(1)
  WRITE(*,'(A,2ES12.4)')'Schwarz timing (s), threaded: first solve, repeated solve = ',times(:,1)
  WRITE(*,'(A,2ES12.4)')'Schwarz timing (s), serial:   first solve, repeated solve = ',times(:,2)
  OPEN(NEWUNIT=io_unit,FILE='ilu_timing.results')
  WRITE(io_unit,*)oft_env%nprocs,omp_get_max_threads(),its(1)
  WRITE(io_unit,*)times(:,1)
  WRITE(io_unit,*)times(:,2)
  CLOSE(io_unit)
END IF
DEALLOCATE(u1)
CALL u%delete; CALL v%delete; CALL b%delete
DEALLOCATE(u,v,b)
CALL lop%delete; CALL mop%delete
DEALLOCATE(lop,mop)
END SUBROUTINE test_schwarz_timing
! Two solves of lop u = b with GMRES and one Schwarz/ILU(0) preconditioner (threaded or
! not); wall times and iterations
SUBROUTINE schwarz_solve(lop,u,v,b,threaded,times,nits)
CLASS(oft_matrix), POINTER, INTENT(in) :: lop
CLASS(oft_vector), INTENT(inout) :: u,v,b
LOGICAL, INTENT(in) :: threaded
REAL(r8), INTENT(out) :: times(2)
INTEGER(i4), INTENT(out) :: nits
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_schwarz_precond), TARGET :: sw
TYPE(oft_ilusolver), TARGET :: sub
TYPE(oft_timer) :: timer
INTEGER(i4) :: i
sub%threaded=threaded
sw%overlap=1
sw%pre=>sub
linv%A=>lop
linv%its=-3
linv%nrits=20
linv%pre=>sw
DO i=1,2
  CALL u%set(0.d0)
  CALL v%add(0.d0,1.d0,b)
  CALL timer%tick
  CALL linv%apply(u,v)
  times(i)=timer%tock()
  IF(i==1)nits=linv%cits
END DO
CALL check('Schwarz: same iterations on repeated solve',linv%cits==nits)
CALL sw%delete
END SUBROUTINE schwarz_solve
! Non-symmetric, diagonally dominant five-point (convection-diffusion like) matrix on an
! m x m grid: ILU(0) drops fill here, a complete LU does not
FUNCTION make_grid(m,scale,no_pattern) RESULT(A)
INTEGER(i4), INTENT(in) :: m
REAL(r8), INTENT(in) :: scale
LOGICAL, OPTIONAL, INTENT(in) :: no_pattern
TYPE(oft_native_matrix), POINTER :: A
INTEGER(i4) :: i,j,k,n,row
INTEGER(i4), ALLOCATABLE :: kr(:),lc(:)
REAL(r8), ALLOCATABLE :: vals(:)
n=m*m
ALLOCATE(kr(n+1),lc(5*n),vals(5*n))
k=0
kr(1)=1
DO j=1,m
  DO i=1,m
    row=(j-1)*m+i
    IF(j>1)THEN
      k=k+1; lc(k)=row-m; vals(k)=-1.3d0
    END IF
    IF(i>1)THEN
      k=k+1; lc(k)=row-1; vals(k)=-1.2d0
    END IF
    k=k+1; lc(k)=row; vals(k)=4.5d0+0.1d0*SIN(REAL(row,r8))
    IF(i<m)THEN
      k=k+1; lc(k)=row+1; vals(k)=-0.8d0
    END IF
    IF(j<m)THEN
      k=k+1; lc(k)=row+m; vals(k)=-0.7d0
    END IF
    kr(row+1)=k+1
  END DO
END DO
A=>make_matrix(n,kr,lc(1:k),scale*vals(1:k),no_pattern)
DEALLOCATE(kr,lc,vals)
END FUNCTION make_grid
! SuperLU ILUTP solver with no dropping (a complete LU for the small matrices used here)
SUBROUTINE super_exact_opts(ilu)
TYPE(oft_ilusolver), INTENT(inout) :: ilu
ilu%package='super'
ilu%drop_tol=0.d0
ilu%fill_factor=1.d3
END SUBROUTINE super_exact_opts
! Apply a new SuperLU ILUTP solver with the options of ilu to A and return the solution
SUBROUTINE super_fresh(ilu,A,g,x,nnz)
TYPE(oft_ilusolver), INTENT(in) :: ilu
TYPE(oft_native_matrix), POINTER, INTENT(in) :: A
TYPE(oft_native_vector), INTENT(inout) :: g
REAL(r8), INTENT(out) :: x(:)
INTEGER(i8), INTENT(out) :: nnz
TYPE(oft_ilusolver) :: fresh
TYPE(oft_native_vector), POINTER :: u
fresh%package='super'
fresh%drop_tol=ilu%drop_tol
fresh%fill_factor=ilu%fill_factor
fresh%milu=ilu%milu
fresh%col_perm=ilu%col_perm
fresh%fill_tol=ilu%fill_tol
fresh%pivot_thresh=ilu%pivot_thresh
u=>make_vec(A%nr)
fresh%A=>A
CALL fresh%apply(u,g)
x=u%v
nnz=fresh%nnz_factor
CALL fresh%delete
CALL free_vec(u)
END SUBROUTINE super_fresh
! With no dropping SuperLU ILUTP is a complete LU (where ILU(0) is not), for one and
! several right hand sides and for matrices without a pattern; the matrix values are
! not modified (they are passed to SuperLU without a copy)
SUBROUTINE test_super_exact()
TYPE(oft_native_matrix), POINTER :: A
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu,ilu0
REAL(r8), ALLOCATABLE :: vsave(:),rhs(:,:)
INTEGER(i4) :: i,k,n
REAL(r8) :: res
A=>make_grid(8,1.d0)
n=A%nr
ALLOCATE(vsave(A%nnz))
vsave=A%M
u=>make_vec(n); g=>make_vec(n)
CALL fill_rhs(g)
!---ILU(0) drops fill for this matrix
ilu0%package='native'
ilu0%A=>A
CALL ilu0%apply(u,g)
CALL check('Super: ILU(0) not exact for grid matrix',rel_resid(A,u%v,g%v)>1.d-6)
CALL ilu0%delete
!---No dropping
CALL super_exact_opts(ilu)
ilu%A=>A
u%v=0.d0
CALL ilu%apply(u,g)
CALL check('Super: exact without dropping',rel_resid(A,u%v,g%v)<1.d-13)
CALL check('Super: no pivots replaced',ilu%npivot_replaced==0)
CALL check('Super: factors have fill',ilu%nnz_factor>A%nnz)
CALL check('Super: matrix values unchanged',same_r(vsave,A%M))
u%v=0.d0
CALL ilu%apply(u,g)
CALL check('Super: reapplied',rel_resid(A,u%v,g%v)<1.d-13)
CALL ilu%delete
!---Several right hand sides solved together
ilu%nrhs=3
ALLOCATE(ilu%sec_rhs(n,2),rhs(n,2))
DO k=1,2
  DO i=1,n
    rhs(i,k)=COS(0.3d0*REAL(k*(n-i),r8))
  END DO
END DO
ilu%sec_rhs=rhs
ilu%A=>A
u%v=0.d0
CALL ilu%apply(u,g)
res=rel_resid(A,u%v,g%v)
DO k=1,2
  res=MAX(res,rel_resid(A,ilu%sec_rhs(:,k),rhs(:,k)))
END DO
CALL check('Super: several right hand sides',res<1.d-13)
CALL ilu%delete
DEALLOCATE(ilu%sec_rhs,rhs)
ilu%nrhs=1
CALL free_mat(A)
!---Matrix without a pattern (keyed on its arrays)
A=>make_grid(8,1.d0,no_pattern=.TRUE.)
ilu%A=>A
u%v=0.d0
CALL ilu%apply(u,g)
CALL check('Super: exact for matrix without pattern',rel_resid(A,u%v,g%v)<1.d-13)
A%M=2.d0*A%M
CALL ilu%update
CALL ilu%apply(u,g)
CALL check('Super: refactor for matrix without pattern',rel_resid(A,u%v,g%v)<1.d-13)
CALL ilu%delete
CALL free_mat(A); CALL free_vec(u); CALL free_vec(g)
DEALLOCATE(vsave)
END SUBROUTINE test_super_exact
! Refactoring: new values with the same pattern reuse the column permutation and
! elimination tree (SamePattern), which must give the same factors as a new
! factorization (bitwise) even with dropping; pattern and size changes are detected
SUBROUTINE test_super_refactor()
TYPE(oft_native_matrix), POINTER :: A1,A2,A3
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu
REAL(r8), ALLOCATABLE :: x(:)
INTEGER(i4) :: i,j
INTEGER(i8) :: nnz
A1=>make_grid(8,1.d0)
u=>make_vec(A1%nr); g=>make_vec(A1%nr)
ALLOCATE(x(A1%nr))
CALL fill_rhs(g)
ilu%package='super'
ilu%drop_tol=5.d-2
ilu%A=>A1
CALL ilu%apply(u,g)
CALL check('Super refactor: dropping is inexact',rel_resid(A1,u%v,g%v)>1.d-8)
!---New values (row dependent, so the dropped entries change): SamePattern refactor
DO i=1,A1%nr
  DO j=A1%kr(i),A1%kr(i+1)-1
    IF(A1%lc(j)/=i)A1%M(j)=A1%M(j)*(1.d0+0.5d0*COS(REAL(3*i+j,r8)))
  END DO
END DO
CALL ilu%update
CALL ilu%apply(u,g)
CALL super_fresh(ilu,A1,g,x,nnz)
CALL check('Super refactor: new values = new factorization',same_r(u%v,x))
CALL check('Super refactor: new values, same fill',ilu%nnz_factor==nnz)
!---No update: factors are kept
A1%M=2.d0*A1%M
CALL ilu%apply(u,g)
CALL check('Super refactor: no update keeps factors',same_r(u%v,x))
CALL free_vec(u); CALL free_vec(g)
!---Larger matrix (new pattern and size), detected without update
CALL super_exact_opts(ilu)
A2=>make_grid(10,1.d0)
u=>make_vec(A2%nr); g=>make_vec(A2%nr)
CALL fill_rhs(g)
ilu%A=>A2
CALL ilu%apply(u,g)
CALL check('Super refactor: larger matrix, new pattern detected',rel_resid(A2,u%v,g%v)<1.d-13)
!---Same size, different pattern object
A3=>make_banded(A2%nr,3,0.5d0)
ilu%A=>A3
CALL ilu%apply(u,g)
CALL check('Super refactor: new pattern object detected',rel_resid(A3,u%v,g%v)<1.d-13)
!---Explicit request for a new symbolic setup
A3%M=3.d0*A3%M
CALL ilu%update(new_pattern=.TRUE.)
CALL ilu%apply(u,g)
CALL check('Super refactor: update(new_pattern)',rel_resid(A3,u%v,g%v)<1.d-13)
CALL free_vec(u); CALL free_vec(g)
!---Back to the smaller matrix
u=>make_vec(A1%nr); g=>make_vec(A1%nr)
CALL fill_rhs(g)
ilu%A=>A1
CALL ilu%apply(u,g)
CALL check('Super refactor: smaller matrix',rel_resid(A1,u%v,g%v)<1.d-13)
CALL ilu%delete
CALL free_mat(A1); CALL free_mat(A2); CALL free_mat(A3)
CALL free_vec(u); CALL free_vec(g)
DEALLOCATE(x)
END SUBROUTINE test_super_refactor
! Option changes: numerical options take effect at the next update, the column
! permutation at the next symbolic setup (triggered by the change); each result
! matches a new solver with the same options bitwise
SUBROUTINE test_super_options()
TYPE(oft_native_matrix), POINTER :: A
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu
REAL(r8), ALLOCATABLE :: x(:)
INTEGER(i8) :: nnz,nnz_ref
INTEGER(i4) :: k
A=>make_grid(10,1.d0)
u=>make_vec(A%nr); g=>make_vec(A%nr)
ALLOCATE(x(A%nr))
CALL fill_rhs(g)
ilu%package='super'
ilu%drop_tol=1.d-2
ilu%A=>A
CALL ilu%apply(u,g)
nnz_ref=ilu%nnz_factor
DO k=1,6
  SELECT CASE(k)
    CASE(1)
      ilu%drop_tol=1.d-1
    CASE(2)
      ilu%fill_factor=2.d0
    CASE(3)
      ilu%milu=2
    CASE(4)
      ilu%pivot_thresh=1.d0
    CASE(5)
      ilu%fill_tol=1.d-4
    CASE(6)
      ilu%col_perm=0
  END SELECT
  IF(k<6)CALL ilu%update
  CALL ilu%apply(u,g)
  CALL super_fresh(ilu,A,g,x,nnz)
  CALL check('Super options: changed option = new solver',same_r(u%v,x).AND.(ilu%nnz_factor==nnz))
  IF(k==1)CALL check('Super options: larger drop tolerance, fewer non-zeros',ilu%nnz_factor<nnz_ref)
  IF(k==2)CALL check('Super options: fill bounded by fill_factor',ilu%nnz_factor<=2*A%nnz+A%nr)
END DO
CALL ilu%delete
CALL free_mat(A); CALL free_vec(u); CALL free_vec(g)
DEALLOCATE(x)
END SUBROUTINE test_super_options
! Zero pivots: SuperLU replaces a pivot only when its whole remaining column is zero, so
! the singular [[1,1],[1,1]] has one replaced pivot (counted, factorization completes);
! the zero diagonal of [[0,1],[1,0]] is never chosen as a pivot, even when diagonal
! pivots are always preferred (pivot_thresh=0)
SUBROUTINE test_super_pivots()
TYPE(oft_native_matrix), POINTER :: A
TYPE(oft_native_vector), POINTER :: u,g
TYPE(oft_ilusolver) :: ilu
A=>make_matrix(2,[1,3,5],[1,2,1,2],[1.d0,1.d0,1.d0,1.d0])
u=>make_vec(2); g=>make_vec(2)
g%v=[1.d0,2.d0]
CALL super_exact_opts(ilu)
ilu%col_perm=0
ilu%A=>A
CALL ilu%apply(u,g)
CALL check('Super pivots: zero pivot replaced',ilu%npivot_replaced==1)
CALL check('Super pivots: solution finite',ALL(ABS(u%v)<HUGE(1.d0)))
CALL ilu%delete
CALL free_mat(A)
A=>make_matrix(2,[1,3,5],[1,2,1,2],[0.d0,1.d0,1.d0,0.d0])
ilu%pivot_thresh=0.d0
ilu%A=>A
CALL ilu%apply(u,g)
CALL check('Super pivots: zero diagonal not used',(ilu%npivot_replaced==0).AND.(rel_resid(A,u%v,g%v)<1.d-14))
CALL ilu%delete
CALL free_mat(A); CALL free_vec(u); CALL free_vec(g)
END SUBROUTINE test_super_pivots
! Copies of a factored solver (as made by block-Jacobi with ALLOCATE(source=)) do not
! share the SuperLU factors: each copy factors anew, and deleting either leaves the
! other intact; deleting twice, or before factoring, is safe
SUBROUTINE test_super_copy()
TYPE(oft_native_matrix), POINTER :: A
TYPE(oft_native_vector), POINTER :: u,g
CLASS(oft_solver), POINTER :: ilu,ilu_copy,ilu_copy2
TYPE(oft_ilusolver) :: tmp
A=>make_grid(8,1.d0)
u=>make_vec(A%nr); g=>make_vec(A%nr)
CALL fill_rhs(g)
CALL super_exact_opts(tmp)
ALLOCATE(ilu, source=tmp)
!---Copy of a solver that was never applied
ALLOCATE(ilu_copy2, source=ilu)
ilu%A=>A
CALL ilu%apply(u,g)
CALL check('Super copy: original',rel_resid(A,u%v,g%v)<1.d-13)
ALLOCATE(ilu_copy, source=ilu)
!---Both in use
u%v=0.d0
CALL ilu_copy%apply(u,g)
CALL check('Super copy: copy factors anew',rel_resid(A,u%v,g%v)<1.d-13)
u%v=0.d0
CALL ilu%apply(u,g)
CALL check('Super copy: original after copy',rel_resid(A,u%v,g%v)<1.d-13)
!---Deleting the original leaves the copy intact
CALL ilu%delete
CALL ilu%delete
DEALLOCATE(ilu)
u%v=0.d0
CALL ilu_copy%apply(u,g)
CALL check('Super copy: copy after original deleted',rel_resid(A,u%v,g%v)<1.d-13)
CALL ilu_copy%delete
DEALLOCATE(ilu_copy)
!---Copy of the template, applied after its original is gone
ilu_copy2%A=>A
u%v=0.d0
CALL ilu_copy2%apply(u,g)
CALL check('Super copy: copy of unapplied solver',rel_resid(A,u%v,g%v)<1.d-13)
CALL ilu_copy2%delete
!---Reuse after delete
ilu_copy2%A=>A
u%v=0.d0
CALL ilu_copy2%apply(u,g)
CALL check('Super copy: apply after delete',rel_resid(A,u%v,g%v)<1.d-13)
CALL ilu_copy2%delete
DEALLOCATE(ilu_copy2)
CALL free_mat(A); CALL free_vec(u); CALL free_vec(g)
END SUBROUTINE test_super_copy
! Native ILU(0) against SuperLU ILUTP for several drop tolerances (and fill bounds and MILU
! variants) on the Laplacian (one process): time of the first application (setup and factorization), refactorization and
! solve (minimum over batches), non-zeros of the factors and GMRES iterations
SUBROUTINE test_compare()
INTEGER(i4), PARAMETER :: ncfg=8,nbatch=3
CHARACTER(LEN=7), PARAMETER :: pkgs(ncfg)=['native','super ','super ','super ','super ','super ','super ','super ']
REAL(r8), PARAMETER :: dtols(ncfg)=[0.d0,1.d-2,1.d-3,1.d-4,1.d-4,1.d-2,1.d-3,1.d-3]
REAL(r8), PARAMETER :: ffacs(ncfg)=[10.d0,10.d0,10.d0,10.d0,30.d0,10.d0,10.d0,10.d0]
INTEGER(i4), PARAMETER :: milus(ncfg)=[0,0,0,0,0,2,2,1]
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_ilusolver), TARGET :: ilu
CLASS(oft_vector), POINTER :: u,v,b
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_native_matrix), POINTER :: A_native
TYPE(oft_timer) :: timer
INTEGER(i4) :: k,i,j,nsolve,nfactor,nits
REAL(r8) :: tfirst,tfactor,tsolve,res
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
CALL ML_oft_lagrange%vec_create(b)
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
IF(.NOT.native_matrix_cast(A_native,lop))CALL oft_abort('Native matrix required','test_compare',__FILE__)
CALL u%set(1.d0)
CALL mop%apply(u,b)
CALL lag_zerob%apply(b)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='ilu_compare.results')
  WRITE(io_unit,*)A_native%nr,A_native%nnz
  WRITE(*,'(A,I10,A,I12)')'ILU comparison: nr = ',A_native%nr,', nnz = ',A_native%nnz
  WRITE(*,'(A)')'ILU comparison: package  drop_tol  fill milu   first (s)  factor (s)   solve (s)  nnz(L+U)/nnz(A)  GMRES its'
END IF
DO k=1,ncfg
  ilu%package=pkgs(k)
  ilu%drop_tol=dtols(k)
  ilu%fill_factor=ffacs(k)
  ilu%milu=milus(k)
  ilu%A=>lop
  CALL v%set(1.d0)
  CALL timer%tick
  CALL ilu%apply(u,v)
  tfirst=timer%tock()
  nsolve=MAX(5,INT(0.2d0/MAX(tfirst,1.d-6)))
  nfactor=MAX(2,nsolve/4)
  tsolve=HUGE(1.d0); tfactor=HUGE(1.d0)
  DO j=1,nbatch
    CALL timer%tick
    DO i=1,nsolve
      CALL ilu%apply(u,v)
    END DO
    tsolve=MIN(tsolve,timer%tock()/REAL(nsolve,r8))
    CALL timer%tick
    DO i=1,nfactor
      CALL ilu%update
      CALL ilu%apply(u,v)
    END DO
    tfactor=MIN(tfactor,timer%tock()/REAL(nfactor,r8))
  END DO
  tfactor=tfactor-tsolve
  !---GMRES with the factored preconditioner
  linv%A=>lop
  linv%its=-3
  linv%nrits=20
  linv%pre=>ilu
  CALL u%set(0.d0)
  CALL v%add(0.d0,1.d0,b)
  CALL linv%apply(u,v)
  nits=linv%cits
  CALL lop%apply(u,v)
  CALL v%add(1.d0,-1.d0,b)
  res=v%dot(v)
  CALL check('Comparison: GMRES converged',SQRT(res)<1.d-6*SQRT(b%dot(b)))
  CALL check('Comparison: iterations',nits>0)
  IF(oft_env%head_proc)THEN
    WRITE(*,'(A,A8,ES9.1,F6.0,I5,3ES12.3,F17.3,I11)')'ILU comparison: ',pkgs(k),dtols(k),ffacs(k),milus(k), &
      tfirst,tfactor,tsolve,REAL(ilu%nnz_factor,r8)/REAL(A_native%nnz,r8),nits
    WRITE(io_unit,*)k,tfirst,tfactor,tsolve,ilu%nnz_factor,nits
  END IF
  NULLIFY(linv%pre)
  CALL linv%delete
  CALL ilu%delete
END DO
IF(oft_env%head_proc)CLOSE(io_unit)
CALL u%delete; CALL v%delete; CALL b%delete
DEALLOCATE(u,v,b)
CALL lop%delete; CALL mop%delete
DEALLOCATE(lop,mop)
END SUBROUTINE test_compare
END PROGRAM test_native_ilu
