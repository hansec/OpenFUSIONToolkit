!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_native_ilu.F90
!
!> Regression tests for the native ILU(0) solver (@ref oft_lu::oft_ilusolver).
!!
!! The current test cases are (selected by `case_id` in `test_ilu_options`):
!! - 1: Local matrices built directly: ILU(0) of banded matrices with a full band
!!   (no fill outside the pattern) is an exact LU factorization, refactoring after
!!   value changes and after switching to matrices with a different pattern or
!!   size (with and without `update(new_pattern=.TRUE.)`), matrices without a
!!   pattern, small-pivot replacement (sign preserving), diagonal shift, deep copies
!!   of a factored solver, and the real and complex ILU(0) kernels called directly.
!! - 2: Poisson equation on a unit cube with GMRES: ILU(0) applied directly to the
!!   operator (single process only) and inside block-Jacobi with `boundary_overlap`
!!   (full local matrix of the distributed operator).
!! - 3: ILU(0) applied directly to a distributed operator (aborts on more than one
!!   process, solves on one).
!! - 4: Timing of the ILU(0) factorization and solve for the Lagrange Laplacian of
!!   order `order` on a single process (written to `ilu_timing.results`).
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
  native_pattern_create, native_matrix_set_pattern, native_matrix_cast
USE oft_native_solvers, ONLY: oft_bjprecond, oft_native_gmres_solver
USE oft_solver_base, ONLY: oft_solver
USE oft_lu, ONLY: oft_ilusolver
USE oft_native_kernels, ONLY: nk_ilu0_factor, nk_ilu0_solve
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_zerob, oft_lag_getlop, oft_lag_getmop
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr,nfail
INTEGER(i4) :: case_id = 1
INTEGER(i4) :: order = 2
INTEGER(i4) :: its(2) = -1
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
TYPE(oft_seam), TARGET :: local_seam
NAMELIST/test_ilu_options/case_id,order
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
CALL ilu%delete
CALL u%delete; CALL g%delete
DEALLOCATE(u,g)
CALL lop%delete
DEALLOCATE(lop)
END SUBROUTINE test_timing
END PROGRAM test_native_ilu
