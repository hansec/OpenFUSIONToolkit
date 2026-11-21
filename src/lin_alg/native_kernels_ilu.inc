!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Module body template for incomplete LU kernels, included by native_kernels.F90
!
! Required macros (undefined by the including file after each instantiation):
!   KT     - Value type (eg. REAL(r8))
!   KZERO  - Zero of type KT
!   KONE   - One of type KT
!
! The including file provides the MODULE and END MODULE statements (so the module
! name is visible to CMake's dependency scanner). Specific procedures are private;
! only the generic names are public, so generics of the same name from each
! instantiation merge in oft_native_kernels.
!
! Factor storage: the factors share the CSR layout (kr, lc) of the matrix. fac(j) is
! the L entry for j < dloc(i) (unit diagonal not stored), the U entry for j > dloc(i),
! and fac(dloc(i)) holds the inverse of the (possibly replaced) pivot of row i.
!
! Threading: the factorization is serial. The solve takes an optional level schedule
! (native_ilu_sched, built by native_ilu_sched_build). Without it, or when
! native_ilu_sched_active is false (single thread, already inside an active parallel
! region, or a schedule that is not expected to pay off), it runs serially and may
! be called from inside a parallel region. With it, it opens one parallel region and
! processes the rows of each level concurrently, with exactly the per-row arithmetic
! of the serial path, so results are bitwise identical for any number of threads.
! The threaded solve reads a copy of the factors reordered to follow the schedules
! (nk_ilu0_sched_values).
!---------------------------------------------------------------------------------
USE omp_lib
USE oft_local
USE oft_native_kernels_base, ONLY: native_ilu_sched, native_ilu_sched_active, native_level_segment
IMPLICIT NONE
PRIVATE
PUBLIC :: nk_ilu0_factor, nk_ilu0_solve, nk_ilu0_sched_values
INTERFACE nk_ilu0_factor
  MODULE PROCEDURE k_ilu0_factor
END INTERFACE nk_ilu0_factor
INTERFACE nk_ilu0_solve
  MODULE PROCEDURE k_ilu0_solve
END INTERFACE nk_ilu0_solve
INTERFACE nk_ilu0_sched_values
  MODULE PROCEDURE k_ilu0_sched_values
END INTERFACE nk_ilu0_sched_values
CONTAINS
!---------------------------------------------------------------------------------
!> Incomplete LU factorization with no fill (ILU(0)) in the CSR layout of the matrix
!!
!! Row-wise (IKJ) elimination restricted to the non-zero pattern, following
!! Y. Saad, "Iterative Methods for Sparse Linear Systems", 2nd ed., SIAM (2003),
!! Sec. 10.3. The factors are written to `fac`, which is
!! aligned with `kr`/`lc`: entries left of the diagonal hold L (unit diagonal not
!! stored), entries right of it hold U, and the diagonal entry `fac(dloc(i))` holds
!! the inverse of the pivot of row `i`.
!!
!! Pivot protection: with \f$ \tau = \f$ `pivot_tol` times the largest magnitude in
!! row `i` of `a`, a pivot \f$ u \f$ with \f$ |u| < \tau \f$ is replaced by
!! \f$ \tau u / |u| \f$, which keeps the sign of a real pivot and the phase of a
!! complex one, or by \f$ +\tau \f$ for an exact zero. Magnitudes are absolute
!! values (moduli for complex values). A zero pivot that is not replaced
!! (`pivot_tol=0` or an all-zero row) stops the factorization with `ierr=i`, a
!! non-finite pivot with `ierr=-i`. With `diag_shift` \f$ = \alpha \f$ the factored
!! matrix is \f$ A + \alpha \mathrm{diag}(A) \f$.
!!
!! @note Requires a square CSR pattern with strictly increasing column indices in
!! each row and a diagonal entry in every row; this is not checked here.
!! @note Serial (no OpenMP directives), unlike most kernels, so it may be called from
!! inside a parallel region (eg. block-Jacobi with several local blocks per process
!! factors its blocks in a parallel loop).
!---------------------------------------------------------------------------------
SUBROUTINE k_ilu0_factor(nr,kr,lc,dloc,a,pivot_tol,diag_shift,fac,nreplaced,ierr)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), CONTIGUOUS, INTENT(in) :: lc(:) !< CSR column indices [nnz]
INTEGER(i4), INTENT(in) :: dloc(nr) !< Location of the diagonal entry in each row [nr]
KT, CONTIGUOUS, INTENT(in) :: a(:) !< Matrix values [nnz]
REAL(r8), INTENT(in) :: pivot_tol !< Relative small-pivot threshold (0 to only catch exact zeros)
REAL(r8), INTENT(in) :: diag_shift !< Relative diagonal shift \f$ \alpha \f$
KT, CONTIGUOUS, INTENT(out) :: fac(:) !< Factor values in the layout of `a` [nnz]
INTEGER(i4), INTENT(out) :: nreplaced !< Number of pivots replaced
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0: success, i>0: zero pivot in row i, i<0: non-finite pivot in row -i)
INTEGER(i4) :: i,j,jj,jd,jrow,jw
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: iw
REAL(r8) :: rowmax,pthresh
KT :: tl,piv
ierr=0
nreplaced=0
!---Location of each column of the current row in fac (0 if not in the row)
ALLOCATE(iw(nr))
iw=0
DO i=1,nr
  jd=dloc(i)
  rowmax=0.d0
  DO j=kr(i),kr(i+1)-1
    rowmax=MAX(rowmax,ABS(a(j)))
    fac(j)=a(j)
    iw(lc(j))=j
  END DO
  fac(jd)=a(jd)*(1.d0+diag_shift)
  !---Eliminate the L part of the row, in increasing column order
  DO j=kr(i),jd-1
    jrow=lc(j)
    tl=fac(j)*fac(dloc(jrow))
    fac(j)=tl
    DO jj=dloc(jrow)+1,kr(jrow+1)-1
      jw=iw(lc(jj))
      IF(jw/=0)fac(jw)=fac(jw)-tl*fac(jj)
    END DO
  END DO
  !---Check and invert pivot
  piv=fac(jd)
  IF(.NOT.(ABS(piv)<=HUGE(rowmax)))THEN
    ierr=-i
    EXIT
  END IF
  pthresh=pivot_tol*rowmax
  IF(ABS(piv)<pthresh)THEN
    IF(piv==KZERO)THEN
      piv=pthresh
    ELSE
      piv=(piv/ABS(piv))*pthresh
    END IF
    nreplaced=nreplaced+1
  ELSE IF(piv==KZERO)THEN
    ierr=i
    EXIT
  END IF
  fac(jd)=KONE/piv
  DO j=kr(i),kr(i+1)-1
    iw(lc(j))=0
  END DO
END DO
DEALLOCATE(iw)
END SUBROUTINE k_ilu0_factor
!---------------------------------------------------------------------------------
!> Approximate solve with ILU(0) factors from @ref k_ilu0_factor, \f$ x = (LU)^{-1} x \f$
!!
!! Forward substitution with L in increasing row order, then backward substitution
!! with U in decreasing row order, in place.
!!
!! Threading: with `sched` and `sfac` (see @ref oft_native_kernels_base::native_ilu_sched_active
!! "native_ilu_sched_active" for when they are used), the rows of each level of the
!! forward and backward schedules are processed concurrently in one parallel region,
!! reading the factors from `sfac` (the copy of `fac` reordered by
!! @ref k_ilu0_sched_values). Each row is computed exactly as in the serial path, with
!! its entries in the same order, from finished values, so the result is bitwise
!! identical for any number of threads.
!!
!! @note Without `sched` the kernel is serial (no OpenMP directives), so it may be
!! called from inside a parallel region (see @ref k_ilu0_factor).
!---------------------------------------------------------------------------------
SUBROUTINE k_ilu0_solve(nr,kr,lc,dloc,fac,x,sched,sfac)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), CONTIGUOUS, INTENT(in) :: lc(:) !< CSR column indices [nnz]
INTEGER(i4), INTENT(in) :: dloc(nr) !< Location of the diagonal entry in each row [nr]
KT, CONTIGUOUS, INTENT(in) :: fac(:) !< Factor values from @ref k_ilu0_factor [nnz]
KT, INTENT(inout) :: x(nr) !< Right-hand side on input, solution on output [nr]
TYPE(native_ilu_sched), OPTIONAL, INTENT(in) :: sched !< Level schedules for threaded execution (optional, with `sfac`)
KT, CONTIGUOUS, OPTIONAL, INTENT(in) :: sfac(:) !< Factor values reordered by @ref k_ilu0_sched_values [nnz] (optional)
INTEGER(i4) :: i,j
KT :: tmp
IF(PRESENT(sched).AND.PRESENT(sfac))THEN
  IF(native_ilu_sched_active(sched,nr))THEN
    CALL ilu0_solve_sched(nr,sfac,x,sched)
    RETURN
  END IF
END IF
!---Forward substitution (unit lower triangle)
DO i=1,nr
  tmp=x(i)
  DO j=kr(i),dloc(i)-1
    tmp=tmp-fac(j)*x(lc(j))
  END DO
  x(i)=tmp
END DO
!---Backward substitution (upper triangle, inverted pivots on the diagonal)
DO i=nr,1,-1
  tmp=x(i)
  DO j=dloc(i)+1,kr(i+1)-1
    tmp=tmp-fac(j)*x(lc(j))
  END DO
  x(i)=fac(dloc(i))*tmp
END DO
END SUBROUTINE k_ilu0_solve
!---------------------------------------------------------------------------------
!> Reorder ILU(0) factors to follow level schedules, for the threaded solve
!!
!! Copies the factors from the CSR layout of the matrix (`fac`, from
!! @ref k_ilu0_factor) to the layout of the reordered values described in
!! @ref oft_native_kernels_base::native_ilu_sched "native_ilu_sched". Threaded when
!! called from serial code on a large matrix.
!---------------------------------------------------------------------------------
SUBROUTINE k_ilu0_sched_values(nr,kr,dloc,fac,sched,sfac)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), INTENT(in) :: dloc(nr) !< Location of the diagonal entry in each row [nr]
KT, CONTIGUOUS, INTENT(in) :: fac(:) !< Factor values from @ref k_ilu0_factor [nnz]
TYPE(native_ilu_sched), INTENT(in) :: sched !< Level schedules
KT, CONTIGUOUS, INTENT(out) :: sfac(:) !< Reordered factor values [nnz]
INTEGER(i4) :: i,j,k,jp
LOGICAL :: threaded
threaded=(nr>OFT_OMP_VTHRESH).AND.(.NOT.omp_in_parallel())
!$omp parallel if(threaded) private(i,j,k,jp)
!$omp do schedule(static)
DO k=1,nr
  i=sched%lower%perm(k)
  jp=sched%lower%kp(k)
  DO j=kr(i),dloc(i)-1
    sfac(jp)=fac(j)
    jp=jp+1
  END DO
END DO
!$omp end do nowait
!$omp do schedule(static)
DO k=1,nr
  i=sched%upper%perm(k)
  jp=sched%upper%kp(k)
  DO j=dloc(i),kr(i+1)-1
    sfac(jp)=fac(j)
    jp=jp+1
  END DO
END DO
!$omp end do
!$omp end parallel
END SUBROUTINE k_ilu0_sched_values
!---------------------------------------------------------------------------------
!> Level-scheduled ILU(0) solve with reordered factors (see @ref k_ilu0_solve)
!!
!! Must be called from serial code. One parallel region covers both sweeps; each
!! execution segment of a schedule ends with a barrier. The forward sweep reads `x`
!! and stores its result by forward-schedule position, the backward sweep stores its
!! result by backward-schedule position, and the solution is finally copied to `x`
!! (see @ref oft_native_kernels_base::native_ilu_sched "native_ilu_sched").
!---------------------------------------------------------------------------------
SUBROUTINE ilu0_solve_sched(nr,sfac,x,sched)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
KT, CONTIGUOUS, INTENT(in) :: sfac(:) !< Reordered factor values [nnz]
KT, INTENT(inout) :: x(nr) !< Right-hand side on input, solution on output [nr]
TYPE(native_ilu_sched), INTENT(in) :: sched !< Level schedules
INTEGER(i4) :: i,s,k1,k2
KT, ALLOCATABLE, DIMENSION(:) :: yl,yu
ALLOCATE(yl(nr),yu(nr))
!$omp parallel private(i,s,k1,k2)
DO s=1,sched%lower%nseg
  CALL native_level_segment(sched%lower,s,k1,k2)
  CALL lsolve_rows(sched%lower%perm,sched%lower%kp,sched%col,k1,k2,sfac,x,yl)
  !$omp barrier
END DO
DO s=1,sched%upper%nseg
  CALL native_level_segment(sched%upper,s,k1,k2)
  CALL usolve_rows(sched%upper%kp,sched%col,k1,k2,sfac,yl,yu)
  !$omp barrier
END DO
!$omp do simd schedule(static)
DO i=1,nr
  x(i)=yu(sched%upper%ipos(i))
END DO
!$omp end parallel
DEALLOCATE(yl,yu)
END SUBROUTINE ilu0_solve_sched
!---------------------------------------------------------------------------------
!> Forward substitution for the rows at positions `k1:k2` of the forward schedule, in
!! order (part of @ref ilu0_solve_sched)
!---------------------------------------------------------------------------------
SUBROUTINE lsolve_rows(perm,kp,col,k1,k2,sfac,x,yl)
INTEGER(i4), CONTIGUOUS, INTENT(in) :: perm(:) !< Rows in execution order [nr]
INTEGER(i4), CONTIGUOUS, INTENT(in) :: kp(:) !< Start of the entries of each position in `sfac` [nr+1]
INTEGER(i4), CONTIGUOUS, INTENT(in) :: col(:) !< Schedule positions of the columns of `sfac` [nnz]
INTEGER(i4), INTENT(in) :: k1 !< First position
INTEGER(i4), INTENT(in) :: k2 !< Last position
KT, CONTIGUOUS, INTENT(in) :: sfac(:) !< Reordered factor values [nnz]
KT, CONTIGUOUS, INTENT(in) :: x(:) !< Right-hand side [nr]
KT, CONTIGUOUS, INTENT(inout) :: yl(:) !< Result by forward-schedule position [nr]
INTEGER(i4) :: j,k
KT :: tmp
DO k=k1,k2
  tmp=x(perm(k))
  DO j=kp(k),kp(k+1)-1
    tmp=tmp-sfac(j)*yl(col(j))
  END DO
  yl(k)=tmp
END DO
END SUBROUTINE lsolve_rows
!---------------------------------------------------------------------------------
!> Backward substitution for the rows at positions `k1:k2` of the backward schedule,
!! in order (part of @ref ilu0_solve_sched)
!---------------------------------------------------------------------------------
SUBROUTINE usolve_rows(kp,col,k1,k2,sfac,yl,yu)
INTEGER(i4), CONTIGUOUS, INTENT(in) :: kp(:) !< Start of the entries (inverted pivot first) of each position in `sfac` [nr+1]
INTEGER(i4), CONTIGUOUS, INTENT(in) :: col(:) !< Schedule positions of the columns of `sfac` [nnz]
INTEGER(i4), INTENT(in) :: k1 !< First position
INTEGER(i4), INTENT(in) :: k2 !< Last position
KT, CONTIGUOUS, INTENT(in) :: sfac(:) !< Reordered factor values [nnz]
KT, CONTIGUOUS, INTENT(in) :: yl(:) !< Result of the forward sweep by forward-schedule position [nr]
KT, CONTIGUOUS, INTENT(inout) :: yu(:) !< Result by backward-schedule position [nr]
INTEGER(i4) :: j,k
KT :: tmp
DO k=k1,k2
  tmp=yl(col(kp(k)))
  DO j=kp(k)+1,kp(k+1)-1
    tmp=tmp-sfac(j)*yu(col(j))
  END DO
  yu(k)=sfac(kp(k))*tmp
END DO
END SUBROUTINE usolve_rows
