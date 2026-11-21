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
!---------------------------------------------------------------------------------
USE oft_local
IMPLICIT NONE
PRIVATE
PUBLIC :: nk_ilu0_factor, nk_ilu0_solve
INTERFACE nk_ilu0_factor
  MODULE PROCEDURE k_ilu0_factor
END INTERFACE nk_ilu0_factor
INTERFACE nk_ilu0_solve
  MODULE PROCEDURE k_ilu0_solve
END INTERFACE nk_ilu0_solve
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
!! @note Serial (no OpenMP directives), so it may be called from inside a parallel
!! region (see @ref k_ilu0_factor). Threaded (level-scheduled) solves are to be added
!! as an optional trailing schedule argument, so existing calls keep the serial path.
!---------------------------------------------------------------------------------
SUBROUTINE k_ilu0_solve(nr,kr,lc,dloc,fac,x)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), CONTIGUOUS, INTENT(in) :: lc(:) !< CSR column indices [nnz]
INTEGER(i4), INTENT(in) :: dloc(nr) !< Location of the diagonal entry in each row [nr]
KT, CONTIGUOUS, INTENT(in) :: fac(:) !< Factor values from @ref k_ilu0_factor [nnz]
KT, INTENT(inout) :: x(nr) !< Right-hand side on input, solution on output [nr]
INTEGER(i4) :: i,j
KT :: tmp
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
