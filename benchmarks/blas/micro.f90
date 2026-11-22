! Micro-benchmark of core BLAS/LAPACK kernels used by OFT (best of nrep)
PROGRAM micro
IMPLICIT NONE
INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15)
INTEGER :: sizes(4) = [500, 1000, 2000, 4000]
INTEGER :: i, n, info, lwork, irep, nrep
REAL(r8), ALLOCATABLE :: A(:,:), B(:,:), C(:,:), A0(:,:), w(:), wi(:), work(:), vl(:,:), vr(:,:)
COMPLEX(r8), ALLOCATABLE :: zA(:,:), zB(:,:), zC(:,:)
INTEGER, ALLOCATABLE :: ipiv(:), iwork(:)
REAL(r8) :: t, tbest, chk
INTEGER(8) :: c0, c1, crate
WRITE(*,'(A)') 'kernel        N      time[s]    GFLOP/s   check'
DO i = 1, SIZE(sizes)
  n = sizes(i)
  nrep = MERGE(3, 1, n <= 2000)
  ALLOCATE(A(n,n), B(n,n), C(n,n), A0(n,n), w(n), wi(n), ipiv(n))
  CALL RANDOM_SEED(PUT=[(42, irep=1,64)])
  CALL RANDOM_NUMBER(A0)
  A0 = A0 + TRANSPOSE(A0)
  DO irep = 1, n
    A0(irep,irep) = A0(irep,irep) + REAL(n, r8)
  END DO
  CALL RANDOM_NUMBER(B)
  ! DGEMM
  tbest = HUGE(1.0_r8)
  DO irep = 1, nrep
    CALL SYSTEM_CLOCK(c0, crate)
    CALL DGEMM('N', 'N', n, n, n, 1.0_r8, A0, n, B, n, 0.0_r8, C, n)
    CALL SYSTEM_CLOCK(c1); tbest = MIN(tbest, REAL(c1-c0, r8)/crate)
  END DO
  chk = ABS(C(n/2,n/3) - DOT_PRODUCT(A0(n/2,:), B(:,n/3)))
  WRITE(*,'(A10,I6,F12.4,F11.2,ES10.2)') 'dgemm', n, tbest, 2.0_r8*REAL(n,r8)**3/tbest/1.E9_r8, chk
  ! DGETRF + DGETRS
  tbest = HUGE(1.0_r8)
  DO irep = 1, nrep
    A = A0
    CALL SYSTEM_CLOCK(c0, crate)
    CALL DGETRF(n, n, A, n, ipiv, info)
    CALL SYSTEM_CLOCK(c1); tbest = MIN(tbest, REAL(c1-c0, r8)/crate)
  END DO
  C(:,1) = B(:,1)
  CALL DGETRS('N', n, 1, A, n, ipiv, C(:,1), n, info)
  chk = MAXVAL(ABS(MATMUL(A0, C(:,1)) - B(:,1)))
  WRITE(*,'(A10,I6,F12.4,F11.2,ES10.2)') 'dgetrf', n, tbest, 2.0_r8/3.0_r8*REAL(n,r8)**3/tbest/1.E9_r8, chk
  ! DGETRI (explicit inverse, used by ThinCurr direct solvers)
  tbest = HUGE(1.0_r8)
  DO irep = 1, nrep
    C = A
    CALL SYSTEM_CLOCK(c0, crate)
    lwork = 64*n
    ALLOCATE(work(lwork))
    CALL DGETRI(n, C, n, ipiv, work, lwork, info)
    DEALLOCATE(work)
    CALL SYSTEM_CLOCK(c1); tbest = MIN(tbest, REAL(c1-c0, r8)/crate)
  END DO
  chk = ABS(DOT_PRODUCT(C(n/2,:), A0(:,n/2)) - 1.0_r8)
  WRITE(*,'(A10,I6,F12.4,F11.2,ES10.2)') 'dgetri', n, tbest, 4.0_r8/3.0_r8*REAL(n,r8)**3/tbest/1.E9_r8, chk
  ! DSYEVD (symmetric eigensolve, vectors)
  IF (n <= 2000) THEN
    A = A0
    ALLOCATE(work(1), iwork(1))
    CALL DSYEVD('V', 'U', n, A, n, w, work, -1, iwork, -1, info)
    lwork = INT(work(1)); irep = iwork(1)
    DEALLOCATE(work, iwork); ALLOCATE(work(lwork), iwork(irep))
    CALL SYSTEM_CLOCK(c0, crate)
    CALL DSYEVD('V', 'U', n, A, n, w, work, lwork, iwork, SIZE(iwork), info)
    CALL SYSTEM_CLOCK(c1); tbest = REAL(c1-c0, r8)/crate
    chk = MAXVAL(ABS(MATMUL(A0, A(:,1)) - w(1)*A(:,1)))
    WRITE(*,'(A10,I6,F12.4,F11.2,ES10.2)') 'dsyevd', n, tbest, 4.0_r8/3.0_r8*REAL(n,r8)**3/tbest/1.E9_r8, chk
    DEALLOCATE(work, iwork)
  END IF
  ! DGEEV (nonsymmetric eigensolve, right vectors; used by ThinCurr direct eigensolver)
  IF (n <= 2000) THEN
    ALLOCATE(vl(1,1), vr(n,n), work(1))
    C = A0; C(1:n/2,:) = C(1:n/2,:) + B(1:n/2,:)  ! make nonsymmetric
    A = C
    CALL DGEEV('N', 'V', n, A, n, w, wi, vl, 1, vr, n, work, -1, info)
    lwork = INT(work(1)); DEALLOCATE(work); ALLOCATE(work(lwork))
    CALL SYSTEM_CLOCK(c0, crate)
    CALL DGEEV('N', 'V', n, A, n, w, wi, vl, 1, vr, n, work, lwork, info)
    CALL SYSTEM_CLOCK(c1); tbest = REAL(c1-c0, r8)/crate
    chk = SUM(w) - SUM([(C(irep,irep), irep=1,n)])  ! trace check
    WRITE(*,'(A10,I6,F12.4,F11.2,ES10.2)') 'dgeev', n, tbest, 0.0_r8, ABS(chk)/ABS(SUM(w))
    DEALLOCATE(vl, vr, work)
  END IF
  ! ZGEMM
  ALLOCATE(zA(n,n), zB(n,n), zC(n,n))
  zA = CMPLX(A0, B, r8); zB = CMPLX(B, A0, r8)
  tbest = HUGE(1.0_r8)
  DO irep = 1, nrep
    CALL SYSTEM_CLOCK(c0, crate)
    CALL ZGEMM('N', 'N', n, n, n, (1.0_r8,0.0_r8), zA, n, zB, n, (0.0_r8,0.0_r8), zC, n)
    CALL SYSTEM_CLOCK(c1); tbest = MIN(tbest, REAL(c1-c0, r8)/crate)
  END DO
  chk = ABS(zC(n/2,n/3) - SUM(zA(n/2,:)*zB(:,n/3)))/ABS(zC(n/2,n/3))
  WRITE(*,'(A10,I6,F12.4,F11.2,ES10.2)') 'zgemm', n, tbest, 8.0_r8*REAL(n,r8)**3/tbest/1.E9_r8, chk
  ! DGEMV (memory bound; ThinCurr time stepping)
  tbest = HUGE(1.0_r8)
  CALL SYSTEM_CLOCK(c0, crate)
  DO irep = 1, 50
    CALL DGEMV('N', n, n, 1.0_r8, A0, n, B(:,1), 1, 0.0_r8, C(:,1), 1)
  END DO
  CALL SYSTEM_CLOCK(c1); tbest = REAL(c1-c0, r8)/crate/50
  WRITE(*,'(A10,I6,F12.4,F11.2,ES10.2)') 'dgemv', n, tbest, 2.0_r8*REAL(n,r8)**2/tbest/1.E9_r8, 0.0_r8
  DEALLOCATE(A, B, C, A0, w, wi, ipiv, zA, zB, zC)
END DO
END PROGRAM micro
