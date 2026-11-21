!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file native_kernels.F90
!
!> Plain-array compute kernels for the native linear algebra backend
!!
!! All inner loops over vector values and matrix entries used by the native vector,
!! matrix and solver classes live here, so each loop exists once and threading,
!! vectorization and (later) offload annotations are applied in one place. Kernels
!! take explicit sizes and contiguous arrays and know nothing about the vector or
!! matrix classes. They contain no MPI communication; reductions return the local
!! (owned) contribution only.
!!
!! Real, complex and mixed real/complex versions are generated from a single source
!! by including a template file once per type combination, each time inside its own
!! module with different type macros:
!! - `native_kernels_vec.inc`: one value type (`nk_set`, `nk_scale`, `nk_maxpy`,
!!   `nk_norm_inf`, `nk_owned_sum`, `nk_owned_norm1`, `nk_gather`, `nk_scatter`,
!!   CSR row/column/diagonal operations and value insertion) for real and complex.
!! - `native_kernels_upd.inc`: destination and source types (`nk_copy`, `nk_axpby`,
!!   `nk_axpbypcz`, `nk_mult`, `nk_div`) for real/real, complex/real and
!!   complex/complex.
!! - `native_kernels_prod.inc`: left (or matrix), right and result types
!!   (`nk_owned_dot`, `nk_owned_mdot`, `nk_csr_matvec`, `nk_csc_matvec`,
!!   `nk_dense_matvec`, `nk_dense_matvect`) for all four real/complex combinations.
!! - `native_kernels_ilu.inc`: one value type, incomplete LU factorization and solve
!!   in the CSR layout of the matrix (`nk_ilu0_factor`, `nk_ilu0_solve`) for real and
!!   complex.
!!
!! Specific procedures in each instantiation are private and share their names; only
!! generic names are public, and the generics of each instantiation merge in
!! @ref oft_native_kernels, which is the module to use. Kernels that only exist for
!! real values (Jacobi smoother sweeps) and the type-independent helpers (including
!! the level schedules of the threaded incomplete LU kernels) are in
!! @ref oft_native_kernels_base.
!!
!! Threading follows the rest of the native backend: each kernel opens its own
!! OpenMP parallel region, conditional on the problem size exceeding
!! `OFT_OMP_VTHRESH`, so kernels must be called from serial code. The exception is
!! the incomplete LU kernels, which are serial unless given a level schedule
!! (@ref oft_native_kernels_base::native_ilu_sched "native_ilu_sched"), and then only
!! open a parallel region when called from serial code, so they may be called from
!! inside a parallel region (block-Jacobi factors and applies several local blocks
!! per process in a parallel loop).
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup doxy_oft_native_la
!---------------------------------------------------------------------------------
#include "local.h"
!---------------------------------------------------------------------------------
!> Type-independent helpers and real-only kernels for the native backend
!---------------------------------------------------------------------------------
MODULE oft_native_kernels_base
USE omp_lib
USE oft_local
IMPLICIT NONE
PRIVATE
!> Offset used to store column indices in 2 bytes (`lc_small = lc - lc_offset`)
INTEGER(i4), PUBLIC, PARAMETER :: lc_offset = HUGE(INT(1,2))
!> Number of rows processed together by multi-vector kernels (`maxpy`, `mdot`), chosen
!! so that one row block of every vector involved stays in cache
INTEGER(i4), PUBLIC, PARAMETER :: native_vec_chunk = 512
!---------------------------------------------------------------------------------
!> Contiguous pointer to real values (for multi-vector kernels)
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: native_rvals_ptr
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: v => NULL() !< Values
END TYPE native_rvals_ptr
!---------------------------------------------------------------------------------
!> Contiguous pointer to complex values (for multi-vector kernels)
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: native_cvals_ptr
  COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: v => NULL() !< Values
END TYPE native_cvals_ptr
!> Default minimum number of rows for a level of a triangular sweep to be executed by
!! all threads of a level-scheduled ILU(0) kernel; runs of consecutive smaller levels
!! are executed by one thread (see @ref native_ilu_sched_build)
INTEGER(i4), PUBLIC, PARAMETER :: native_ilu_level_rows = 64
!> Minimum average number of non-zeros per synchronization (barrier) for the
!! level-scheduled ILU(0) solve (see @ref native_ilu_sched_build)
INTEGER(i4), PUBLIC, PARAMETER :: native_ilu_solve_nnz = 4096
!> Minimum average number of rows per synchronization (barrier) for the level-scheduled
!! ILU(0) factorization (see @ref native_ilu_sched_build)
INTEGER(i4), PUBLIC, PARAMETER :: native_ilu_factor_rows = 32
!---------------------------------------------------------------------------------
!> Execution schedule for one triangular sweep, rows grouped by dependency level
!!
!! Level 1 holds the rows that depend on no other row of the sweep, and level
!! \f$ k+1 \f$ the rows whose dependencies are all in levels \f$ \le k \f$, so the
!! rows of one level may be processed concurrently once the previous levels are done.
!! Consecutive levels are combined into execution segments: a parallel segment is a
!! single level whose rows are shared between threads, a serial segment is a run of
!! small consecutive levels processed in order by one thread. Each segment ends with
!! a barrier.
!!
!! `kp` locates the entries of each row, taken in execution order, in the factor
!! values reordered for the threaded solve (see @ref native_ilu_sched), so a sweep
!! reads them sequentially, and `ipos` is the inverse of `perm`.
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: native_level_sched
  INTEGER(i4) :: nlevels = 0 !< Number of levels
  INTEGER(i4) :: maxrows = 0 !< Largest number of rows in one level
  INTEGER(i4) :: nseg = 0 !< Number of execution segments
  INTEGER(i4) :: npar = 0 !< Number of rows in parallel segments
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: perm !< Rows grouped by level, increasing row order within a level [nr]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lptr !< Start of each level in `perm` [nlevels+1]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: sptr !< Start of each segment in `perm` [nseg+1]
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: spar !< Segment is shared between threads (single level) [nseg]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: kp !< Start of the entries of row `perm(k)` in the reordered values [nr+1]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ipos !< Position of each row in `perm` [nr]
END TYPE native_level_sched
!---------------------------------------------------------------------------------
!> Level schedules for threaded ILU(0) factorization and solves on a square CSR pattern
!!
!! Built from the structure alone by @ref native_ilu_sched_build. The forward (L)
!! schedule also orders the factorization, since row `i` of ILU(0) only depends on
!! the rows referenced by its L part. `use_factor` and `use_solve` record whether the
!! threaded kernels are expected to be faster than the serial ones; the kernels only
!! take the scheduled path when the corresponding flag is set (see
!! @ref native_ilu_sched_active).
!!
!! The threaded solve works in the execution order of the schedules, so that each
!! thread reads and writes contiguous blocks of memory:
!! - It reads the factors from a copy reordered to follow the schedules (made by
!!   `nk_ilu0_sched_values`), since processing rows level by level directly in the
!!   CSR layout of the matrix loses the locality of the serial sweeps (about 2x
!!   slower on one thread for FE matrices). The copy holds, for the rows of the
!!   forward schedule in execution order, the entries left of the diagonal (entries 1
!!   to `nnzl`), followed, for the rows of the backward schedule in execution order,
!!   by the inverted pivot and the entries right of the diagonal.
!! - Intermediate results are stored by position in the schedule rather than by row,
!!   since rows of one level are spread through the vector and writing them in place
!!   makes threads share cache lines (no speedup at all on 2 to 4 threads). `col` maps
!!   each reordered entry to the position of its column in the schedule of its sweep;
!!   the entry of each pivot instead holds the position of the row in the forward
!!   schedule, where the backward sweep gets its starting value.
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: native_ilu_sched
  INTEGER(i4) :: nr = 0 !< Number of rows
  INTEGER(i4) :: nnz = 0 !< Number of non-zeros
  INTEGER(i4) :: nnzl = 0 !< Number of non-zeros left of the diagonal
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: col !< Schedule positions for the reordered factor values [nnz]
  LOGICAL :: use_factor = .FALSE. !< Use the scheduled factorization
  LOGICAL :: use_solve = .FALSE. !< Use the scheduled solve
  TYPE(native_level_sched) :: lower !< Forward (L) sweep and factorization
  TYPE(native_level_sched) :: upper !< Backward (U) sweep
END TYPE native_ilu_sched
PUBLIC native_row_split, native_thread_rows, nk_jacobi_start, nk_jacobi_sweep, nk_jacobi_restart
PUBLIC native_ilu_sched_build, native_ilu_sched_active, native_ilu_sched_free, native_ilu_sched_stats
PUBLIC native_level_segment
CONTAINS
!---------------------------------------------------------------------------------
!> Get the first row of one part of a non-zero balanced, contiguous row partition
!!
!! Rows are split into `nparts` contiguous ranges that each contain approximately
!! the same number of non-zeros. Part `ipart` covers rows
!! `native_row_split(kr,nr,ipart,nparts)` to `native_row_split(kr,nr,ipart+1,nparts)-1`,
!! so the parts are disjoint and together cover all rows.
!---------------------------------------------------------------------------------
PURE FUNCTION native_row_split(kr,nr,ipart,nparts) RESULT(irow)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), INTENT(in) :: ipart !< Part index (0-based, 0 <= ipart <= nparts)
INTEGER(i4), INTENT(in) :: nparts !< Number of parts
INTEGER(i4) :: irow !< First row in part (`nr+1` when `ipart>=nparts`)
INTEGER(i4) :: ilow,ihigh,imid
INTEGER(i8) :: nnz_target
IF(ipart<=0)THEN
  irow=1
  RETURN
ELSE IF(ipart>=nparts)THEN
  irow=nr+1
  RETURN
END IF
nnz_target=INT(kr(1),i8)+(INT(kr(nr+1)-kr(1),i8)*INT(ipart,i8))/INT(nparts,i8)
!---First row whose entries start at or after the target
ilow=1
ihigh=nr+1
DO WHILE(ilow<ihigh)
  imid=(ilow+ihigh)/2
  IF(INT(kr(imid),i8)<nnz_target)THEN
    ilow=imid+1
  ELSE
    ihigh=imid
  END IF
END DO
irow=ilow
END FUNCTION native_row_split
!---------------------------------------------------------------------------------
!> Get the range of rows to be processed by the calling OpenMP thread
!!
!! Uses @ref native_row_split with one part per thread in the current team, so each
!! thread processes a contiguous block of rows with a similar number of non-zeros.
!! Outside of an active parallel region all rows are returned. The same splitting
!! applies to columns when a CSC column pointer is passed instead.
!---------------------------------------------------------------------------------
SUBROUTINE native_thread_rows(kr,nr,ifirst,ilast)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), INTENT(out) :: ifirst !< First row for this thread
INTEGER(i4), INTENT(out) :: ilast !< Last row for this thread (`ilast<ifirst` if none)
INTEGER(i4) :: ithread,nthreads
ithread=omp_get_thread_num()
nthreads=omp_get_num_threads()
ifirst=native_row_split(kr,nr,ithread,nthreads)
ilast=native_row_split(kr,nr,ithread+1,nthreads)-1
END SUBROUTINE native_thread_rows
!---------------------------------------------------------------------------------
!> Start a damped Jacobi iteration from a zero guess
!!
!! \f$ u = 0 \f$, \f$ p = \omega D^{-1} g \f$
!---------------------------------------------------------------------------------
SUBROUTINE nk_jacobi_start(n,df,d,g,u,p)
INTEGER(i4), INTENT(in) :: n !< Number of entries
REAL(r8), INTENT(in) :: df !< Damping factor \f$ \omega \f$
REAL(r8), INTENT(in) :: d(n) !< Diagonal \f$ D \f$ [n]
REAL(r8), INTENT(in) :: g(n) !< Residual [n]
REAL(r8), INTENT(out) :: u(n) !< Solution [n]
REAL(r8), INTENT(out) :: p(n) !< Correction [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH/2)
DO i=1,n
  u(i)=0.d0
  p(i)=df*g(i)/d(i)
END DO
END SUBROUTINE nk_jacobi_start
!---------------------------------------------------------------------------------
!> One fused damped Jacobi sweep, given \f$ q = A p \f$ for the current correction
!!
!! \f$ u = u + p \f$ (only if `u` is present), \f$ g = g - q \f$,
!! \f$ p = \omega D^{-1} g \f$, and optionally saves the updated `u` and `g`.
!---------------------------------------------------------------------------------
SUBROUTINE nk_jacobi_sweep(n,df,d,q,g,p,u,usave,gsave)
INTEGER(i4), INTENT(in) :: n !< Number of entries
REAL(r8), INTENT(in) :: df !< Damping factor \f$ \omega \f$
REAL(r8), INTENT(in) :: d(n) !< Diagonal \f$ D \f$ [n]
REAL(r8), INTENT(in) :: q(n) !< Product of matrix with correction [n]
REAL(r8), INTENT(inout) :: g(n) !< Residual [n]
REAL(r8), INTENT(inout) :: p(n) !< Correction [n]
REAL(r8), OPTIONAL, INTENT(inout) :: u(n) !< Solution [n] (optional)
REAL(r8), OPTIONAL, INTENT(out) :: usave(n) !< Copy of updated solution [n] (optional, requires `u`)
REAL(r8), OPTIONAL, INTENT(out) :: gsave(n) !< Copy of updated residual [n] (optional, with `usave`)
INTEGER(i4) :: i
IF(PRESENT(u))THEN
  IF(PRESENT(usave))THEN
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH/2)
    DO i=1,n
      u(i)=u(i)+p(i)
      g(i)=g(i)-q(i)
      p(i)=df*g(i)/d(i)
      usave(i)=u(i)
      gsave(i)=g(i)
    END DO
  ELSE
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH/2)
    DO i=1,n
      u(i)=u(i)+p(i)
      g(i)=g(i)-q(i)
      p(i)=df*g(i)/d(i)
    END DO
  END IF
ELSE
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH/2)
  DO i=1,n
    g(i)=g(i)-q(i)
    p(i)=df*g(i)/d(i)
  END DO
END IF
END SUBROUTINE nk_jacobi_sweep
!---------------------------------------------------------------------------------
!> Restart a Jacobi iteration from a saved state
!!
!! \f$ p = u - u_{save} \f$, \f$ g = g_{save} \f$
!---------------------------------------------------------------------------------
SUBROUTINE nk_jacobi_restart(n,u,usave,gsave,g,p)
INTEGER(i4), INTENT(in) :: n !< Number of entries
REAL(r8), INTENT(in) :: u(n) !< Solution [n]
REAL(r8), INTENT(in) :: usave(n) !< Saved solution [n]
REAL(r8), INTENT(in) :: gsave(n) !< Saved residual [n]
REAL(r8), INTENT(out) :: g(n) !< Residual [n]
REAL(r8), INTENT(out) :: p(n) !< Correction [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  p(i)=u(i)-usave(i)
  g(i)=gsave(i)
END DO
END SUBROUTINE nk_jacobi_restart
!---------------------------------------------------------------------------------
!> Build the level schedules of the ILU(0) factorization and solve for a CSR pattern
!!
!! The dependencies of row `i` are the columns `lc(j)<i` of its entries for the
!! forward sweep and the factorization, and the columns `lc(j)>i` for the backward
!! sweep, so only `kr` and `lc` are needed. Levels are found in one pass over the
!! entries (in increasing row order for the forward sweep, decreasing for the
!! backward one), and rows are sorted into levels with a stable counting sort.
!!
!! Execution segments (see @ref native_level_sched): a level with at least
!! `min_rows` rows is a parallel segment; runs of consecutive smaller levels form
!! serial segments. Each segment costs one barrier, and only rows in parallel levels
!! are shared between threads, so the schedules are only used when most rows are in
!! parallel levels (at least half, in each schedule used) and segments are large on
!! average:
!! - Solve: `use_solve` when `nr>OFT_OMP_VTHRESH` and the two sweeps together
!!   average at least @ref native_ilu_solve_nnz non-zeros per segment.
!! - Factorization: `use_factor` when `nr>OFT_OMP_VTHRESH/4` and the forward
!!   schedule averages at least @ref native_ilu_factor_rows rows per segment (a row of
!!   the factorization costs far more than a row of a sweep).
!!
!! These thresholds come from timings of Lagrange Laplacians (orders 1-3, 559 to
!! 68705 rows, natural ordering) on 2 and 4 threads: the threaded solve was 1.3-1.4x
!! (2 threads) and 1.7-2.0x (4 threads) faster at 5700 or more non-zeros per segment
!! and no faster at 2100 or fewer, and the threaded factorization 1.4-1.7x (2 threads)
!! and 2.4-2.5x (4 threads) faster at 9009 rows or more (over 90% of rows in parallel
!! levels), with no gain for matrices with 5-40% of rows in parallel levels.
!!
!! The flags may be changed after the build (eg. to force the scheduled kernels in
!! tests); results do not depend on them.
!!
!! @note The layout of the reordered factor values assumes the pattern requirements
!! of the ILU(0) kernels (sorted columns and a diagonal entry in every row). Without
!! exactly one diagonal entry in every row it is not built and the schedules are
!! never used (@ref native_ilu_sched_active is false).
!---------------------------------------------------------------------------------
SUBROUTINE native_ilu_sched_build(nr,kr,lc,sched,min_rows)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), CONTIGUOUS, INTENT(in) :: lc(:) !< CSR column indices [nnz]
TYPE(native_ilu_sched), INTENT(inout) :: sched !< Schedules (previous contents are released)
INTEGER(i4), OPTIONAL, INTENT(in) :: min_rows !< Minimum rows in a parallel level (default @ref native_ilu_level_rows)
INTEGER(i4) :: i,j,k,lev,minr,jp
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: level
CALL native_ilu_sched_free(sched)
minr=native_ilu_level_rows
IF(PRESENT(min_rows))minr=min_rows
sched%nr=nr
sched%nnz=kr(nr+1)-kr(1)
ALLOCATE(level(nr))
!---Forward sweep: rows referenced left of the diagonal precede the row
DO i=1,nr
  lev=0
  DO j=kr(i),kr(i+1)-1
    IF(lc(j)<i)lev=MAX(lev,level(lc(j)))
  END DO
  level(i)=lev+1
END DO
CALL level_sort(level,minr,sched%lower)
!---Backward sweep: rows referenced right of the diagonal precede the row
DO i=nr,1,-1
  lev=0
  DO j=kr(i),kr(i+1)-1
    IF(lc(j)>i)lev=MAX(lev,level(lc(j)))
  END DO
  level(i)=lev+1
END DO
CALL level_sort(level,minr,sched%upper)
DEALLOCATE(level)
!---Reordered layout needs exactly one diagonal entry per row (otherwise not threaded)
jp=0
DO i=1,nr
  DO j=kr(i),kr(i+1)-1
    IF(lc(j)==i)jp=jp+1
  END DO
END DO
IF(jp/=nr)RETURN
!---Layout of the reordered factor values (L entries, then pivot and U entries)
ALLOCATE(sched%col(sched%nnz))
ALLOCATE(sched%lower%kp(nr+1),sched%upper%kp(nr+1))
jp=1
DO k=1,nr
  i=sched%lower%perm(k)
  sched%lower%kp(k)=jp
  DO j=kr(i),kr(i+1)-1
    IF(lc(j)<i)THEN
      sched%col(jp)=sched%lower%ipos(lc(j))
      jp=jp+1
    END IF
  END DO
END DO
sched%lower%kp(nr+1)=jp
sched%nnzl=jp-1
DO k=1,nr
  i=sched%upper%perm(k)
  sched%upper%kp(k)=jp
  sched%col(jp)=sched%lower%ipos(i)
  jp=jp+1
  DO j=kr(i),kr(i+1)-1
    IF(lc(j)>i)THEN
      sched%col(jp)=sched%upper%ipos(lc(j))
      jp=jp+1
    END IF
  END DO
END DO
sched%upper%kp(nr+1)=jp
!---Profitability (see description)
sched%use_solve=(nr>OFT_OMP_VTHRESH) &
  .AND.(2*sched%lower%npar>=nr).AND.(2*sched%upper%npar>=nr) &
  .AND.(INT(sched%nnz,i8)>=INT(native_ilu_solve_nnz,i8)*INT(sched%lower%nseg+sched%upper%nseg,i8))
sched%use_factor=(nr>OFT_OMP_VTHRESH/4).AND.(2*sched%lower%npar>=nr) &
  .AND.(nr>=native_ilu_factor_rows*sched%lower%nseg)
CONTAINS
! Group rows by level (stable) and split levels into execution segments
SUBROUTINE level_sort(level,minr,ls)
INTEGER(i4), INTENT(in) :: level(:)
INTEGER(i4), INTENT(in) :: minr
TYPE(native_level_sched), INTENT(inout) :: ls
INTEGER(i4) :: i,k,lev,nlev,nrows
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: knext
LOGICAL :: in_serial
nlev=0
IF(nr>0)nlev=MAXVAL(level)
ls%nlevels=nlev
ALLOCATE(ls%lptr(nlev+1),ls%perm(nr))
ls%lptr=0
DO i=1,nr
  ls%lptr(level(i)+1)=ls%lptr(level(i)+1)+1
END DO
ls%lptr(1)=1
DO k=1,nlev
  ls%lptr(k+1)=ls%lptr(k+1)+ls%lptr(k)
END DO
ALLOCATE(knext(nlev),ls%ipos(nr))
knext=ls%lptr(1:nlev)
DO i=1,nr
  ls%perm(knext(level(i)))=i
  ls%ipos(i)=knext(level(i))
  knext(level(i))=knext(level(i))+1
END DO
DEALLOCATE(knext)
!---Segments: count, then fill
ls%maxrows=0
DO k=1,2
  ls%nseg=0
  ls%npar=0
  in_serial=.FALSE.
  DO lev=1,nlev
    nrows=ls%lptr(lev+1)-ls%lptr(lev)
    IF(k==1)ls%maxrows=MAX(ls%maxrows,nrows)
    IF(nrows>=minr)THEN
      ls%nseg=ls%nseg+1
      ls%npar=ls%npar+nrows
      in_serial=.FALSE.
      IF(k==2)THEN
        ls%sptr(ls%nseg)=ls%lptr(lev)
        ls%spar(ls%nseg)=.TRUE.
      END IF
    ELSE IF(.NOT.in_serial)THEN
      ls%nseg=ls%nseg+1
      in_serial=.TRUE.
      IF(k==2)THEN
        ls%sptr(ls%nseg)=ls%lptr(lev)
        ls%spar(ls%nseg)=.FALSE.
      END IF
    END IF
  END DO
  IF(k==1)ALLOCATE(ls%sptr(ls%nseg+1),ls%spar(ls%nseg))
END DO
ls%sptr(ls%nseg+1)=nr+1
END SUBROUTINE level_sort
END SUBROUTINE native_ilu_sched_build
!---------------------------------------------------------------------------------
!> Release the storage of ILU(0) level schedules
!---------------------------------------------------------------------------------
SUBROUTINE native_ilu_sched_free(sched)
TYPE(native_ilu_sched), INTENT(inout) :: sched !< Schedules
CALL free_level(sched%lower)
CALL free_level(sched%upper)
IF(ALLOCATED(sched%col))DEALLOCATE(sched%col)
sched%nr=0
sched%nnz=0
sched%nnzl=0
sched%use_factor=.FALSE.
sched%use_solve=.FALSE.
CONTAINS
SUBROUTINE free_level(ls)
TYPE(native_level_sched), INTENT(inout) :: ls
IF(ALLOCATED(ls%perm))DEALLOCATE(ls%perm)
IF(ALLOCATED(ls%lptr))DEALLOCATE(ls%lptr)
IF(ALLOCATED(ls%sptr))DEALLOCATE(ls%sptr)
IF(ALLOCATED(ls%spar))DEALLOCATE(ls%spar)
IF(ALLOCATED(ls%kp))DEALLOCATE(ls%kp)
IF(ALLOCATED(ls%ipos))DEALLOCATE(ls%ipos)
ls%nlevels=0
ls%maxrows=0
ls%nseg=0
ls%npar=0
END SUBROUTINE free_level
END SUBROUTINE native_ilu_sched_free
!---------------------------------------------------------------------------------
!> Check whether a level-scheduled ILU(0) kernel would run threaded
!!
!! True when the schedule is built for `nr` rows, the corresponding profitability
!! flag is set, more than one thread is available, and the caller is not already
!! inside an active parallel region (eg. block-Jacobi with several local blocks per
!! process, which factors and applies its blocks in a parallel loop).
!!
!! @result Use the scheduled path
!---------------------------------------------------------------------------------
FUNCTION native_ilu_sched_active(sched,nr,factor) RESULT(active)
TYPE(native_ilu_sched), INTENT(in) :: sched !< Schedules
INTEGER(i4), INTENT(in) :: nr !< Number of rows of the factored matrix
LOGICAL, INTENT(in) :: factor !< Check the factorization (otherwise the solve)
LOGICAL :: active
active=.FALSE.
IF(sched%nr/=nr.OR.(.NOT.ALLOCATED(sched%col)))RETURN
IF(factor)THEN
  IF(.NOT.sched%use_factor)RETURN
ELSE
  IF(.NOT.sched%use_solve)RETURN
END IF
IF(omp_in_parallel())RETURN
active=(omp_get_max_threads()>1)
END FUNCTION native_ilu_sched_active
!---------------------------------------------------------------------------------
!> Write a summary of ILU(0) level schedules (levels, rows per level and segments)
!---------------------------------------------------------------------------------
SUBROUTINE native_ilu_sched_stats(sched,prefix)
TYPE(native_ilu_sched), INTENT(in) :: sched !< Schedules
CHARACTER(LEN=*), INTENT(in) :: prefix !< Leading text of each line
CALL level_stats('L',sched%lower)
CALL level_stats('U',sched%upper)
WRITE(*,'(2A,L1,A,L1)')prefix,'threaded factor = ',sched%use_factor,', threaded solve = ',sched%use_solve
CONTAINS
SUBROUTINE level_stats(name,ls)
CHARACTER(LEN=*), INTENT(in) :: name
TYPE(native_level_sched), INTENT(in) :: ls
REAL(r8) :: avg_lev,avg_seg
avg_lev=REAL(sched%nr,r8)/REAL(MAX(ls%nlevels,1),r8)
avg_seg=REAL(sched%nr,r8)/REAL(MAX(ls%nseg,1),r8)
WRITE(*,'(4A,I0,A,F0.1,A,I0,A,I0,A,F0.1,A,F0.1,A)')prefix,'sweep ',name,': ',ls%nlevels, &
  ' levels (',avg_lev,' rows/level, max ',ls%maxrows,'), ',ls%nseg,' segments (', &
  avg_seg,' rows/segment, ',100.d0*REAL(ls%npar,r8)/REAL(MAX(sched%nr,1),r8),'% in parallel levels)'
END SUBROUTINE level_stats
END SUBROUTINE native_ilu_sched_stats
!---------------------------------------------------------------------------------
!> Positions in `perm` processed by the calling thread for one execution segment
!!
!! Rows of a parallel segment are split into contiguous blocks, one per thread, with
!! nearly equal numbers of reordered factor entries plus rows (from `kp`), since the
!! cost of a row varies widely (eg. vertex and edge rows of quadratic elements). A
!! serial segment is processed entirely by thread 0 (`k2<k1` for others). Must be
!! called from inside the parallel region executing the schedule.
!---------------------------------------------------------------------------------
SUBROUTINE native_level_segment(ls,s,k1,k2)
TYPE(native_level_sched), INTENT(in) :: ls !< Level schedule
INTEGER(i4), INTENT(in) :: s !< Segment
INTEGER(i4), INTENT(out) :: k1 !< First position in `perm`
INTEGER(i4), INTENT(out) :: k2 !< Last position in `perm`
INTEGER(i4) :: ithread,nthreads
ithread=omp_get_thread_num()
nthreads=omp_get_num_threads()
IF(ls%spar(s))THEN
  k1=split_point(ithread)
  k2=split_point(ithread+1)-1
ELSE IF(ithread==0)THEN
  k1=ls%sptr(s)
  k2=ls%sptr(s+1)-1
ELSE
  k1=1
  k2=0
END IF
CONTAINS
! First position of part ipart: the first position whose cumulative cost (entries plus
! rows before it in the segment) reaches the fraction ipart/nthreads of the total
FUNCTION split_point(ipart) RESULT(kpos)
INTEGER(i4), INTENT(in) :: ipart
INTEGER(i4) :: kpos,ks,ke,ilow,ihigh,imid
INTEGER(i8) :: target
ks=ls%sptr(s)
ke=ls%sptr(s+1)
IF(ipart<=0)THEN
  kpos=ks
  RETURN
ELSE IF(ipart>=nthreads)THEN
  kpos=ke
  RETURN
END IF
target=(INT(ls%kp(ke)-ls%kp(ks)+ke-ks,i8)*INT(ipart,i8))/INT(nthreads,i8)
ilow=ks
ihigh=ke
DO WHILE(ilow<ihigh)
  imid=(ilow+ihigh)/2
  IF(INT(ls%kp(imid)-ls%kp(ks)+imid-ks,i8)<target)THEN
    ilow=imid+1
  ELSE
    ihigh=imid
  END IF
END DO
kpos=ilow
END FUNCTION split_point
END SUBROUTINE native_level_segment
END MODULE oft_native_kernels_base
!---------------------------------------------------------------------------------
! Single-type kernels
!---------------------------------------------------------------------------------
#define KT REAL(r8)
#define KZERO 0.d0
#define KONE 1.d0
#define KPTR native_rvals_ptr
MODULE oft_native_kernels_vec_r8
#include "native_kernels_vec.inc"
END MODULE oft_native_kernels_vec_r8
#undef KT
#undef KZERO
#undef KONE
#undef KPTR
#define KT COMPLEX(c8)
#define KZERO (0.d0,0.d0)
#define KONE (1.d0,0.d0)
#define KPTR native_cvals_ptr
MODULE oft_native_kernels_vec_c8
#include "native_kernels_vec.inc"
END MODULE oft_native_kernels_vec_c8
#undef KT
#undef KZERO
#undef KONE
#undef KPTR
!---------------------------------------------------------------------------------
! Incomplete LU kernels
!---------------------------------------------------------------------------------
#define KT REAL(r8)
#define KZERO 0.d0
#define KONE 1.d0
MODULE oft_native_kernels_ilu_r8
#include "native_kernels_ilu.inc"
END MODULE oft_native_kernels_ilu_r8
#undef KT
#undef KZERO
#undef KONE
#define KT COMPLEX(c8)
#define KZERO (0.d0,0.d0)
#define KONE (1.d0,0.d0)
MODULE oft_native_kernels_ilu_c8
#include "native_kernels_ilu.inc"
END MODULE oft_native_kernels_ilu_c8
#undef KT
#undef KZERO
#undef KONE
!---------------------------------------------------------------------------------
! Update kernels (destination/source)
!---------------------------------------------------------------------------------
#define KY REAL(r8)
#define KX REAL(r8)
#define KYZERO 0.d0
#define KTOY(a) a
MODULE oft_native_kernels_upd_rr
#include "native_kernels_upd.inc"
END MODULE oft_native_kernels_upd_rr
#undef KY
#undef KX
#undef KYZERO
#undef KTOY
#define KY COMPLEX(c8)
#define KX REAL(r8)
#define KYZERO (0.d0,0.d0)
#define KTOY(a) CMPLX(a,0.d0,c8)
MODULE oft_native_kernels_upd_cr
#include "native_kernels_upd.inc"
END MODULE oft_native_kernels_upd_cr
#undef KY
#undef KX
#undef KYZERO
#undef KTOY
#define KY COMPLEX(c8)
#define KX COMPLEX(c8)
#define KYZERO (0.d0,0.d0)
#define KTOY(a) a
MODULE oft_native_kernels_upd_cc
#include "native_kernels_upd.inc"
END MODULE oft_native_kernels_upd_cc
#undef KY
#undef KX
#undef KYZERO
#undef KTOY
!---------------------------------------------------------------------------------
! Product kernels (left or matrix/right/result)
!---------------------------------------------------------------------------------
#define KA REAL(r8)
#define KB REAL(r8)
#define KR REAL(r8)
#define KRZERO 0.d0
#define KBPTR native_rvals_ptr
#define KCONJ(a) a
#define KREAL
MODULE oft_native_kernels_prod_rr
#include "native_kernels_prod.inc"
END MODULE oft_native_kernels_prod_rr
#undef KA
#undef KB
#undef KR
#undef KRZERO
#undef KBPTR
#undef KCONJ
#undef KREAL
#define KA REAL(r8)
#define KB COMPLEX(c8)
#define KR COMPLEX(c8)
#define KRZERO (0.d0,0.d0)
#define KBPTR native_cvals_ptr
#define KCONJ(a) a
MODULE oft_native_kernels_prod_rc
#include "native_kernels_prod.inc"
END MODULE oft_native_kernels_prod_rc
#undef KA
#undef KB
#undef KR
#undef KRZERO
#undef KBPTR
#undef KCONJ
#define KA COMPLEX(c8)
#define KB REAL(r8)
#define KR COMPLEX(c8)
#define KRZERO (0.d0,0.d0)
#define KBPTR native_rvals_ptr
#define KCONJ(a) CONJG(a)
MODULE oft_native_kernels_prod_cr
#include "native_kernels_prod.inc"
END MODULE oft_native_kernels_prod_cr
#undef KA
#undef KB
#undef KR
#undef KRZERO
#undef KBPTR
#undef KCONJ
#define KA COMPLEX(c8)
#define KB COMPLEX(c8)
#define KR COMPLEX(c8)
#define KRZERO (0.d0,0.d0)
#define KBPTR native_cvals_ptr
#define KCONJ(a) CONJG(a)
MODULE oft_native_kernels_prod_cc
#include "native_kernels_prod.inc"
END MODULE oft_native_kernels_prod_cc
#undef KA
#undef KB
#undef KR
#undef KRZERO
#undef KBPTR
#undef KCONJ
!---------------------------------------------------------------------------------
!> Native backend kernels for all supported value types
!!
!! Collects the generic kernel interfaces of every instantiation (see the file
!! description of native_kernels.F90) together with @ref oft_native_kernels_base.
!! Generic kernels:
!! - Single type: `nk_set`, `nk_scale`, `nk_maxpy`, `nk_norm_inf`, `nk_owned_sum`,
!!   `nk_owned_norm1`, `nk_gather`, `nk_scatter`, `nk_csr_zero_rows`,
!!   `nk_csr_zero_cols`, `nk_csr_identity_rows`, `nk_csr_get_diag`, `nk_set_at`,
!!   `nk_add_at`, `nk_atomic_add_at` (real and complex)
!! - Destination/source: `nk_copy`, `nk_axpby`, `nk_axpbypcz`, `nk_mult`, `nk_div`
!!   (real/real, complex/real, complex/complex)
!! - Products: `nk_owned_dot`, `nk_owned_mdot`, `nk_csr_matvec` (4- and 2-byte
!!   column indices), `nk_csc_matvec`, `nk_dense_matvec`, `nk_dense_matvect`
!!   (all real/complex combinations; the left operand of dot products is conjugated)
!! - Incomplete LU: `nk_ilu0_factor`, `nk_ilu0_solve` (real and complex; serial, or
!!   level-scheduled with OpenMP when given schedules from `native_ilu_sched_build`)
!! - Real only: `nk_jacobi_start`, `nk_jacobi_sweep`, `nk_jacobi_restart`
!---------------------------------------------------------------------------------
MODULE oft_native_kernels
USE oft_native_kernels_base
USE oft_native_kernels_vec_r8
USE oft_native_kernels_vec_c8
USE oft_native_kernels_ilu_r8
USE oft_native_kernels_ilu_c8
USE oft_native_kernels_upd_rr
USE oft_native_kernels_upd_cr
USE oft_native_kernels_upd_cc
USE oft_native_kernels_prod_rr
USE oft_native_kernels_prod_rc
USE oft_native_kernels_prod_cr
USE oft_native_kernels_prod_cc
IMPLICIT NONE
PUBLIC
END MODULE oft_native_kernels
