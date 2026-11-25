!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!> @file oft_patch_solvers.F90
!
!> Overlapping patch (Vanka/star-type) smoothers and preconditioners
!!
!! Patches are lists of local matrix rows attached to the matrix by the code
!! that builds it (see @ref oft_la_base::oft_patch_info "oft_patch_info"), for
!! example vertex stars built from the FE topology by
!! @ref fem_utils::comp_fem_vertex_patches "comp_fem_vertex_patches". Each patch
!! is solved exactly with a dense LU factorization of the corresponding block
!! of the assembled matrix.
!!
!! Supported modes:
!! - Additive: average (partition of unity) of the patch corrections
!! - Multiplicative: Gauss-Seidel-like sweep over patches, ordered by color,
!!   with patches of one color processed in parallel (OpenMP)
!! - Symmetric: forward followed by backward multiplicative sweep
!!
!! Distributed meshes: patches whose defining entity is shared between
!! processes ("seam" patches) are completed across processes. Every process
!! holding part of a seam patch sends its rows of the patch block, built from its
!! unassembled local matrix, to the owning process, which sums them into the
!! complete patch block and factors it. During application the owner gathers
!! the residual, solves and returns the correction to all processes holding part
!! of the patch, so the overlap of the patches is the same as in serial. In
!! the multiplicative modes the interior patches of all processes are swept
!! first (they do not interact across processes), followed by an additive
!! correction from the seam patches using the updated global residual.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup doxy_oft_lin_alg
!------------------------------------------------------------------------------
MODULE oft_patch_solvers
USE oft_local
USE oft_base
USE oft_sort, ONLY: sort_array, search_array
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_patch_info
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_setup_full
USE oft_solver_base, ONLY: oft_solver, solver_setup, solver_delete
IMPLICIT NONE
#include "local.h"
PRIVATE
!------------------------------------------------------------------------------
!> Dense LU factorization of a single patch
!------------------------------------------------------------------------------
TYPE :: oft_patch_factor
  INTEGER(i4) :: ptype = 0 !< Patch type (0: interior, 1: owned seam, 2: contributed seam)
  INTEGER(i4) :: n = 0 !< Number of rows in patch
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: ind => NULL() !< Local row for each patch row (0 if not held locally)
  INTEGER(i8), POINTER, CONTIGUOUS, DIMENSION(:) :: gid => NULL() !< Global row index for each patch row (seam patches)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: piv => NULL() !< LU pivots
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: lu => NULL() !< LU factors (double precision)
  REAL(r4), POINTER, CONTIGUOUS, DIMENSION(:,:) :: lu4 => NULL() !< LU factors (single precision)
END TYPE oft_patch_factor
!------------------------------------------------------------------------------
!> Seam patch communication with one neighboring process
!------------------------------------------------------------------------------
TYPE :: oft_patch_nbr
  INTEGER(i4) :: rank = -1 !< Neighbor process
  INTEGER(i4) :: nsend = 0 !< Number of local patches contributed to neighbor (owned by neighbor)
  INTEGER(i4) :: nrecv = 0 !< Number of owned patches the neighbor contributes to
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lsend => NULL() !< Contributed patches (ordered by global key)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lrecv => NULL() !< Owned patches (ordered by global key)
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: pos => NULL() !< Positions of the neighbor's rows in each owned patch
  INTEGER(i4) :: nsbuf = 0 !< Size of send buffer for vectors
  INTEGER(i4) :: nrbuf = 0 !< Size of receive buffer for vectors
END TYPE oft_patch_nbr
!------------------------------------------------------------------------------
!> Overlapping patch smoother/preconditioner
!!
!! When `its<=0` the object acts as a preconditioner (\f$ u = df B g \f$).
!! When `its>0` it performs `its` damped Richardson iterations
!! \f$ u \leftarrow u + df B (g - A u) \f$ starting from the input `u` and
!! returns the final residual in `g`, so it can be used directly as a
!! multigrid smoother.
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_patch_precond
  INTEGER(i4) :: mode = 1 !< Patch combination (1: additive, 2: multiplicative, 3: symmetric multiplicative)
  REAL(r8) :: df = 1.d0 !< Damping factor
  LOGICAL :: single = .FALSE. !< Store factors in single precision?
  LOGICAL :: refactor = .TRUE. !< Refactor patches on next application?
  LOGICAL :: new_structure = .TRUE. !< Rebuild index structures on next application?
  LOGICAL :: seams = .FALSE. !< Are there seam patches on any process?
  INTEGER(i4) :: npatch = 0 !< Number of patches
  INTEGER(i4) :: ncolors = 0 !< Number of patch colors
  INTEGER(i4) :: maxn = 0 !< Maximum patch size
  INTEGER(i4) :: nnbr = 0 !< Number of neighbor processes for seam patches
  INTEGER(i8) :: nfactor = 0 !< Total number of stored factor entries
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kcolor => NULL() !< Pointer into color patch list
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lcolor => NULL() !< Interior patch list ordered by color
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lseam => NULL() !< Seam patches owned locally
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: rep => NULL() !< Representative row for each local row
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kcopy => NULL() !< Pointer into copy list for each representative row
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lcopy => NULL() !< Local copies of each representative row
  INTEGER(i8), POINTER, CONTIGUOUS, DIMENSION(:) :: glob => NULL() !< Global index of each local row
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: isrep => NULL() !< Is row a representative (non-redundant)?
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: iso => NULL() !< Is row decoupled (diagonal only)?
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: diag => NULL() !< Diagonal values for decoupled rows
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: wts => NULL() !< Combination weights (additive: all patches, otherwise seam patches)
  TYPE(oft_patch_factor), POINTER, DIMENSION(:) :: pf => NULL() !< Patch factorizations
  TYPE(oft_patch_nbr), POINTER, DIMENSION(:) :: nbr => NULL() !< Seam patch neighbors
  CLASS(oft_vector), POINTER :: rv => NULL() !< Residual work vector
  CLASS(oft_vector), POINTER :: duv => NULL() !< Correction work vector
CONTAINS
  !> Apply smoother/preconditioner
  PROCEDURE :: apply => patch_precond_apply
  !> Mark for refactorization after matrix update
  PROCEDURE :: update => patch_precond_update
  !> Setup from XML node
  PROCEDURE :: setup_from_xml => patch_precond_setup_xml
  !> Print solver information
  PROCEDURE :: view => patch_precond_view
  !> Check thread safety
  PROCEDURE :: check_thread => patch_precond_check_thread
  !> Clean-up internal storage
  PROCEDURE :: delete => patch_precond_delete
END TYPE oft_patch_precond
CONTAINS
!------------------------------------------------------------------------------
!> Exchange real buffers with neighbor processes (sizes known on both sides)
!------------------------------------------------------------------------------
SUBROUTINE nbr_exchange(self,sbuf,rbuf)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
TYPE(oft_1d_real), INTENT(inout) :: sbuf(:) !< Send buffers [nnbr]
TYPE(oft_1d_real), INTENT(inout) :: rbuf(:) !< Receive buffers [nnbr]
#ifdef HAVE_MPI
INTEGER(i4) :: i,ierr
#ifdef OFT_MPI_F08
TYPE(mpi_request) :: reqs(2*self%nnbr)
#else
INTEGER(i4) :: reqs(2*self%nnbr)
#endif
IF(self%nnbr==0)RETURN
DO i=1,self%nnbr
  CALL MPI_IRECV(rbuf(i)%v,rbuf(i)%n,OFT_MPI_R8,self%nbr(i)%rank,31,oft_env%comm,reqs(i),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','nbr_exchange',__FILE__)
END DO
DO i=1,self%nnbr
  CALL MPI_ISEND(sbuf(i)%v,sbuf(i)%n,OFT_MPI_R8,self%nbr(i)%rank,31,oft_env%comm,reqs(self%nnbr+i),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','nbr_exchange',__FILE__)
END DO
CALL oft_mpi_waitall(2*self%nnbr,reqs,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','nbr_exchange',__FILE__)
#else
IF(self%nnbr>0)CALL oft_abort('Seam patches require MPI','nbr_exchange',__FILE__)
#endif
END SUBROUTINE nbr_exchange
!------------------------------------------------------------------------------
!> Exchange real buffers with neighbor processes, communicating sizes first
!------------------------------------------------------------------------------
SUBROUTINE nbr_exchange_sized(self,sbuf,rbuf)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
TYPE(oft_1d_real), INTENT(inout) :: sbuf(:) !< Send buffers [nnbr]
TYPE(oft_1d_real), INTENT(inout) :: rbuf(:) !< Receive buffers [nnbr] (allocated here)
TYPE(oft_1d_real), ALLOCATABLE, DIMENSION(:) :: ssz,rsz
INTEGER(i4) :: i
ALLOCATE(ssz(self%nnbr),rsz(self%nnbr))
DO i=1,self%nnbr
  ssz(i)%n=1
  ALLOCATE(ssz(i)%v(1))
  ssz(i)%v(1)=REAL(sbuf(i)%n,8)
  rsz(i)%n=1
  ALLOCATE(rsz(i)%v(1))
END DO
CALL nbr_exchange(self,ssz,rsz)
DO i=1,self%nnbr
  rbuf(i)%n=NINT(rsz(i)%v(1))
  ALLOCATE(rbuf(i)%v(rbuf(i)%n))
  DEALLOCATE(ssz(i)%v,rsz(i)%v)
END DO
DEALLOCATE(ssz,rsz)
CALL nbr_exchange(self,sbuf,rbuf)
END SUBROUTINE nbr_exchange_sized
!------------------------------------------------------------------------------
!> Free a set of communication buffers
!------------------------------------------------------------------------------
SUBROUTINE free_bufs(bufs)
TYPE(oft_1d_real), ALLOCATABLE, INTENT(inout) :: bufs(:) !< Buffers
INTEGER(i4) :: i
IF(.NOT.ALLOCATED(bufs))RETURN
DO i=1,SIZE(bufs)
  IF(ASSOCIATED(bufs(i)%v))DEALLOCATE(bufs(i)%v)
END DO
DEALLOCATE(bufs)
END SUBROUTINE free_bufs
!------------------------------------------------------------------------------
!> Build representative row map, decoupled row flags, patch lists, seam
!! communication structures and combination weights
!------------------------------------------------------------------------------
SUBROUTINE patch_precond_structure(self,A)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
CLASS(oft_native_matrix), INTENT(inout) :: A !< Matrix (with patches)
TYPE(oft_patch_info), POINTER :: pinfo
TYPE(oft_1d_real), ALLOCATABLE, DIMENSION(:) :: sbuf,rbuf
INTEGER(i4) :: i,j,k,ip,nr,n,ic,maxn,jn,off,nloc
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: isort,cnt,nbr_map,gsort
INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: lge,glist,gkeys
LOGICAL, ALLOCATABLE, DIMENSION(:) :: isred,mark
REAL(r8), POINTER, DIMENSION(:) :: vals
DEBUG_STACK_PUSH
nr=A%nr
pinfo=>A%patch_info
IF(.NOT.ASSOCIATED(pinfo))CALL oft_abort("Matrix has no patches defined", &
  "patch_precond_structure",__FILE__)
CALL free_structure(self)
!---Representative rows (local periodic copies share a global index)
ALLOCATE(self%rep(nr),self%isrep(nr),self%iso(nr),self%diag(nr),self%wts(nr),self%glob(nr))
ALLOCATE(isred(nr),lge(nr),isort(nr))
isred=.FALSE.
IF(ASSOCIATED(A%redind))THEN
  DO i=1,A%nred
    isred(A%redind(i))=.TRUE.
  END DO
END IF
DO i=1,A%ni
  DO j=1,A%i_map(i)%n
    lge(A%i_map(i)%offset+j)=A%i_map(i)%offsetg+A%i_map(i)%lge(j)
  END DO
END DO
self%glob=lge
DO i=1,nr
  isort(i)=i
END DO
CALL sort_array(lge,isort,nr)
i=1
DO WHILE(i<=nr)
  j=i
  DO WHILE(j<nr)
    IF(lge(j+1)/=lge(i))EXIT
    j=j+1
  END DO
  !---Group i:j shares a global index, pick first non-redundant row
  k=-1
  DO ip=i,j
    IF(.NOT.isred(isort(ip)))THEN
      k=isort(ip)
      EXIT
    END IF
  END DO
  IF(k<0)k=isort(i)
  DO ip=i,j
    self%rep(isort(ip))=k
  END DO
  i=j+1
END DO
DO i=1,nr
  self%isrep(i)=(self%rep(i)==i)
END DO
!---Copies of each representative row
ALLOCATE(self%kcopy(nr+1),self%lcopy(nr))
self%kcopy=0
DO i=1,nr
  self%kcopy(self%rep(i)+1)=self%kcopy(self%rep(i)+1)+1
END DO
self%kcopy(1)=1
DO i=1,nr
  self%kcopy(i+1)=self%kcopy(i+1)+self%kcopy(i)
END DO
isort=0
DO i=1,nr
  j=self%rep(i)
  self%lcopy(self%kcopy(j)+isort(j))=i
  isort(j)=isort(j)+1
END DO
DEALLOCATE(isred,lge,isort)
!---Decoupled rows: no non-zero off-diagonal entries in representative columns
!   on any process (combined over processes with a sum of local indicators)
NULLIFY(vals)
ALLOCATE(vals(nr))
!$omp parallel do private(k)
DO i=1,nr
  vals(i)=0.d0
  self%diag(i)=0.d0
  DO k=A%kr(i),A%kr(i+1)-1
    IF(A%lc(k)==i)THEN
      self%diag(i)=A%Mfull(k)
    ELSE IF(self%isrep(A%lc(k)).AND.A%Mfull(k)/=0.d0)THEN
      vals(i)=1.d0
    END IF
  END DO
  IF(.NOT.self%isrep(i))vals(i)=0.d0
END DO
CALL self%rv%set(0.d0)
CALL self%rv%restore_local(vals,add=.TRUE.)
CALL self%rv%get_local(vals)
DO i=1,nr
  self%iso(i)=(vals(self%rep(i))<0.5d0).AND.(self%diag(self%rep(i))/=0.d0)
END DO
!---Build patch row lists (representatives only, decoupled rows removed)
self%npatch=pinfo%npatch
ALLOCATE(self%pf(self%npatch))
maxn=0
!$omp parallel private(j,k,n,isort) reduction(max:maxn)
ALLOCATE(isort(nr))
!$omp do schedule(dynamic,16)
DO ip=1,self%npatch
  n=0
  DO j=1,pinfo%patches(ip)%n
    k=self%rep(pinfo%patches(ip)%v(j))
    IF(self%iso(k))CYCLE
    n=n+1
    isort(n)=k
  END DO
  CALL sort_unique(isort,n)
  self%pf(ip)%n=n
  ALLOCATE(self%pf(ip)%ind(n))
  self%pf(ip)%ind=isort(1:n)
  IF(pinfo%share(ip)%n==0)THEN
    self%pf(ip)%ptype=0
  ELSE IF(pinfo%owner(ip)==oft_env%rank)THEN
    self%pf(ip)%ptype=1
  ELSE
    self%pf(ip)%ptype=2
  END IF
  maxn=MAX(maxn,n)
END DO
DEALLOCATE(isort)
!$omp end parallel
!---Seam patches: neighbors, ordered send/receive lists and global row lists
self%seams=ANY(self%pf(:)%ptype>0)
#ifdef HAVE_MPI
IF(oft_env%nprocs>1)self%seams=oft_mpi_max(MERGE(1,0,self%seams))>0
#endif
self%nnbr=0
IF(self%seams)THEN
  !---Neighbor list (processes sharing a seam patch)
  ALLOCATE(nbr_map(0:oft_env%nprocs-1))
  nbr_map=0
  DO ip=1,self%npatch
    DO j=1,pinfo%share(ip)%n
      nbr_map(pinfo%share(ip)%v(j))=1
    END DO
  END DO
  self%nnbr=SUM(nbr_map)
  ALLOCATE(self%nbr(self%nnbr))
  j=0
  DO i=0,oft_env%nprocs-1
    IF(nbr_map(i)>0)THEN
      j=j+1
      self%nbr(j)%rank=i
      nbr_map(i)=j
    END IF
  END DO
  !---Patches exchanged with each neighbor, ordered by global key
  ALLOCATE(gkeys(self%npatch),gsort(self%npatch))
  DO ip=1,self%npatch
    gkeys(ip)=pinfo%gid(ip)
    gsort(ip)=ip
  END DO
  CALL sort_array(gkeys,gsort,self%npatch)
  DO i=1,self%nnbr
    ALLOCATE(self%nbr(i)%lsend(self%npatch),self%nbr(i)%lrecv(self%npatch))
  END DO
  DO k=1,self%npatch
    ip=gsort(k)
    IF(self%pf(ip)%ptype==2)THEN
      i=nbr_map(pinfo%owner(ip))
      self%nbr(i)%nsend=self%nbr(i)%nsend+1
      self%nbr(i)%lsend(self%nbr(i)%nsend)=ip
    ELSE IF(self%pf(ip)%ptype==1)THEN
      DO j=1,pinfo%share(ip)%n
        i=nbr_map(pinfo%share(ip)%v(j))
        self%nbr(i)%nrecv=self%nbr(i)%nrecv+1
        self%nbr(i)%lrecv(self%nbr(i)%nrecv)=ip
      END DO
    END IF
  END DO
  DEALLOCATE(gkeys,gsort,nbr_map)
  !---Send global row lists of contributed patches to owners
  ALLOCATE(sbuf(self%nnbr),rbuf(self%nnbr))
  DO i=1,self%nnbr
    n=0
    DO k=1,self%nbr(i)%nsend
      n=n+1+self%pf(self%nbr(i)%lsend(k))%n
    END DO
    sbuf(i)%n=n
    ALLOCATE(sbuf(i)%v(n))
    n=0
    DO k=1,self%nbr(i)%nsend
      ip=self%nbr(i)%lsend(k)
      n=n+1
      sbuf(i)%v(n)=REAL(self%pf(ip)%n,8)
      DO j=1,self%pf(ip)%n
        sbuf(i)%v(n+j)=REAL(self%glob(self%pf(ip)%ind(j)),8)
      END DO
      n=n+self%pf(ip)%n
    END DO
  END DO
  CALL nbr_exchange_sized(self,sbuf,rbuf)
  !---Owned seam patches: union of global rows over contributors
  ALLOCATE(cnt(self%npatch))
  cnt=0
  DO ip=1,self%npatch
    IF(self%pf(ip)%ptype==1)cnt(ip)=self%pf(ip)%n
  END DO
  DO i=1,self%nnbr
    n=0
    DO k=1,self%nbr(i)%nrecv
      ip=self%nbr(i)%lrecv(k)
      jn=NINT(rbuf(i)%v(n+1))
      cnt(ip)=cnt(ip)+jn
      n=n+1+jn
    END DO
  END DO
  !---Collect global ids into per patch lists
  DO ip=1,self%npatch
    IF(self%pf(ip)%ptype/=1)CYCLE
    ALLOCATE(self%pf(ip)%gid(cnt(ip)))
    DO j=1,self%pf(ip)%n
      self%pf(ip)%gid(j)=self%glob(self%pf(ip)%ind(j))
    END DO
    cnt(ip)=self%pf(ip)%n
  END DO
  DO i=1,self%nnbr
    n=0
    DO k=1,self%nbr(i)%nrecv
      ip=self%nbr(i)%lrecv(k)
      jn=NINT(rbuf(i)%v(n+1))
      DO j=1,jn
        self%pf(ip)%gid(cnt(ip)+j)=INT(rbuf(i)%v(n+1+j),8)
      END DO
      cnt(ip)=cnt(ip)+jn
      n=n+1+jn
    END DO
  END DO
  !---Sort, remove duplicates and map to local rows (0 for rows not held locally)
  ALLOCATE(glist(nr),gsort(nr))
  nloc=0
  DO i=1,nr
    IF(self%isrep(i))THEN
      nloc=nloc+1
      glist(nloc)=self%glob(i)
      gsort(nloc)=i
    END IF
  END DO
  CALL sort_array(glist,gsort,nloc)
  DO ip=1,self%npatch
    IF(self%pf(ip)%ptype/=1)CYCLE
    n=cnt(ip)
    CALL sort_unique8(self%pf(ip)%gid,n)
    DEALLOCATE(self%pf(ip)%ind)
    ALLOCATE(self%pf(ip)%ind(n))
    self%pf(ip)%n=n
    DO j=1,n
      k=INT(search_array(self%pf(ip)%gid(j),glist,INT(nloc,8)),4)
      self%pf(ip)%ind(j)=0
      IF(k>0)self%pf(ip)%ind(j)=gsort(k)
    END DO
    maxn=MAX(maxn,n)
  END DO
  DEALLOCATE(glist,gsort,cnt)
  !---Positions of each contributor's rows in the owned patches
  DO i=1,self%nnbr
    ALLOCATE(self%nbr(i)%pos(self%nbr(i)%nrecv))
    n=0
    self%nbr(i)%nrbuf=0
    DO k=1,self%nbr(i)%nrecv
      ip=self%nbr(i)%lrecv(k)
      jn=NINT(rbuf(i)%v(n+1))
      self%nbr(i)%pos(k)%n=jn
      ALLOCATE(self%nbr(i)%pos(k)%v(jn))
      DO j=1,jn
        self%nbr(i)%pos(k)%v(j)=INT(search_array(INT(rbuf(i)%v(n+1+j),8),self%pf(ip)%gid,INT(self%pf(ip)%n,8)),4)
        IF(self%nbr(i)%pos(k)%v(j)<=0)CALL oft_abort("Seam patch row not found","patch_precond_structure",__FILE__)
      END DO
      self%nbr(i)%nrbuf=self%nbr(i)%nrbuf+jn
      n=n+1+jn
    END DO
    self%nbr(i)%nsbuf=0
    DO k=1,self%nbr(i)%nsend
      self%nbr(i)%nsbuf=self%nbr(i)%nsbuf+self%pf(self%nbr(i)%lsend(k))%n
    END DO
  END DO
  CALL free_bufs(sbuf)
  CALL free_bufs(rbuf)
END IF
self%maxn=maxn
!---Owned seam patch list
n=COUNT(self%pf(:)%ptype==1)
ALLOCATE(self%lseam(n))
n=0
DO ip=1,self%npatch
  IF(self%pf(ip)%ptype==1)THEN
    n=n+1
    self%lseam(n)=ip
  END IF
END DO
!---Interior patches by color
self%ncolors=MAX(1,pinfo%ncolors)
ALLOCATE(self%kcolor(self%ncolors+1),self%lcolor(self%npatch))
self%kcolor=0
DO ip=1,self%npatch
  IF(self%pf(ip)%ptype/=0)CYCLE
  ic=pinfo%color(ip)
  self%kcolor(ic+1)=self%kcolor(ic+1)+1
END DO
self%kcolor(1)=1
DO ic=1,self%ncolors
  self%kcolor(ic+1)=self%kcolor(ic+1)+self%kcolor(ic)
END DO
ALLOCATE(cnt(self%ncolors))
cnt=0
DO ip=1,self%npatch
  IF(self%pf(ip)%ptype/=0)CYCLE
  ic=pinfo%color(ip)
  self%lcolor(self%kcolor(ic)+cnt(ic))=ip
  cnt(ic)=cnt(ic)+1
END DO
DEALLOCATE(cnt)
!---Check that interior patches of one color do not overlap (required for threading)
ALLOCATE(mark(nr))
mark=.FALSE.
DO ic=1,self%ncolors
  DO k=self%kcolor(ic),self%kcolor(ic+1)-1
    ip=self%lcolor(k)
    DO j=1,self%pf(ip)%n
      IF(mark(self%pf(ip)%ind(j)))CALL oft_abort("Patches of one color overlap", &
        "patch_precond_structure",__FILE__)
      mark(self%pf(ip)%ind(j))=.TRUE.
    END DO
  END DO
  DO k=self%kcolor(ic),self%kcolor(ic+1)-1
    ip=self%lcolor(k)
    mark(self%pf(ip)%ind)=.FALSE.
  END DO
END DO
DEALLOCATE(mark)
!---Combination weights: 1/(number of patches containing the row), counting all
!   patches in additive mode and only seam patches otherwise. Every process
!   holding a row sees all patches containing it, so counts are consistent.
vals=0.d0
DO ip=1,self%npatch
  IF((self%mode/=1).AND.(self%pf(ip)%ptype==0))CYCLE
  DO j=1,self%pf(ip)%n
    IF(self%pf(ip)%ind(j)>0)vals(self%pf(ip)%ind(j))=vals(self%pf(ip)%ind(j))+1.d0
  END DO
END DO
DO i=1,nr
  self%wts(i)=0.d0
  IF(vals(i)>0.5d0)self%wts(i)=1.d0/vals(i)
  IF(self%isrep(i).AND.(.NOT.self%iso(i)).AND.self%mode==1.AND.vals(i)<0.5d0) &
    CALL oft_abort("Row not contained in any patch","patch_precond_structure",__FILE__)
END DO
DEALLOCATE(vals)
self%new_structure=.FALSE.
DEBUG_STACK_POP
CONTAINS
!---Sort list and remove duplicates
SUBROUTINE sort_unique(list,nl)
INTEGER(i4), INTENT(inout) :: list(:),nl
INTEGER(i4) :: jl,nu
IF(nl==0)RETURN
CALL sort_array(list,nl)
nu=1
DO jl=2,nl
  IF(list(jl)/=list(nu))THEN
    nu=nu+1
    list(nu)=list(jl)
  END IF
END DO
nl=nu
END SUBROUTINE sort_unique
!---Sort list and remove duplicates (8-byte)
SUBROUTINE sort_unique8(list,nl)
INTEGER(i8), INTENT(inout) :: list(:)
INTEGER(i4), INTENT(inout) :: nl
INTEGER(i4) :: jl,nu
IF(nl==0)RETURN
CALL sort_array(list,INT(nl,8))
nu=1
DO jl=2,nl
  IF(list(jl)/=list(nu))THEN
    nu=nu+1
    list(nu)=list(jl)
  END IF
END DO
nl=nu
END SUBROUTINE sort_unique8
END SUBROUTINE patch_precond_structure
!------------------------------------------------------------------------------
!> Free index structures and patch factorizations
!------------------------------------------------------------------------------
SUBROUTINE free_structure(self)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
INTEGER(i4) :: ip,i,k
IF(ASSOCIATED(self%pf))THEN
  DO ip=1,SIZE(self%pf)
    IF(ASSOCIATED(self%pf(ip)%ind))DEALLOCATE(self%pf(ip)%ind)
    IF(ASSOCIATED(self%pf(ip)%gid))DEALLOCATE(self%pf(ip)%gid)
    IF(ASSOCIATED(self%pf(ip)%piv))DEALLOCATE(self%pf(ip)%piv)
    IF(ASSOCIATED(self%pf(ip)%lu))DEALLOCATE(self%pf(ip)%lu)
    IF(ASSOCIATED(self%pf(ip)%lu4))DEALLOCATE(self%pf(ip)%lu4)
  END DO
  DEALLOCATE(self%pf)
END IF
IF(ASSOCIATED(self%nbr))THEN
  DO i=1,SIZE(self%nbr)
    IF(ASSOCIATED(self%nbr(i)%lsend))DEALLOCATE(self%nbr(i)%lsend)
    IF(ASSOCIATED(self%nbr(i)%lrecv))DEALLOCATE(self%nbr(i)%lrecv)
    IF(ASSOCIATED(self%nbr(i)%pos))THEN
      DO k=1,SIZE(self%nbr(i)%pos)
        IF(ASSOCIATED(self%nbr(i)%pos(k)%v))DEALLOCATE(self%nbr(i)%pos(k)%v)
      END DO
      DEALLOCATE(self%nbr(i)%pos)
    END IF
  END DO
  DEALLOCATE(self%nbr)
END IF
self%nnbr=0
IF(ASSOCIATED(self%kcolor))DEALLOCATE(self%kcolor,self%lcolor)
IF(ASSOCIATED(self%lseam))DEALLOCATE(self%lseam)
IF(ASSOCIATED(self%rep))DEALLOCATE(self%rep,self%isrep,self%iso,self%diag,self%wts,self%glob)
IF(ASSOCIATED(self%kcopy))DEALLOCATE(self%kcopy,self%lcopy)
END SUBROUTINE free_structure
!------------------------------------------------------------------------------
!> Extract and factor dense patch matrices
!!
!! Interior patches use the fully assembled local matrix. Seam patch blocks are
!! summed from the unassembled local matrices of all processes holding part of
!! the patch (every cell coupling two rows of a vertex star contains the vertex,
!! so the sum is the assembled block).
!------------------------------------------------------------------------------
SUBROUTINE patch_precond_factor(self,A)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
CLASS(oft_native_matrix), INTENT(inout) :: A !< Matrix
TYPE(oft_1d_real), ALLOCATABLE, DIMENSION(:) :: sbuf,rbuf
INTEGER(i4) :: ip,i,j,k,n,info,nfail,m,jj,kk,off
INTEGER(i8) :: nfactor
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: atmp
DEBUG_STACK_PUSH
nfail=0
nfactor=0
!---Allocate factors and assemble owned seam patches from local contributions
DO ip=1,self%npatch
  n=self%pf(ip)%n
  IF(self%pf(ip)%ptype==2.OR.n==0)CYCLE
  IF(.NOT.ASSOCIATED(self%pf(ip)%piv))ALLOCATE(self%pf(ip)%piv(n))
  IF(.NOT.ASSOCIATED(self%pf(ip)%lu))ALLOCATE(self%pf(ip)%lu(n,n))
END DO
IF(self%nnbr>0)THEN
  ALLOCATE(sbuf(self%nnbr),rbuf(self%nnbr))
  DO i=1,self%nnbr
    n=0
    DO k=1,self%nbr(i)%nsend
      n=n+self%pf(self%nbr(i)%lsend(k))%n**2
    END DO
    sbuf(i)%n=n
    ALLOCATE(sbuf(i)%v(n))
    n=0
    DO k=1,self%nbr(i)%nrecv
      n=n+self%nbr(i)%pos(k)%n**2
    END DO
    rbuf(i)%n=n
    ALLOCATE(rbuf(i)%v(n))
  END DO
END IF
!$omp parallel private(ip,i,j,k,n,m,off,info,loc,atmp) reduction(+:nfail) reduction(+:nfactor)
ALLOCATE(loc(A%nr))
loc=0
ALLOCATE(atmp(self%maxn,self%maxn))
!---Contributed seam patches: partial blocks from the unassembled local matrix
!$omp do schedule(dynamic,4)
DO i=1,self%nnbr
  off=0
  DO k=1,self%nbr(i)%nsend
    ip=self%nbr(i)%lsend(k)
    n=self%pf(ip)%n
    CALL partial_block(self%pf(ip)%ind,n,loc,atmp(1:n,1:n))
    DO j=1,n
      sbuf(i)%v(off+(j-1)*n+1:off+j*n)=atmp(1:n,j)
    END DO
    off=off+n*n
  END DO
END DO
!$omp do schedule(dynamic,4)
DO ip=1,self%npatch
  n=self%pf(ip)%n
  IF(n==0)CYCLE
  SELECT CASE(self%pf(ip)%ptype)
  CASE(0)
    CALL fill_patch(self%pf(ip)%lu,self%pf(ip)%ind,n,loc)
  CASE(1)
    !---Local part of owned seam patch (rows not held locally are skipped)
    CALL partial_block_owned(self%pf(ip)%ind,n,loc,self%pf(ip)%lu)
  END SELECT
END DO
DEALLOCATE(loc,atmp)
!$omp end parallel
!---Sum contributions into owned seam patches
IF(self%nnbr>0)THEN
  CALL nbr_exchange(self,sbuf,rbuf)
  DO i=1,self%nnbr
    off=0
    DO k=1,self%nbr(i)%nrecv
      ip=self%nbr(i)%lrecv(k)
      m=self%nbr(i)%pos(k)%n
      DO jj=1,m
        DO kk=1,m
          self%pf(ip)%lu(self%nbr(i)%pos(k)%v(kk),self%nbr(i)%pos(k)%v(jj)) = &
            self%pf(ip)%lu(self%nbr(i)%pos(k)%v(kk),self%nbr(i)%pos(k)%v(jj)) + rbuf(i)%v(off+(jj-1)*m+kk)
        END DO
      END DO
      off=off+m*m
    END DO
  END DO
  CALL free_bufs(sbuf)
  CALL free_bufs(rbuf)
END IF
!---Factor interior and owned seam patches
!$omp parallel do private(n,info) reduction(+:nfail) reduction(+:nfactor) schedule(dynamic,4)
DO ip=1,self%npatch
  n=self%pf(ip)%n
  IF(self%pf(ip)%ptype==2.OR.n==0)CYCLE
  IF(self%single)THEN
    IF(.NOT.ASSOCIATED(self%pf(ip)%lu4))ALLOCATE(self%pf(ip)%lu4(n,n))
    self%pf(ip)%lu4=REAL(self%pf(ip)%lu,4)
    DEALLOCATE(self%pf(ip)%lu)
    CALL sgetrf(n,n,self%pf(ip)%lu4,n,self%pf(ip)%piv,info)
  ELSE
    CALL dgetrf(n,n,self%pf(ip)%lu,n,self%pf(ip)%piv,info)
  END IF
  IF(info/=0)nfail=nfail+1
  nfactor=nfactor+INT(n,8)*INT(n,8)
END DO
self%nfactor=nfactor
IF(nfail>0)CALL oft_abort("Singular patch matrix encountered","patch_precond_factor",__FILE__)
self%refactor=.FALSE.
DEBUG_STACK_POP
CONTAINS
!---Copy dense patch block from fully assembled local matrix
SUBROUTINE fill_patch(amat,ind,np,lmap)
REAL(r8), INTENT(out) :: amat(:,:)
INTEGER(i4), INTENT(in) :: ind(:),np
INTEGER(i4), INTENT(inout) :: lmap(:)
INTEGER(i4) :: jl,kl,il,cl
DO jl=1,np
  lmap(ind(jl))=jl
END DO
amat=0.d0
DO jl=1,np
  il=ind(jl)
  DO kl=A%kr(il),A%kr(il+1)-1
    cl=lmap(A%lc(kl))
    IF(cl>0)amat(jl,cl)=A%Mfull(kl)
  END DO
END DO
DO jl=1,np
  lmap(ind(jl))=0
END DO
END SUBROUTINE fill_patch
!---Partial block from unassembled local matrix: rows of all local copies of each
!   representative row, columns mapped to representatives
SUBROUTINE partial_block(ind,np,lmap,amat)
INTEGER(i4), INTENT(in) :: ind(:),np
INTEGER(i4), INTENT(inout) :: lmap(:)
REAL(r8), INTENT(out) :: amat(:,:)
INTEGER(i4) :: jl,kl,il,cl,ml
DO jl=1,np
  lmap(ind(jl))=jl
END DO
amat=0.d0
DO jl=1,np
  DO ml=self%kcopy(ind(jl)),self%kcopy(ind(jl)+1)-1
    il=self%lcopy(ml)
    DO kl=A%kr(il),A%kr(il+1)-1
      cl=lmap(self%rep(A%lc(kl)))
      IF(cl>0)amat(jl,cl)=amat(jl,cl)+A%M(kl)
    END DO
  END DO
END DO
DO jl=1,np
  lmap(ind(jl))=0
END DO
END SUBROUTINE partial_block
!---Local partial block of an owned seam patch, placed at patch positions
SUBROUTINE partial_block_owned(ind,np,lmap,amat)
INTEGER(i4), INTENT(in) :: ind(:),np
INTEGER(i4), INTENT(inout) :: lmap(:)
REAL(r8), INTENT(out) :: amat(:,:)
INTEGER(i4) :: jl,kl,il,cl,ml
DO jl=1,np
  IF(ind(jl)>0)lmap(ind(jl))=jl
END DO
amat=0.d0
DO jl=1,np
  IF(ind(jl)<=0)CYCLE
  DO ml=self%kcopy(ind(jl)),self%kcopy(ind(jl)+1)-1
    il=self%lcopy(ml)
    DO kl=A%kr(il),A%kr(il+1)-1
      cl=lmap(self%rep(A%lc(kl)))
      IF(cl>0)amat(jl,cl)=amat(jl,cl)+A%M(kl)
    END DO
  END DO
END DO
DO jl=1,np
  IF(ind(jl)>0)lmap(ind(jl))=0
END DO
END SUBROUTINE partial_block_owned
END SUBROUTINE patch_precond_factor
!------------------------------------------------------------------------------
!> Solve with one dense patch factorization (in place)
!------------------------------------------------------------------------------
SUBROUTINE patch_solve(self,ip,b)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
INTEGER(i4), INTENT(in) :: ip !< Patch index
REAL(r8), INTENT(inout) :: b(:) !< RHS (input), solution (output)
INTEGER(i4) :: n,info
REAL(r4) :: b4(self%pf(ip)%n)
n=self%pf(ip)%n
IF(self%single)THEN
  b4=REAL(b(1:n),4)
  CALL sgetrs('N',n,1,self%pf(ip)%lu4,n,self%pf(ip)%piv,b4,n,info)
  b(1:n)=REAL(b4,8)
ELSE
  CALL dgetrs('N',n,1,self%pf(ip)%lu,n,self%pf(ip)%piv,b,n,info)
END IF
END SUBROUTINE patch_solve
!------------------------------------------------------------------------------
!> Additive correction from seam patches for residual `r`, accumulated into `acc`
!!
!! Residuals are gathered by the owning processes, which solve and return the
!! correction to every process holding part of each patch.
!------------------------------------------------------------------------------
SUBROUTINE patch_precond_seam(self,r,acc)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
REAL(r8), INTENT(in) :: r(:) !< Residual (local values)
REAL(r8), INTENT(inout) :: acc(:) !< Accumulated correction (local values)
TYPE(oft_1d_real), ALLOCATABLE, DIMENSION(:) :: sbuf,rbuf
TYPE(oft_1d_real), ALLOCATABLE, DIMENSION(:) :: xs
INTEGER(i4) :: i,j,k,ip,n,off,is
DEBUG_STACK_PUSH
!---Send residual values of contributed patches to owners
ALLOCATE(sbuf(self%nnbr),rbuf(self%nnbr))
DO i=1,self%nnbr
  sbuf(i)%n=self%nbr(i)%nsbuf
  rbuf(i)%n=self%nbr(i)%nrbuf
  ALLOCATE(sbuf(i)%v(sbuf(i)%n),rbuf(i)%v(rbuf(i)%n))
  off=0
  DO k=1,self%nbr(i)%nsend
    ip=self%nbr(i)%lsend(k)
    DO j=1,self%pf(ip)%n
      sbuf(i)%v(off+j)=r(self%pf(ip)%ind(j))
    END DO
    off=off+self%pf(ip)%n
  END DO
END DO
CALL nbr_exchange(self,sbuf,rbuf)
!---Assemble RHS of owned seam patches
ALLOCATE(xs(SIZE(self%lseam)))
DO is=1,SIZE(self%lseam)
  ip=self%lseam(is)
  n=self%pf(ip)%n
  xs(is)%n=n
  ALLOCATE(xs(is)%v(n))
  xs(is)%v=0.d0
  DO j=1,n
    IF(self%pf(ip)%ind(j)>0)xs(is)%v(j)=r(self%pf(ip)%ind(j))
  END DO
END DO
DO i=1,self%nnbr
  off=0
  DO k=1,self%nbr(i)%nrecv
    ip=self%nbr(i)%lrecv(k)
    is=seam_index(ip)
    DO j=1,self%nbr(i)%pos(k)%n
      xs(is)%v(self%nbr(i)%pos(k)%v(j))=rbuf(i)%v(off+j)
    END DO
    off=off+self%nbr(i)%pos(k)%n
  END DO
END DO
!---Solve owned seam patches
!$omp parallel do schedule(dynamic,4)
DO is=1,SIZE(self%lseam)
  IF(xs(is)%n>0)CALL patch_solve(self,self%lseam(is),xs(is)%v)
END DO
!---Return corrections to contributors and accumulate local parts
DO i=1,self%nnbr
  off=0
  DO k=1,self%nbr(i)%nrecv
    ip=self%nbr(i)%lrecv(k)
    is=seam_index(ip)
    DO j=1,self%nbr(i)%pos(k)%n
      rbuf(i)%v(off+j)=xs(is)%v(self%nbr(i)%pos(k)%v(j))
    END DO
    off=off+self%nbr(i)%pos(k)%n
  END DO
END DO
DO is=1,SIZE(self%lseam)
  ip=self%lseam(is)
  DO j=1,self%pf(ip)%n
    IF(self%pf(ip)%ind(j)>0)acc(self%pf(ip)%ind(j))=acc(self%pf(ip)%ind(j))+xs(is)%v(j)
  END DO
END DO
CALL nbr_exchange(self,rbuf,sbuf)
DO i=1,self%nnbr
  off=0
  DO k=1,self%nbr(i)%nsend
    ip=self%nbr(i)%lsend(k)
    DO j=1,self%pf(ip)%n
      acc(self%pf(ip)%ind(j))=acc(self%pf(ip)%ind(j))+sbuf(i)%v(off+j)
    END DO
    off=off+self%pf(ip)%n
  END DO
END DO
CALL free_bufs(sbuf)
CALL free_bufs(rbuf)
CALL free_bufs(xs)
DEBUG_STACK_POP
CONTAINS
!---Position of owned seam patch in list
FUNCTION seam_index(jp) RESULT(js)
INTEGER(i4), INTENT(in) :: jp
INTEGER(i4) :: js
DO js=1,SIZE(self%lseam)
  IF(self%lseam(js)==jp)RETURN
END DO
CALL oft_abort("Seam patch not found","seam_index",__FILE__)
END FUNCTION seam_index
END SUBROUTINE patch_precond_seam
!------------------------------------------------------------------------------
!> Compute patch correction \f$ du = B r \f$ (local arrays, consistent on output)
!------------------------------------------------------------------------------
SUBROUTINE patch_precond_sweep(self,A,r,du)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
CLASS(oft_native_matrix), INTENT(inout) :: A !< Matrix
REAL(r8), INTENT(in) :: r(:) !< Residual (local values)
REAL(r8), INTENT(inout) :: du(:) !< Correction (local values)
INTEGER(i4) :: i,ic,k,ip
REAL(r8), POINTER, DIMENSION(:) :: acc,r2
DEBUG_STACK_PUSH
!---Decoupled rows
!$omp parallel do
DO i=1,A%nr
  du(i)=0.d0
  IF(self%isrep(i).AND.self%iso(i))du(i)=r(i)/self%diag(i)
END DO
IF(self%mode==1)THEN
  !---Additive: average of all patch corrections
  ALLOCATE(acc(A%nr))
  acc=0.d0
  DO ic=1,self%ncolors
    !$omp parallel do private(ip) schedule(dynamic,4)
    DO k=self%kcolor(ic),self%kcolor(ic+1)-1
      ip=self%lcolor(k)
      IF(self%pf(ip)%n==0)CYCLE
      CALL solve_patch(ip,r,acc,.FALSE.)
    END DO
  END DO
  IF(self%seams)CALL patch_precond_seam(self,r,acc)
  !$omp parallel do
  DO i=1,A%nr
    IF(self%isrep(i).AND.(.NOT.self%iso(i)))du(i)=du(i)+self%wts(i)*acc(i)
  END DO
  DEALLOCATE(acc)
ELSE
  !---Forward sweep over interior patches
  CALL interior_sweep(.FALSE.)
  !---Seam patches with updated global residual
  IF(self%seams)THEN
    NULLIFY(r2)
    CALL copy_reps()
    CALL self%duv%restore_local(du)
    CALL A%apply(self%duv,self%rv)
    CALL self%rv%get_local(r2)
    !$omp parallel do
    DO i=1,A%nr
      r2(i)=r(i)-r2(i)
    END DO
    ALLOCATE(acc(A%nr))
    acc=0.d0
    CALL patch_precond_seam(self,r2,acc)
    !$omp parallel do
    DO i=1,A%nr
      IF(self%isrep(i).AND.(.NOT.self%iso(i)))du(i)=du(i)+self%wts(i)*acc(i)
    END DO
    DEALLOCATE(acc,r2)
  END IF
  !---Backward sweep over interior patches
  IF(self%mode==3)CALL interior_sweep(.TRUE.)
END IF
CALL copy_reps()
DEBUG_STACK_POP
CONTAINS
!---Copy representative values to redundant copies
SUBROUTINE copy_reps()
INTEGER(i4) :: j
!$omp parallel do
DO j=1,A%nr
  IF(.NOT.self%isrep(j))du(j)=du(self%rep(j))
END DO
END SUBROUTINE copy_reps
!---Multiplicative sweep over interior patches by color
SUBROUTINE interior_sweep(backward)
LOGICAL, INTENT(in) :: backward
INTEGER(i4) :: jc,jcc,kp,jp
DO jcc=1,self%ncolors
  jc=jcc
  IF(backward)jc=self%ncolors-jcc+1
  !$omp parallel do private(jp) schedule(dynamic,4)
  DO kp=self%kcolor(jc),self%kcolor(jc+1)-1
    jp=self%lcolor(kp)
    IF(self%pf(jp)%n==0)CYCLE
    CALL solve_patch(jp,r,du,.TRUE.)
  END DO
END DO
END SUBROUTINE interior_sweep
!---Solve a single interior patch and update the correction
SUBROUTINE solve_patch(jp,rin,uacc,resid)
INTEGER(i4), INTENT(in) :: jp
REAL(r8), INTENT(in) :: rin(:)
REAL(r8), INTENT(inout) :: uacc(:)
LOGICAL, INTENT(in) :: resid
INTEGER(i4) :: jj,kk,ii,n
REAL(r8) :: b(self%pf(jp)%n)
n=self%pf(jp)%n
DO jj=1,n
  ii=self%pf(jp)%ind(jj)
  b(jj)=rin(ii)
  !---Local residual with current correction (representative columns only)
  IF(resid)THEN
    DO kk=A%kr(ii),A%kr(ii+1)-1
      IF(self%isrep(A%lc(kk)))b(jj)=b(jj)-A%Mfull(kk)*uacc(A%lc(kk))
    END DO
  END IF
END DO
CALL patch_solve(self,jp,b)
DO jj=1,n
  ii=self%pf(jp)%ind(jj)
  uacc(ii)=uacc(ii)+b(jj)
END DO
END SUBROUTINE solve_patch
END SUBROUTINE patch_precond_sweep
!------------------------------------------------------------------------------
!> Apply patch preconditioner (`its<=0`) or smoother (`its>0`)
!------------------------------------------------------------------------------
SUBROUTINE patch_precond_apply(self,u,g)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
CLASS(oft_vector), INTENT(inout) :: u !< Guess (input), Solution (output)
CLASS(oft_vector), INTENT(inout) :: g !< RHS (input), Residual (output)
INTEGER(i4) :: i,k
REAL(r8), POINTER, DIMENSION(:) :: rloc,duloc
DEBUG_STACK_PUSH
SELECT TYPE(A=>self%A)
CLASS IS(oft_native_matrix)
  !---Setup
  IF(.NOT.self%initialized)THEN
    CALL native_matrix_setup_full(A,u)
    CALL u%new(self%rv)
    CALL u%new(self%duv)
    CALL solver_setup(self)
  END IF
  IF(self%refactor.OR.self%new_structure)THEN
    CALL A%update_slice
    CALL patch_precond_structure(self,A)
    CALL patch_precond_factor(self,A)
    IF(oft_debug_print(1))CALL self%view
  END IF
  NULLIFY(rloc,duloc)
  ALLOCATE(duloc(A%nr))
  IF(self%its<=0)THEN
    !---Preconditioner: u = df * B g
    CALL g%get_local(rloc)
    CALL patch_precond_sweep(self,A,rloc,duloc)
    duloc=self%df*duloc
    CALL u%restore_local(duloc)
  ELSE
    !---Damped Richardson iteration
    CALL self%A%apply(u,self%rv)
    CALL self%rv%add(-1.d0,1.d0,g)
    DO k=1,self%its
      CALL self%rv%get_local(rloc)
      CALL patch_precond_sweep(self,A,rloc,duloc)
      duloc=self%df*duloc
      CALL self%duv%restore_local(duloc)
      IF(ASSOCIATED(self%bc))CALL self%bc%apply(self%duv)
      CALL u%add(1.d0,1.d0,self%duv)
      CALL self%A%apply(u,self%rv)
      CALL self%rv%add(-1.d0,1.d0,g)
    END DO
    CALL g%add(0.d0,1.d0,self%rv)
  END IF
  DEALLOCATE(rloc,duloc)
CLASS DEFAULT
  CALL oft_abort("Patch smoother requires a native matrix","patch_precond_apply",__FILE__)
END SELECT
DEBUG_STACK_POP
END SUBROUTINE patch_precond_apply
!------------------------------------------------------------------------------
!> Mark patch factorizations for update after a matrix change
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE patch_precond_update(self,new_pattern)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
LOGICAL, OPTIONAL, INTENT(in) :: new_pattern !< Update matrix non-zero pattern? (optional)
self%refactor=.TRUE.
IF(PRESENT(new_pattern))THEN
  IF(new_pattern)self%new_structure=.TRUE.
END IF
END SUBROUTINE patch_precond_update
!------------------------------------------------------------------------------
!> Setup patch smoother from XML definition
!!
!! Supported elements:
!! - `its`: Number of Richardson iterations, `<=0` for preconditioner use (per-level list allowed)
!! - `df`: Damping factor (per-level list allowed)
!! - `mode`: `additive`, `multiplicative` or `symmetric`
!! - `single`: Store patch factors in single precision (`T`/`F`)
!------------------------------------------------------------------------------
SUBROUTINE patch_precond_setup_xml(self,solver_node,level)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
TYPE(xml_node) :: current_node
INTEGER(i4) :: val_level,ierr
INTEGER(i4), POINTER, DIMENSION(:) :: read_int
REAL(r8), POINTER, DIMENSION(:) :: read_real
CHARACTER(LEN=:), ALLOCATABLE :: mode_str
DEBUG_STACK_PUSH
val_level=1
IF(PRESENT(level))val_level=level
!---Iterations
CALL xml_get_element(solver_node,"its",current_node,ierr)
IF(ierr==0)THEN
  NULLIFY(read_int)
  CALL xml_read_content(current_node,read_int,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `its` node","patch_precond_setup_xml",__FILE__)
  IF(ASSOCIATED(read_int))THEN
    IF(SIZE(read_int)>1)THEN
      IF(SIZE(read_int)<val_level)CALL oft_abort("Not enough `its` values specified", &
        "patch_precond_setup_xml",__FILE__)
      self%its=read_int(val_level)
    ELSE
      self%its=read_int(1)
    END IF
    DEALLOCATE(read_int)
  END IF
END IF
!---Damping factor
CALL xml_get_element(solver_node,"df",current_node,ierr)
IF(ierr==0)THEN
  NULLIFY(read_real)
  CALL xml_read_content(current_node,read_real,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `df` node","patch_precond_setup_xml",__FILE__)
  IF(ASSOCIATED(read_real))THEN
    IF(SIZE(read_real)>1)THEN
      IF(SIZE(read_real)<val_level)CALL oft_abort("Not enough `df` values specified", &
        "patch_precond_setup_xml",__FILE__)
      self%df=read_real(val_level)
    ELSE
      self%df=read_real(1)
    END IF
    DEALLOCATE(read_real)
  END IF
END IF
!---Combination mode
CALL xml_get_element(solver_node,"mode",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,mode_str,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `mode` node","patch_precond_setup_xml",__FILE__)
  SELECT CASE(TRIM(ADJUSTL(mode_str)))
  CASE("additive")
    self%mode=1
  CASE("multiplicative")
    self%mode=2
  CASE("symmetric")
    self%mode=3
  CASE DEFAULT
    CALL oft_abort("Unknown patch smoother mode '"//TRIM(mode_str)//"'","patch_precond_setup_xml",__FILE__)
  END SELECT
  DEALLOCATE(mode_str)
END IF
!---Single precision factors
CALL xml_get_element(solver_node,"single",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%single,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `single` node","patch_precond_setup_xml",__FILE__)
END IF
DEBUG_STACK_POP
END SUBROUTINE patch_precond_setup_xml
!------------------------------------------------------------------------------
!> Print patch smoother information
!------------------------------------------------------------------------------
SUBROUTINE patch_precond_view(self)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
CHARACTER(LEN=14), PARAMETER :: mode_names(3) = (/'additive      ','multiplicative','symmetric     '/)
REAL(r8) :: mem
mem=REAL(self%nfactor,8)*8.d0/1.d6
IF(self%single)mem=mem/2.d0
WRITE(*,'(2X,A)')'Patch smoother:'
WRITE(*,'(4X,A,A)')    '- Mode:      ',TRIM(mode_names(self%mode))
WRITE(*,'(4X,A,I8)')   '- Its:       ',self%its
WRITE(*,'(4X,A,ES11.3)')'- Damping:   ',self%df
WRITE(*,'(4X,A,I8)')   '- Patches:   ',self%npatch
WRITE(*,'(4X,A,I8)')   '- Colors:    ',self%ncolors
WRITE(*,'(4X,A,I8)')   '- Owned seam:',SIZE(self%lseam)
WRITE(*,'(4X,A,I8)')   '- Neighbors: ',self%nnbr
WRITE(*,'(4X,A,I8)')   '- Max size:  ',self%maxn
WRITE(*,'(4X,A,F11.1)')'- Mem (MB):  ',mem
END SUBROUTINE patch_precond_view
!------------------------------------------------------------------------------
!> Patch smoother is not thread safe (it uses OpenMP internally)
!------------------------------------------------------------------------------
FUNCTION patch_precond_check_thread(self) RESULT(thread_safe)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
LOGICAL :: thread_safe
thread_safe=.FALSE.
END FUNCTION patch_precond_check_thread
!------------------------------------------------------------------------------
!> Destroy patch smoother and deallocate all internal storage
!------------------------------------------------------------------------------
SUBROUTINE patch_precond_delete(self,propogate)
CLASS(oft_patch_precond), INTENT(inout) :: self !< Solver object
LOGICAL, OPTIONAL, INTENT(in) :: propogate !< Delete child solvers? (optional)
DEBUG_STACK_PUSH
CALL free_structure(self)
IF(ASSOCIATED(self%rv))THEN
  CALL self%rv%delete
  DEALLOCATE(self%rv)
END IF
IF(ASSOCIATED(self%duv))THEN
  CALL self%duv%delete
  DEALLOCATE(self%duv)
END IF
self%npatch=0
self%ncolors=0
self%refactor=.TRUE.
self%new_structure=.TRUE.
self%initialized=.FALSE.
CALL solver_delete(self,propogate)
DEBUG_STACK_POP
END SUBROUTINE patch_precond_delete
END MODULE oft_patch_solvers
