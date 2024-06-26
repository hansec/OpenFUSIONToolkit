!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file taylor.f90
!
!> @defgroup doxy_oft_physics Physics
!! Physics modules for the Open FUSION Toolkit
!
!> Subroutines and fields for Taylor state calculations using
!! mimetic operators.
!! - Force-Free eigenmodes
!! - Vacuum fields for geometries with cut planes
!! - Inhomogeneous force free states from vacuum fields
!!
!! @author Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
module taylor
USE oft_base
USE oft_io
USE multigrid, ONLY: mg_mesh
!---
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_cg_solver, oft_native_cg_eigsolver, &
  oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---
USE fem_base, ONLY: oft_fem_type
USE fem_utils, ONLY: fem_interp
USE oft_lag_basis, ONLY: oft_lagrange_level, oft_lagrange_lin_level, &
oft_lagrange_nlevels, oft_lag_set_level
USE oft_lag_operators, ONLY: lag_zerob, lag_getlop_pre, oft_lag_getlop
USE oft_h0_basis, ONLY: oft_h0, oft_h0_nlevels, oft_h0_geval_all, oft_h0_fem
USE oft_h0_operators, ONLY: h0_zerogrnd, h0_getlop_pre, oft_h0_getlop
USE oft_hcurl_basis, ONLY: oft_hcurl, oft_hcurl_level, oft_hcurl_nlevels, &
oft_hcurl_blevel, oft_hcurl_eval_all, oft_hcurl_ceval_all, oft_hcurl_set_level, &
oft_hcurl_get_cgops, oft_hcurl_fem
USE oft_hcurl_fields, ONLY: oft_hcurl_create!, oft_hcurl_save
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, oft_hcurl_orthog, &
oft_hcurl_divout, hcurl_getwop_pre, hcurl_zerob, oft_hcurl_getmop, oft_hcurl_getkop, &
oft_hcurl_getwop, hcurl_interp, oft_hcurl_getjmlb, hcurl_getjmlb_pre
USE oft_h1_basis, ONLY: oft_h1, oft_h1_level, oft_h1_nlevels, oft_h1_set_level
USE oft_h1_fields, ONLY: oft_h1_create!, oft_h1_save
USE oft_h1_operators, ONLY: oft_h1_divout, h1_getmop, h1_mc
!---
USE diagnostic, ONLY: tfluxfun
implicit none
#include "local.h"
!---------------------------------------------------------------------------
! CLASS oft_taylor_interp
!---------------------------------------------------------------------------
!> Interpolate a Taylor state field
!!
!! Taylor state fields consist of a gradient component, the vacuum field, and
!! a curl component, the plasma field.
!---------------------------------------------------------------------------
type, extends(fem_interp) :: oft_taylor_rinterp
  class(oft_vector), pointer :: ua => NULL() !< Plasma vector potential
  class(oft_vector), pointer :: uvac => NULL() !< Vacuum magnectic field
  real(r8), pointer, dimension(:) :: vac_grad => NULL() !< Local vacuum field (gradient)
  real(r8), pointer, dimension(:) :: vac_curl => NULL() !< Local vacuum field (curl)
  real(r8), pointer, dimension(:) :: acurl => NULL() !< Local vector potential
  class(oft_h0_fem), pointer :: hgrad_rep => NULL() !< H1(Grad) FE representation
  class(oft_hcurl_fem), pointer :: hcurl_rep => NULL() !< H1(Curl) FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => taylor_rinterp_setup
  !> Reconstruct field
  procedure :: interp => taylor_rinterp
end type oft_taylor_rinterp
!---Force-free eigenmode variables
integer(i4) :: taylor_minlev=-1 !< Lowest FE level for MG solvers
integer(i4) :: taylor_nm=0 !< Number of force-free fields to be computed
real(r8), pointer, dimension(:,:) :: taylor_hlam => NULL() !< Homogeneous force-free lambdas
real(r8), pointer, dimension(:,:) :: taylor_htor => NULL() !< Homogeneous force-free toroidal fluxes
type(oft_vector_ptr), pointer, dimension(:,:) :: taylor_hffa => NULL() !< Homogeneous force-free fields
!---Inhomogeneous fields
integer(i4) :: taylor_nh=2 !< Number of jump planes in current geometry
real(r8), pointer, dimension(:,:)  :: taylor_hcpc => NULL() !< Center points of jump planes
real(r8), pointer, dimension(:,:)  :: taylor_hcpv => NULL() !< Normal vectors for jump planes
real(r8)  :: taylor_jtol = 1.d-6 !< Tolerance for identifying edges on jump plane
integer(i4), parameter :: taylor_tag_size = 4 !< Size of  jump planes character tags
character(LEN=taylor_tag_size), pointer, dimension(:) :: taylor_htag => NULL() !< Injector names
type(oft_vector_ptr), pointer, dimension(:,:) :: taylor_hvac => NULL() !< Vacuum magnetic fields
type(oft_vector_ptr), pointer, dimension(:,:) :: taylor_hcur => NULL() !< Inhomogeneous source fields
type(oft_vector_ptr), pointer, dimension(:,:) :: taylor_gffa => NULL() !< Inhomogeneous force-free fields
!---General
logical :: taylor_rst=.TRUE. !< Save solutions to data files
contains
!---------------------------------------------------------------------------
! SUBROUTINE: taylor_hmodes
!---------------------------------------------------------------------------
!> Compute 'taylor_nm' Force-Free eignemodes.
!!
!! @note When `taylor_rst=.TRUE.` one restart files will be generated for
!! each computed mode on each MG level `hffa_*.rst`.
!!
!! @param[in] nm Number of modes to compute (optional: 1)
!---------------------------------------------------------------------------
subroutine taylor_hmodes(nm)
integer(i4), optional, intent(in) :: nm
class(oft_vector), pointer :: u,tmp
!--- Taylor eigenvalue solver
TYPE(oft_native_cg_eigsolver) :: eigsolver
type(oft_hcurl_orthog), target :: orthog
!--- Divergence cleaner
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_hcurl_divout) :: divout
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_wop => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_lop => NULL()
INTEGER(i4) :: nlevels_wop,nlevels_lop
!--- Matrix pointers
class(oft_matrix), pointer :: kop => NULL()
class(oft_matrix), pointer :: wop => NULL()
class(oft_matrix), pointer :: lop => NULL()
!--- Local variables
type(oft_hcurl_cinterp) :: Bfield
real(r8) :: alam,elapsed_time
integer(i4) :: i,j,k,ierr
character(2) :: pnum,mnum
character(40) :: filename
logical :: rst
type(oft_timer) :: mytimer
!---
DEBUG_STACK_PUSH
IF(PRESENT(nm))THEN
  IF(nm<=0.OR.nm>20)CALL oft_abort("Invalid number of modes requested.", "taylor_hmodes", __FILE__)
  taylor_nm=nm
ELSE
  taylor_nm=1
END IF
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'============================'
  WRITE(*,'(A)')'Starting calculation of Taylor states'
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
  CALL mytimer%tick
END IF
IF(taylor_minlev<0)taylor_minlev=oft_hcurl_nlevels
!---Allocate storage
ALLOCATE(taylor_hffa(taylor_nm,oft_hcurl_nlevels))
ALLOCATE(taylor_hlam(taylor_nm,oft_hcurl_nlevels))
ALLOCATE(taylor_htor(taylor_nm,oft_hcurl_nlevels))
!---Setup orthogonalization
orthog%orthog=>taylor_hffa
!---------------------------------------------------------------------------
! Create ML Matrices
!---------------------------------------------------------------------------
nlevels_wop=oft_hcurl_nlevels-taylor_minlev+1
nlevels_lop=oft_lagrange_lin_level-taylor_minlev+1
ALLOCATE(ml_wop(nlevels_wop),ml_lop(nlevels_lop))
DO i=1,nlevels_wop
  CALL oft_hcurl_set_level(taylor_minlev+i-1)
  NULLIFY(ml_wop(i)%M)
  CALL oft_hcurl_getwop(ml_wop(i)%M,'zerob')
  IF(oft_hcurl_level<=oft_lagrange_lin_level)THEN
    CALL oft_lag_set_level(oft_hcurl_level)
    NULLIFY(ml_lop(i)%M)
    CALL oft_lag_getlop(ml_lop(i)%M,"zerob")
  END IF
END DO
CALL oft_hcurl_set_level(oft_hcurl_nlevels)
!---------------------------------------------------------------------------
! Loop over desired number of modes
!---------------------------------------------------------------------------
do k=1,taylor_nm
  orthog%nm=k-1
  call oft_hcurl_set_level(taylor_minlev)
  !---Loop over levels for each mode
  do i=taylor_minlev,oft_hcurl_nlevels
    !---Build level field
    call oft_hcurl_create(taylor_hffa(k,i)%f,i)
    !---Alias to general field
    u=>taylor_hffa(k,i)%f
!---------------------------------------------------------------------------
! Initialize guess
!---------------------------------------------------------------------------
    if(i==taylor_minlev)then
      call u%set(1.d0,random=.TRUE.)
    else
      call hcurl_interp(taylor_hffa(k,i-1)%f,u)
      if(k/=1)call u%set(1.d0,random=.TRUE.)
      if(oft_hcurl_level==oft_hcurl_blevel+1)cycle
    end if
    if(taylor_rst)then
      write(pnum,'(I2.2)')oft_hcurl%order
      write(mnum,'(I2.2)')k
      filename='hffa_r'//mg_mesh%rlevel//'_p'//pnum//'_m'//mnum//'.rst'
      if(oft_file_exist(filename))then
        CALL oft_hcurl%vec_load(u,filename,'hffa')
      end if
    end if
    !---Apply BC
    CALL hcurl_zerob(u)
    !---Setup Solver
    wop=>ml_wop(i-taylor_minlev+1)%M
    orthog%wop=>wop
    NULLIFY(kop)
    CALL oft_hcurl_getkop(kop,"zerob")
!---------------------------------------------------------------------------
! Create eigenvalue solver
!---------------------------------------------------------------------------
    eigsolver%A=>wop
    eigsolver%M=>kop
    eigsolver%its=-2
    eigsolver%bc=>hcurl_zerob
    eigsolver%orthog=>orthog
    IF(k>1)THEN
      eigsolver%nrestarts=4
      IF(i==taylor_minlev)THEN
        eigsolver%ninner=1000
      ELSE
        eigsolver%ninner=100
      END IF
    END IF
!---------------------------------------------------------------------------
! Setup preconditioner
!---------------------------------------------------------------------------
    if(i==taylor_minlev)then ! Lowest level uses diag precond
      !---Setup Preconditioner
      CALL create_diag_pre(eigsolver%pre)
    else ! Higher levels use MG
      CALL hcurl_getwop_pre(eigsolver%pre,ml_wop,nlevels=i-taylor_minlev+1)
    end if
    !---Solve
    CALL eigsolver%apply(u,taylor_hlam(k,i))
    CALL eigsolver%pre%delete
    DEALLOCATE(eigsolver%pre)
    CALL eigsolver%delete
    CALL kop%delete
    DEALLOCATE(kop)
    !---
    Bfield%u=>u
    CALL Bfield%setup
    taylor_htor(k,i) = tfluxfun(Bfield,oft_hcurl%quad%order)
    CALL Bfield%delete
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
    CALL oft_lag_set_level(MIN(oft_hcurl_level,oft_lagrange_lin_level))
    IF(nlevels_lop<=0)THEN
      NULLIFY(lop)
      CALL oft_lag_getlop(lop,"zerob")
    ELSE
      lop=>ml_lop(oft_lagrange_level-taylor_minlev+1)%M
    END IF
    if(oft_lagrange_level<=taylor_minlev)then ! Lowest level uses diag precond
      CALL create_cg_solver(linv)
      CALL create_diag_pre(linv%pre)
    else ! Higher levels use MG
      CALL create_cg_solver(linv, force_native=.TRUE.)
      CALL lag_getlop_pre(linv%pre,ml_lop,nlevels=oft_lagrange_level-taylor_minlev+1)
    end if
    linv%A=>lop
    linv%its=-2
    divout%solver=>linv
    divout%bc=>lag_zerob
    divout%pm=oft_env%pm
    CALL divout%apply(u)
    CALL linv%pre%delete
    DEALLOCATE(linv%pre)
    CALL linv%delete
    DEALLOCATE(linv)
    IF(nlevels_lop<=0)CALL lop%delete
!---------------------------------------------------------------------------
! Write restart files
!---------------------------------------------------------------------------
    if(taylor_rst)then
      write(pnum,'(I2.2)')oft_hcurl%order
      write(mnum,'(I2.2)')k
      filename='hffa_r'//mg_mesh%rlevel//'_p'//pnum//'_m'//mnum//'.rst'
      if(.NOT.oft_file_exist(filename))then
        CALL oft_mpi_barrier(ierr)
        CALL oft_hcurl%vec_save(u,filename,'hffa')
      end if
    end if
  end do ! End level loop
end do ! End mode loop
if(oft_env%head_proc)then
  DO i=taylor_minlev,oft_hcurl_nlevels
    WRITE(*,'(2X,A,I3)')'Level =',i
    DO k=1,taylor_nm
      WRITE(*,'(4X,A,I4,A,ES14.6)')'Mode =',k,'   Lambda = ',taylor_hlam(k,i)
    END DO
  END DO
  elapsed_time=mytimer%tock()
  WRITE(*,*)
  WRITE(*,'(2X,A,F12.3)')'Time Elapsed = ',elapsed_time
end if
!---------------------------------------------------------------------------
! Deallocate ML matrices
!---------------------------------------------------------------------------
DO i=1,nlevels_wop
  CALL ml_wop(i)%M%delete
  DEALLOCATE(ml_wop(i)%M)
  IF(i<nlevels_lop)THEN
    CALL ml_lop(i)%M%delete
    DEALLOCATE(ml_lop(i)%M)
  END IF
END DO
DEALLOCATE(ml_wop,ml_lop)
CALL oft_lag_set_level(oft_lagrange_nlevels)
DEBUG_STACK_POP
end subroutine taylor_hmodes
!---------------------------------------------------------------------------
! SUBROUTINE: taylor_vacuum
!---------------------------------------------------------------------------
!> Generate vacuum fields for a geometry with cut planes.
!!
!! @note When `taylor_rst=.TRUE.` two restart files will be generated for
!! each jump plane `hvac_*.rst` and `hcur_*.rst`
!!
!! @param[in] nh Number of jump planes
!! @param[in] hcpc Jump plane center possitions [3,nh]
!! @param[in] hcpv Jump plane normal vectors [3,nh]
!! @param[in] htags Names for each jump plane [LEN=taylor_tag_size,nh] (optional)
!---------------------------------------------------------------------------
subroutine taylor_vacuum(nh,hcpc,hcpv,htags,energy)
integer(i4), intent(in) :: nh
real(r8), intent(in) :: hcpc(3,nh)
real(r8), intent(in) :: hcpv(3,nh)
character(LEN=taylor_tag_size), optional, intent(in) :: htags(nh)
real(r8), optional, intent(out) :: energy(nh)
!---H1 Divergence cleaner
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_h1_divout) :: divout
!---WOP solver
CLASS(oft_solver), POINTER :: winv => NULL()
!---H1(Curl) Divergence cleaner
CLASS(oft_solver), POINTER :: linv_lag => NULL()
TYPE(oft_hcurl_divout), TARGET :: hcurl_divout
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_int => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_lop => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_wop => NULL()
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: levels,nu
REAL(r8), ALLOCATABLE, DIMENSION(:) :: df
!---
class(oft_matrix), pointer :: lop => NULL()
class(oft_matrix), pointer :: lop_lag => NULL()
class(oft_matrix), pointer :: mop => NULL()
class(oft_matrix), pointer :: mop_hcurl => NULL()
class(oft_matrix), pointer :: wop => NULL()
!---
class(oft_vector), pointer :: u,b,tmp
real(r8), pointer, dimension(:) :: vals => NULL()
real(r8) :: alam,venergy
integer(i4) :: i,j,k,nlevels,ierr
logical :: rst=.FALSE.
character(2) :: pnum,mnum
character(40) :: filename
DEBUG_STACK_PUSH
NULLIFY(lop,lop_lag,mop,mop_hcurl,wop)
NULLIFY(ml_int,ml_lop,ml_wop)
!---Create taylor module variables
taylor_nh=nh
ALLOCATE(taylor_hcpc(3,taylor_nh),taylor_hcpv(3,taylor_nh))
taylor_hcpc=hcpc
taylor_hcpv=hcpv
ALLOCATE(taylor_htag(taylor_nh))
IF(PRESENT(htags))THEN
  taylor_htag=htags
ELSE
  DO i=1,taylor_nh
    WRITE(taylor_htag(i),'(A3,I1)')'INJ',i
  END DO
END IF
!---Loop over cut planes
IF(taylor_rst)THEN
  rst=.TRUE.
  DO i=1,taylor_nh
    WRITE(pnum,'(I2.2)')oft_hcurl%order
    WRITE(mnum,'(I2.2)')i
    filename='hvac_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    rst=rst.AND.oft_file_exist(filename)
    filename='hcur_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    rst=rst.AND.oft_file_exist(filename)
  END DO
ELSE
  rst=.FALSE.
END IF
IF(taylor_minlev<0)taylor_minlev=oft_hcurl_nlevels
IF(.NOT.rst)THEN
!---------------------------------------------------------------------------
! Setup H0::LOP preconditioner
!---------------------------------------------------------------------------
  if(taylor_minlev==oft_h0_nlevels-1)then ! Lowest level uses diag precond
    CALL oft_h0_getlop(lop,'grnd')
    CALL create_cg_solver(linv)
    CALL create_diag_pre(linv%pre)
  else ! Nested levels use MG
    CALL create_cg_solver(linv, force_native=.TRUE.)
    CALL h0_getlop_pre(linv%pre,ml_lop,nlevels=oft_h0_nlevels-taylor_minlev+1)
    lop=>ml_lop(oft_h0_nlevels-taylor_minlev+1)%M
  end if
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
  linv%A=>lop
  linv%its=-2
  divout%solver=>linv
  divout%bc=>h0_zerogrnd
ELSE
  CALL create_cg_solver(linv)
  CALL create_diag_pre(linv%pre)
END IF
!---------------------------------------------------------------------------
! Setup H1(Curl)::WOP preconditioner
!---------------------------------------------------------------------------
CALL create_cg_solver(winv, force_native=.TRUE.)
IF((taylor_minlev==oft_hcurl_nlevels).OR.rst)THEN ! Lowest level uses diag precond
  NULLIFY(wop)
  CALL oft_hcurl_getwop(wop,'zerob')
  CALL create_diag_pre(winv%pre)
ELSE ! Nested levels use MG
  CALL hcurl_getwop_pre(winv%pre,ml_wop,nlevels=oft_hcurl_nlevels-taylor_minlev+1)
  wop=>ml_wop(oft_hcurl_nlevels-taylor_minlev+1)%M
END IF
!---------------------------------------------------------------------------
! Create H1(Curl)::WOP solver
!---------------------------------------------------------------------------
winv%A=>wop
winv%its=-3
winv%atol=1.d-9
SELECT TYPE(this=>winv)
CLASS IS(oft_native_cg_solver)
    this%cleaner=>hcurl_divout
  CLASS DEFAULT
    CALL oft_abort('Error allocating winv solver', 'taylor_vacuum', __FILE__)
END SELECT
!---------------------------------------------------------------------------
! Create H1(Curl) divergence cleaner
!---------------------------------------------------------------------------
CALL oft_lag_set_level(MIN(oft_hcurl_level,oft_lagrange_lin_level))
CALL oft_lag_getlop(lop_lag,"zerob")
CALL create_cg_solver(linv_lag)
linv_lag%A=>lop_lag
!linv_lag%its=-2
CALL create_diag_pre(linv_lag%pre)
linv_lag%its=40
hcurl_divout%solver=>linv_lag
hcurl_divout%bc=>lag_zerob
hcurl_divout%app_freq=2
CALL oft_hcurl_getmop(mop_hcurl,'zerob')
hcurl_divout%mop=>mop_hcurl
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
NULLIFY(tmp)
!---Get H1 mass matrix
CALL h1_getmop(mop,'none')
!---Allocate vacuum and current field containers
ALLOCATE(taylor_hvac(taylor_nh,oft_h1_nlevels))
ALLOCATE(taylor_hcur(taylor_nh,oft_h1_nlevels))
!---Create temporary H1(Curl) vector
CALL oft_hcurl_create(b)
!---Loop over cut planes
DO i=1,taylor_nh
!---------------------------------------------------------------------------
! Compute vacuum fields
!---------------------------------------------------------------------------
  !---Setup level fields
  CALL oft_h1_create(taylor_hvac(i,oft_h1_level)%f)
  u=>taylor_hvac(i,oft_h1_level)%f
  IF(.NOT.ASSOCIATED(tmp))CALL u%new(tmp)
  rst=.FALSE.
  IF(taylor_rst)THEN
    WRITE(pnum,'(I2.2)')oft_hcurl%order
    WRITE(mnum,'(I2.2)')i
    filename='hvac_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    IF(oft_file_exist(filename))THEN
      CALL oft_h1%vec_load(u,filename,'hvac')
      rst=.TRUE.
    END IF
  END IF
  !---Compute jump and solve
  IF(.NOT.rst)THEN
    CALL h1_mc(u,taylor_hcpc(:,i),taylor_hcpv(:,i),taylor_jtol)
    venergy=u%dot(u)
    IF(venergy<1.d-12)CALL oft_abort('Plane does not intersect mesh.','taylor_vacuum',__FILE__)
    divout%pm=oft_env%pm
    CALL divout%apply(u)
  END IF
!---------------------------------------------------------------------------
! Write restart file
!---------------------------------------------------------------------------
  IF(taylor_rst)THEN
    WRITE(pnum,'(I2.2)')oft_hcurl%order
    WRITE(mnum,'(I2.2)')i
    filename='hvac_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    IF(.NOT.oft_file_exist(filename))THEN
      CALL oft_mpi_barrier(ierr)
      CALL oft_h1%vec_save(u,filename,'hvac')
    END IF
  END IF
  !---Compute field energy
  CALL mop%apply(u,tmp)
  venergy=u%dot(tmp)
  CALL u%scale(1.d0/venergy)
  IF(oft_env%head_proc)THEN
    WRITE(*,*)'Injector      = ',taylor_htag(i)
    WRITE(*,*)'Vacuum Energy = ',1.d0/venergy
  END IF
  IF(PRESENT(energy))energy(i)=1.d0/venergy
!---------------------------------------------------------------------------
! Compute current field
!---------------------------------------------------------------------------
  !---Use vacuum field as source term
  call tmp%scale(1.d0/venergy)
  call oft_hcurl_create(taylor_hcur(i,oft_hcurl_level)%f)
  u=>taylor_hcur(i,oft_hcurl_level)%f
  rst=.FALSE.
  IF(taylor_rst)THEN
    WRITE(pnum,'(I2.2)')oft_hcurl%order
    WRITE(mnum,'(I2.2)')i
    filename='hcur_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    IF(oft_file_exist(filename))THEN
      CALL oft_hcurl%vec_load(u,filename,'hcur')
      rst=.TRUE.
    END IF
  END IF
  IF(.NOT.rst)THEN
    !---Copy H1(Curl) subpace into new vector
    CALL b%set(0.d0)
    NULLIFY(vals)
    CALL tmp%get_slice(vals,1)
    CALL b%restore_slice(vals)
    !---Compute current field
    hcurl_divout%pm=.FALSE.
    hcurl_divout%mop=>mop_hcurl
    CALL hcurl_zerob(b)
    CALL winv%apply(u,b)
  END IF
  !---Clean divergence
  NULLIFY(hcurl_divout%mop)
  hcurl_divout%app_freq=1
  hcurl_divout%pm=.TRUE.
  CALL hcurl_divout%apply(u)
  !---
  DO j=1,taylor_nm
    CALL wop%apply(taylor_hffa(j,oft_hcurl_level)%f,b)
    venergy = u%dot(b)
    IF(oft_env%head_proc)WRITE(*,'(A,I3,A,E10.3)')'Mode ',j,' Coupling = ',venergy
  END DO
!---------------------------------------------------------------------------
! Write restart file
!---------------------------------------------------------------------------
  IF(taylor_rst)THEN
    WRITE(pnum,'(I2.2)')oft_hcurl%order
    WRITE(mnum,'(I2.2)')i
    filename='hcur_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    IF(.NOT.oft_file_exist(filename))THEN
      CALL oft_mpi_barrier(ierr)
      CALL oft_hcurl%vec_save(u,filename,'hcur')
    END IF
  END IF
END DO
!---
CALL tmp%delete
CALL b%delete
NULLIFY(tmp,b)
!---
CALL mop%delete
CALL mop_hcurl%delete
CALL lop_lag%delete
IF(ASSOCIATED(ml_wop))THEN
  DO i=1,SIZE(ml_wop)
    CALL ml_wop(i)%M%delete
    DEALLOCATE(ml_wop(i)%M)
  END DO
  DEALLOCATE(ml_wop)
END IF
IF(ASSOCIATED(ml_lop))THEN
  DO i=1,SIZE(ml_lop)
    CALL ml_lop(i)%M%delete
    DEALLOCATE(ml_lop(i)%M)
  END DO
  DEALLOCATE(ml_lop)
END IF
CALL linv%pre%delete
DEALLOCATE(linv%pre)
CALL linv%delete
DEALLOCATE(linv)
CALL winv%pre%delete
DEALLOCATE(winv%pre)
CALL winv%delete
DEALLOCATE(winv)
CALL linv_lag%pre%delete
DEALLOCATE(linv_lag%pre)
CALL linv_lag%delete
DEALLOCATE(linv_lag)
CALL oft_lag_set_level(oft_lagrange_nlevels)
DEBUG_STACK_POP
end subroutine taylor_vacuum
!---------------------------------------------------------------------------
! SUBROUTINE: taylor_injectors
!---------------------------------------------------------------------------
!> Compute force-free plasma response to external fields generated by @ref
!! taylor::taylor_injectors "taylor_injectors"
!!
!! @note When `taylor_rst=.TRUE.` one restart files will be generated for
!! each jump plane `gffa_*.rst`
!!
!! @param[in] lambda Desired lambda for force-free state
!---------------------------------------------------------------------------
subroutine taylor_injectors(lambda)
real(r8), intent(in) :: lambda
CLASS(oft_matrix), POINTER :: lop_lag => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_jmlb => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_wop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: kop => NULL()
CLASS(oft_matrix), POINTER :: wop => NULL()
CLASS(oft_matrix), POINTER :: jmlb_mat => NULL()
!---JMLB solver
TYPE(oft_native_gmres_solver), TARGET :: jmlb_inv
!---H1(Curl) Divergence cleaner
CLASS(oft_solver), POINTER :: linv_lag => NULL()
TYPE(oft_hcurl_divout), TARGET :: hcurl_divout
!---
class(oft_vector), pointer :: u,b,tmp
type(oft_hcurl_orthog), target :: orthog
real(r8) :: venergy,lam_file
integer(i4) :: i,k,ierr
logical :: do_orthog,rst
character(2) :: pnum,mnum
character(40) :: filename
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(taylor_hcpc))CALL oft_abort("Vacuum fields not available", &
"taylor_injectors", __FILE__)
IF(taylor_rst)THEN
  rst=.TRUE.
  DO i=1,taylor_nh
    WRITE(pnum,'(I2.2)')oft_hcurl%order
    WRITE(mnum,'(I2.2)')i
    filename='gffa_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    rst=rst.AND.oft_file_exist(filename)
  END DO
ELSE
  rst=.FALSE.
END IF
IF(taylor_minlev<0)taylor_minlev=oft_hcurl_nlevels
NULLIFY(mop,kop,wop,jmlb_mat,ml_jmlb)
IF(.NOT.rst)THEN
  !---Orthogonalize (if within 5% of Taylor state)
  do_orthog=.FALSE.
  IF(taylor_nm>0)THEN
    IF(ABS((lambda-taylor_hlam(1,oft_hcurl_level))/taylor_hlam(1,oft_hcurl_level))<5.d-2)THEN
      CALL oft_hcurl_getwop(wop,'zerob')
      orthog%orthog=>taylor_hffa
      orthog%nm=1
      orthog%wop=>wop
      do_orthog=.TRUE.
    END IF
  END IF
!---------------------------------------------------------------------------
! Setup H1(Curl)::WOP preconditioner
!---------------------------------------------------------------------------
  IF(taylor_minlev==oft_hcurl_nlevels)THEN ! Lowest level uses diag precond
    CALL oft_hcurl_getjmlb(jmlb_mat,lambda,'zerob')
    CALL create_diag_pre(jmlb_inv%pre)
  ELSE ! Nested levels use MG
    CALL hcurl_getjmlb_pre(jmlb_inv%pre,ml_jmlb,lambda,nlevels=oft_hcurl_nlevels-taylor_minlev+1)
    jmlb_mat=>ml_jmlb(oft_hcurl_nlevels-taylor_minlev+1)%M
  END IF
!---------------------------------------------------------------------------
! Create H1(Curl)::WOP solver
!---------------------------------------------------------------------------
  jmlb_inv%A=>jmlb_mat
  jmlb_inv%atol=1.d-8
  jmlb_inv%its=-3
  jmlb_inv%nrits=20
  jmlb_inv%itplot=1
  !jmlb_inv%bc=>hcurl_zerob
END IF
!---------------------------------------------------------------------------
! Create H1(Curl) divergence cleaner
!---------------------------------------------------------------------------
CALL oft_lag_set_level(MIN(oft_hcurl_level,oft_lagrange_lin_level))
CALL oft_lag_getlop(lop_lag,"zerob")
CALL create_cg_solver(linv_lag)
linv_lag%A=>lop_lag
linv_lag%its=-2
CALL create_diag_pre(linv_lag%pre)
linv_lag%its=40
hcurl_divout%solver=>linv_lag
hcurl_divout%bc=>lag_zerob
hcurl_divout%app_freq=10
CALL oft_hcurl_getmop(mop,'zerob')
hcurl_divout%mop=>mop
!jmlb_inv%cleaner=>hcurl_divout
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
call oft_h1_set_level(oft_h1_nlevels)
!---
call oft_hcurl_create(tmp)
allocate(taylor_gffa(taylor_nh,oft_h1_nlevels))
IF(.NOT.rst)CALL oft_hcurl_getkop(kop,'zerob')
!---
do i=1,taylor_nh
  !---
  CALL oft_hcurl_create(taylor_gffa(i,oft_h1_level)%f)
  u=>taylor_gffa(i,oft_h1_level)%f
  rst=.FALSE.
  IF(taylor_rst)THEN
    WRITE(pnum,'(I2.2)')oft_hcurl%order
    WRITE(mnum,'(I2.2)')i
    filename='gffa_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    IF(oft_file_exist(filename))THEN
      CALL hdf5_read(lam_file,filename,'lambda')
      IF(ABS(lam_file-lambda)<1.d-5)THEN
        CALL oft_hcurl%vec_load(u,filename,'gffa')
        rst=.TRUE.
      END IF
    END IF
  END IF
  IF(.NOT.rst)THEN
    CALL u%add(0.d0,1.d0,taylor_hcur(i,oft_h1_level)%f)
    IF(do_orthog)CALL orthog%apply(u)
    CALL kop%apply(u,tmp)
    !---Solve
    b=>tmp
    CALL hcurl_zerob(b)
    CALL u%set(0.d0)
    CALL jmlb_inv%apply(u,b)
  END IF
  !---Clean divergence
  NULLIFY(hcurl_divout%mop)
  hcurl_divout%app_freq=1
  hcurl_divout%pm=.TRUE.
  CALL hcurl_divout%apply(u)
!---------------------------------------------------------------------------
! Write restart file
!---------------------------------------------------------------------------
  IF(taylor_rst)THEN
    WRITE(pnum,'(I2.2)')oft_hcurl%order
    WRITE(mnum,'(I2.2)')i
    filename='gffa_r'//mg_mesh%rlevel//'_p'//pnum//'_h'//mnum//'.rst'
    IF(.NOT.oft_file_exist(filename))THEN
      CALL oft_mpi_barrier(ierr)
      CALL oft_hcurl%vec_save(u,filename,'gffa')
      IF(oft_env%head_proc)CALL hdf5_write(lambda,filename,'lambda')
    END IF
  END IF
  !---
  CALL u%add(lambda**2,lambda,taylor_hcur(i,oft_h1_level)%f)
end do
!---
CALL tmp%delete
NULLIFY(tmp,b)
!---
CALL mop%delete
IF(ASSOCIATED(kop))THEN
  CALL kop%delete
  DEALLOCATE(kop)
END IF
CALL lop_lag%delete
IF(ASSOCIATED(ml_jmlb))THEN
  DO i=1,SIZE(ml_jmlb)
    CALL ml_jmlb(i)%M%delete
    DEALLOCATE(ml_jmlb(i)%M)
  END DO
  DEALLOCATE(ml_jmlb)
ELSE
  IF(ASSOCIATED(ml_jmlb))THEN
    CALL jmlb_mat%delete
    DEALLOCATE(jmlb_mat)
  END IF
  IF(ASSOCIATED(wop))THEN
    CALL wop%delete
    DEALLOCATE(wop)
  END IF
END IF
!---
IF(.NOT.rst)CALL jmlb_inv%pre%delete
CALL oft_lag_set_level(oft_lagrange_nlevels)
DEBUG_STACK_POP
end subroutine taylor_injectors
!---------------------------------------------------------------------------
! SUBROUTINE: taylor_injector_single
!---------------------------------------------------------------------------
!> Compute force-free plasma response to external fields generated by @ref
!! taylor::taylor_vacuum "taylor_vacuum"
!!
!! @param[in] lambda Desired lambda for force-free state
!! @param[in] fluxes Flux for each injector
!! @param[in] gffa Plasma component (non-vacuum) of injector field
!---------------------------------------------------------------------------
subroutine taylor_injector_single(lambda,fluxes,gffa)
real(r8), intent(in) :: lambda
real(r8), intent(in) :: fluxes(:)
class(oft_vector), pointer, intent(inout) :: gffa
CLASS(oft_matrix), POINTER :: lop_lag => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_jmlb => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: kop => NULL()
CLASS(oft_matrix), POINTER :: wop => NULL()
CLASS(oft_matrix), POINTER :: jmlb_mat => NULL()
!---JMLB solver
TYPE(oft_native_gmres_solver), TARGET :: jmlb_inv
!---H1(Curl) Divergence cleaner
CLASS(oft_solver), POINTER :: linv_lag => NULL()
TYPE(oft_hcurl_divout), TARGET :: hcurl_divout
!---
class(oft_vector), pointer :: tmp
type(oft_hcurl_orthog), target :: orthog
real(r8) :: venergy,lam_file
integer(i4) :: i,k
character(2) :: pnum,mnum
character(40) :: filename
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(taylor_hcpc))CALL oft_abort("Vacuum fields not available", &
"taylor_injector_single", __FILE__)
NULLIFY(mop,kop,wop,ml_jmlb)
IF(taylor_minlev<0)taylor_minlev=oft_hcurl_nlevels
!---------------------------------------------------------------------------
! Setup H1(Curl)::WOP preconditioner
!---------------------------------------------------------------------------
IF(taylor_minlev==oft_hcurl_nlevels)THEN ! Lowest level uses diag precond
  CALL oft_hcurl_getjmlb(jmlb_mat,lambda,'zerob')
  CALL create_diag_pre(jmlb_inv%pre)
ELSE ! Nested levels use MG
  CALL hcurl_getjmlb_pre(jmlb_inv%pre,ml_jmlb,lambda,nlevels=1)
  jmlb_mat=>ml_jmlb(1)%M
END IF
!---------------------------------------------------------------------------
! Create H1(Curl)::WOP solver
!---------------------------------------------------------------------------
jmlb_inv%A=>jmlb_mat
jmlb_inv%atol=1.d-7
jmlb_inv%its=-3
jmlb_inv%nrits=20
jmlb_inv%itplot=1
!jmlb_inv%bc=>hcurl_zerob
!---------------------------------------------------------------------------
! Create H1(Curl) divergence cleaner
!---------------------------------------------------------------------------
CALL oft_lag_set_level(MIN(oft_hcurl_level,oft_lagrange_lin_level))
CALL oft_lag_getlop(lop_lag,"zerob")
CALL create_cg_solver(linv_lag)
linv_lag%A=>lop_lag
linv_lag%its=-2
CALL create_diag_pre(linv_lag%pre)
linv_lag%its=40
hcurl_divout%solver=>linv_lag
hcurl_divout%bc=>lag_zerob
hcurl_divout%app_freq=10
CALL oft_hcurl_getmop(mop,'zerob')
hcurl_divout%mop=>mop
!jmlb_inv%cleaner=>hcurl_divout
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
call oft_h1_set_level(oft_h1_nlevels)
!---
call oft_hcurl_create(gffa)
call oft_hcurl_create(tmp)
!---
CALL gffa%set(0.d0)
do i=1,taylor_nh
  CALL gffa%add(1.d0,fluxes(i),taylor_hcur(i,oft_h1_level)%f)
end do
!---Orthogonalize (if within 5% of Taylor state)
IF(taylor_nm>0)THEN
  IF(ABS((lambda-taylor_hlam(1,oft_hcurl_level))/taylor_hlam(1,oft_hcurl_level))<5.d-2)THEN
    CALL oft_hcurl_getwop(wop,'zerob')
    orthog%orthog=>taylor_hffa
    orthog%nm=1
    orthog%wop=>wop
    CALL orthog%apply(gffa)
    CALL wop%delete
    DEALLOCATE(wop)
  END IF
END IF
!---Get H1(Curl) helicity matrix
CALL oft_hcurl_getkop(kop,'zerob')
CALL kop%apply(gffa,tmp)
CALL kop%delete
DEALLOCATE(kop)
!---Solve
CALL hcurl_zerob(tmp)
CALL gffa%set(0.d0)
CALL jmlb_inv%apply(gffa,tmp)
!---Clean divergence
NULLIFY(hcurl_divout%mop)
hcurl_divout%app_freq=1
hcurl_divout%pm=.TRUE.
CALL hcurl_divout%apply(gffa)
!---
CALL gffa%scale(lambda**2)
do i=1,taylor_nh
  CALL gffa%add(1.d0,fluxes(i)*lambda,taylor_hcur(i,oft_h1_level)%f)
end do
!---
CALL tmp%delete
DEALLOCATE(tmp)
!---
CALL lop_lag%delete
IF(ASSOCIATED(ml_jmlb))THEN
  DO i=1,SIZE(ml_jmlb)
    CALL ml_jmlb(i)%M%delete
    DEALLOCATE(ml_jmlb(i)%M)
  END DO
  DEALLOCATE(ml_jmlb)
ELSE
  CALL jmlb_mat%delete
  DEALLOCATE(jmlb_mat)
END IF
!---
CALL jmlb_inv%pre%delete
CALL oft_lag_set_level(oft_lagrange_nlevels)
DEBUG_STACK_POP
end subroutine taylor_injector_single
!---------------------------------------------------------------------------
! SUBROUTINE: taylor_rinterp_setup
!---------------------------------------------------------------------------
!> Setup interpolator for composite Taylor state fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_taylor_rinterp or children
!---------------------------------------------------------------------------
subroutine taylor_rinterp_setup(self)
class(oft_taylor_rinterp), intent(inout) :: self
DEBUG_STACK_PUSH
!---Get local slice
CALL self%ua%get_local(self%acurl)
CALL self%uvac%get_local(self%vac_curl,1)
CALL self%uvac%get_local(self%vac_grad,2)
self%hgrad_rep=>oft_h0
self%hcurl_rep=>oft_hcurl
DEBUG_STACK_POP
end subroutine taylor_rinterp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: taylor_rinterp
!---------------------------------------------------------------------------
!> Reconstruct a composite Taylor state field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [3]
!---------------------------------------------------------------------------
subroutine taylor_rinterp(self,cell,f,gop,val)
class(oft_taylor_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: cgop(3,6)
real(r8), allocatable :: rop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vac_grad))CALL oft_abort('Setup has not been called!','taylor_rinterp',__FILE__)
val=0.d0
!---Get curl dofs
allocate(j(self%hcurl_rep%nce),rop(3,self%hcurl_rep%nce))
call self%hcurl_rep%ncdofs(cell,j) ! get DOFs
CALL oft_hcurl_get_cgops(gop,cgop)
!---Reconstruct field
call oft_hcurl_eval_all(self%hcurl_rep,cell,f,rop,gop)
do jc=1,self%hcurl_rep%nce
  val=val+self%vac_curl(j(jc))*rop(:,jc)
end do
call oft_hcurl_ceval_all(self%hcurl_rep,cell,f,rop,cgop)
do jc=1,self%hcurl_rep%nce
  val=val+self%acurl(j(jc))*rop(:,jc)
end do
deallocate(j,rop)
!---Get curl dofs
allocate(j(self%hgrad_rep%nce),rop(3,self%hgrad_rep%nce))
call self%hgrad_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
call oft_h0_geval_all(self%hgrad_rep,cell,f,rop,gop)
do jc=1,self%hgrad_rep%nce
  val=val+self%vac_grad(j(jc))*rop(:,jc)
end do
deallocate(j,rop)
DEBUG_STACK_POP
end subroutine taylor_rinterp
end module taylor
