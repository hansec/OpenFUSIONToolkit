!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_mapping_jac.F90
!
!> Regression tests for tetmesh Jacobian transformations. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! -
!!
!! @authors Chris Hansen
!! @date August 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_mapping_jac
USE oft_base
USE oft_quadrature
USE oft_mesh_type, ONLY: mesh
USE oft_mesh_sphere, ONLY: mesh_sphere_id
USE multigrid_build, ONLY: multigrid_construct
!---
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange, oft_lag_npos, oft_lag_geval, &
oft_lag_d2eval
IMPLICIT NONE
INTEGER(i4) :: xi,xj,ierr,nfail,i,ntests,io_unit
INTEGER(i4) :: order
REAL(r8) :: check_vec(6),tol=1.d-6
NAMELIST/test_mapping_options/order
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_mapping_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct
IF(mesh%cad_type/=mesh_sphere_id)CALL oft_abort('Wrong mesh type, test for SPHERE only.','main',__FILE__)
IF(oft_env%nprocs>1)CALL oft_abort('Test is for serial meshes only.','main',__FILE__)
!---Setup FEM
CALL oft_lag_setup(order)
!---Run test cases
nfail=0
OPEN(NEWUNIT=io_unit,FILE='mapping_jac.tests')
READ(io_unit,*)ntests
DO i=1,ntests
  !---Evaluate errors
  READ(io_unit,'(2I12,6E25.17)')xi,xj,check_vec
  IF((xj==0).OR.(order/mesh%order>=2))nfail=nfail+check_jac2(xi,xj,check_vec,order,tol)
END DO
CLOSE(io_unit)
OPEN(NEWUNIT=io_unit,FILE='mapping_jac.results')
WRITE(io_unit,*)nfail
CLOSE(io_unit)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------
! FUNCTION check_jac2
!---------------------------------------------------------------------------
!> Validate 2nd order spatial derivatives by comparing the computed derivatives
!! to a known result. Currently, setup to test all quadratic functions in 3
!! dimensions.
!---------------------------------------------------------------------------
FUNCTION check_jac2(xi,xj,check_vec,order,tol) RESULT(fail_count)
INTEGER(i4), INTENT(in) :: xi,xj,order
REAL(r8), INTENT(in) :: check_vec(6),tol
!---
INTEGER(i4) :: i,j,k,fail_count
REAL(r8) :: f(4),rop(3),gop(3,4),v,val(6),vloc(6)
REAL(r8), ALLOCATABLE :: pt_loc(:,:),g2op(:,:),Kmat(:,:)
TYPE(oft_quad_type), POINTER :: quad
CHARACTER(LEN=1), PARAMETER :: coords(0:3)=(/'s','x','y','z'/)
WRITE(*,*)'Testing ',coords(xi),coords(xj)
!---
fail_count=0
quad=>oft_lagrange%quad
!$omp parallel private(k,j,f,pt_loc,v,vloc,val,rop,gop,g2op,Kmat) reduction(+:fail_count)
ALLOCATE(pt_loc(0:3,oft_lagrange%nce))
IF(mesh%type==3)THEN
  ALLOCATE(g2op(6,6),Kmat(6,3))
ELSE
  ALLOCATE(g2op(6,10),Kmat(10,3))
END IF
pt_loc=1.d0
!$omp do
DO i=1,mesh%nc
  DO k=1,oft_lagrange%nce
    CALL oft_lag_npos(oft_lagrange,i,k,f)
    pt_loc(1:3,k)=mesh%log2phys(i,f)
  END DO
  !---
  DO j=1,quad%np ! Loop over quadrature points
    !---
    CALL mesh%jacobian(i,quad%pts(:,j),gop,v)
    CALL mesh%hessian(i,quad%pts(:,j),g2op,Kmat)
    vloc=0.d0
    DO k=1,oft_lagrange%nce
      CALL oft_lag_geval(oft_lagrange,i,k,quad%pts(:,j),rop,gop)
      CALL oft_lag_d2eval(oft_lagrange,i,k,quad%pts(:,j),val,g2op)
      vloc=vloc+pt_loc(xi,k)*pt_loc(xj,k)*(val-MATMUL(g2op,MATMUL(Kmat,rop)))
    END DO
    IF(SQRT(SUM((vloc-check_vec)**2))>tol)fail_count=fail_count+1
 END DO
END DO
DEALLOCATE(pt_loc,g2op,Kmat)
!$omp end parallel
END FUNCTION check_jac2
END PROGRAM test_mapping_jac
