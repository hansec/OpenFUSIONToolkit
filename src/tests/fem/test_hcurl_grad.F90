!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_hcurl_grad.F90
!
!> Regression tests for vector H(Curl) + Grad(H^1) finite elements. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the equation \f$ M B = M \hat{I} \f$
!! - Check boundary condition handling in the mass matrix (`bc_test`)
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
program test_hcurl_grad
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
! USE oft_io
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_native_mlpre, create_cg_solver, create_diag_pre
!---
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
!---Grad(H^1) FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_getmop, hcurl_grad_setup_interp, hcurl_grad_getmop_pre, hcurl_grad_mloptions, &
  oft_hcurl_grad_rinterp, oft_hcurl_grad_gzerop, oft_hcurl_grad_zerob
IMPLICIT NONE
INTEGER(i4) :: order,ierr,io_unit,minlev
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_h1,ML_h1grad
TYPE(oft_ml_fem_type), TARGET :: ML_oft_hcurl
TYPE(oft_ml_fem_comp_type), TARGET :: ML_hcurl_grad
TYPE(oft_hcurl_grad_gzerop), TARGET :: hcurl_grad_gzerop
TYPE(oft_hcurl_grad_zerob), TARGET :: hcurl_grad_zerob
LOGICAL :: mg_test = .FALSE.
LOGICAL :: bc_test = .FALSE.
NAMELIST/test_hcurl_grad_options/order,mg_test,bc_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_hcurl_grad_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!------------------------------------------------------------------------------
! Build FE structures
!------------------------------------------------------------------------------
minlev=2
IF(oft_env%nprocs>1)minlev=mg_mesh%nbase+1
IF(mg_mesh%mesh%type==3)minlev=mg_mesh%mgdim
!--- Grad(H^1) subspace
CALL oft_h1_setup(mg_mesh,order+1,ML_oft_h1,minlev=minlev)
IF(mg_test)CALL h1_setup_interp(ML_oft_h1)
!--- H(Curl) subspace
CALL oft_hcurl_setup(mg_mesh,order,ML_oft_hcurl,minlev=minlev)
IF(mg_test)CALL hcurl_setup_interp(ML_oft_hcurl)
!--- Full H(Curl) space
CALL oft_hcurl_grad_setup(ML_oft_hcurl,ML_oft_h1,ML_hcurl_grad,ML_h1grad,minlev=minlev)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>ML_hcurl_grad
hcurl_grad_zerob%ML_hcurl_grad_rep=>ML_hcurl_grad
IF(mg_test)THEN
  CALL hcurl_grad_setup_interp(ML_hcurl_grad,ML_oft_h1,create_full=.TRUE.)
  CALL hcurl_grad_mloptions
END IF
!---Run tests
oft_env%pm=.FALSE.
IF(bc_test)THEN
  CALL test_mop_bc
ELSE IF(mg_test)THEN
  CALL test_mopmg
ELSE
  CALL test_mop
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H(Curl) elements.
!---------------------------------------------------------------------------------
SUBROUTINE test_mop
!---Create solver objects
CLASS(oft_solver), POINTER :: winv => NULL()
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_hcurl_grad%vec_create(u)
CALL ML_hcurl_grad%vec_create(v)
!---Get FE operators
CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop,'none')
!---Setup matrix solver
CALL create_cg_solver(winv)
winv%A=>mop
winv%its=-3
winv%rtol=1.d-10
CALL create_diag_pre(winv%pre)
!---Solve
CALL u%set(1.d0)
CALL hcurl_grad_gzerop%apply(u)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu/(u%ng-mg_mesh%mesh%global%np)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
END SUBROUTINE test_mop
!---------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H(Curl) elements.
!---------------------------------------------------------------------------------
SUBROUTINE test_mopmg
!---Create solver objects
! TYPE(oft_native_gmres_solver) :: winv
CLASS(oft_solver), POINTER :: winv => NULL()
CLASS(oft_solver), POINTER :: pre
TYPE(oft_matrix_ptr), POINTER :: mats(:) => NULL()
!---Local variables
INTEGER(i4) :: i,nlevels
REAL(r8) :: uu,r,rsq
REAL(r8), POINTER :: bvout(:,:)
REAL(r8), POINTER :: vals(:) => NULL()
CLASS(oft_vector), POINTER :: u,v,x,y
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(oft_hcurl_grad_rinterp) :: field
!---Set FE level
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_hcurl_grad%vec_create(u)
CALL ML_hcurl_grad%vec_create(v)
!---Setup matrix solver
nlevels=ML_hcurl_grad%nlevels-ML_hcurl_grad%minlev
CALL create_cg_solver(winv,force_native=.TRUE.)
winv%its=-3
winv%rtol=1.d-10
CALL hcurl_grad_getmop_pre(ML_hcurl_grad,winv%pre,mats,nlevels=nlevels)
mop=>mats(nlevels)%M
winv%A=>mop
!---Solve
CALL u%set(1.d0)
CALL hcurl_grad_gzerop%apply(u)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu/(u%ng-mg_mesh%mesh%global%np)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
END SUBROUTINE test_mopmg
!---------------------------------------------------------------------------------
! Check boundary condition handling in the mass matrix: for bc = "none" and "zerob"
! eliminated rows must be identity rows and all other rows must match the
! unmodified matrix. The maximum deviation for each case is written to file.
!---------------------------------------------------------------------------------
SUBROUTINE test_mop_bc
INTEGER(i4) :: i,k
REAL(r8) :: errs(2,2)
REAL(r8), POINTER, DIMENSION(:) :: uv,mv,v0v,v1v
CLASS(oft_vector), POINTER :: u,mask,v0,v1
CLASS(oft_matrix), POINTER :: mop0,mop1
CHARACTER(LEN=5), PARAMETER :: bcs(2) = (/'none ','zerob'/)
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
CALL ML_hcurl_grad%vec_create(u)
CALL ML_hcurl_grad%vec_create(mask)
CALL ML_hcurl_grad%vec_create(v0)
CALL ML_hcurl_grad%vec_create(v1)
!---Reference matrix without boundary conditions
NULLIFY(mop0,uv,mv,v0v,v1v)
CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop0,'full')
DO k=1,2
  NULLIFY(mop1)
  CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop1,TRIM(bcs(k)))
  !---Mask of retained (1) and eliminated (0) rows
  CALL mask%set(1.d0)
  CALL hcurl_grad_gzerop%apply(mask)
  IF(k==2)CALL hcurl_grad_zerob%apply(mask)
  !---Compare matrix-vector products for a non-trivial vector
  CALL u%get_local(uv)
  DO i=1,u%n
    uv(i)=SIN(1.3d0*i+oft_env%rank)+0.5d0*COS(0.37d0*i*i)
  END DO
  CALL u%restore_local(uv)
  CALL mop0%apply(u,v0)
  CALL mop1%apply(u,v1)
  CALL u%get_local(uv)
  CALL mask%get_local(mv)
  CALL v0%get_local(v0v)
  CALL v1%get_local(v1v)
  errs(:,k)=0.d0
  DO i=1,u%n
    IF(mv(i)>0.5d0)THEN
      errs(2,k)=MAX(errs(2,k),ABS(v1v(i)-v0v(i)))
    ELSE
      errs(1,k)=MAX(errs(1,k),ABS(v1v(i)-uv(i)))
    END IF
  END DO
  CALL mop1%delete
  DEALLOCATE(mop1)
END DO
errs(:,1)=oft_mpi_max(errs(:,1),2)
errs(:,2)=oft_mpi_max(errs(:,2),2)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad_bc.results')
  WRITE(io_unit,*)MAXVAL(errs(:,1))
  WRITE(io_unit,*)MAXVAL(errs(:,2))
  CLOSE(io_unit)
END IF
!---Cleanup
CALL mop0%delete
DEALLOCATE(mop0,uv,mv,v0v,v1v)
CALL u%delete
CALL mask%delete
CALL v0%delete
CALL v1%delete
DEALLOCATE(u,mask,v0,v1)
END SUBROUTINE test_mop_bc
end program test_hcurl_grad
