!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_solver_xml.F90
!
!> Regression tests for scalar Lagrange finite elements. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$
!! - Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$, with MG
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_solver_xml
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast
USE oft_la_utils, ONLY: combine_matrices
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_solver_xml
!
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions, &
  oft_lag_zerob, df_lop, nu_lop, oft_lag_getlop, oft_lag_getmop, lag_getlop_pre
IMPLICIT NONE
INTEGER(i4), PARAMETER :: order=4
INTEGER(i4) :: io_unit,ierr
INTEGER(i4) :: test_type = 1
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
NAMELIST/test_xml_options/test_type
!---Initialize enviroment
CALL oft_init
!---Read optional test options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_xml_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
IF(test_type==2)THEN
  CALL oft_lag_setup(mg_mesh,2,ML_oft_lagrange,minlev=-1)
  lag_zerob%ML_lag_rep=>ML_oft_lagrange
  oft_env%pm=.FALSE.
  CALL test_comp_bj
ELSE
  CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
  lag_zerob%ML_lag_rep=>ML_oft_lagrange
  oft_env%pm=.FALSE.
  CALL test_lap
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
!> Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$ and output
!! required iterataions and final field energy.
!---------------------------------------------------------------------------------
SUBROUTINE test_lap
!---Create solver objects
CLASS(oft_solver), POINTER :: linv
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---
integer(i4) :: nnodes
TYPE(xml_node) :: solver_node
!---Set FE level
CALL ML_oft_lagrange%set_level(ML_oft_lagrange%nlevels)
!---Create solver fields
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
!---Get FE operators
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
!---Setup matrix solver
CALL xml_get_element(oft_env%xml,"solver",solver_node,ierr)
IF(ierr==0)THEN
  CALL create_solver_xml(linv,solver_node)
ELSE
  CALL oft_abort('Could not find XML node.','test_lap',__FILE__)
END IF
linv%A=>lop
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_zerob%apply(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='xml.results')
  WRITE(io_unit,*)linv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrices
CALL lop%delete
CALL mop%delete
DEALLOCATE(lop,mop)
!---Destory preconditioner
IF(ASSOCIATED(linv%pre))CALL linv%pre%delete
!---Destory solver
CALL linv%delete
END SUBROUTINE test_lap
!---------------------------------------------------------------------------------
! Test 3-field coupled block matrix with XML block-Jacobi slice groups
!---------------------------------------------------------------------------------
SUBROUTINE test_comp_bj
CLASS(oft_solver), POINTER :: linv
REAL(r8) :: uu
INTEGER(i4) :: i
CLASS(oft_vector), POINTER :: u,v,u1,v1
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: lop12 => NULL()
CLASS(oft_matrix), POINTER :: lop23 => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: cmat => NULL()
CLASS(oft_native_matrix), POINTER :: nmat
TYPE(oft_fem_comp_type) :: comp_rep
INTEGER(i4) :: mask(3,3)
TYPE(oft_matrix_ptr) :: mats(3,3)
TYPE(xml_node) :: solver_node
REAL(r8), POINTER :: vtmp(:)
!---Set FE level
CALL ML_oft_lagrange%set_level(ML_oft_lagrange%nlevels)
!---Setup 3-field composite FE representation
comp_rep%nfields=3
ALLOCATE(comp_rep%fields(3),comp_rep%field_tags(3))
DO i=1,3
  comp_rep%fields(i)%fe=>ML_oft_lagrange%current_level
END DO
comp_rep%field_tags=(/'u1','u2','u3'/)
!---Build single-field operators
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop12,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,lop23,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
IF(native_matrix_cast(nmat,lop12))nmat%M=0.8d0*nmat%M
IF(native_matrix_cast(nmat,lop23))nmat%M=5.d0*nmat%M
!---Build 3x3 coupled composite matrix: strong (1,2) coupling, moderate (2,3) coupling
mask=0
DO i=1,3
  mask(i,i)=1
  mats(i,i)%m=>lop
END DO
mask(1,2)=1; mask(2,1)=1
mask(2,3)=1; mask(3,2)=1
mats(1,2)%m=>lop12
mats(2,1)%m=>lop12
mats(2,3)%m=>lop23
mats(3,2)%m=>lop23
CALL comp_rep%vec_create(u)
CALL comp_rep%vec_create(v)
CALL comp_rep%mat_create(cmat,mask=mask)
CALL combine_matrices(mats,3,3,cmat)
CALL cmat%assemble
!---Setup RHS
CALL ML_oft_lagrange%vec_create(u1)
CALL ML_oft_lagrange%vec_create(v1)
CALL u1%set(1.d0)
CALL mop%apply(u1,v1)
CALL lag_zerob%apply(v1)
NULLIFY(vtmp)
CALL v1%get_local(vtmp)
CALL v%restore_local(vtmp,iblock=1)
CALL v%restore_local(vtmp,iblock=2)
CALL v%restore_local(vtmp,iblock=3)
DEALLOCATE(vtmp)
!---Setup XML solver
CALL xml_get_element(oft_env%xml,"solver",solver_node,ierr)
IF(ierr==0)THEN
  CALL create_solver_xml(linv,solver_node)
ELSE
  CALL oft_abort('Could not find XML node.','test_comp_bj',__FILE__)
END IF
linv%A=>cmat
!---Solve
CALL u%set(0.d0)
CALL linv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='xml.results')
  WRITE(io_unit,*)linv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Cleanup
CALL u%delete
CALL v%delete
CALL u1%delete
CALL v1%delete
DEALLOCATE(u,v,u1,v1)
CALL lop%delete
CALL lop12%delete
CALL lop23%delete
CALL mop%delete
CALL cmat%delete
DEALLOCATE(lop,lop12,lop23,mop,cmat)
IF(ASSOCIATED(linv%pre))CALL linv%pre%delete
CALL linv%delete
CALL comp_rep%delete
END SUBROUTINE test_comp_bj
END PROGRAM test_solver_xml
