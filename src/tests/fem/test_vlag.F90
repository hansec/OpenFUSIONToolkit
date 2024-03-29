!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_vlag.F90
!
!> Regression tests for vector Lagrange finite elements. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Invert the mass matrix
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_vlag
USE oft_base
USE oft_mesh_type, ONLY: mesh
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid_build, ONLY: multigrid_construct
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels, oft_lag_set_level, &
  oft_lagrange, oft_lagrange_ops
USE oft_lag_fields, ONLY: oft_lag_vcreate
USE oft_lag_operators, ONLY: oft_lag_vgetmop, lag_vinterp, lag_vinject, lag_vzerob, &
  df_lop, nu_lop, lag_setup_interp, lag_mloptions, oft_lag_getlop
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr, oft_graph_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_la_utils, ONLY: create_matrix, combine_matrices
USE oft_solver_utils, ONLY: create_mlpre, create_cg_solver, create_diag_pre
IMPLICIT NONE
INTEGER(i4), PARAMETER :: minlev=2
INTEGER(i4) :: order,ierr,io_unit
LOGICAL :: mg_test=.FALSE.
NAMELIST/test_lag_options/order,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_lag_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct
IF(mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
CALL oft_lag_setup(order,minlev)
IF(mg_test)THEN
  CALL lag_setup_interp(.TRUE.)
  CALL lag_mloptions
END IF
!---Run tests
oft_env%pm=.FALSE.
IF(mg_test)THEN
  CALL test_mopmg
ELSE
  CALL test_mop
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE: test_mop
!------------------------------------------------------------------------------
!> Invert the mass matrix and output required iterataions and final field energy.
!------------------------------------------------------------------------------
SUBROUTINE test_mop
!---Create solver objects
CLASS(oft_solver), POINTER :: minv => NULL()
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL oft_lag_set_level(oft_lagrange_nlevels)
!---Create solver fields
CALL oft_lag_vcreate(u)
CALL oft_lag_vcreate(v)
!---Get FE operators
CALL oft_lag_vgetmop(mop,'none')
!---Setup matrix solver
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
CALL create_diag_pre(minv%pre)
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='lagrange.results')
  WRITE(io_unit,*)minv%cits
  WRITE(io_unit,*)uu/u%ng
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrix
CALL mop%delete
DEALLOCATE(mop)
!---Destory preconditioner
CALL minv%pre%delete
!---Destory solver
CALL minv%delete
END SUBROUTINE test_mop
!------------------------------------------------------------------------------
! SUBROUTINE: test_mopmg
!------------------------------------------------------------------------------
!> Invert the mass matrix and output required iterataions and final field energy.
!------------------------------------------------------------------------------
SUBROUTINE test_mopmg
!---Create solver objects
CLASS(oft_solver), POINTER :: linv => NULL()
!--- ML structures for MG-preconditioner
type(oft_matrix_ptr), pointer :: ml_int(:) => NULL()
type(oft_matrix_ptr), pointer :: ml_vlop(:) => NULL()
type(oft_matrix_ptr), pointer :: ml_lop(:) => NULL()
INTEGER(i4), ALLOCATABLE :: levels(:)
REAL(r8), ALLOCATABLE :: df(:)
INTEGER(i4), ALLOCATABLE :: nu(:)
!---Local variables
CLASS(oft_vector), pointer :: fvec,cvec
TYPE(oft_graph_ptr) :: graphs(3,3)
TYPE(oft_matrix_ptr) :: mats(3,3)
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
REAL(r8) :: uu
INTEGER(i4) :: i,nlevels
!---------------------------------------------------------------------------
! Create ML Matrices
!---------------------------------------------------------------------------
nlevels=oft_lagrange_nlevels-minlev+1
!---Get FE operators
ALLOCATE(ml_int(nlevels-1),ml_vlop(nlevels),ml_lop(nlevels))
ALLOCATE(levels(nlevels),df(nlevels),nu(nlevels))
DO i=1,nlevels
  CALL oft_lag_set_level(minlev+i-1)
  levels(i)=minlev+(i-1)
  df(i)=df_lop(levels(i))
  nu(i)=nu_lop(levels(i))
  !---
  NULLIFY(ml_lop(i)%M)
  CALL oft_lag_getlop(ml_lop(i)%M,'zerob')
!---------------------------------------------------------------------------
! Create composite matrix
!---------------------------------------------------------------------------
  ALLOCATE(graphs(1,1)%g)
  graphs(1,1)%g%nr=oft_lagrange%ne
  graphs(1,1)%g%nrg=oft_lagrange%global%ne
  graphs(1,1)%g%nc=oft_lagrange%ne
  graphs(1,1)%g%ncg=oft_lagrange%global%ne
  graphs(1,1)%g%nnz=oft_lagrange%nee
  graphs(1,1)%g%kr=>oft_lagrange%kee
  graphs(1,1)%g%lc=>oft_lagrange%lee
  !---
  graphs(2,2)%g=>graphs(1,1)%g
  graphs(3,3)%g=>graphs(1,1)%g
  !---Get coarse and fine vectors
  CALL oft_lag_vcreate(cvec)
  CALL oft_lag_vcreate(fvec)
  NULLIFY(ml_vlop(i)%M)
  !---Construct matrix
  CALL create_matrix(ml_vlop(i)%M,graphs,fvec,cvec)
  DEALLOCATE(graphs(1,1)%g)
!---------------------------------------------------------------------------
! Combine child matrices into composite matrix
!---------------------------------------------------------------------------
  !---Specify child matrices
  mats(1,1)%m=>ml_lop(i)%M
  mats(2,2)%m=>ml_lop(i)%M
  mats(3,3)%m=>ml_lop(i)%M
  !---Combine matrices
  CALL combine_matrices(mats,3,3,ml_vlop(i)%M)
  CALL ml_vlop(i)%M%assemble(fvec)
  !---Delete temporaries
  CALL cvec%delete
  CALL fvec%delete
  DEALLOCATE(cvec,fvec)
  !---
  IF(i>1)ml_int(i-1)%M=>oft_lagrange_ops%vinterp
END DO
CALL oft_lag_set_level(oft_lagrange_nlevels)
!---Create solver fields
CALL oft_lag_vcreate(u)
CALL oft_lag_vcreate(v)
!---Get FE operators
CALL oft_lag_vgetmop(mop,'none')
!---Setup matrix solver
CALL create_cg_solver(linv,force_native=.TRUE.)
linv%A=>ml_vlop(nlevels)%M
linv%its=-3
linv%itplot=1
CALL create_mlpre(linv%pre,ml_vlop(1:nlevels),levels, &
  nlevels=nlevels,create_vec=oft_lag_vcreate, &
  interp=lag_vinterp,inject=lag_vinject, &
  bc=lag_vzerob,stype=1,df=df,nu=nu)
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_vzerob(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='lagrange.results')
  WRITE(io_unit,*)linv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrix
CALL mop%delete
DEALLOCATE(mop)
!---Destory preconditioner
CALL linv%pre%delete
DEALLOCATE(linv%pre)
!---Destory solver
CALL linv%delete
END SUBROUTINE test_mopmg
END PROGRAM test_vlag
