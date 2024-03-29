!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_hcurl.F90
!
!> Regression tests for vector H1(Curl) finite elements. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
program test_hcurl
USE oft_base
USE oft_mesh_type, ONLY: mesh
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: mg_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_set_level, oft_hcurl_nlevels
USE oft_hcurl_fields, ONLY: oft_hcurl_create
USE oft_hcurl_operators, ONLY: oft_hcurl_getkop, oft_hcurl_getwop, hcurl_zerob, &
  hcurl_setup_interp, hcurl_mloptions, hcurl_getwop_pre
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
IMPLICIT NONE
INTEGER(i4) :: minlev
INTEGER(i4) :: order,ierr,io_unit
LOGICAL :: mg_test
NAMELIST/test_hcurl_options/order,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_hcurl_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct
IF(mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
minlev=2
IF(mesh%type==3)minlev=mg_mesh%mgmax
CALL oft_hcurl_setup(order,minlev)
IF(mg_test)THEN
  CALL hcurl_setup_interp
  CALL hcurl_mloptions
END IF
!---Run tests
oft_env%pm=.FALSE.
IF(mg_test)THEN
  CALL test_wopmg
ELSE
  CALL test_wop
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE: test_wop
!------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H1(Curl) elements.
!------------------------------------------------------------------------------
SUBROUTINE test_wop
!---Create solver objects
CLASS(oft_solver), POINTER :: winv => NULL()
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: wop => NULL()
CLASS(oft_matrix), POINTER :: kop => NULL()
!---Set FE level
CALL oft_hcurl_set_level(oft_hcurl_nlevels)
!---Create solver fields
CALL oft_hcurl_create(u)
CALL oft_hcurl_create(v)
!---Get FE operators
CALL oft_hcurl_getwop(wop,'zerob')
CALL oft_hcurl_getkop(kop,'none')
!---Setup matrix solver
CALL create_cg_solver(winv)
winv%A=>wop
winv%its=-3
CALL create_diag_pre(winv%pre)
!---Solve
CALL u%set(1.d0)
CALL kop%apply(u,v)
CALL hcurl_zerob(v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
CALL kop%apply(u,v)
uu=u%dot(v)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrices
CALL kop%delete
CALL wop%delete
DEALLOCATE(kop,wop)
!---Destory preconditioner
CALL winv%pre%delete
!---Destory solver
CALL winv%delete
END SUBROUTINE test_wop
!------------------------------------------------------------------------------
! SUBROUTINE: test_wopmg
!------------------------------------------------------------------------------
!> Same as \ref test_hcurl::test_wop "test_wop" but use MG preconditioning.
!------------------------------------------------------------------------------
SUBROUTINE test_wopmg
!---Solver object
CLASS(oft_solver), POINTER :: winv => NULL()
!---Local variables
REAL(r8) :: uu
INTEGER(i4) :: i,nlevels
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: wop => NULL()
CLASS(oft_matrix), POINTER :: kop => NULL()
TYPE(oft_matrix_ptr), POINTER :: ml_wop(:) => NULL()
!---------------------------------------------------------------------------
! Create ML Matrices
!---------------------------------------------------------------------------
nlevels=oft_hcurl_nlevels-minlev+1
!---Create solver fields
CALL oft_hcurl_create(u)
CALL oft_hcurl_create(v)
!---Get FE operators
CALL oft_hcurl_getkop(kop,'none')
!---------------------------------------------------------------------------
! Setup matrix solver
!---------------------------------------------------------------------------
CALL create_cg_solver(winv,force_native=.TRUE.)
winv%its=-3
!---Setup MG preconditioner
CALL hcurl_getwop_pre(winv%pre,ml_wop,nlevels=nlevels)
wop=>ml_wop(nlevels)%M
winv%A=>wop
!---------------------------------------------------------------------------
! Solve system
!---------------------------------------------------------------------------
CALL u%set(1.d0)
CALL kop%apply(u,v)
CALL hcurl_zerob(v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
CALL kop%apply(u,v)
uu=u%dot(v)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrices
CALL wop%delete
CALL kop%delete
DEALLOCATE(wop,kop)
!---Destory preconditioner
CALL winv%pre%delete
DEALLOCATE(winv%pre)
!---Destory solver
CALL winv%delete
END SUBROUTINE test_wopmg
end program test_hcurl
