!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_taylor.F90
!
!> Regression tests for computing Taylor states in a spherical geometry.
!! States are computed at different polynomial orders and the resulting
!! eigenvalues and toroidal currents are compared to reference cases.
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_taylor
USE oft_base
USE oft_mesh_type, ONLY: mesh
USE oft_mesh_sphere, ONLY: mesh_sphere_id
USE multigrid, ONLY: mg_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_nlevels
USE oft_hcurl_fields, ONLY: oft_hcurl_create
USE oft_hcurl_operators, ONLY: hcurl_setup_interp, hcurl_mloptions
USE taylor, ONLY: taylor_minlev, taylor_hmodes, taylor_hlam, &
  taylor_htor
implicit none
INTEGER(i4) :: order=1,nm=1,ierr,io_unit
LOGICAL :: mg_test=.FALSE.
NAMELIST/test_taylor_options/order,nm,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_taylor_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct
IF(mesh%cad_type/=mesh_sphere_id)CALL oft_abort('Wrong mesh type, test for SPHERE only.','main',__FILE__)
IF(mg_test)THEN
  taylor_minlev=2
  IF(mesh%type==3)taylor_minlev=mg_mesh%mgmax
ELSE
  taylor_minlev=mg_mesh%mgmax+order-1
END IF
!---
CALL oft_hcurl_setup(order,taylor_minlev)
IF(mg_test)THEN
  CALL hcurl_setup_interp
  CALL hcurl_mloptions
END IF
!---
CALL oft_lag_setup(order,taylor_minlev)
IF(mg_test)THEN
  CALL lag_setup_interp
  CALL lag_mloptions
END IF
!---Run tests
oft_env%pm=.FALSE.
CALL taylor_hmodes(nm)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='taylor.results')
  WRITE(io_unit,*)taylor_hlam(:,oft_hcurl_nlevels)
  WRITE(io_unit,*)taylor_htor(:,oft_hcurl_nlevels)
  CLOSE(io_unit)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_taylor
