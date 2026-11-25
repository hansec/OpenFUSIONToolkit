!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Regression test for the overlapping patch (vertex-star) smoother in
! oft_patch_solvers, applied to a quadratic Lagrange Laplacian on a cube mesh.
! The smoother is used as a GMRES preconditioner (additive, multiplicative and
! symmetric modes, or a single patch covering all rows, which is an exact solve)
! and as a stationary (Richardson) iteration.
!---------------------------------------------------------------------------------
PROGRAM test_native_patch
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_native_la, ONLY: oft_native_matrix
USE oft_native_solvers, ONLY: oft_native_gmres_solver
USE oft_patch_solvers, ONLY: oft_patch_precond
!
USE fem_base, ONLY: oft_fem_type, oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: comp_fem_vertex_patches
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_zerob, oft_lag_getlop, oft_lag_getmop
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
INTEGER(i4), PARAMETER :: order=2
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
INTEGER(i4) :: mode = 1
INTEGER(i4) :: rich_its = 0
LOGICAL :: whole = .FALSE.
LOGICAL :: single = .FALSE.
NAMELIST/test_patch_options/mode,rich_its,whole,single
IF(use_petsc)THEN
  WRITE(*,*)'SKIP TEST'
  STOP
END IF
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_patch_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---Setup FE representation
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
lag_zerob%ML_lag_rep=>ML_oft_lagrange
CALL test_lap()
!---Finalize enviroment
CALL oft_finalize
CONTAINS
SUBROUTINE test_lap()
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_patch_precond), TARGET :: patch_pre
TYPE(oft_fem_comp_type) :: comp_rep
INTEGER(i4) :: i,k
REAL(r8) :: uu,res0,res
CLASS(oft_vector), POINTER :: u,v,r
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Build operators and vectors
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
CALL ML_oft_lagrange%vec_create(r)
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
!---Vertex-star patches through a single-field composite representation
comp_rep%nfields=1
ALLOCATE(comp_rep%fields(1))
comp_rep%fields(1)%fe=>ML_oft_lagrange%current_level
IF(whole)THEN
  ALLOCATE(lop%patch_info)
  lop%patch_info%npatch=1
  lop%patch_info%ncolors=1
  ALLOCATE(lop%patch_info%patches(1),lop%patch_info%color(1),lop%patch_info%gid(1))
  ALLOCATE(lop%patch_info%owner(1),lop%patch_info%share(1))
  lop%patch_info%patches(1)%n=lop%nr
  ALLOCATE(lop%patch_info%patches(1)%v(lop%nr))
  DO i=1,lop%nr
    lop%patch_info%patches(1)%v(i)=i
  END DO
  lop%patch_info%color=1
  lop%patch_info%gid=1
  lop%patch_info%owner=oft_env%rank
  lop%patch_info%share(1)%n=0
  ALLOCATE(lop%patch_info%share(1)%v(0))
ELSE
  CALL comp_fem_vertex_patches(comp_rep,lop%patch_info)
END IF
!---Setup patch smoother
patch_pre%mode=mode
patch_pre%single=single
!---Right hand side
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_zerob%apply(v)
CALL u%set(0.d0)
IF(rich_its>0)THEN
  !---Stationary iteration: report residual reduction factor
  patch_pre%A=>lop
  patch_pre%its=rich_its
  patch_pre%bc=>lag_zerob
  res0=SQRT(v%dot(v))
  CALL r%add(0.d0,1.d0,v)
  CALL patch_pre%apply(u,r)
  res=SQRT(r%dot(r))
  k=rich_its
  uu=(res/res0)**(1.d0/REAL(rich_its,8))
ELSE
  !---Preconditioned GMRES
  linv%A=>lop
  linv%its=-3
  linv%itplot=1
  linv%nrits=20
  linv%pre=>patch_pre
  CALL linv%apply(u,v)
  k=linv%cits
  uu=SQRT(u%dot(u))
END IF
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='patch.results')
  WRITE(io_unit,*)k
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Cleanup
CALL u%delete
CALL v%delete
CALL r%delete
DEALLOCATE(u,v,r)
CALL lop%delete
CALL mop%delete
DEALLOCATE(lop,mop)
CALL patch_pre%delete
IF(rich_its<=0)CALL linv%delete
END SUBROUTINE test_lap
END PROGRAM test_native_patch
