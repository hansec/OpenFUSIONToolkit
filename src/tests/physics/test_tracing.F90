!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_tracing.F90
!
!> Regression tests for 3D field line tracing
!!
!! Traces field lines of an analytic helical field, B = (-y, x, Bz), through a cube
!! mesh (distributed over processors for parallel runs) and compares the results to
!! the exact solution. Covers tracer types 1 (LSODE) and 4 (Runge-Kutta with steps limited
!! at cell boundaries), including transfers of the integrator state between domains.
!!
!! On periodic meshes field lines are continued across the periodic boundaries. Three
!! analytic fields (`field_type`), each invariant under the periodic maps, are available:
!! - (1) Helical field above on a cube mesh that is periodic in z
!! - (2) B = (1, 1, 1/3 + 0.1 cos(pi x)) on a cube mesh that is periodic in x and y (and
!!   optionally z). Lines starting at x=y=0 cross the x and y boundaries at the same time
!!   (edges of the periodic domain). On the triply periodic mesh the line starting at the
!!   origin also crosses the z boundaries there (corners of the periodic domain) every third
!!   edge crossing
!! - (3) B = (Bx, -z, y), rotation about the x-axis, on a native mesh of a wedge about the
!!   x-axis with revolved periodicity, so that lines repeatedly cross the periodic planes
!!
!! A time variable is traced as a fourth equation, so the end point can be compared to the
!! exact solution, mapped into the periodic domain. For field 1 a Poincare mode is also
!! available, which checks crossings of the plane x=0 and that the jumps across the periodic
!! boundaries are not counted as crossings of the plane z=0.
!!
!! Three modes are available on non-periodic meshes:
!! - Single field line (`tracing_line`): the end point is compared to the exact solution
!! - Poincare section (`tracing_poincare`): several field lines are traced in parallel
!!   and the crossings of the plane x=0 are compared to the exact helices
!! - Field with jumps between cells (`field_jump>0`): the helical field is scaled by a factor
!!   that varies from cell to cell, so the field jumps at every cell face as for finite element
!!   fields while its normal component keeps its sign. The field line still follows the exact
!!   helix, which is checked here, while the time along it (traced as a fourth equation) is
!!   compared by the driver to a reference trace with a tight tolerance
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!--------------------------------------------------------------------------------
MODULE test_tracing_helpers
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
INTEGER(i4) :: field_type = 1 !< Analytic field (1: helical, 2: planar drift, 3: rotation about x)
REAL(r8) :: swirl_bz = 0.3d0 !< Axial field (pitch) of the helical field
REAL(r8) :: swirl_jump = 0.d0 !< Amplitude of the variation of the field strength between cells
REAL(r8), PARAMETER :: drift_v(3) = [1.d0,1.d0,1.d0/3.d0] !< Mean field of the planar drift field
REAL(r8), PARAMETER :: drift_eps = 0.1d0 !< Amplitude of the variation of the planar drift field with x
REAL(r8), PARAMETER :: rot_bx = 0.1d0 !< Axial field of the rotation about the x-axis
!> Analytic helical field
TYPE, EXTENDS(fem_interp) :: swirl_interp
CONTAINS
  PROCEDURE :: interp => swirl_interp_apply
END TYPE swirl_interp
CONTAINS
!> Evaluate field at a logical position in a cell
SUBROUTINE swirl_interp_apply(self,cell,f,gop,val)
CLASS(swirl_interp), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cell
REAL(r8), INTENT(in) :: f(:)
REAL(r8), INTENT(in) :: gop(3,4)
REAL(r8), INTENT(out) :: val(:)
REAL(r8) :: pt(3),pc(3)
pt=self%mesh%log2phys(cell,f)
SELECT CASE(field_type)
  CASE(2)
    val(1:3)=drift_v
    val(3)=val(3)+drift_eps*COS(pi*pt(1))
    RETURN
  CASE(3)
    val(1:3)=[rot_bx,-pt(3),pt(2)]
    RETURN
END SELECT
val(1)=-pt(2)
val(2)=pt(1)
val(3)=swirl_bz
IF(swirl_jump>0.d0)THEN
  !---Term set by the cell center, so it is independent of the cell numbering and decomposition
  IF(self%mesh%type==3)THEN
    pc=self%mesh%log2phys(cell,[0.5d0,0.5d0,0.5d0,0.d0])
  ELSE
    pc=self%mesh%log2phys(cell,[0.25d0,0.25d0,0.25d0,0.25d0])
  END IF
  val(1:3)=val(1:3)*(1.d0+swirl_jump*SIN(37.d0*pc(1)+53.d0*pc(2)+71.d0*pc(3)))
END IF
END SUBROUTINE swirl_interp_apply
!> Field line ODE with time as a fourth equation
SUBROUTINE time_ydot(t,y,B,n,ydot)
REAL(r8), INTENT(in) :: t
REAL(r8), INTENT(in) :: y(n)
REAL(r8), INTENT(in) :: B(3)
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(out) :: ydot(n)
ydot(1:3)=B
ydot(4)=1.d0
END SUBROUTINE time_ydot
!> Exact field line position at time t, starting from pt0 at t=0 (not mapped into the periodic domain)
FUNCTION exact_point(pt0,t) RESULT(pt)
REAL(r8), INTENT(in) :: pt0(3)
REAL(r8), INTENT(in) :: t
REAL(r8) :: pt(3)
SELECT CASE(field_type)
  CASE(2)
    pt=pt0+drift_v*t
    pt(3)=pt(3)+drift_eps/(pi*drift_v(1))*(SIN(pi*pt(1))-SIN(pi*pt0(1)))
  CASE(3)
    pt(1)=pt0(1)+rot_bx*t
    pt(2)=pt0(2)*COS(t)-pt0(3)*SIN(t)
    pt(3)=pt0(2)*SIN(t)+pt0(3)*COS(t)
  CASE DEFAULT
    pt(1)=pt0(1)*COS(t)-pt0(2)*SIN(t)
    pt(2)=pt0(1)*SIN(t)+pt0(2)*COS(t)
    pt(3)=pt0(3)+swirl_bz*t
END SELECT
END FUNCTION exact_point
!> Distance between a traced point and the exact point, using the periodic image of the exact
!! point that is closest (wrapped into the domain, and possibly mapped once more across a boundary
!! as the two points can lie on opposite sides of a periodic boundary)
FUNCTION periodic_error(mesh,pt,pt_exact) RESULT(err)
USE oft_mesh_type, ONLY: oft_mesh, mesh_periodic_wrap, periodic_map_inverse
CLASS(oft_mesh), INTENT(in) :: mesh
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(in) :: pt_exact(3)
REAL(r8) :: err
INTEGER(i4) :: p
REAL(r8) :: pe(3),amap(3,4)
LOGICAL :: mapped
pe=pt_exact
mapped=mesh_periodic_wrap(mesh,pe,amap)
err=magnitude(pt-pe)
DO p=1,mesh%periodic%nper
  amap=mesh%periodic%pmap(:,:,p)
  err=MIN(err,magnitude(pt-MATMUL(amap(:,1:3),pe)-amap(:,4)))
  amap=periodic_map_inverse(amap)
  err=MIN(err,magnitude(pt-MATMUL(amap(:,1:3),pe)-amap(:,4)))
END DO
END FUNCTION periodic_error
END MODULE test_tracing_helpers
!
PROGRAM test_tracing
USE oft_base
USE oft_io, ONLY: hdf5_create_file
USE oft_mesh_cube, ONLY: mesh_cube_id
USE oft_mesh_native, ONLY: mesh_native_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE tracing, ONLY: oft_tracer, create_tracer, tracing_line, tracing_poincare
USE test_tracing_helpers, ONLY: swirl_interp, swirl_bz, swirl_jump, time_ydot, field_type, &
  exact_point, periodic_error
IMPLICIT NONE
INTEGER(i4) :: iounit,ierr
REAL(r8) :: pt(3),phase,err
TYPE(multigrid_mesh), TARGET :: mg_mesh
TYPE(swirl_interp), TARGET :: Bfield
CLASS(oft_tracer), POINTER :: tracer
INTEGER(i4) :: tracer_type=1
INTEGER(i4) :: max_steps=100000
LOGICAL :: test_poincare=.FALSE.
REAL(r8) :: field_jump=0.d0
REAL(r8) :: tol=1.d-10
NAMELIST/test_tracing_options/tracer_type,test_poincare,field_jump,tol,field_type,max_steps
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=iounit,FILE=oft_env%ifile)
READ(iounit,test_tracing_options,IOSTAT=ierr)
CLOSE(iounit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(field_type==3)THEN
  IF(mg_mesh%mesh%cad_type/=mesh_native_id)CALL oft_abort('Wrong mesh type, test for NATIVE only.','test_tracing',__FILE__)
ELSE
  IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','test_tracing',__FILE__)
END IF
!---Setup tracer
Bfield%mesh=>mg_mesh%mesh
Bfield%dim=3
CALL create_tracer(tracer,tracer_type)
tracer%tol=tol
tracer%maxsteps=max_steps
tracer%maxtrans=1000
tracer%B=>Bfield
tracer%ml_mesh=>mg_mesh
IF(mg_mesh%mesh%periodic%nper>0)THEN
  IF(test_poincare)THEN
    CALL run_periodic_poincare_test
  ELSE
    CALL run_periodic_test
  END IF
  CALL tracer%delete()
  CALL oft_finalize
END IF
IF(test_poincare)THEN
  CALL run_poincare_test
  CALL tracer%delete()
  CALL oft_finalize
END IF
IF(field_jump>0.d0)THEN
  CALL run_jump_test
  CALL tracer%delete()
  CALL oft_finalize
END IF
!---Trace from (0.5,0,-0.9): circle of radius 0.5 advancing in z until exiting at z=1
pt=[0.5d0,0.d0,-0.9d0]
IF(oft_env%head_proc)CALL hdf5_create_file('trace_test')
CALL oft_mpi_barrier(ierr)
CALL tracing_line(tracer,pt,'trace_test')
!---Check result against exact solution x=0.5*cos(t), y=0.5*sin(t), z=z0+0.3*t
IF(tracer%status/=-1)CALL oft_abort('Tracer did not exit mesh as expected','test_tracing',__FILE__)
IF(tracer%nsteps<10)CALL oft_abort('Too few tracer steps','test_tracing',__FILE__)
IF(oft_env%nprocs>1.AND.tracer%ntrans<1)CALL oft_abort('No domain transfers occurred','test_tracing',__FILE__)
IF(tracer%y(3)<0.9d0)CALL oft_abort('Tracer exited early','test_tracing',__FILE__)
phase=(tracer%y(3)-pt(3))/0.3d0
err=SQRT((tracer%y(1)-0.5d0*COS(phase))**2+(tracer%y(2)-0.5d0*SIN(phase))**2)
IF(oft_env%head_proc)WRITE(*,'(A,ES11.3,2I8)')'  Final error, steps, transfers = ',err,tracer%nsteps,tracer%ntrans
IF(err>1.d-6)CALL oft_abort('Tracer end point differs from exact solution','test_tracing',__FILE__)
CALL tracer%delete()
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---Poincare section of several field lines, compared to the exact helices
SUBROUTINE run_poincare_test
INTEGER(i4), PARAMETER :: nlines=3
INTEGER(i4), PARAMETER :: ncross=6
REAL(r8), PARAMETER :: z0=-0.9d0
INTEGER(i4) :: i,k,iline,ncount(nlines),nfound(ncross,nlines),io_stat
REAL(r8) :: pts(3,nlines),radius(nlines),pt_cross(3),z_exact,err,err_max
LOGICAL :: have_omp
!---Poincare tracing requires OpenMP
have_omp=.FALSE.
!$ have_omp=.TRUE.
IF(.NOT.have_omp)THEN
  IF(oft_env%head_proc)WRITE(*,*)'SKIP TEST'
  CALL oft_finalize
END IF
!---Lower pitch so each line makes several turns before exiting at z=1
swirl_bz=0.05d0
radius=[0.25d0,0.5d0,0.75d0]
DO i=1,nlines
  pts(:,i)=[radius(i),0.d0,z0]
END DO
CALL tracing_poincare(tracer,pts,nlines,'poincare_test.dat')
CALL oft_mpi_barrier(ierr)
!---Exact crossings of x=0 (x increasing): theta=3*pi/2+2*pi*k -> (0,-r,z0+Bz*theta)
!---Lines exit at z=1 (theta=38), so k=0,...,5 are crossed and the next at z=1.22 is not
IF(oft_env%head_proc)THEN
  ncount=0
  nfound=0
  err_max=0.d0
  OPEN(NEWUNIT=iounit,FILE='poincare_test.dat',STATUS='OLD')
  DO
    READ(iounit,*,IOSTAT=io_stat)iline,pt_cross
    IF(io_stat/=0)EXIT
    IF(iline<1.OR.iline>nlines)CALL oft_abort('Invalid field line index','test_tracing',__FILE__)
    ncount(iline)=ncount(iline)+1
    k=NINT(((pt_cross(3)-z0)/swirl_bz-1.5d0*pi)/(2.d0*pi))
    IF(k<0.OR.k>=ncross)CALL oft_abort('Unexpected crossing','test_tracing',__FILE__)
    nfound(k+1,iline)=nfound(k+1,iline)+1
    z_exact=z0+swirl_bz*(1.5d0*pi+2.d0*pi*k)
    err=SQRT(pt_cross(1)**2+(pt_cross(2)+radius(iline))**2+(pt_cross(3)-z_exact)**2)
    err_max=MAX(err_max,err)
  END DO
  CLOSE(iounit)
  WRITE(*,'(A,ES11.3,3I4)')'  Max crossing error, crossings per line = ',err_max,ncount
  IF(ANY(nfound/=1))CALL oft_abort('Wrong number of crossings','test_tracing',__FILE__)
  IF(err_max>1.d-3)CALL oft_abort('Crossings differ from exact solution','test_tracing',__FILE__)
END IF
END SUBROUTINE run_poincare_test
!---Single field line with a field that jumps at every cell face
SUBROUTINE run_jump_test
REAL(r8) :: pt4(4)
swirl_jump=field_jump
tracer%neq=4
tracer%ydot=>time_ydot
pt4=[0.5d0,0.d0,-0.9d0,0.d0]
IF(oft_env%head_proc)CALL hdf5_create_file('trace_test')
CALL oft_mpi_barrier(ierr)
CALL tracing_line(tracer,pt4,'trace_test')
IF(tracer%status/=-1)CALL oft_abort('Tracer did not exit mesh as expected','test_tracing',__FILE__)
IF(oft_env%nprocs>1.AND.tracer%ntrans<1)CALL oft_abort('No domain transfers occurred','test_tracing',__FILE__)
IF(tracer%y(3)<0.9d0)CALL oft_abort('Tracer exited early','test_tracing',__FILE__)
!---The jumps only change the field strength, so the field line is still the exact helix
phase=(tracer%y(3)-pt4(3))/swirl_bz
err=SQRT((tracer%y(1)-0.5d0*COS(phase))**2+(tracer%y(2)-0.5d0*SIN(phase))**2)
IF(oft_env%head_proc)WRITE(*,'(A,ES11.3,2I8)')'  Final error, steps, transfers = ',err,tracer%nsteps,tracer%ntrans
IF(err>MAX(1.d2*tol,1.d-8))CALL oft_abort('Tracer end point is not on the exact helix','test_tracing',__FILE__)
END SUBROUTINE run_jump_test
!---Single field line on a periodic mesh, continued across the periodic boundaries. The time is
!---traced as a fourth equation, so the end point can be compared to the exact solution
SUBROUTINE run_periodic_test
INTEGER(i4) :: expected_status
REAL(r8) :: pt4(4),pt_exact(3),t_end,ncross
tracer%neq=4
tracer%ydot=>time_ydot
expected_status=-3 ! Line stays on the mesh until the step limit
SELECT CASE(field_type)
  CASE(2)
    pt4=[0.d0,0.d0,0.d0,0.d0]
    IF(mg_mesh%mesh%periodic%nper<3)THEN
      pt4(3)=-0.9d0
      expected_status=-1 ! Line leaves the mesh at z=1
    END IF
  CASE(3)
    pt4=[0.05d0,0.75d0,0.d0,0.d0]
    expected_status=-1 ! Line leaves the mesh at x=1
  CASE DEFAULT
    pt4=[0.5d0,0.d0,-0.9d0,0.d0]
END SELECT
IF(oft_env%head_proc)CALL hdf5_create_file('trace_test')
CALL oft_mpi_barrier(ierr)
CALL tracing_line(tracer,pt4,'trace_test')
t_end=tracer%y(4)
pt_exact=exact_point(pt4(1:3),t_end)
err=periodic_error(mg_mesh%mesh,tracer%y(1:3),pt_exact)
!---Number of periodic boundaries crossed by the exact solution (along the first periodic direction)
IF(mg_mesh%mesh%periodic%revolved)THEN
  ncross=t_end/(2.d0*mg_mesh%mesh%periodic%prange(2,1))
ELSE
  ncross=ABS(pt_exact(mg_mesh%mesh%periodic%paxis(1))-pt4(mg_mesh%mesh%periodic%paxis(1))) &
    /(mg_mesh%mesh%periodic%prange(2,1)-mg_mesh%mesh%periodic%prange(1,1))
END IF
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,ES11.3,2I8)')'  Final error, steps, transfers = ',err,tracer%nsteps,tracer%ntrans
  WRITE(*,'(A,2F10.4,I4)')'  Final time, periodic crossings, status = ',t_end,ncross,tracer%status
END IF
IF(tracer%status/=expected_status)CALL oft_abort('Unexpected tracer status','test_tracing',__FILE__)
IF(oft_env%nprocs>1.AND.tracer%ntrans<1)CALL oft_abort('No domain transfers occurred','test_tracing',__FILE__)
IF(ncross<2.d0)CALL oft_abort('Too few periodic crossings','test_tracing',__FILE__)
IF(field_type==2.AND.expected_status==-1.AND.tracer%y(3)<0.99d0)CALL oft_abort('Tracer exited early','test_tracing',__FILE__)
IF(field_type==3.AND.tracer%y(1)<0.99d0)CALL oft_abort('Tracer exited early','test_tracing',__FILE__)
IF(err>1.d-6)CALL oft_abort('Tracer end point differs from exact solution','test_tracing',__FILE__)
END SUBROUTINE run_periodic_test
!---Poincare sections of several field lines descending through a mesh that is periodic in z.
!---Crossings of x=0 are compared to the exact helices, with z mapped into the periodic domain,
!---and the plane z=0 is never crossed upwards, so the jumps from z=-1 to z=1 at the periodic
!---boundary must not give crossings
SUBROUTINE run_periodic_poincare_test
INTEGER(i4), PARAMETER :: nlines=3
INTEGER(i4), PARAMETER :: kmax=200
REAL(r8), PARAMETER :: z0=0.9d0
INTEGER(i4) :: i,k,kbest,iline,ncount(nlines),nfound(0:kmax,nlines),io_stat
REAL(r8) :: pts(3,nlines),radius(nlines),pt_cross(3),dz,dz_min,err,err_max
LOGICAL :: have_omp
!---Poincare tracing requires OpenMP
have_omp=.FALSE.
!$ have_omp=.TRUE.
IF(.NOT.have_omp)THEN
  IF(oft_env%head_proc)WRITE(*,*)'SKIP TEST'
  CALL oft_finalize
END IF
swirl_bz=-0.2d0
radius=[0.25d0,0.5d0,0.75d0]
DO i=1,nlines
  pts(:,i)=[radius(i),0.d0,z0]
END DO
CALL tracing_poincare(tracer,pts,nlines,'poincare_test.dat')
CALL oft_mpi_barrier(ierr)
CALL tracing_poincare(tracer,pts,nlines,'poincare_z.dat',offset=0.d0,pcoord=3)
CALL oft_mpi_barrier(ierr)
IF(oft_env%head_proc)THEN
  !---Exact crossings of x=0 (x increasing): theta=3*pi/2+2*pi*k -> (0,-r,z0+Bz*theta), mapped into [-1,1)
  ncount=0
  nfound=0
  err_max=0.d0
  OPEN(NEWUNIT=iounit,FILE='poincare_test.dat',STATUS='OLD')
  DO
    READ(iounit,*,IOSTAT=io_stat)iline,pt_cross
    IF(io_stat/=0)EXIT
    IF(iline<1.OR.iline>nlines)CALL oft_abort('Invalid field line index','test_tracing',__FILE__)
    ncount(iline)=ncount(iline)+1
    !---Find the crossing with the closest periodic image
    dz_min=HUGE(1.d0)
    kbest=-1
    DO k=0,kmax
      dz=pt_cross(3)-(z0+swirl_bz*(1.5d0*pi+2.d0*pi*k))
      dz=dz-2.d0*ANINT(dz/2.d0)
      IF(ABS(dz)<dz_min)THEN
        dz_min=ABS(dz)
        kbest=k
      END IF
    END DO
    nfound(kbest,iline)=nfound(kbest,iline)+1
    err=SQRT(pt_cross(1)**2+(pt_cross(2)+radius(iline))**2+dz_min**2)
    err_max=MAX(err_max,err)
  END DO
  CLOSE(iounit)
  WRITE(*,'(A,ES11.3,3I4)')'  Max crossing error, crossings per line = ',err_max,ncount
  !---Each crossing up to the last one is found exactly once (the lines run until the step limit)
  DO iline=1,nlines
    IF(ncount(iline)<3)CALL oft_abort('Too few crossings','test_tracing',__FILE__)
    IF(ncount(iline)>kmax)CALL oft_abort('Too many crossings','test_tracing',__FILE__)
    IF(ANY(nfound(0:ncount(iline)-1,iline)/=1))CALL oft_abort('Wrong crossings','test_tracing',__FILE__)
  END DO
  IF(err_max>1.d-3)CALL oft_abort('Crossings differ from exact solution','test_tracing',__FILE__)
  !---No crossings of z=0
  k=0
  OPEN(NEWUNIT=iounit,FILE='poincare_z.dat',STATUS='OLD')
  DO
    READ(iounit,*,IOSTAT=io_stat)iline,pt_cross
    IF(io_stat/=0)EXIT
    k=k+1
  END DO
  CLOSE(iounit)
  WRITE(*,'(A,I4)')'  Crossings of z=0 = ',k
  IF(k/=0)CALL oft_abort('Periodic boundary crossing counted as section crossing','test_tracing',__FILE__)
END IF
END SUBROUTINE run_periodic_poincare_test
END PROGRAM test_tracing
