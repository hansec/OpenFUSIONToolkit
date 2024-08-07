!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thincurr_td.F90
!
!> @defgroup doxy_thincurr ThinCurr
!! ThinCurr drivers for 3D thin-wall modeling
!! @ingroup doxy_oft_bin
!
!> Run time-dependent thin wall simulations using ThinCurr
!!
!! **Option group:** `thincurr_td_options`
!! |  Option                 |  Description  | Type [dim] |
!! |-------------------------|---------------|------------|
!! |  `mesh_file="none"`     |  Surface mesh filename (Cubit) | str |
!! |  `dt=1.E-4`             |  Time step for time dependent run | float |
!! |  `nsteps=400`           |  Number of steps for time dependent run | int |
!! |  `timestep_cn=T`        |  Use Crank-Nicolson timestep | bool |
!! |  `cg_tol=1.E-6`         |  Convergence tolerance for `direct=F` | float |
!! |  `nplot=10`             |  Restart save frequency for time dependent run | int |
!! |  `plot_run=F`           |  Produce plot files from stored restart files | bool |
!! |  `direct=T`             |  Use direct solver | bool |
!! |  `save_mat=T`           |  Store inverted matrix for later use | bool |
!! |  `compute_b=F`          |  Compute magnetic fields on cell centers | bool |
!!
!! @authors Chris Hansen
!! @date Feb 2022
!! @ingroup doxy_thincurr
!---------------------------------------------------------------------------
PROGRAM thincurr_td
USE oft_base
USE oft_io, ONLY: hdf5_create_timestep, oft_bin_file
USE oft_mesh_type, ONLY: smesh
USE oft_mesh_native, ONLY: native_read_nodesets, native_read_sidesets
#ifdef HAVE_NCDF
USE oft_mesh_cubit, ONLY: cubit_read_nodesets, cubit_read_sidesets
#endif
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_lu, ONLY: lapack_cholesky
USE oft_native_la, ONLY: oft_native_dense_matrix
USE oft_deriv_matrices, ONLY: oft_sum_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE mhd_utils, ONLY: mu0
USE thin_wall
USE thin_wall_hodlr
USE thin_wall_solvers, ONLY: run_td_sim
IMPLICIT NONE
#include "local.h"
TYPE(tw_type), TARGET :: tw_sim
TYPE(tw_sensors) :: sensors
!
LOGICAL :: exists
INTEGER(4) :: i,j,k,n,ierr,io_unit,ncols,ntimes
REAL(8), ALLOCATABLE :: curr_ic(:)
REAL(8), POINTER, DIMENSION(:,:) :: curr_waveform,volt_waveform
TYPE(oft_timer) :: mytimer
CLASS(oft_vector), POINTER :: uio
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_ssets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL()
TYPE(oft_tw_hodlr_op), TARGET :: tw_hodlr
!
REAL(8) :: dt = 1.d-4
REAL(8) :: cg_tol=1.d-6
INTEGER(4) :: nsteps = 400
INTEGER(4) :: nstatus = 10
INTEGER(4) :: nplot = 10
INTEGER(4) :: jumper_start = -1
LOGICAL :: timestep_cn=.TRUE.
LOGICAL :: direct = .TRUE.
LOGICAL :: save_L = .FALSE.
LOGICAL :: save_Mcoil = .FALSE.
LOGICAL :: save_Msen = .FALSE.
LOGICAL :: plot_run = .FALSE.
LOGICAL :: compute_B = .FALSE.
CHARACTER(LEN=OFT_PATH_SLEN) :: curr_file="none"
CHARACTER(LEN=OFT_PATH_SLEN) :: volt_file="none"
NAMELIST/thincurr_td_options/curr_file,volt_file,dt,nsteps,nstatus,nplot,direct,save_L,save_Mcoil,save_Msen,  &
  plot_run,compute_B,timestep_cn,cg_tol,jumper_start
!---
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit,thincurr_td_options, IOSTAT=ierr)
CLOSE(io_unit)
if(ierr<0)call oft_abort('No thin-wall options found in input file.','thincurr_td',__FILE__)
if(ierr>0)call oft_abort('Error parsing thin-wall options in input file.','thincurr_td',__FILE__)
!---Setup mesh
CALL multigrid_construct_surf
! ALLOCATE(mg_mesh)
! mg_mesh%mgmax=1
! mg_mesh%nbase=1
! oft_env%nbase=1
! mg_mesh%mgdim=mg_mesh%mgmax
! CALL smesh_cubit_load
SELECT CASE(smesh%cad_type)
CASE(0)
  CALL native_read_nodesets(mesh_nsets)
  CALL native_read_sidesets(mesh_ssets)
CASE(2)
#ifdef HAVE_NCDF
  CALL cubit_read_nodesets(mesh_nsets)
  CALL cubit_read_sidesets(mesh_ssets)
#endif
CASE DEFAULT
  CALL oft_abort("Unsupported mesh type","thincurr_td",__FILE__)
END SELECT
IF(ASSOCIATED(mesh_ssets))THEN
  IF(mesh_ssets(1)%n>0)THEN
    tw_sim%nclosures=mesh_ssets(1)%n
    ALLOCATE(tw_sim%closures(tw_sim%nclosures))
    tw_sim%closures=mesh_ssets(1)%v
  END IF
END IF
tw_sim%mesh=>smesh
IF(jumper_start>0)THEN
  n=SIZE(mesh_nsets)
  hole_nsets=>mesh_nsets(1:jumper_start-1)
  jumper_nsets=>mesh_nsets(jumper_start:n)
ELSE
  hole_nsets=>mesh_nsets
END IF
CALL tw_sim%setup(hole_nsets)
IF((TRIM(curr_file)=="none").AND.(tw_sim%n_icoils>0))CALL oft_abort("No waveform filename specified", &
  "thincurr_td",__FILE__)
!---Setup I/0
CALL smesh%setup_io(1)
IF(oft_debug_print(1))CALL tw_sim%save_debug()
!---------------------------------------------------------------------------
! Time-dependent run
!---------------------------------------------------------------------------
!---Load drivers and sensors
CALL tw_load_sensors('floops.loc',tw_sim,sensors,jumper_nsets)
!---Compute inductances
WRITE(*,*)
IF(.NOT.plot_run)THEN
  IF(save_Mcoil)THEN
    CALL tw_compute_Ael2dr(tw_sim,'Mcoil.save')
  ELSE
    CALL tw_compute_Ael2dr(tw_sim)
  END IF
END IF
IF(save_Msen)THEN
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops,'Msen.save')
ELSE
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
END IF
!---------------------------------------------------------------------------
! Load or build element to element mutual matrix
!---------------------------------------------------------------------------
tw_hodlr%tw_obj=>tw_sim
CALL tw_hodlr%setup(.FALSE.)
IF(.NOT.plot_run)THEN
  IF(tw_hodlr%L_svd_tol>0.d0)THEN
    IF(direct)CALL oft_abort('HODLR compression does not support "direct=T"','thincurr_td',__FILE__)
    IF(save_L)CALL oft_abort('HODLR compression does not support "save_L=T"','thincurr_td',__FILE__)
    CALL tw_hodlr%compute_L()
  ELSE
    IF(save_L)THEN
      CALL tw_compute_LmatDirect(tw_sim,tw_sim%Lmat,save_file='Lmat.save')
    ELSE
      CALL tw_compute_LmatDirect(tw_sim,tw_sim%Lmat)
    END IF
  END IF
END IF
!---------------------------------------------------------------------------
! Run main calculation or plots
!---------------------------------------------------------------------------
IF(plot_run)THEN
  CALL plot_td_sim(tw_sim,dt,nsteps,nplot,sensors)
ELSE
  !---Setup resistivity matrix
  CALL tw_compute_Rmat(tw_sim,.FALSE.)
  !---Load I-coil waveform
  IF(TRIM(curr_file)/="none")THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(curr_file))
    READ(io_unit,*)ncols,ntimes
    ALLOCATE(curr_waveform(ntimes,ncols))
    DO i=1,ntimes
      READ(io_unit,*)curr_waveform(i,:)
    END DO
    CLOSE(io_unit)
  ELSE
    NULLIFY(curr_waveform)
  END IF
  !---Load V-coil waveform
  IF(TRIM(volt_file)/="none")THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(volt_file))
    READ(io_unit,*)ncols,ntimes
    ALLOCATE(volt_waveform(ntimes,ncols))
    DO i=1,ntimes
      READ(io_unit,*)volt_waveform(i,:)
    END DO
    CLOSE(io_unit)
  ELSE
    NULLIFY(volt_waveform)
  END IF
  !---Run time-dependent simulation
  ALLOCATE(curr_ic(tw_sim%nelems))
  curr_ic=0.d0
  oft_env%pm=.FALSE.
  IF(tw_hodlr%L_svd_tol>0.d0)THEN
    CALL run_td_sim(tw_sim,dt,nsteps,curr_ic,direct,cg_tol,timestep_cn,nstatus, &
      nplot,sensors,curr_waveform,volt_waveform,tw_hodlr)
  ELSE
    CALL run_td_sim(tw_sim,dt,nsteps,curr_ic,direct,cg_tol,timestep_cn,nstatus, &
      nplot,sensors,curr_waveform,volt_waveform)
  END IF
END IF
!---
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE plot_sim
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE plot_td_sim(self,dt,nsteps,nplot,sensors)
TYPE(tw_type), TARGET, INTENT(inout) :: self
REAL(8), INTENT(in) :: dt
INTEGER(4), INTENT(in) :: nsteps
INTEGER(4), INTENT(in) :: nplot
TYPE(tw_sensors), INTENT(in) :: sensors
!---
INTEGER(4) :: i,j,jj,k,ntimes_curr,ncols,itime,io_unit,face,ind1,ind2
REAL(8) :: uu,t,tmp,area,tmp2,val_prev
REAL(8), ALLOCATABLE, DIMENSION(:) :: coil_vec,senout,jumpout
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: cc_vals,coil_waveform
REAL(8), POINTER, DIMENSION(:) :: vals,vtmp
CLASS(oft_vector), POINTER :: u,Bx,By,Bz
TYPE(oft_bin_file) :: floop_hist,jumper_hist
LOGICAL :: exists
CHARACTER(LEN=4) :: pltnum
WRITE(*,*)'Post-processing simulation'
!---Load coil waveform
IF(TRIM(curr_file)/="none")THEN
  OPEN(NEWUNIT=io_unit,FILE=TRIM(curr_file))
  READ(io_unit,*)ncols,ntimes_curr
  IF(ncols-1/=self%n_icoils)CALL oft_abort('# of drivers in waveform does not match coil definitions', &
    'run_sim',__FILE__)
  ALLOCATE(coil_waveform(ntimes_curr,ncols))
  ALLOCATE(coil_vec(ncols-1))
  DO i=1,ntimes_curr
    READ(io_unit,*)coil_waveform(i,:)
    coil_waveform(i,2:ncols)=coil_waveform(i,2:ncols)*mu0 ! Convert to magnetic units
  END DO
  CLOSE(io_unit)
ELSE
  ntimes_curr=0
END IF
CALL self%Uloc%new(u)
!---
ALLOCATE(vals(self%nelems))
t=0.d0
IF(ntimes_curr>0)THEN
  DO j=1,self%n_icoils
    coil_vec(j)=linterp(coil_waveform(:,1),coil_waveform(:,j+1),ntimes_curr,t,1)
  END DO
END IF
!
CALL u%get_local(vals)
IF(sensors%nfloops>0)THEN
  ALLOCATE(senout(sensors%nfloops+1))
  senout=0.d0
  IF(ntimes_curr>0)CALL dgemv('N',sensors%nfloops,self%n_icoils,1.d0,self%Adr2sen, &
    sensors%nfloops,coil_vec,1,0.d0,senout(2),1)
  !---Setup history file
  IF(oft_env%head_proc)THEN
    floop_hist%filedesc = 'ThinCurr flux loop history file'
    CALL floop_hist%setup('floops.hist')
    CALL floop_hist%add_field('time', 'r8', desc="Simulation time [s]")
    DO i=1,sensors%nfloops
      CALL floop_hist%add_field(sensors%floops(i)%name, 'r8')
    END DO
    CALL floop_hist%write_header
    CALL floop_hist%open
    senout(1)=t
    CALL floop_hist%write(data_r8=senout)
  END IF
END IF
IF(sensors%njumpers>0)THEN
  ALLOCATE(jumpout(sensors%njumpers+1))
  DO j=1,sensors%njumpers
    tmp=0.d0
    val_prev=0.d0
    ind1=self%pmap(sensors%jumpers(j)%points(1))
    IF(ind1>0)val_prev=vals(ind1)
    DO k=1,sensors%jumpers(j)%np-1
      ind1=self%pmap(sensors%jumpers(j)%points(k+1))
      IF(ind1>0)THEN
        tmp=tmp+vals(ind1)-val_prev
        val_prev=vals(ind1)
      ELSE
        tmp=tmp-val_prev
        val_prev=0.d0
      END IF
    END DO
    DO k=1,self%nholes
      tmp=tmp+vals(self%np_active+k)*sensors%jumpers(j)%hole_facs(k)
    END DO
    jumpout(j+1)=tmp/mu0
  END DO
  !---Setup history file
  IF(oft_env%head_proc)THEN
    jumper_hist%filedesc = 'ThinCurr current jumper history file'
    CALL jumper_hist%setup('jumpers.hist')
    CALL jumper_hist%add_field('time', 'r8', desc="Simulation time [s]")
    DO i=1,sensors%njumpers
      CALL jumper_hist%add_field(sensors%jumpers(i)%name, 'r8')
    END DO
    CALL jumper_hist%write_header
    CALL jumper_hist%open
    jumpout(1)=t
    CALL jumper_hist%write(data_r8=jumpout)
  END IF
END IF
IF(compute_B)THEN
  IF(.NOT.ALLOCATED(cc_vals))ALLOCATE(cc_vals(3,self%mesh%np))
  IF(tw_hodlr%B_svd_tol>0.d0)THEN
    CALL tw_hodlr%compute_B()
    CALL self%Uloc_pts%new(Bx)
    CALL self%Uloc_pts%new(By)
    CALL self%Uloc_pts%new(Bz)
  ELSE
    CALL tw_compute_Bops(self)
  END IF
END IF
DO i=1,nsteps
  t=t+dt
  IF(MOD(i,nplot)==0)THEN
    IF(ntimes_curr>0)THEN
      DO j=1,self%n_icoils
        coil_vec(j)=linterp(coil_waveform(:,1),coil_waveform(:,j+1),ntimes_curr,t,1)
      END DO
    END IF
    !
    WRITE(pltnum,'(I4.4)')i
    CALL tw_rst_load(u,'pThinCurr_'//pltnum//'.rst','U')
    CALL u%get_local(vals)
    CALL hdf5_create_timestep(t)
    CALL tw_save_pfield(self,vals,'J')
    !
    IF(compute_B)THEN
      IF(tw_hodlr%B_svd_tol>0.d0)THEN
        CALL tw_hodlr%apply_bop(u,Bx,By,Bz)
        NULLIFY(vtmp)
        CALL Bx%get_local(vtmp)
        cc_vals(1,:)=vtmp
        CALL By%get_local(vtmp)
        cc_vals(2,:)=vtmp
        CALL Bz%get_local(vtmp)
        cc_vals(3,:)=vtmp
        IF(ntimes_curr>0)THEN
          !$omp parallel do private(k,tmp) collapse(2)
          DO j=1,self%n_icoils
            DO jj=1,3
              tmp=0.d0
              !$omp simd reduction(+:tmp)
              DO k=1,smesh%np
                tmp=tmp+tw_hodlr%Icoil_Bmat(k,j,jj)
              END DO
              cc_vals(jj,k)=cc_vals(jj,k)+tmp*coil_vec(j)
            END DO
          END DO
        END IF
      ELSE
        !$omp parallel do private(j,jj,tmp)
        DO k=1,smesh%np
          DO jj=1,3
            tmp=0.d0
            !$omp simd reduction(+:tmp)
            DO j=1,self%nelems
              tmp=tmp+vals(j)*self%Bel(j,k,jj)
            END DO
            cc_vals(jj,k)=tmp
          END DO
        END DO
        IF(ntimes_curr>0)THEN
          !$omp parallel do private(j,jj,tmp)
          DO k=1,smesh%np
            DO jj=1,3
              tmp=0.d0
              !$omp simd reduction(+:tmp)
              DO j=1,self%n_icoils
                tmp=tmp+coil_vec(j)*self%Bdr(j,k,jj)
              END DO
              cc_vals(jj,k)=cc_vals(jj,k)+tmp
            END DO
          END DO
        END IF
      END IF
      CALL self%mesh%save_vertex_vector(cc_vals,'B_v')
    END IF
    !
    IF(sensors%nfloops>0)THEN
      IF(ntimes_curr>0)CALL dgemv('N',sensors%nfloops,self%nelems,1.d0,self%Ael2sen, &
        sensors%nfloops,vals,1,0.d0,senout(2),1)
      CALL dgemv('N',sensors%nfloops,self%n_icoils,1.d0,self%Adr2sen,sensors%nfloops, &
        coil_vec,1,1.d0,senout(2),1)
      senout(1)=t
      CALL floop_hist%write(data_r8=senout)
    END IF
    IF(sensors%njumpers>0)THEN
      DO j=1,sensors%njumpers
        tmp=0.d0
        val_prev=0.d0
        ind1=self%pmap(sensors%jumpers(j)%points(1))
        IF(ind1>0)val_prev=vals(ind1)
        DO k=1,sensors%jumpers(j)%np-1
          ind1=self%pmap(sensors%jumpers(j)%points(k+1))
          IF(ind1>0)THEN
            tmp=tmp+vals(ind1)-val_prev
            val_prev=vals(ind1)
          ELSE
            tmp=tmp-val_prev
            val_prev=0.d0
          END IF
        END DO
        DO k=1,self%nholes
          tmp=tmp+vals(self%np_active+k)*sensors%jumpers(j)%hole_facs(k)
        END DO
        jumpout(j+1)=tmp/mu0
      END DO
      jumpout(1)=t
      CALL jumper_hist%write(data_r8=jumpout)
    END IF
  END IF
END DO
!---Cleanup
IF(sensors%nfloops>0)THEN
  CALL floop_hist%close()
  DEALLOCATE(senout)
END IF
IF(sensors%njumpers>0)THEN
  CALL jumper_hist%close()
  DEALLOCATE(jumpout)
END IF
CALL u%delete
DEALLOCATE(u)
DEALLOCATE(vals,vtmp)
IF(ntimes_curr>0)DEALLOCATE(coil_waveform,coil_vec)
END SUBROUTINE plot_td_sim
END PROGRAM thincurr_td
