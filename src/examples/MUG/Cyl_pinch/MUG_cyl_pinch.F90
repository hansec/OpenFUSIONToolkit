!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Cylindrical pinch kink {#doc_mug_cyl_kink}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of \ref doc_mhd_main "MUG" to compute the linear growth rate of
!! an ideal external kink mode in a periodic screw pinch ("straight tokamak") surrounded by a conducting wall.
!!
!!\section doc_mug_cyl_kink_theory Background
!!
!! The plasma column, with radius \f$ a \f$, carries a uniform (flat) axial current density in a uniform
!! axial magnetic field \f$ B_z \f$ and is separated from a perfectly conducting wall at \f$ r = b \f$ by a vacuum region.
!! The axial direction is periodic with length \f$ L = 2 \pi R_0 \f$, so that axial mode \f$ n \f$ has
!! wavenumber \f$ k = n/R_0 \f$ and the safety factor at the plasma edge is
!! \f$ q_a = a B_z / (R_0 B_{\theta}(a)) \f$, which is constant across the plasma for a flat current.
!!
!! [Shafranov](https://ui.adsabs.harvard.edu/abs/1970SPTP...15..175S/abstract) showed that, for large aspect ratio
!! (\f$ a/R_0 \ll 1 \f$) and low \f$ \beta \f$, external kink modes with poloidal mode number \f$ m \ge 2 \f$
!! are unstable for \f$ m - 1 + (a/b)^{2m} < n q_a < m \f$ with growth rate
!!
!! \f[ (\gamma \tau_A)^2 = \frac{2}{q_a^2} (m - n q_a) \left[ 1 - \frac{m - n q_a}{1 - (a/b)^{2m}} \right] \f]
!!
!! where \f$ \tau_A = R_0 / V_A \f$ and \f$ V_A = B_z / \sqrt{\mu_0 \rho} \f$. This result is also derived in
!! J. P. Freidberg, *Ideal MHD* (Cambridge University Press, 2014) and is used as a verification test for
!! other 3D MHD codes, for example the SpeCyl code in
!! [Galeazzi et al. (2026)](https://arxiv.org/abs/2610.13428).
!!
!! In this example we use \f$ a = 0.1 \f$ m, \f$ b = 2 a \f$, \f$ R_0 = 0.5 \f$ m (\f$ a/R_0 = 0.2 \f$), \f$ B_z = 0.8 \f$ T,
!! a deuterium density of \f$ 10^{20} \f$ m\f$ ^{-3} \f$ and \f$ q_a = 1.6 \f$. For these parameters the \f$ m/n = 2/1 \f$
!! mode is unstable, with \f$ \gamma \tau_A = 0.423 \f$ (\f$ \gamma = 1.04 \times 10^6 \f$ s\f$ ^{-1} \f$), even with the wall present.
!!
!!\section doc_mug_cyl_kink_code Code Walk Through
!!
!! The code consists of a helper module, which defines the initial conditions, diagnostics and vacuum
!! resistivity profile, and a main program that constructs the initial conditions and runs the simulation.
!!
!!\subsection doc_mug_cyl_kink_code_helper Helper Module
!!
!! The helper module defines three objects used by the main program:
!!  - `cyl_pinch_init`: an interpolation object that evaluates the equilibrium density, temperature and magnetic field, loaded
!!    from a radial profile file, and the initial magnetic perturbation.
!!  - `cyl_probe`: a probe object that records the magnetic field and velocity at points around the plasma
!!    edge during post-processing (`plot_run=T`).
!!  - `vac_diff_shape`: a spatial shaping function for the resistivity, used to model the vacuum region.
! START SOURCE
MODULE cyl_pinch
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh
USE fem_base, ONLY: oft_fem_type
USE fem_utils, ONLY: fem_interp
USE oft_lag_operators, ONLY: oft_lag_vrinterp
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_rinterp
USE diagnostic, ONLY: field_probe
USE mhd_utils, ONLY: mu0
USE xmhd, ONLY: xmhd_sub_fields, oft_xmhd_probe, xmhd_ML_lagrange, xmhd_ML_hcurl_grad
IMPLICIT NONE
!------------------------------------------------------------------------------
!> Equilibrium and perturbation fields for the cylindrical pinch
!!
!! The equilibrium is read from a file with the number of points on the first line,
!! followed by one line per radial point containing \f$ r \f$ [m], \f$ B_{\theta} \f$ [T],
!! \f$ B_z \f$ [T], \f$ n \f$ [m^-3] and \f$ T \f$ [eV]. The perturbation is the curl-free field of a
!! helical \f$ (m,n) \f$ current sheet at the plasma edge with zero normal field at the wall.
!------------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: cyl_pinch_init
  CHARACTER(LEN=80) :: filename = 'none' !< Equilibrium profile file
  CHARACTER(LEN=1) :: field = 'n' !< Field component to evaluate ('n', 't', 'b' or 'd')
  INTEGER(i4) :: npts = 0 !< Number of points in profile
  INTEGER(i4) :: m = 2 !< Poloidal mode number of perturbation
  INTEGER(i4) :: n = 1 !< Axial mode number of perturbation
  REAL(r8) :: L = 1.d0 !< Periodic length of cylinder [m]
  REAL(r8), POINTER, DIMENSION(:,:) :: vals => NULL() !< Profile values [5,npts]
CONTAINS
  !> Load equilibrium profile
  PROCEDURE :: setup => cyl_pinch_setup
  !> Evaluate field
  PROCEDURE :: interp => cyl_pinch_apply
END TYPE cyl_pinch_init
!------------------------------------------------------------------------------
!> Synthetic probes at the plasma edge for post-processing
!------------------------------------------------------------------------------
TYPE, EXTENDS(oft_xmhd_probe) :: cyl_probe
  LOGICAL :: initialized = .FALSE. !< Flag to indicate setup has been called
  INTEGER(i4) :: npts = 180 !< Number of probes (evenly spaced in \f$ \theta \f$ at \f$ r = a \f$, \f$ z = 0 \f$)
  TYPE(oft_hcurl_grad_rinterp), POINTER :: Bfield => NULL() !< Magnetic field interpolation object
  TYPE(oft_lag_vrinterp), POINTER :: Vfield => NULL() !< Velocity field interpolation object
  TYPE(field_probe) :: b_probes !< Magnetic field probes
  TYPE(field_probe) :: v_probes !< Velocity probes
CONTAINS
  !> Sample probe signals
  PROCEDURE :: apply => cyl_probe_apply
END TYPE cyl_probe
REAL(r8) :: r_plasma = 0.1d0 !< Plasma radius [m]
REAL(r8) :: r_wall = 0.2d0 !< Wall radius [m] (set from mesh)
REAL(r8) :: eta_vac = 1.d4 !< Resistivity enhancement factor in vacuum region
REAL(r8) :: vac_start = 0.005d0 !< Distance outside the plasma radius where resistivity ramp starts [m]
REAL(r8) :: vac_width = 0.01d0 !< Width of resistivity ramp [m]
CONTAINS
!------------------------------------------------------------------------------
!> Load equilibrium profile and set geometric parameters from the mesh
!------------------------------------------------------------------------------
SUBROUTINE cyl_pinch_setup(self,mesh)
CLASS(cyl_pinch_init), INTENT(inout) :: self
CLASS(oft_mesh), TARGET, INTENT(in) :: mesh !< Mesh
INTEGER(i4) :: i,io_unit
REAL(r8) :: zmin,zmax
self%mesh=>mesh
!---Get periodic length and wall radius from mesh
zmin=-oft_mpi_max(-MINVAL(mesh%r(3,:)))
zmax=oft_mpi_max(MAXVAL(mesh%r(3,:)))
r_wall=oft_mpi_max(MAXVAL(SQRT(mesh%r(1,:)**2+mesh%r(2,:)**2)))
self%L=zmax-zmin
!---Load equilibrium profile
OPEN(NEWUNIT=io_unit,FILE=TRIM(self%filename))
READ(io_unit,*)self%npts
ALLOCATE(self%vals(5,self%npts))
DO i=1,self%npts
  READ(io_unit,*)self%vals(:,i)
END DO
CLOSE(io_unit)
IF(oft_env%head_proc)THEN
  WRITE(*,'(2A)')'Cylindrical pinch equilibrium loaded: ',TRIM(self%filename)
  WRITE(*,'(A,ES11.3)')'  L       = ',self%L
  WRITE(*,'(A,ES11.3)')'  a       = ',r_plasma
  WRITE(*,'(A,ES11.3)')'  b       = ',r_wall
  WRITE(*,'(A,ES11.3)')'  B_z     = ',self%vals(3,1)
  WRITE(*,'(A,ES11.3)')'  I_p     = ',self%vals(2,self%npts)*2.d0*pi*self%vals(1,self%npts)/mu0
  WRITE(*,'(A,I8)')    '  # of points: ',self%npts
END IF
END SUBROUTINE cyl_pinch_setup
!------------------------------------------------------------------------------
!> Evaluate the requested field at a given point
!!
!! For the perturbation ('d') the radial and poloidal components of
!! \f$ \nabla \times (\hat{z} \psi) \f$ are used, where
!! \f$ \psi \propto \hat{r}^m \f$ for \f$ \hat{r} = r/a < 1 \f$ and
!! \f$ \psi \propto (\hat{r}_w^{2m} \hat{r}^{-m} - \hat{r}^m)/(\hat{r}_w^{2m} - 1) \f$ outside,
!! with \f$ \hat{r}_w = b/a \f$.
!------------------------------------------------------------------------------
SUBROUTINE cyl_pinch_apply(self,cell,f,gop,val)
CLASS(cyl_pinch_init), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cell !< Cell for interpolation
REAL(r8), INTENT(in) :: f(:) !< Position in cell in logical coord [4]
REAL(r8), INTENT(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
REAL(r8), INTENT(out) :: val(:) !< Field value at point
REAL(r8) :: pt(3),rad,theta,Bp,Bz,rhat(3),that(3),phase,rw,fr,ft
pt=self%mesh%log2phys(cell,f)
rad=SQRT(pt(1)**2+pt(2)**2)
theta=ATAN2(pt(2),pt(1))
rhat=[COS(theta),SIN(theta),0.d0]
that=[-SIN(theta),COS(theta),0.d0]
SELECT CASE(self%field)
  CASE('n') ! Density
    val=linterp(self%vals(1,:),self%vals(4,:),self%npts,rad,1)
  CASE('t') ! Temperature
    val=linterp(self%vals(1,:),self%vals(5,:),self%npts,rad,1)
  CASE('b') ! Equilibrium magnetic field
    Bp=linterp(self%vals(1,:),self%vals(2,:),self%npts,rad,1)
    Bz=linterp(self%vals(1,:),self%vals(3,:),self%npts,rad,1)
    val=Bp*that+[0.d0,0.d0,Bz]
  CASE('d') ! Magnetic perturbation
    rad=rad/r_plasma
    rw=r_wall/r_plasma
    phase=self%m*theta-self%n*2.d0*pi*pt(3)/self%L
    IF(rad<=1.d0)THEN
      fr=rad**(self%m-1)
      ft=rad**(self%m-1)
    ELSE
      fr=(rw**(2*self%m)*rad**(-self%m-1)-rad**(self%m-1))/(rw**(2*self%m)-1.d0)
      ft=-(rw**(2*self%m)*rad**(-self%m-1)+rad**(self%m-1))/(rw**(2*self%m)-1.d0)
    END IF
    val=fr*COS(phase)*rhat-ft*SIN(phase)*that
  CASE DEFAULT
    CALL oft_abort('Unknown field component','cyl_pinch_apply',__FILE__)
END SELECT
END SUBROUTINE cyl_pinch_apply
!------------------------------------------------------------------------------
!> Sample magnetic field and velocity at the plasma edge
!!
!! Signals are saved to `b_probes.hist` and `v_probes.hist`
!------------------------------------------------------------------------------
SUBROUTINE cyl_probe_apply(self,sub_fields,t)
CLASS(cyl_probe), INTENT(inout) :: self
TYPE(xmhd_sub_fields), INTENT(inout) :: sub_fields !< MHD fields at current time
REAL(r8), INTENT(in) :: t !< Current time [s]
INTEGER(i4) :: i
REAL(r8) :: theta
REAL(r8), ALLOCATABLE :: pts(:,:)
CHARACTER(LEN=80) :: file_tmp
!---Setup probes on first call
IF(.NOT.self%initialized)THEN
  ALLOCATE(pts(3,self%npts))
  DO i=1,self%npts
    theta=(i-1)*2.d0*pi/REAL(self%npts,8)
    pts(:,i)=[r_plasma*COS(theta),r_plasma*SIN(theta),0.d0]
  END DO
  SELECT TYPE(this=>xmhd_ML_lagrange%current_level)
    CLASS IS(oft_fem_type)
      self%v_probes%mesh=>this%mesh
      self%b_probes%mesh=>this%mesh
  END SELECT
  !---Setup velocity probes
  ALLOCATE(self%Vfield)
  CALL self%v_probes%setup(self%npts,pts)
  file_tmp='v_probes.hist'
  CALL self%v_probes%setup_save(file_tmp,'V')
  self%v_probes%B=>self%Vfield
  !---Setup magnetic field probes
  ALLOCATE(self%Bfield)
  CALL self%b_probes%setup(self%npts,pts)
  file_tmp='b_probes.hist'
  CALL self%b_probes%setup_save(file_tmp,'B')
  self%b_probes%B=>self%Bfield
  self%initialized=.TRUE.
  DEALLOCATE(pts)
END IF
!---Update interpolators and sample
self%Bfield%u=>sub_fields%B
CALL self%Bfield%setup(xmhd_ML_hcurl_grad%current_level)
self%Vfield%u=>sub_fields%V
CALL self%Vfield%setup(xmhd_ML_lagrange%current_level)
CALL self%b_probes%save(t)
CALL self%v_probes%save(t)
END SUBROUTINE cyl_probe_apply
!------------------------------------------------------------------------------
!> Resistivity shaping function for the vacuum region
!!
!! The resistivity is increased by a factor of `eta_vac` outside the plasma,
!! using a linear ramp that starts `vac_start` outside the plasma radius and
!! has a width of `vac_width`.
!!
!! @result Resistivity scale factor at point
!------------------------------------------------------------------------------
FUNCTION vac_diff_shape(pt) RESULT(val)
REAL(r8), INTENT(in) :: pt(3) !< Location of evaluation
REAL(r8) :: val,rad
rad=SQRT(pt(1)**2+pt(2)**2)
val=1.d0+(eta_vac-1.d0)*MAX(MIN((rad-r_plasma-vac_start)/vac_width,1.d0),0.d0)
END FUNCTION vac_diff_shape
END MODULE cyl_pinch
!!\subsection doc_mug_cyl_kink_code_inc Main Program: Module Includes
!! The main program begins by including the required OFT modules and the helper module above.
PROGRAM MUG_cyl_pinch
!---Runtime
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre, create_bjacobi_pre, create_ilu_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_getmop, oft_lag_project, &
  oft_lag_vgetmop, oft_lag_vproject
!---H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, oft_hcurl_grad_rinterp, &
  oft_hcurl_grad_project, hcurl_grad_getmop, oft_hcurl_grad_gzerop
!---Physics
USE xmhd, ONLY: xmhd_lin_run, xmhd_plot, xmhd_minlev, xmhd_taxis, xmhd_sub_fields, &
  xmhd_ML_hcurl, xmhd_ML_H1, xmhd_ML_hcurl_grad, xmhd_ML_H1grad, xmhd_ML_lagrange, xmhd_ML_vlagrange, &
  res_shaping
!---Self
USE cyl_pinch, ONLY: cyl_pinch_init, cyl_probe, vac_diff_shape, r_plasma, eta_vac, vac_start, vac_width
IMPLICIT NONE
!!\subsection doc_mug_cyl_kink_code_vars Local Variables
!! Next we define the local variables needed to initialize our case and
!! run the time-dependent solve and post-processing.
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
TYPE(multigrid_mesh) :: mg_mesh
TYPE(xdmf_plot_file) :: plot_file
TYPE(xmhd_sub_fields) :: ic_fields,pert_fields
TYPE(oft_hcurl_grad_rinterp), TARGET :: bfield
TYPE(oft_hcurl_grad_gzerop) :: hcurl_grad_gzerop
TYPE(cyl_pinch_init), TARGET :: init_field
TYPE(cyl_probe) :: probes
REAL(r8), POINTER :: vals(:) => NULL()
REAL(r8), POINTER :: vec_vals(:,:) => NULL()
INTEGER(i4) :: i,ierr,io_unit
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = -1
INTEGER(i4) :: mm = 2
INTEGER(i4) :: nn = 1
REAL(r8) :: a = 0.1d0
REAL(r8) :: db = 1.d-5
LOGICAL :: pm = .FALSE.
LOGICAL :: plot_run = .FALSE.
CHARACTER(LEN=80) :: eq_file = 'none'
NAMELIST/cyl_pinch_options/order,minlev,mm,nn,a,db,eq_file,eta_vac,vac_start,vac_width,plot_run,pm
!!\subsection doc_mug_cyl_kink_code_setup OFT Initialization
!! See \ref doc_api_ex1 for a detailed description of calls in this section.
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,cyl_pinch_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)CALL oft_abort('Error reading "cyl_pinch_options" group in input file','main',__FILE__)
IF(TRIM(eq_file)=='none')CALL oft_abort('No equilibrium file specified ("eq_file")','main',__FILE__)
r_plasma=a
!------------------------------------------------------------------------------
! Setup grid
!------------------------------------------------------------------------------
CALL multigrid_construct(mg_mesh,[0.d0,1.d0,0.d0])
!------------------------------------------------------------------------------
! Build FE structures
!------------------------------------------------------------------------------
!--- Lagrange
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
!--- Grad(H^1) subspace
ALLOCATE(xmhd_ML_H1)
CALL oft_h1_setup(mg_mesh,order+1,xmhd_ML_H1,minlev=minlev)
CALL h1_setup_interp(xmhd_ML_H1)
!--- H(Curl) subspace
ALLOCATE(xmhd_ML_hcurl)
CALL oft_hcurl_setup(mg_mesh,order,xmhd_ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(xmhd_ML_hcurl)
!--- Full H(Curl) space
ALLOCATE(xmhd_ML_hcurl_grad,xmhd_ML_H1grad)
CALL oft_hcurl_grad_setup(xmhd_ML_hcurl,xmhd_ML_H1,xmhd_ML_hcurl_grad,xmhd_ML_H1grad,minlev=minlev)
CALL hcurl_grad_setup_interp(xmhd_ML_hcurl_grad,xmhd_ML_H1)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>xmhd_ML_hcurl_grad
!---Load equilibrium profile
init_field%filename=eq_file
init_field%m=mm
init_field%n=nn
CALL init_field%setup(mg_mesh%mesh)
!!\subsection doc_mug_cyl_kink_code_plot Perform post-processing
!!
!! To visualize the solution fields once a simulation has completed, the \ref xmhd::xmhd_plot
!! "xmhd_plot" subroutine is used. This subroutine steps through the restart files
!! produced by \ref xmhd::xmhd_lin_run "xmhd_lin_run" and generates plot files, and
!! samples the `cyl_probe` signals at evenly spaced points in time as specified in the
!! input file, see \ref doc_mug_cyl_kink_post.
IF(plot_run)THEN
  CALL xmhd_plot(probes)
  CALL oft_finalize
END IF
!!\subsection doc_mug_cyl_kink_code_ic Computing Initial Conditions
!!
!! The equilibrium density, temperature and magnetic field, and the magnetic perturbation, are defined by
!! `cyl_pinch_init` and are mapped onto the finite element spaces by \f$ L_2 \f$ projection,
!! which requires the solution of a mass matrix system for each field. The initial fields are also
!! saved to the plot file for visualization.
CALL plot_file%setup("cyl_pinch")
CALL mg_mesh%mesh%setup_io(plot_file,order)
!---------------------------------------------------------------------------
! Equilibrium density and temperature (scalar Lagrange)
!---------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_getmop(xmhd_ML_lagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
CALL create_diag_pre(minv%pre)
minv%A=>mop
minv%its=-2
CALL xmhd_ML_lagrange%vec_create(u)
CALL xmhd_ML_lagrange%vec_create(v)
init_field%field='n'
CALL oft_lag_project(xmhd_ML_lagrange%current_level,init_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL xmhd_ML_lagrange%vec_create(ic_fields%Ne)
CALL ic_fields%Ne%add(0.d0,1.d0,u)
CALL ic_fields%Ne%get_local(vals)
CALL mg_mesh%mesh%save_vertex_scalar(vals,plot_file,'N0')
init_field%field='t'
CALL oft_lag_project(xmhd_ML_lagrange%current_level,init_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL xmhd_ML_lagrange%vec_create(ic_fields%Ti)
CALL ic_fields%Ti%add(0.d0,1.d0,u)
CALL ic_fields%Ti%get_local(vals)
CALL mg_mesh%mesh%save_vertex_scalar(vals,plot_file,'T0')
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%delete(.TRUE.)
DEALLOCATE(minv)
!---------------------------------------------------------------------------
! Equilibrium and perturbed magnetic fields (H(Curl) + Grad(H^1))
!---------------------------------------------------------------------------
NULLIFY(mop)
CALL hcurl_grad_getmop(xmhd_ML_hcurl_grad%current_level,mop,"none")
CALL create_cg_solver(minv)
CALL create_bjacobi_pre(minv%pre,-1)
DEALLOCATE(minv%pre%pre)
CALL create_ilu_pre(minv%pre%pre)
minv%A=>mop
minv%its=-2
CALL xmhd_ML_hcurl_grad%vec_create(u)
CALL xmhd_ML_hcurl_grad%vec_create(v)
!---Equilibrium field
init_field%field='b'
CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,init_field,v)
CALL hcurl_grad_gzerop%apply(v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL xmhd_ML_hcurl_grad%vec_create(ic_fields%B)
CALL ic_fields%B%add(0.d0,1.d0,u)
!---Perturbation field
init_field%field='d'
CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,init_field,v)
CALL hcurl_grad_gzerop%apply(v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL xmhd_ML_hcurl_grad%vec_create(pert_fields%B)
CALL pert_fields%B%add(0.d0,db,u)
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%delete(.TRUE.)
DEALLOCATE(minv)
!---------------------------------------------------------------------------
! Project magnetic fields to vector Lagrange basis for plotting
!---------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
CALL create_diag_pre(minv%pre)
minv%A=>mop
minv%its=-2
CALL xmhd_ML_vlagrange%vec_create(u)
CALL xmhd_ML_vlagrange%vec_create(v)
ALLOCATE(vec_vals(3,xmhd_ML_lagrange%current_level%ne))
!---Equilibrium field
bfield%u=>ic_fields%B
CALL bfield%setup(xmhd_ML_hcurl_grad%current_level)
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,bfield,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
DO i=1,3
  vals=>vec_vals(i,:)
  CALL u%get_local(vals,i)
END DO
CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'B0')
!---Perturbation field
bfield%u=>pert_fields%B
CALL bfield%setup(xmhd_ML_hcurl_grad%current_level)
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,bfield,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
DO i=1,3
  vals=>vec_vals(i,:)
  CALL u%get_local(vals,i)
END DO
CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'dB')
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop,vec_vals)
CALL minv%delete(.TRUE.)
DEALLOCATE(minv)
!!\subsection doc_mug_cyl_kink_code_run Run Simulation
!!
!! Finally, the remaining equilibrium (zero flow) and perturbation fields are created and the
!! linear simulation is run using \ref xmhd::xmhd_lin_run "xmhd_lin_run". Before calling the
!! solver the vacuum resistivity profile is set by associating the `res_shaping` function pointer in the \ref xmhd "xmhd"
!! module with `vac_diff_shape`.
!!
!! The vacuum region between the plasma and the wall is modeled in two parts. Just outside the plasma, a
!! low-density fluid region with a resistivity enhanced by `eta_vac` allows the plasma boundary to move freely
!! and the vacuum field to diffuse rapidly compared to the growth time of the mode. Outside this layer, mesh
!! region 2 is treated as a solid region (see \ref doc_mug_cyl_kink_input_solver), where no flow is allowed and
!! the same enhanced resistivity is applied.
CALL xmhd_ML_vlagrange%vec_create(ic_fields%V)
CALL xmhd_ML_vlagrange%vec_create(pert_fields%V)
CALL xmhd_ML_lagrange%vec_create(pert_fields%Ne)
CALL xmhd_ML_lagrange%vec_create(pert_fields%Ti)
xmhd_minlev=minlev
xmhd_taxis=3
res_shaping=>vac_diff_shape
oft_env%pm=pm
CALL xmhd_lin_run(ic_fields,pert_fields)
!---Finalize enviroment
CALL oft_finalize
END PROGRAM MUG_cyl_pinch
! STOP SOURCE
!!
!!\section doc_mug_cyl_kink_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment.
!! Some annotation of the options is provided inline below, for more information on the
!! available options in the \c xmhd_options group see \ref xmhd::xmhd_lin_run "xmhd_lin_run".
!! As the density and temperature are advanced (`xmhd_adv_den=T` and `xmhd_adv_temp=T`) the perturbed pressure
!! is computed self-consistently, which is required to correctly linearize the equilibrium force balance
!! \f$ \nabla p = \textbf{J} \times \textbf{B} \f$, which is not small compared to the drive for this instability.
!!
! OFT_DOC_INCLUDE: oft.in
!!
!!\subsection doc_mug_cyl_kink_input_solver Solver and region specification
!!
!! The `oft_in.xml` file below defines mesh region 2 as a solid region (`type=2`) with the same base resistivity
!! as the plasma, and specifies the preconditioner for the time-dependent solve.
!!
!! The time step used here (\f$ \Delta t \approx \tau_A / 4 \f$) is much larger than the fast magnetosonic
!! CFL limit, especially in the low-density vacuum region, which couples the solution strongly over long distances along the
!! cylinder within each time step. For this case, block-Jacobi based preconditioners (with or without polynomial multigrid)
!! converge poorly, so a distributed direct solve using [SuperLU-DIST](https://github.com/xiaoyeli/superlu_dist) is used
!! as the preconditioner instead. As the case is linear and the time step is fixed, the factorization is only computed once
!! (`xmhd_prefreq` is larger than `nsteps`) and each time step requires only one or two GMRES iterations.
!!
!! @note This example requires OFT to be built with SuperLU-DIST (eg. `--build_superlu_dist=1` with `build_libs.py`).
!! Using 2 MPI processes (with `OMP_NUM_THREADS=1`), the factorization requires approximately 30 GB of memory and the
!! simulation takes approximately 25 minutes, including approximately 8 minutes for the factorization.
!!
!! The simulation can be run by specifying both the FORTRAN input and XML input files to the executable as below.
!!
!!\verbatim
!!~$ mpirun -np 2 ./MUG_cyl_pinch oft.in oft_in.xml
!!\endverbatim
!!
! OFT_DOC_INCLUDE: oft_in.xml
!!
!!\section doc_mug_cyl_kink_eq Equilibrium Creation
!!
!! The equilibrium profile file `kink_eq.dat` is provided with this example and was created using the python script
!! `create_eq.py`, which also prints the growth rate predicted by the flat-current model for the same parameters.
!! The current and density profiles are flat inside the plasma and smoothed at the edge with a
!! \f$ \tanh \f$ of width 4 mm. The density outside the plasma is 1% of the plasma density. The total current is chosen to give
!! the desired \f$ q_a \f$, so the field outside the plasma matches the flat-current model. The temperature profile
!! (\f$ T_e = T_i \f$) is set to satisfy radial force balance, giving a central temperature of 250 eV and \f$ \beta \approx 3\% \f$,
!! with a temperature of 5 eV outside the plasma.
!! Equilibria for other parameters can be generated using the command line options of the script (see `--help`).
!!
!!\verbatim
!!~$ python create_eq.py --qa=1.6
!!\endverbatim
!!
!!\section doc_mug_cyl_kink_post Post-Processing options
!!
!!\subsection doc_mug_cyl_kink_post_growth Growth rate
!!
!! During the simulation the evolution of some global quantities, including the kinetic and perturbed magnetic energies,
!! is written to the history file `xmhd.hist`. The python script `plot_growth.py` reads this file, fits the growth rate
!! of the kinetic energy (\f$ E_K \propto e^{2 \gamma t} \f$) over the second half of the simulation and compares it to the
!! growth rate from the flat-current model.
!!
!!\verbatim
!!~$ python plot_growth.py --hist_file=xmhd.hist
!!\endverbatim
!!
!! After a short initial transient the perturbation grows exponentially with
!! \f$ \gamma = 9.28 \times 10^5 \f$ s\f$ ^{-1} \f$ (\f$ \gamma \tau_A = 0.376 \f$), which is 11% below the flat-current model.
!! The remaining difference is expected, as the MUG model includes physics and numerical choices that are not part of the
!! flat-current model, most notably a resistive, finite-density vacuum region in place of a true vacuum (increasing the vacuum density
!! from 1% to 5% of the plasma density reduces the growth rate by approximately 8%), a smoothed current edge and finite aspect ratio.
!!
!! \image html MUG_cyl_kink_ex-Growth.png "Time evolution of kinetic and perturbed magnetic energy (top) and instantaneous growth rate (bottom), compared to the flat-current model."
!!
!!\subsection doc_mug_cyl_kink_post_plot Creating plot files
!!
!! To generate 3D plots and sample the probe signals, a plot run can be performed by setting `plot_run=T` in the
!! `cyl_pinch_options` input group, which calls \ref xmhd::xmhd_plot "xmhd_plot". With this option additional run time
!! options are available in the `xmhd_plot_options` group that control how restart files are sampled for plotting (see the input
!! file above, where `rst_end` must be set to the index of the last restart file). The `cyl_probe` object records the perturbed magnetic field and velocity at 180 points around the plasma edge
!! (\f$ r = a \f$, \f$ z = 0 \f$) to the files `b_probes.hist` and `v_probes.hist`.
!!
!!\verbatim
!!~$ mpirun -np 2 ./MUG_cyl_pinch oft.in oft_in.xml
!!\endverbatim
!!
!! Once the post-processing run is complete `bin/OFT_build_xdmf.py` can be used to generate `*.xmf` files that can be loaded by
!! [VisIt](https://visit-dav.github.io/visit-website/index.html), [ParaView](https://www.paraview.org/), or other visualization programs.
!! The python script `plot_mode.py` uses the plot files and probe signals to show the structure of the mode in the
!! \f$ z = 0 \f$ plane and around the plasma edge, which has the expected \f$ m = 2 \f$ structure.
!!
!!\verbatim
!!~$ python plot_mode.py
!!\endverbatim
!!
!! \image html MUG_cyl_kink_ex-Mode.png "Perturbed radial magnetic field (color) and velocity inside the plasma (arrows) in the z=0 plane (left) and perturbed radial field around the plasma edge (right) at the end of the simulation."
!!
!!\section doc_mug_cyl_kink_mesh Mesh Creation
!! A mesh file `MUG-cyl_pinch.h5` is provided with this example. It was created using the following
!! [CUBIT](https://cubit.sandia.gov/) script, which generates one half of the periodic cylinder
!! (\f$ 0 \le z \le L/2 \f$) using hexahedral elements.
!!
! OFT_DOC_INCLUDE: cyl_pinch.jou
!!
!! Once complete, the mesh should be converted into the native mesh format using the `OFT_convert_cubit.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python OFT_convert_cubit.py --in_file=MUG-cyl_pinch.g --periodic_nodeset=1
!!\endverbatim
!!
!! @note The option `--periodic_nodeset=1` indicates that the surface corresponding to the first, and only, nodeset
!! should be treated as periodic following reflection. The full periodic mesh is created when the mesh is loaded
!! by reflecting this mesh across the \f$ z = 0 \f$ plane (`reflect=T` and `ref_periodic=T` in the
!! `native_mesh_options` group).
