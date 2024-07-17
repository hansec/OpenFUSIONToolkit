!!MUG Example: Accretion disk    {#doc_mug_accretion_ex}
!!============================
!!
!![TOC]
!!
!!\section doc_mug_accretion_ex_code_helper Helper module
!!
!! Need docs
! START SOURCE
MODULE accretion_helpers
USE oft_base
USE oft_io, ONLY: oft_bin_file
USE oft_mesh_type, ONLY: mesh
USE fem_utils, ONLY: fem_interp
USE mhd_utils, ONLY: mu0
USE oft_la_base, ONLY: oft_vector
USE oft_lag_fields, ONLY: oft_lag_vcreate, oft_lag_create
USE oft_lag_operators, ONLY: oft_lag_rinterp
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: h1curl_zerob, oft_h1_rinterp
USE xmhd, ONLY: oft_xmhd_driver, xmhd_sub_fields
IMPLICIT NONE
!---------------------------------------------------------------------------
! Field interpolation object for intial conditions
!---------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: accretion_interp
  REAL(8) :: rmin = 0.5d0
  REAL(8) :: rmax = 1.d0
  CHARACTER(LEN=2) :: field = 'n' ! Field component to initialize
CONTAINS
  PROCEDURE :: interp => accretion_interp_apply ! Reconstruct field
END TYPE accretion_interp
CONTAINS
!!\subsection doc_mug_accretion_ex_code_helper_interp Initial condition field interpolator
!!
!! Need docs
SUBROUTINE accretion_interp_apply(self,cell,f,gop,val)
CLASS(accretion_interp), INTENT(inout) :: self ! Needs docs
INTEGER(4), INTENT(in) :: cell ! Needs docs
REAL(8), INTENT(in) :: f(:) ! Needs docs
REAL(8), INTENT(in) :: gop(3,4) ! Needs docs
REAL(8), INTENT(out) :: val(:) ! Needs docs
REAL(8) :: pt(3),theta,rad
! Map logical positionto physical coordinates
pt = mesh%log2phys(cell,f)
theta = ATAN2(pt(2),pt(1))
rad = SQRT(pt(1)**2 + pt(2)**2)
! Return requested field evaluated at "(cell,f) -> pt"
SELECT CASE(TRIM(self%field))
  CASE('n0') ! Density
    val = 1.d0
  CASE('b0') ! Equilibrium magnetic field
    val = [0.d0,0.d0,1.d0]
  CASE('V0') ! Equilibrium velocity field
    val = [-SIN(theta),COS(theta),0.d0]*(1.d0/rad)
  ! CASE('db') ! Perturbed magnetic field
  !   Bper(1) = -pi*COS(2.d0*pi*pt(1)/self%Lx)*SIN(pi*pt(3)/self%Lz)/self%Lz
  !   Bper(2) = 0.d0
  !   Bper(3) = 2.d0*pi*SIN(2.d0*pi*pt(1)/self%Lx)*COS(pi*pt(3)/self%Lz)/self%Lx
  !   val = Bper
  CASE DEFAULT
    CALL oft_abort('Unknown field component','accretion_interp_apply',__FILE__)
END SELECT
END SUBROUTINE accretion_interp_apply
END MODULE accretion_helpers
!!\section doc_ex6_code_driver Driver program
!!
!! Need docs
PROGRAM MUG_slab_recon
USE oft_base
!--Grid
USE oft_mesh_type, ONLY: mesh, rgrnd
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre, create_bjacobi_pre, &
  create_ilu_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_fields, ONLY: oft_lag_vcreate, oft_lag_create
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, &
  oft_lag_getmop, oft_lag_project
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_setup_interp
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: h1_setup_interp, h1_getmop, oft_h1_project, h1grad_zerop, &
  oft_h1_rinterp, oft_h1_cinterp
!---Physics
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, temp_floor, den_floor, den_scale, &
  xmhd_sub_fields, xmhd_lin_run
!---Self
USE accretion_helpers, ONLY: accretion_interp
IMPLICIT NONE
!---Fixed scaling parameters
REAL(8), PARAMETER :: N0 = 1.d18
REAL(8), PARAMETER :: T0 = 10.d0
REAL(8), PARAMETER :: B0 = 0.05d0
REAL(8), PARAMETER :: V0 = 1.d3
!---Mass matrix solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
!---Local variables
TYPE(xmhd_sub_fields) :: ic_fields,pert_fields
TYPE(oft_h1_rinterp), TARGET :: bfield
TYPE(oft_h1_cinterp), TARGET :: jfield
TYPE(accretion_interp), TARGET :: accretion_field
REAL(8), POINTER :: vals(:),vec_vals(:,:)
INTEGER(4) :: i,ierr,io_unit
!---Input file options
INTEGER(4) :: order = 2
INTEGER(4) :: minlev = 1
REAL(8) :: db = 1.d-1
LOGICAL :: pm = .FALSE.
LOGICAL :: linear = .FALSE.
LOGICAL :: plot_run = .FALSE.
LOGICAL :: view_ic = .FALSE.
NAMELIST/accretion_options/order,minlev,linear,db,plot_run,view_ic,pm
!!\subsection doc_ex6_code_driver_init Grid and FE setup
!!
!! Need docs
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,accretion_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
rgrnd=(/0.d0,0.d0,1.d0/)
CALL multigrid_construct
!---Setup I/0
IF(view_ic.OR.plot_run)CALL mesh%setup_io(order)
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!---Lagrange
CALL oft_lag_setup(order,minlev)
CALL lag_setup_interp
!---H1(Curl) subspace
CALL oft_hcurl_setup(order,minlev)
CALL hcurl_setup_interp
!---H1(Grad) subspace
CALL oft_h0_setup(order+1,minlev)
CALL h0_setup_interp
!---H1 full space
CALL oft_h1_setup(order,minlev)
CALL h1_setup_interp
!---------------------------------------------------------------------------
! If plot run, just run and finish
!---------------------------------------------------------------------------
IF(plot_run)THEN
  CALL xmhd_plot()
  CALL oft_finalize
END IF
!!\subsection doc_ex6_code_driver_ic Initial conditions
!!
!! Need docs
!---Set constant initial temperature
CALL oft_lag_create(ic_fields%Ti)
CALL ic_fields%Ti%set(T0)
!---------------------------------------------------------------------------
! Set intial velocity from analytic definition
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_lag_vgetmop(mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL oft_lag_vcreate(u)
CALL oft_lag_vcreate(v)
!---Project onto scalar Lagrange basis
accretion_field%field='V0'
CALL oft_lag_vproject(accretion_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
!---Create density field and set values
CALL oft_lag_vcreate(ic_fields%V)
CALL ic_fields%V%add(0.d0,1.d0,u)
CALL ic_fields%V%scale(V0) ! Scale to desired value
!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)
!---------------------------------------------------------------------------
! Set intial density from analytic definition
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_lag_getmop(mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL oft_lag_create(u)
CALL oft_lag_create(v)
!---Project onto scalar Lagrange basis
accretion_field%field='n0'
CALL oft_lag_project(accretion_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
!---Create density field and set values
CALL oft_lag_create(ic_fields%Ne)
CALL ic_fields%Ne%add(0.d0,1.d0,u)
CALL ic_fields%Ne%scale(N0) ! Scale to desired value
!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)
!---------------------------------------------------------------------------
! Set intial magnetic field from analytic definition
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL h1_getmop(mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
! CALL create_diag_pre(minv%pre) ! Setup Preconditioner
CALL create_bjacobi_pre(minv%pre,-1)
DEALLOCATE(minv%pre%pre)
CALL create_ilu_pre(minv%pre%pre)
!---Create fields for solver
CALL oft_h1_create(u)
CALL oft_h1_create(v)
!---Project onto vector H(Curl) basis
accretion_field%field='b0'
CALL oft_h1_project(accretion_field,v)
CALL h1grad_zerop(v) ! Zero out redundant vertex degrees of freedom
CALL u%set(0.d0)
CALL minv%apply(u,v)
!---Create magnetic field and set values
CALL oft_h1_create(ic_fields%B)
CALL ic_fields%B%add(0.d0,1.d0,u)
CALL ic_fields%B%scale(B0) ! Scale to desired value
!---Compute perturbed magnetic field
CALL oft_h1_create(pert_fields%B)
! accretion_field%field='db'
! CALL oft_h1_project(accretion_field,v)
! CALL h1grad_zerop(v) ! Zero out redundant vertex degrees of freedom
! CALL u%set(0.d0)
! CALL minv%apply(u,v)
! CALL pert_fields%B%add(0.d0,1.d0,u)
! CALL pert_fields%B%scale(B0*db) ! Scale to desired value
!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre%pre)
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)
!---------------------------------------------------------------------------
! Save initial conditions if desired
!---------------------------------------------------------------------------
IF(view_ic)THEN
  !---Set up output arrays for plotting
  NULLIFY(vals)
  ALLOCATE(vec_vals(3,ic_fields%Ne%n))
  !---Plot density
  CALL ic_fields%Ne%get_local(vals) ! Fetch local values
  CALL mesh%save_vertex_scalar(vals,'N0') ! Add field to plotting file
  !---Plot temperature
  CALL ic_fields%Ti%get_local(vals) ! Fetch local values
  CALL mesh%save_vertex_scalar(vals,'T0') ! Add field to plotting file
  !---Plot velocity
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL ic_fields%V%get_local(vals,i)
  END DO
  CALL mesh%save_vertex_vector(vec_vals,'V0') ! Add field to plotting file
  !---------------------------------------------------------------------------
  ! Project B and J for plotting
  !---------------------------------------------------------------------------
  !---Generate mass matrix
  NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
  CALL oft_lag_vgetmop(mop,"none") ! Construct mass matrix with "none" BC
  !---Setup linear solver
  CALL create_cg_solver(minv)
  minv%A=>mop ! Set matrix to be solved
  minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
  CALL create_diag_pre(minv%pre) ! Setup Preconditioner
  !---Create fields for solver
  CALL oft_lag_vcreate(u)
  CALL oft_lag_vcreate(v)
  !---Project B onto vector Lagrange basis
  bfield%u=>ic_fields%B
  CALL bfield%setup
  CALL oft_lag_vproject(bfield,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  !---Retrieve and save projected magnetic field
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL u%get_local(vals,i)
  END DO
  CALL mesh%save_vertex_vector(vec_vals,'B0') ! Add field to plotting file
  !---Project B onto vector Lagrange basis
  bfield%u=>pert_fields%B
  CALL bfield%setup
  CALL oft_lag_vproject(bfield,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  !---Retrieve and save projected magnetic field
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL u%get_local(vals,i)
  END DO
  CALL mesh%save_vertex_vector(vec_vals,'dB') ! Add field to plotting file
  !---Project J onto vector Lagrange basis
  jfield%u=>ic_fields%B
  CALL jfield%setup
  CALL oft_lag_vproject(jfield,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  !---Retrieve and save projected current density
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL u%get_local(vals,i)
  END DO
  CALL mesh%save_vertex_vector(vec_vals,'J0') ! Add field to plotting file
  !---Cleanup objects used for projection
  CALL u%delete ! Destroy LHS vector
  CALL v%delete ! Destroy RHS vector
  CALL mop%delete ! Destroy mass matrix
  DEALLOCATE(u,v,mop) ! Deallocate objects
  CALL minv%pre%delete ! Destroy preconditioner
  DEALLOCATE(minv%pre)
  CALL minv%delete ! Destroy solver
  DEALLOCATE(minv)
  !---Finalize enviroment
  CALL oft_finalize
END IF
!!\subsection doc_ex6_code_driver_run Run simulation
!!
!! Need docs
xmhd_minlev=minlev  ! Set minimum level for multigrid preconditioning
den_scale=N0        ! Set density scale
oft_env%pm=pm       ! Show linear iteration progress?
IF(linear)THEN
  CALL oft_lag_vcreate(pert_fields%V)
  CALL oft_lag_create(pert_fields%Ti)
  CALL oft_lag_create(pert_fields%Ne)
  CALL xmhd_lin_run(ic_fields,pert_fields)
ELSE
  CALL ic_fields%B%add(1.d0,1.d0,pert_fields%B)
  !---Run simulation
  temp_floor=T0*1.d-2 ! Set temperature floor
  den_floor=N0*1.d-2  ! Set density floor
  CALL xmhd_run(ic_fields)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM MUG_slab_recon
! STOP SOURCE
!!
!!\section doc_mug_accretion_ex_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment.
!! As with \ref doc_ex6 "Example 6" this example should only be run with multiple processes.
!! Some annotation of the options is provided inline below, for more information on the
!! available options in the \c xmhd_options group see \ref xmhd::xmhd_plot "xmhd_plot".
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='test'
!! cad_type=0
!! nlevels=2
!! nbase=1
!! grid_order=2
!! fix_boundary=T
!!/
!!
!!&native_mesh_options
!! filename='cyl_accretion.h5'
!! reflect=T
!!/
!!
!!&accretion_options
!! order=2           ! FE order
!! minlev=2          ! Minimum level for MG preconditioning
!! view_ic=F         ! View initial conditions but do not run simulation
!! plot_run=F        ! Run plotting instead of simulation
!!/
!!
!!&xmhd_options
!! mu_ion=1.         ! Ion mass (atomic units)
!! xmhd_ohmic=T      ! Include Ohmic heating
!! xmhd_visc_heat=T  ! Include viscous heating
!! bbc='bc'          ! Perfectly-conducting BC for B-field
!! vbc='all'         ! Zero-flow BC for velocity
!! nbc='n'           ! Neumann BC for density
!! tbc='n'           ! Neumann BC for temperature
!! dt=4.e-7          ! Maximum time step
!! eta=978.7         ! Constant resistivity
!! visc_type='iso'   ! Use isotropic viscosity tensor
!! nu_par=9877.0     ! Fluid viscosity
!! d_dens=10.        ! Density diffusion
!! kappa_par=3914.9  ! Parallel thermal conduction (fixed)
!! kappa_perp=3914.9 ! Perpendicular thermal conduction (fixed)
!! nsteps=2000       ! Number of time steps to take
!! rst_ind=0         ! Index of file to restart from (0 -> use subroutine arguments)
!! rst_freq=10       ! Restart file frequency
!! xmhd_mfnk=T       ! Use matrix-free method
!! lin_tol=1.E-8     ! Linear solver tolerance
!! nl_tol=1.E-6      ! Non-linear solver tolerance
!! nu_xmhd=0,20,5    ! Number of smoother iterations for default preconditioner
!! ittarget=30       ! Target for # of linear iterations per time step
!! xmhd_prefreq=20   ! Preconditioner update frequency
!!/
!!\endverbatim
!!
!!\subsection doc_mug_accretion_ex_input_solver Solver specification
!!
!! Time dependent MHD solvers are accelerated significantly by the use of
!! a more sophisticated preconditioner than the default method. Below is
!! an example `oft_in.xml` file that constructs an appropriate ILU(0) preconditioner.
!! Currently, this preconditioner method is the suggest starting preconditioner for all
!! time-dependent MHD solves.
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ ./MUG_accretion oft.in oft_in.xml
!!\endverbatim
!!
!!```xml
!!<oft>
!!  <xmhd>
!!    <pre type="gmres">
!!      <its>8</its>
!!      <nrits>8</nrits>
!!      <pre type="block_jacobi">
!!        <nlocal>-1</nlocal>
!!        <solver type="ilu"></solver>
!!      </pre>
!!    </pre>
!!  </xmhd>
!!</oft>
!!```
!!
!!\subsection doc_mug_accretion_ex_input_plot Post-Processing options
!!
!! When running the code for post-processing additional run time options are available.
!!
!!\verbatim
!!&xmhd_plot_options
!! t0=1.E-8
!! dt=1.E-6
!! rst_start=0
!! rst_end=1000
!!/
!!\endverbatim
!!
!! \image html example_gem-result.png "Resulting current distribution for the first eigenmode"
!!
!!\section doc_mug_accretion_ex_mesh Mesh Creation
!! A mesh file `cyl_accretion.h5` is provided with this example. Instructions to generate your
!! own mesh for the geometry using [CUBIT](https://cubit.sandia.gov/).
!!
!!\subsection doc_mug_accretion_ex_cubit Meshing with CUBIT
!!
!! A suitable mesh for this example, with radius of 1m and height of 2m, can be created using
!! the CUBIT script below.
!!
!!\verbatim
!!reset
!!
!!create Cylinder height 1 radius 1
!!create Cylinder height 2 radius 0.5
!!subtract volume 2 from volume 1
!!move Volume 1  x 0 y 0 z 0.5 include_merged
!!
!!nodeset 1 add surface 9
!!
!!volume 1 scheme Tetmesh
!!set tetmesher interior points on
!!set tetmesher optimize level 3 optimize overconstrained  off sliver  off
!!set tetmesher boundary recovery  off
!!volume 1 size .1
!!mesh volume 1
!!
!!set duplicate block elements off
!!block 1 add volume 1 
!!block 1 element type tetra10
!!
!!set large exodus file on
!!export Genesis  "cyl_accretion.g" overwrite block 1
!!\endverbatim
!!
!! Once complete the mesh should be converted into the native mesh format using the `convert_cubit.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python convert_cubit.py --in_file=cyl_accretion.g
!!\endverbatim
