!!MUG Example 2: Spheromak Heating    {#doc_mhd_ex2}
!!============================
!!
!![TOC]
!!
!!This example demonstrates the use of the \ref xmhd "extended MHD" module in the Open FUSION Toolkit (OFT). In
!!this example self-heating of a spheromak in a unit cylinder will be simulated. This process
!!provides a simple test case illustrating the basic ascpects of an extended MHD simulation,
!!including: 1) Temperature dependent resistivity, 2) Anisotropic thermal conduction and 3)
!!Ohmic and Viscous heating.
!!
!!The dynamics in this test case will be prdominetly limited to heating. However, if the
!!simulation is run long enough evolution of the equilibrium profile will be observed and
!!eventual instability due to current peaking will occur.
!!
!!\section doc_ex6_cubit Mesh Creation with CUBIT
!!
!!A suitable mesh for this example, with radius of 1m and height of 1m, can be created using
!!the CUBIT script below.
!!
!!\verbatim
!!reset
!!create Cylinder height 1 radius 1
!!volume 1 scheme Tetmesh
!!set tetmesher interior points on
!!set tetmesher optimize level 3 optimize overconstrained  off sliver  off
!!set tetmesher boundary recovery  off
!!volume 1 size .15
!!mesh volume 1
!!refine parallel fileroot 'cyl' overwrite no_execute
!!\endverbatim
!!
!!\section doc_ex6_code Code Walk Through
!!
!!The code consists of three basic sections, required imports and variable definitions,
!!finite element setup, and system creation and solution.
!!
!!\subsection doc_ex6_code_inc Module Includes
! START SOURCE
PROGRAM xmhd_cyl
!---Runtime
USE oft_base
!---Grid
USE oft_mesh_type, ONLY: mesh, rgrnd
USE multigrid_build, ONLY: multigrid_construct, multigrid_add_quad
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels, oft_lag_set_level
USE oft_lag_fields, ONLY: oft_lag_vcreate, oft_lag_create
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_level, oft_hcurl_nlevels
USE oft_hcurl_operators, ONLY: hcurl_setup_interp, hcurl_mloptions, hcurl_zerob
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_setup_interp, oft_h0_getlop, h0_zerogrnd
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: oft_h1_divout, h1_zeroi, h1_mc, h1curl_zerob, h1_setup_interp
!---Physics
USE taylor, ONLY: taylor_hmodes, taylor_minlev, taylor_hffa, taylor_hlam
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, xmhd_taxis, vel_scale, den_scale, &
  den_floor, temp_floor, xmhd_sub_fields
IMPLICIT NONE
!!\section doc_ex6_code_vars Local Variables
!---H1 divergence cleaner
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_h1_divout) :: divout
CLASS(oft_matrix), POINTER :: lop => NULL()
!---Local variables
INTEGER(i4) :: ierr,io_unit
REAL(r8), POINTER, DIMENSION(:) :: tmp => NULL()
TYPE(xmhd_sub_fields) :: ic_fields
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: b0_scale = 1.E-1_r8
REAL(r8) :: n0 = 1.d19
REAL(r8) :: t0 = 6.d0
LOGICAL :: plot_run=.FALSE.
NAMELIST/cyl_options/order,minlev,b0_scale,n0,t0,plot_run
!!\section doc_ex6_code_setup OFT Initialization
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,cyl_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
rgrnd=(/2.d0,0.d0,0.d0/)
CALL multigrid_construct
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!---Lagrange
CALL oft_lag_setup(order)
CALL lag_setup_interp
CALL lag_mloptions
!---H1(Curl) subspace
CALL oft_hcurl_setup(order)
CALL hcurl_setup_interp
CALL hcurl_mloptions
!---H1(Grad) subspace
CALL oft_h0_setup(order+1)
CALL h0_setup_interp
!---H1 full space
CALL oft_h1_setup(order)
CALL h1_setup_interp
!!\section doc_ex6_code_taylor Computing Initial Conditions
!!
!! For this simulation we only need the spheromak mode, which is the lowest
!! force-free eignstate in this geometry. As a result the initial condition
!! is stable to all types of mode activity.
taylor_minlev=minlev
CALL taylor_hmodes(1)
CALL oft_lag_set_level(oft_lagrange_nlevels)
!!\subsection doc_ex6_code_taylor_gauge Setting Magnetic Boundary Conditions
!!
!! As in \ref doc_ex5 "Example 5" we must transform the gauge of the Taylor
!! state solution to the appropriate magnetic field BCs. For more information
!! on this see the description in \ref doc_ex5_code_taylor_gauge "Setting Magnetic Boundary Conditions"
!! of Example 5.
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
NULLIFY(lop)
CALL oft_h0_getlop(lop,"grnd")
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-2
CALL create_diag_pre(linv%pre) ! Setup Preconditioner
divout%solver=>linv
divout%bc=>h0_zerogrnd
divout%keep_boundary=.TRUE.
!---------------------------------------------------------------------------
! Setup initial conditions
!---------------------------------------------------------------------------
CALL oft_h1_create(ic_fields%B)
CALL taylor_hffa(1,oft_hcurl_level)%f%get_local(tmp)
CALL ic_fields%B%restore_local(tmp,1)
CALL divout%apply(ic_fields%B)
!!\subsection doc_ex6_code_ic Set Initial Conditions
!!
!! Now we set initial conditions for the simulation using the computed taylor
!! state, flat temperature and density profiles and zero initial velocity. The
!! lower bound for density (\ref xmhd::den_floor "den_floor") and temperature
!! (\ref xmhd::temp_floor "temp_floor") in the simulation are also set. These
!! variables act to prevent very low density and negative temperature regions
!! in the simulation. The scale factors for the velocity (\ref xmhd::vel_scale
!! "vel_scale") and density (\ref xmhd::den_scale "den_scale") evolution
!! equations are also set. These variables are used to scale the corresponding
!! rows in the non-linear and linear operators to provide even weighting in
!! the residual calculations. In general these scale factors should be set to
!! the order of magnitude expected for the corresponding variables, \f$ km/s \f$ and
!! \f$ 10^{19} m^{-3} \f$ in this simulation.
CALL ic_fields%B%scale(b0_scale*taylor_hlam(1,oft_hcurl_level))
!---Clean up temporary matrices and fields
CALL lop%delete
DEALLOCATE(tmp,lop)
!---Create velocity field
CALL oft_lag_vcreate(ic_fields%V)
vel_scale = 1.d3
!---Create density field
CALL oft_lag_create(ic_fields%Ne)
CALL ic_fields%Ne%set(n0)
den_scale = n0
den_floor = n0*1.d-2
!---Create temperature field
CALL oft_lag_create(ic_fields%Ti)
CALL ic_fields%Ti%set(t0)
temp_floor = t0*1.d-2
!!\section doc_ex6_code_run Run Simulation
!!
!! Finally, the simulation can be run using the driver routine for non-linear
!! extended MHD (\ref xmhd::xmhd_run "xmhd_run"). This routine advances the
!! solution in time with the physics specified in the input file, see the
!! documentation for \ref xmhd::xmhd_run "xmhd_run", and produces restart files
!! that contain the solution at different times.
!!
!! By default a MG preconditioner is used with the coarsest level specified by
!! \ref xmhd::xmhd_minlev "xmhd_minlev". Several
!! quantities are also recorded to a history file \c "xmhd.hist" during the simulation,
!! including the toroidal current (where the symmetry axis is specified by \ref xmhd::xmhd_taxis
!! "xmhd_taxis"). The data in the history file may be plotted using the script
!! \c "src/utilities/scripts/plot_mhd_hist.py"
!!
!! \note OFT plotting scripts require the python packages NUMPY and MATPLOTLIB as well
!! as path access to the python modules provided in "src/utilities".
!!
!! To visualize the solution fields once a simulation has completed the \ref xmhd::xmhd_plot
!! "xmhd_plot" subroutine is used. This subroutine steps through the restart files
!! produced by \ref xmhd::xmhd_run "xmhd_run" and generates plot files, and optionally
!! probe signals at evenly spaced points in time as specified in the input file, see
!! \ref xmhd::xmhd_plot "xmhd_plot".
xmhd_minlev=minlev
xmhd_taxis=3
oft_env%pm=.FALSE.
IF(plot_run)THEN
  !---Setup I/0
  CALL mesh%setup_io(order)
  !---Run post-processing routine
  CALL xmhd_plot
ELSE
  !---Run simulation
  CALL xmhd_run(ic_fields)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM xmhd_cyl
! STOP SOURCE
!!
!!\section doc_ex6_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment.
!! As with \ref doc_ex5 "Example 5" this example should only be run with multiple processes.
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
!! cad_type=2
!! nlevels=2
!! nbase=1
!! grid_order=2
!! fix_boundary=T
!!/
!!
!!&cubit_options
!! filename='cyl.in.e'
!! inpname='cyl.3dm'
!! lf_file=T
!!/
!!
!!&hcurl_op_options
!! df_wop=0.,.65,.372,.324
!! nu_wop=0,64,2,1
!!/
!!
!!&lag_op_options
!! df_lop=0.,.9,0.86,0.64
!! nu_lop=0,64,2,1
!!/
!!
!!&cyl_options
!! order=3
!! minlev=2
!! b0_scale=1.e-1
!! n0=1.e19
!! t0=6.
!! plot_run=F
!!/
!!
!!&xmhd_options
!! xmhd_ohmic=T      ! Include Ohmic heating
!! xmhd_visc_heat=T  ! Include viscous heating
!! bbc='bc'          ! Perfectly-conducting BC for B-field
!! vbc='all'         ! Zero-flow BC for velocity
!! nbc='d'           ! Dirichlet BC for density
!! tbc='d'           ! Dirichlet BC for temperature
!! dt=6.e-8          ! Maximum time step
!! eta=25.           ! Resistivity at reference temperature (Spitzer-like)
!! eta_temp=6.       ! Reference temperature for resistivity
!! nu_par=400.       ! Fluid viscosity
!! d_dens=10.        ! Density diffusion
!! kappa_par=1.E4    ! Parallel thermal conduction (fixed)
!! kappa_perp=1.E2   ! Perpendicular thermal conduction (fixed)
!! nsteps=2000       ! Number of time steps to take
!! rst_freq=10       ! Restart file frequency
!! lin_tol=1.E-10    ! Linear solver tolerance
!! nl_tol=1.E-5      ! Non-linear solver tolerance
!! nu_xmhd=0,1,8,2   ! Number of smoother iterations for default preconditioner
!! rst_ind=0         ! Index of file to restart from (0 -> use subroutine arguments)
!! ittarget=40       ! Target for # of linear iterations per time step
!! mu_ion=2.          ! Ion mass (atomic units)
!! xmhd_prefreq=20   ! Preconditioner update frequency
!!/
!!\endverbatim
!!
!!\subsection doc_ex6_input_plot Post-Processing options
!!
!! When running the code for post-processing additional run time options are available.
!!
!!\verbatim
!!&xmhd_plot_options
!! t0=1.E-8
!! dt=1.E-6
!! rst_start=0
!! rst_end=2000
!!/
!!\endverbatim
!!
!!\subsection doc_ex6_input_solver Solver specification
!!
!! Time dependent MHD solvers are accelerated significantly by the use of
!! a more sophisticated preconditioner than the default method. Below is
!! an example `oft_in.xml` file that constructs an appropriate MG preconditioner.
!! Currently, this preconditioner method is the suggest preconditioner for all
!! time-dependent MHD solves.
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ ./example6 oft.in oft_in.xml
!!\endverbatim
!!
!! \warning Use of this preconditioner requires OFT be built with the PETSc and
!! FOX libraries.
!!
!!\verbatim
!!<oft>
!!  <xmhd>
!!    <pre type="mg">
!!      <smoother direction="both">
!!        <solver type="gmres">
!!          <its>0,0,2,2</its>
!!          <nrits>0,0,2,2</nrits>
!!          <pre type="add_schwarz">
!!            <nlocal>0,0,-1,-1</nlocal>
!!            <solver type="lu">
!!              <type>lu</type>
!!              <package>superd</package>
!!            </solver>
!!          </pre>
!!        </solver>
!!      </smoother>
!!      <coarse>
!!        <solver type="gmres">
!!          <its>12</its>
!!          <nrits>12</nrits>
!!          <pre type="add_schwarz">
!!            <nlocal>1</nlocal>
!!            <solver type="lu">
!!              <type>lu</type>
!!              <package>superd</package>
!!            </solver>
!!          </pre>
!!        </solver>
!!      </coarse>
!!    </pre>
!!  </xmhd>
!!</oft>
!!\endverbatim
