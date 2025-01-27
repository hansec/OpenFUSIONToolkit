PROGRAM tDiff_circle
!---Runtime
USE oft_base
!---Grid
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector
!
USE oft_lag_basis, ONLY: oft_blagrange
USE oft_blag_operators, ONLY: blag_zerob
USE thermal_diffusion_2d
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_tdiff_sim) :: tDiff_sim
CLASS(oft_vector), POINTER :: vec_tmp
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 100
INTEGER(i4) :: rst_freq = 10
REAL(r8) :: ti0 = 1.d0
REAL(r8) :: te0 = 4.d0
REAL(r8) :: kappa_i = 1.d0
REAL(r8) :: kappa_e = 1.d0
REAL(r8) :: tau_eq = -1.d0
REAL(r8) :: dt = 1.d-3
LOGICAL :: pm=.FALSE.
LOGICAL :: use_mfnk=.FALSE.
NAMELIST/tdiff_options/order,ti0,te0,kappa_i,kappa_e,tau_eq, &
dt,nsteps,rst_freq,use_mfnk,pm
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,tdiff_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf
!
CALL tDiff_sim%setup(order)
!
NULLIFY(vec_tmp,vec_vals)
CALL oft_blagrange%vec_create(vec_tmp)
CALL vec_tmp%set(ti0)
CALL blag_zerob(vec_tmp)
CALL vec_tmp%get_local(vec_vals)
CALL tDiff_sim%u%restore_local(vec_vals,1)
CALL vec_tmp%set(te0)
CALL blag_zerob(vec_tmp)
CALL vec_tmp%get_local(vec_vals)
CALL tDiff_sim%u%restore_local(vec_vals,2)
!
tDiff_sim%kappa_i=kappa_i
tDiff_sim%kappa_e=kappa_e
tDiff_sim%tau_eq=tau_eq
tDiff_sim%dt=dt
tDiff_sim%nsteps=nsteps
tDiff_sim%rst_freq=rst_freq
tDiff_sim%mfnk=use_mfnk
oft_env%pm=pm
CALL tDiff_sim%run_simulation()
!---Finalize enviroment
CALL oft_finalize
END PROGRAM tDiff_circle