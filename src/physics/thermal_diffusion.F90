!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thermal_diffusion.F90
!
!> Example module for modeling thermal diffusion of two species, with equilibration, in 2D
!!
!! The equation to be solved is a simple 
!!
!! \f[ \frac{\partial T_i}{\partial t} = kappa_i T_i^(5/2) \nabla^2 T_i  + \tau_eq (T_e - T_i) + Q_i \f]
!!
!! \f[ \frac{\partial T_e}{\partial t} = kappa_e T_e^(5/2) \nabla^2 T_e  + \tau_eq (T_i - T_e) + Q_e \f]
!!
!! @authors Chris Hansen
!! @date January 2025
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE thermal_diffusion
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat
USE oft_deriv_matrices, ONLY: oft_noop_matrix
USE fem_composite, ONLY: oft_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup_bmesh, oft_scalar_bfem, oft_blagrange, &
  oft_lag_setup
IMPLICIT NONE
#include "local.h"
PRIVATE
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
type, extends(oft_noop_matrix) :: tdiff_nlfun
  real(r8) :: dt = -1.d0 !< Time step
!   real(r8) :: diag_vals(7) = 0.d0 !< Diagnostic values
contains
  !> Apply the matrix
  procedure :: apply_real => nlfun_apply
end type tdiff_nlfun
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
type, extends(oft_noop_matrix) :: tdiff_massfun
!   real(r8) :: diag_vals(7) = 0.d0 !< Diagnostic values
contains
  !> Apply the matrix
  procedure :: apply_real => massfun_apply
end type tdiff_massfun
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
type, public :: oft_tdiff_sim
  real(r8) :: dt = -1.d0 !< Needs docs
  real(r8) :: t = 0.d0 !< Needs docs
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: Ti_bc => NULL() !< Ti BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: Te_bc => NULL() !< Te BC flag
  INTEGER(i4)L, CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  type(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  class(oft_vector), pointer :: u => NULL() !< Needs docs
  class(oft_matrix), pointer :: jacobian => NULL() !< Needs docs
  type(tdiff_nlfun) :: nlfun !< Needs docs
  type(tdiff_massfun) :: massfun !< Needs docs
contains
  !> Setup
  procedure :: setup => setup
  !> Run simulation
  procedure :: run_simulation => run_simulation
  !> Save restart file
  procedure :: rst_save => rst_save
  !> Load restart file
  procedure :: rst_load => rst_load
end type oft_tdiff_sim
CONTAINS
!---------------------------------------------------------------------------
!> Main driver subroutine for extended MHD time advance
!!
!! Runtime options are set in the main input file using the group
!! \c xmhd_options group, see \ref xmhd_read_settings.
!---------------------------------------------------------------------------
subroutine run_simulation(self)
class(oft_tdiff_sim), intent(inout) :: self
type(oft_nksolver) :: nksolver
!---Solver objects
class(oft_vector), pointer :: u,v,up
type(oft_native_gmres_solver), target :: solver
type(oft_timer) :: mytimer
!---History file fields
TYPE(oft_bin_file) :: hist_file
integer(i4) :: hist_i4(3)
real(4) :: hist_r4(4)
!---
real(r8) :: t,dt
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)

t=0.d0
CALL u%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
WRITE(rst_char,104)0
CALL self%rst_save(u, t, dt, 'tDiff_'//rst_char//'.rst', 'U')

!---
CALL build_approx_jacobian(self,u)

!---------------------------------------------------------------------------
! Setup linear solver
!---------------------------------------------------------------------------
NULLIFY(solver%pre)
IF(ASSOCIATED(jacobian_pre_node))THEN
  CALL create_solver_xml(solver%pre,jacobian_pre_node)
ELSE
  CALL create_diag_pre(solver%pre)
END IF
solver%pre%A=>self%jacobian

!---------------------------------------------------------------------------
! Setup nonlinear solver
!---------------------------------------------------------------------------
nksolver%A=>self%nlfun
nksolver%J_inv=>solver
nksolver%its=30
nksolver%atol=nl_tol
nksolver%rtol=1.d-20 ! Disable relative tolerance
IF(xmhd_mfnk)THEN
  nksolver%J_update=>xmhd_mfnk_update
  nksolver%up_freq=1
ELSE
  nksolver%J_update=>xmhd_set_ops
  nksolver%up_freq=4
END IF

!---------------------------------------------------------------------------
! Setup history file
!---------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('oft_tdiff.hist', desc="History file for non-linear thermal diffusion run")
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('ti',   'r4', desc="Average ion temperature [eV]")
  CALL hist_file%add_field('te',   'r4', desc="Average electron temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
END IF

!---------------------------------------------------------------------------
! Begin time stepping
!---------------------------------------------------------------------------
DO i=1,nsteps
  CALL up%add(0.d0,1.d0,u)
  CALL self%massfun%apply(u,v)
  ! IF(ASSOCIATED(self%driver))CALL self%driver%apply(up,u,v,t,dt)
  CALL nksolver%apply(u,v)
  !---------------------------------------------------------------------------
  ! Write out initial solution progress
  !---------------------------------------------------------------------------
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    hist_i4=[rst_ind+i-1,nksolver%lits,nksolver%nlits]
    hist_r4=REAL([t,Ti_avg,Te_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,2X,I4,F12.3,ES12.2)
    WRITE(*,103)rst_ind+i,t,nksolver%lits,nksolver%nlits,elapsed_time,dt
    IF(oft_debug_print(1))WRITE(*,*)
    CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
  END IF
  !---------------------------------------------------------------------------
  ! Update timestep and save solution
  !---------------------------------------------------------------------------
  t=t+dt
  IF(MOD(i,rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Create restart file
    WRITE(rst_char,104)rst_ind+i
    READ(rst_char,104,IOSTAT=io_stat)rst_tmp
    IF((io_stat/=0).OR.(rst_tmp/=rst_ind+i))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
    CALL self%rst_save(u, t, dt, 'tDiff_'//rst_char//'.rst', 'U')
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
      CALL hist_file%flush
    END IF
  END IF
END DO

end subroutine run_simulation
!---------------------------------------------------------------------------
!> Compute the NL error function, where we are solving F(x) = 0
!!
!! b = F(a)
!---------------------------------------------------------------------------
subroutine nlfun_apply(self,a,b)
class(tdiff_nlfun), intent(inout) :: self !< NL function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
!$omp parallel private(curved,cell_dofs)
!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  Ti_weights_loc = Ti_weights(cell_dofs)
  Te_weights_loc = Te_weights(cell_dofs)
!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    CALL oft_lag_eval_all(oft_blagrange,i,quad%pts(:,m),basis_vals)
    CALL oft_lag_geval_all(oft_blagrange,i,quad%pts(:,m),basis_grads,jac_mat)
    !---Reconstruct values of solution fields
    T_i = 0.d0; dT_i = 0.d0; T_e = 0.d0; dT_e = 0.d0
    DO jr=1,oft_blagrange%nce
      T_i = T_i + Ti_weights_loc(jr)*basis_vals(jr)
      T_e = T_e + Te_weights_loc(jr)*basis_vals(jr)
      dT_i = dT_i + Ti_weights_loc(jr)*basis_grads(:,jr)
      dT_e = dT_e + Te_weights_loc(jr)*basis_grads(:,jr)
    END DO
    !---Compute local function contributions
    DO jr=1,oft_blagrange%nce
      res_loc(jr,1) = res_loc(jr,1) &
        + basis_vals(jr)*T_i*jac_det*quad%wts(m) &
        + self%dt*kappa_i*(T_i**2.5d0)*DOT_PRODUCT(basis_grads(:,jr),dT_i)*jac_det*quad%wts(m) &
        + self%dt*tau_eq*(T_e-T_i)*jac_det*quad%wts(m)
      res_loc(jr,2) = res_loc(jr,2) &
        + basis_vals(jr)*T_e*jac_det*quad%wts(m) &
        + self%dt*kappa_e*(T_e**2.5d0)*DOT_PRODUCT(basis_grads(:,jr),dT_e)*jac_det*quad%wts(m) &
        + self%dt*tau_eq*(T_i-T_e)*jac_det*quad%wts(m)
    END DO
  END DO
  !$omp critical (tdiff_nlfun)
  DO jr=1,oft_blagrange%nce
    res_full(cell_dofs(jr),:) = res_full(cell_dofs(jr),:) + res_loc(jr,:)
  END DO
  !$omp end critical (tdiff_nlfun)
END DO
!$omp end parallel
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,Ti_weights,res_full(:,1),self%Ti_bc)
CALL fem_dirichlet_vec(oft_blagrange,Te_weights,res_full(:,2),self%Te_bc)
!---Put results into full vector
CALL b%restore_local(res_full(:,1),1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(res_full(:,2),2,add=.TRUE.)
end subroutine nlfun_apply
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_approx_jacobian(self,a)
class(oft_tdiff_sim), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Solution for computing jacobian
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
CALL self%jacobian%zero
!--Setup thread locks
ALLOCATE(tlocks(xmhd_rep%nfields))
DO i=1,xmhd_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
!$omp parallel private(curved,cell_dofs)
allocate(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
  Ti_weights_loc = Ti_weights(cell_dofs)
  Te_weights_loc = Te_weights(cell_dofs)
!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    CALL oft_lag_eval_all(oft_blagrange,i,quad%pts(:,m),basis_vals)
    CALL oft_lag_geval_all(oft_blagrange,i,quad%pts(:,m),basis_grads,jac_mat)
    !---Reconstruct values of solution fields
    T_i = 0.d0; dT_i = 0.d0; T_e = 0.d0; dT_e = 0.d0
    DO jr=1,oft_blagrange%nce
      T_i = T_i + Ti_weights_loc(jr)*basis_vals(jr)
      T_e = T_e + Te_weights_loc(jr)*basis_vals(jr)
      dT_i = dT_i + Ti_weights_loc(jr)*basis_grads(:,jr)
      dT_e = dT_e + Te_weights_loc(jr)*basis_grads(:,jr)
    END DO
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        !---Ion rows
        jac_loc(1,1)%m = jac_loc(1,1)%m &
          + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)
          + self%dt*kappa_i*(T_i**2.5d0)*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m) &
          + self%dt*kappa_i*(T_i**1.5d0)*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr),dT_i)*jac_det*quad%wts(m) &
          + self%dt*tau_eq*(-basis_vals(jc))*jac_det*quad%wts(m)
        jac_loc(1,2)%m = jac_loc(1,2)%m &
          + self%dt*tau_eq*(basis_vals(jc))*jac_det*quad%wts(m)
        !---Electron rows
        jac_loc(2,2)%m = jac_loc(2,2)%m &
          + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)
          + self%dt*kappa_e*(T_e**2.5d0)*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m) &
          + self%dt*kappa_e*(T_e**1.5d0)*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr),dT_e)*jac_det*quad%wts(m) &
          + self%dt*tau_eq*(-basis_vals(jc))*jac_det*quad%wts(m)
        jac_loc(2,1)%m = jac_loc(2,1)%m &
          + self%dt*tau_eq*(basis_vals(jc))*jac_det*quad%wts(m)
      END DO
    END DO
  END DO
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%Ti_bc(cell_dofs),1)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%Te_bc(cell_dofs),2)
  CALL self%fe_rep%mat_add_local(self%jacobian,jac_loc,iloc,tlocks)
END DO
!$omp end parallel
!--Destroy thread locks
DO i=1,xmhd_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
IF(oft_debug_print(2))write(*,'(4X,A)')'Setting BCs'
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%Ti_bc,1)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%Te_bc,2)
!
call self%fe_rep%vec_create(tmp)
call self%jacobian%assemble(tmp)
call tmp%delete
DEALLOCATE(tmp)
end subroutine build_approx_jacobian
!---------------------------------------------------------------------------
!> Setup composite FE representation and ML environment
!---------------------------------------------------------------------------
subroutine setup(self,order)
class(oft_tdiff_sim), intent(inout) :: self
integer(i4), intent(in) :: order
IF(ASSOCIATED(self%fe_rep))CALL oft_abort("Setup can only be called once","setup",__FILE__)
IF(.NOT.ASSOCAITED(oft_blagrange))THEN
  IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building lagrange FE space'
  smesh%tess_order=order
  CALL oft_blag_setup(order, -1)
END IF
!
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Creating FE type'
ALLOCATE(self%fe_rep)
self%fe_rep%nfields=2
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields))
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange
self%fe_rep%field_tags(1)='Ti'
self%fe_rep%fields(2)%fe=>oft_blagrange
self%fe_rep%field_tags(2)='Te'
end subroutine setup
!---------------------------------------------------------------------------
!> Save xMHD solution state to a restart file
!---------------------------------------------------------------------------
subroutine rst_save(self,u,t,dt,filename,path)
class(oft_tdiff_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to save
real(r8), intent(in) :: t !< Current solution time
real(r8), intent(in) :: dt !< Current timestep
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_save(u,filename,path)
IF(oft_env%head_proc)THEN
  CALL hdf5_write(t,filename,'t')
  CALL hdf5_write(dt,filename,'dt')
END IF
DEBUG_STACK_POP
end subroutine rst_save
!---------------------------------------------------------------------------
!> Load xMHD solution state from a restart file
!---------------------------------------------------------------------------
subroutine rst_load(self,u,filename,path,t,dt)
class(oft_tdiff_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to load
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
real(r8), optional, intent(out) :: t !< Time of loaded solution
real(r8), optional, intent(out) :: dt !< Timestep at loaded time
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_load(u,filename,path)
IF(PRESENT(t))CALL hdf5_read(t,filename,'t')
IF(PRESENT(dt))CALL hdf5_read(dt,filename,'dt')
DEBUG_STACK_POP
end subroutine rst_load
END MODULE thermal_diffusion