!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file diag_hcurl_grad_cond.F90
!
!> Diagnostic driver (not a regression test): dumps the H(Curl) + Grad(H^1) mass
!! matrix and the block partitions used by MUG's local smoothers, for offline
!! conditioning and block-Jacobi analysis.
!!
!! Written to the current directory (serial runs only):
!! - `mass.txt`: header `nr nnz`, then `i j value` (1-based) for every stored entry
!! - `dofs.txt`: header `nfields np ne nf nc`, then for each field `ne g1 g2 g3 g4`
!! - `color_entity.txt`: block of each DoF from the entity-aligned FE partition
!!   (as used by MUG with `xmhd_nparts`)
!! - `color_metis.txt`: block of each DoF from METIS on the scalar DoF graph (as used by
!!   block-Jacobi without a coloring)
!! - `cell_blocks_<tag>.txt`: cell-overlap blocks (`nblocks`, then for each block
!!   `n`, the DoFs, `na`, the assigned DoFs) for overlap 0 and 1 (point/face adjacency)
!! - `curlcurl.txt`: curl-curl matrix of the H(Curl) field alone, same format as `mass.txt`
!!   (its DoFs are the leading rows/columns of the H(Curl) + Grad(H^1) numbering)
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!---------------------------------------------------------------------------------
PROGRAM diag_hcurl_grad_cond
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_mesh_local_util, ONLY: mesh_local_orient, oriented_edges, oriented_faces
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_graph, oft_overlap_dofs, overlap_dofs_delete, &
  OFT_CELL_ADJ_POINT, OFT_CELL_ADJ_FACE
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast, partition_graph
USE fem_base, ONLY: oft_ml_fem_type, oft_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE fem_cell_overlap, ONLY: oft_fem_cell_overlap, fem_cell_overlap_create
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: oft_hcurl_getwop
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_getmop
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_vgetmop
IMPLICIT NONE
INTEGER(i4) :: order,nparts,ierr,io_unit,i,j,minlev
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: color
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_h1,ML_h1grad,ML_hcurl
TYPE(oft_ml_fem_comp_type), TARGET :: ML_hcurl_grad,ML_vlag
TYPE(oft_ml_fem_type), TARGET :: ML_lag
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_native_matrix), POINTER :: mop_native
CLASS(oft_vector), POINTER :: u
TYPE(oft_fem_cell_overlap) :: cell_ovl
TYPE(oft_graph) :: graph
NAMELIST/diag_hcurl_grad_options/order,nparts
!---Initialize enviroment
CALL oft_init
order=1
nparts=8
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,diag_hcurl_grad_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(oft_env%nprocs>1)CALL oft_abort('Serial runs only','diag_hcurl_grad_cond',__FILE__)
!---Setup grid and FE spaces (finest level only)
CALL multigrid_construct(mg_mesh)
minlev=mg_mesh%mgdim
CALL oft_h1_setup(mg_mesh,order+1,ML_h1,minlev=minlev)
CALL oft_hcurl_setup(mg_mesh,order,ML_hcurl,minlev=minlev)
CALL oft_hcurl_grad_setup(ML_hcurl,ML_h1,ML_hcurl_grad,ML_h1grad,minlev=minlev)
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---Mass matrix
CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop,'none')
CALL dump_all('',ML_hcurl_grad)
!---Curl-curl matrix of the H(Curl) field
NULLIFY(mop)
CALL oft_hcurl_getwop(ML_hcurl%current_level,mop,'none')
IF(.NOT.native_matrix_cast(mop_native,mop))CALL oft_abort('Native matrix required','diag_hcurl_grad_cond',__FILE__)
OPEN(NEWUNIT=io_unit,FILE='curlcurl.txt')
WRITE(io_unit,'(2I12)')mop_native%nr,mop_native%nnz
DO i=1,mop_native%nr
  DO j=mop_native%kr(i),mop_native%kr(i+1)-1
    WRITE(io_unit,'(2I10,ES26.17E3)')i,mop_native%lc(j),mop_native%M(j)
  END DO
END DO
CLOSE(io_unit)
!---Mesh and local orientation (tetrahedra only)
IF(mg_mesh%mesh%type==1)CALL dump_mesh
!---Vector Lagrange (same order) baseline
CALL oft_lag_setup(mg_mesh,order,ML_lag,ML_vlag_obj=ML_vlag,minlev=minlev)
CALL ML_vlag%set_level(ML_vlag%nlevels,propogate=.TRUE.)
NULLIFY(mop)
CALL oft_lag_vgetmop(ML_vlag%current_level,mop,'none')
CALL dump_all('lag_',ML_vlag)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!> Dump matrix, DoF layout, partitions and cell overlap blocks for one representation
SUBROUTINE dump_all(prefix,ML_rep)
CHARACTER(LEN=*), INTENT(in) :: prefix
TYPE(oft_ml_fem_comp_type), INTENT(inout) :: ML_rep
IF(.NOT.native_matrix_cast(mop_native,mop))CALL oft_abort('Native matrix required','diag_hcurl_grad_cond',__FILE__)
NULLIFY(u)
CALL ML_rep%vec_create(u)
OPEN(NEWUNIT=io_unit,FILE=prefix//'mass.txt')
WRITE(io_unit,'(2I12)')mop_native%nr,mop_native%nnz
DO i=1,mop_native%nr
  DO j=mop_native%kr(i),mop_native%kr(i+1)-1
    WRITE(io_unit,'(2I10,ES26.17E3)')i,mop_native%lc(j),mop_native%M(j)
  END DO
END DO
CLOSE(io_unit)
!---DoF layout
OPEN(NEWUNIT=io_unit,FILE=prefix//'dofs.txt')
WRITE(io_unit,'(5I10)')ML_rep%current_level%nfields,mg_mesh%mesh%np,mg_mesh%mesh%ne, &
  mg_mesh%mesh%nf,mg_mesh%mesh%nc
DO i=1,ML_rep%current_level%nfields
  SELECT TYPE(this=>ML_rep%current_level%fields(i)%fe)
  CLASS IS(oft_fem_type)
    WRITE(io_unit,'(5I10)')this%ne,this%gstruct
  END SELECT
END DO
CLOSE(io_unit)
!---Entity-aligned partition (as MUG with xmhd_nparts) and cell overlap blocks
ALLOCATE(color(mop_native%nr))
CALL fem_cell_overlap_create(cell_ovl,mg_mesh,ML_rep%current_level,nparts,color)
CALL write_ints(prefix//'color_entity.txt',color)
CALL dump_cell_blocks(prefix//'cell_blocks_ov0.txt',0,OFT_CELL_ADJ_POINT)
CALL dump_cell_blocks(prefix//'cell_blocks_ov1p.txt',1,OFT_CELL_ADJ_POINT)
CALL dump_cell_blocks(prefix//'cell_blocks_ov1f.txt',1,OFT_CELL_ADJ_FACE)
!---METIS on the scalar DoF graph (block-Jacobi without a coloring)
graph%nr=mop_native%nr
graph%nnz=mop_native%nnz
graph%kr=>mop_native%kr
graph%lc=>mop_native%lc
color=0
CALL partition_graph(graph,nparts,color)
CALL write_ints(prefix//'color_metis.txt',color)
DEALLOCATE(color)
CALL cell_ovl%delete
END SUBROUTINE dump_all
!> Write points and, per cell: vertices, edges, oriented edges, faces, oriented faces (1-based)
SUBROUTINE dump_mesh
INTEGER(i4) :: k,iunit
OPEN(NEWUNIT=iunit,FILE='mesh.txt')
WRITE(iunit,'(2I10)')mg_mesh%mesh%np,mg_mesh%mesh%nc
DO k=1,mg_mesh%mesh%np
  WRITE(iunit,'(3ES26.17E3)')mg_mesh%mesh%r(:,k)
END DO
DO k=1,mg_mesh%mesh%nc
  CALL mesh_local_orient(mg_mesh%mesh,k)
  WRITE(iunit,'(4I10,6I10,12I3,4I10,12I3)')mg_mesh%mesh%lc(:,k),ABS(mg_mesh%mesh%lce(:,k)), &
    oriented_edges,ABS(mg_mesh%mesh%lcf(:,k)),oriented_faces
END DO
CLOSE(iunit)
END SUBROUTINE dump_mesh
!> Write an integer array, one value per line, with the count as header
SUBROUTINE write_ints(filename,vals)
CHARACTER(LEN=*), INTENT(in) :: filename
INTEGER(i4), INTENT(in) :: vals(:)
INTEGER(i4) :: k,iunit
OPEN(NEWUNIT=iunit,FILE=filename)
WRITE(iunit,'(I12)')SIZE(vals)
DO k=1,SIZE(vals)
  WRITE(iunit,'(I8)')vals(k)
END DO
CLOSE(iunit)
END SUBROUTINE write_ints
!> Build and write cell overlap blocks (global row ids, 1-based)
SUBROUTINE dump_cell_blocks(filename,layers,adj)
CHARACTER(LEN=*), INTENT(in) :: filename
INTEGER(i4), INTENT(in) :: layers,adj
TYPE(oft_overlap_dofs) :: dofs
INTEGER(i4) :: b,k,iunit
CALL cell_ovl%build(mop,u,layers,adj,-1,dofs)
OPEN(NEWUNIT=iunit,FILE=filename)
WRITE(iunit,'(I8)')dofs%nblocks
DO b=1,dofs%nblocks
  WRITE(iunit,'(I10)')dofs%blk(b)%n
  WRITE(iunit,'(10I10)')(INT(dofs%gdof(dofs%blk(b)%v(k)),4),k=1,dofs%blk(b)%n)
  WRITE(iunit,'(I10)')dofs%asg(b)%n
  WRITE(iunit,'(10I10)')(INT(dofs%gdof(dofs%asg(b)%v(k)),4),k=1,dofs%asg(b)%n)
END DO
CLOSE(iunit)
CALL overlap_dofs_delete(dofs)
END SUBROUTINE dump_cell_blocks
END PROGRAM diag_hcurl_grad_cond
