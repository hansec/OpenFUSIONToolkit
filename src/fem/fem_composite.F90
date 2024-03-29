!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file fem_composite.F90
!
!> Classes and infrastructure for composite FE representations
!!
!! @authors Chris Hansen
!! @date April 2014
!! @ingroup doxy_oft_fem
!---------------------------------------------------------------------------
MODULE fem_composite
USE oft_base
USE oft_stitching, ONLY: oft_seam, seam_list
USE oft_io, ONLY: hdf5_rst, hdf5_write, hdf5_read, hdf5_rst_destroy, hdf5_create_file
! USE oft_io, ONLY: hdf5_rst, hdf5_rst_write, hdf5_rst_read, hdf5_rst_destroy, &
!   hdf5_rst_read_version, hdf5_rst_write_version
!---
USE oft_la_base, ONLY: oft_vector, oft_map, map_list, oft_graph_ptr, &
  oft_matrix, oft_matrix_ptr, oft_local_mat
USE oft_native_la, ONLY: oft_native_vector, native_vector_cast, &
  native_vector_slice_push, native_vector_slice_pop
USE oft_la_utils, ONLY: create_vector, combine_matrices, create_matrix, &
  create_identity_graph
!---
USE fem_base, ONLY: oft_fem_ptr, oft_ml_fem_ptr, fem_max_levels, &
  fem_common_linkage, fem_gen_legacy, fem_idx_ver, fem_idx_path
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
! TYPE oft_fem_type
!---------------------------------------------------------------------------
!> Composite FE type
!---------------------------------------------------------------------------
TYPE, PUBLIC :: oft_fem_comp_type
  INTEGER(i4) :: nfields = 0 !< Number of fields in composite representation
  CHARACTER(LEN=4), POINTER, DIMENSION(:) :: field_tags => NULL() !< Character tags for fields
  TYPE(oft_fem_ptr), POINTER, DIMENSION(:) :: fields => NULL() !< Individual representations
  TYPE(oft_seam), POINTER :: linkage => NULL() !< Global linkage information
  TYPE(oft_map), POINTER :: map(:) => NULL() !< Linear algebra mapping
  CLASS(oft_vector), POINTER :: cache_PETSc => NULL() !< PETSc vector cache
  CLASS(oft_vector), POINTER :: cache_native => NULL() !< Native vector cache
CONTAINS
  !> Create vector for FE representation
  PROCEDURE :: vec_create => fem_vec_create
  !> Save vector to HDF5 file
  PROCEDURE :: vec_save => fem_vec_save
  !> Load vector from HDF5 file
  PROCEDURE :: vec_load => fem_vec_load
  !> Create matrix for FE representation
  PROCEDURE :: mat_create => fem_mat_create
  !> Create matrix for FE representation
  PROCEDURE :: mat_setup_local => fem_mat_setup_local
  !> Create matrix for FE representation
  PROCEDURE :: mat_destroy_local => fem_mat_destroy_local
  !> Create matrix for FE representation
  PROCEDURE :: mat_zero_local => fem_mat_zero_local
  !> Create matrix for FE representation
  PROCEDURE :: mat_zero_local_rows => fem_mat_zero_local_rows
  !> Create matrix for FE representation
  PROCEDURE :: mat_add_local => fem_mat_add_local
END TYPE oft_fem_comp_type
!---------------------------------------------------------------------------
! TYPE oft_fem_comp_ptr
!---------------------------------------------------------------------------
!> Base FE type
!---------------------------------------------------------------------------
TYPE, PUBLIC :: oft_fem_comp_ptr
  TYPE(oft_fem_comp_type), POINTER :: fe => NULL() !< Finite element object
END TYPE oft_fem_comp_ptr
!---------------------------------------------------------------------------
! TYPE oft_ml_fem_comp_type
!---------------------------------------------------------------------------
!> Base FE type
!---------------------------------------------------------------------------
TYPE, PUBLIC :: oft_ml_fem_comp_type
  INTEGER(i4) :: nfields = 0 !< Number of fields in composite representation
  INTEGER(i4) :: nlevels = 0 !< Number of FE levels
  INTEGER(i4) :: level = 0 !< Current FE level
  INTEGER(i4) :: abs_level = 0 !< Asoblute FE refinement level
  INTEGER(i4) :: blevel = 0 !< FE base level
  TYPE(oft_fem_comp_ptr) :: levels(fem_max_levels)
  TYPE(oft_fem_comp_type), POINTER :: current_level => NULL()
  TYPE(oft_ml_fem_ptr), POINTER, DIMENSION(:) :: ml_fields => NULL()
  CHARACTER(LEN=4), POINTER, DIMENSION(:) :: field_tags => NULL() !< Character tags for fields
  TYPE(oft_matrix_ptr) :: interp_matrices(fem_max_levels)
CONTAINS
  !> Setup composite FE structure from sub-field FE definitions
  PROCEDURE :: setup => ml_fem_setup
  !> Create vector for FE representation
  PROCEDURE :: vec_create => ml_fem_vec_create
  !> Build interpolation operators from sub-fields
  PROCEDURE :: build_interp => ml_fem_build_interp
  !> Set level in ML framework if available
  PROCEDURE :: set_level => ml_fem_set_level
END TYPE oft_ml_fem_comp_type
CONTAINS
!---------------------------------------------------------------------------
! SUBROUTINE: fem_vec_create
!---------------------------------------------------------------------------
!> Create weight vector for FE representation
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine fem_vec_create(self,new,cache,native)
class(oft_fem_comp_type), intent(inout) :: self
class(oft_vector), pointer, intent(out) :: new
logical, optional, intent(in) :: cache
logical, optional, intent(in) :: native
INTEGER(i4) :: i
TYPE(map_list), POINTER :: maps(:)
TYPE(seam_list), POINTER :: stitches(:)
logical :: do_cache,force_native
DEBUG_STACK_PUSH
do_cache=.TRUE.
force_native=.FALSE.
IF(PRESENT(cache))do_cache=cache
IF(PRESENT(native))force_native=native
!---
IF(use_petsc.AND.(.NOT.force_native))THEN
  IF(ASSOCIATED(self%cache_PETSc))THEN
    CALL self%cache_PETSc%new(new)
  ELSE
    !---Create vector
    ALLOCATE(stitches(self%nfields),maps(self%nfields))
    DO i=1,self%nfields
      stitches(i)%s=>self%fields(i)%fe%linkage
      maps(i)%m=>self%fields(i)%fe%map
    END DO
    CALL create_vector(new,stitches,maps)
    DEALLOCATE(stitches,maps)
    !---
    IF(do_cache)THEN
      CALL new%new(self%cache_PETSc)
      IF(.NOT.ASSOCIATED(self%linkage))self%linkage=>self%cache_PETSc%stitch_info
      IF(.NOT.ASSOCIATED(self%map))self%map=>self%cache_PETSc%map
    END IF
  END IF
ELSE
  IF(ASSOCIATED(self%cache_native))THEN
    CALL self%cache_native%new(new)
  ELSE
    !---Create vector
    ALLOCATE(stitches(self%nfields),maps(self%nfields))
    DO i=1,self%nfields
      stitches(i)%s=>self%fields(i)%fe%linkage
      maps(i)%m=>self%fields(i)%fe%map
    END DO
    CALL create_vector(new,stitches,maps,native=native)
    DEALLOCATE(stitches,maps)
    !---
    IF(do_cache)THEN
      CALL new%new(self%cache_native)
      IF(.NOT.ASSOCIATED(self%linkage))self%linkage=>self%cache_native%stitch_info
      IF(.NOT.ASSOCIATED(self%map))self%map=>self%cache_native%map
    END IF
  END IF
END IF
DEBUG_STACK_POP
end subroutine fem_vec_create
!---------------------------------------------------------------------------
! SUBROUTINE: fem_vec_save
!---------------------------------------------------------------------------
!> Save a Lagrange scalar field to a HDF5 restart file.
!!
!! @param[in] source Source field
!! @param[in] filen Name of destination file
!! @param[in] tag Field label in file
!---------------------------------------------------------------------------
subroutine fem_vec_save(self,source,filename,path,append)
class(oft_fem_comp_type), intent(inout) :: self
class(oft_vector), target, intent(inout) :: source
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, optional, intent(in) :: append
INTEGER(i4) :: i
real(r8), pointer, dimension(:) :: valtmp
class(oft_vector), pointer :: outfield
class(oft_native_vector), pointer :: outvec
type(hdf5_rst) :: rst_info
logical :: do_append
DEBUG_STACK_PUSH
do_append=.FALSE.
IF(PRESENT(append))do_append=append
IF(.NOT.do_append)THEN
  IF(oft_env%head_proc)THEN
    CALL hdf5_create_file(filename)
    CALL hdf5_write(fem_idx_ver,filename,fem_idx_path)
  END IF
END IF
NULLIFY(valtmp)
!---
IF(oft_debug_print(1))WRITE(*,'(6A)')oft_indent,'Writing "',TRIM(path), &
  '" to file "',TRIM(filename),'"'
!---
DO i=1,self%nfields
  IF(oft_debug_print(2))WRITE(*,'(4X,3A)')'Field -> "', &
  TRIM(path)//'_'//TRIM(self%field_tags(i)),'"'
  CALL self%fields(i)%fe%vec_create(outfield,native=.TRUE.)
  IF(native_vector_cast(outvec,outfield)<0)CALL oft_abort('Failed to create "outfield".', &
    'fem_vec_save',__FILE__)
  !---
  CALL source%get_local(valtmp,i)
  CALL outfield%restore_local(valtmp)
  CALL native_vector_slice_push(outvec,self%fields(i)%fe%global%le,rst_info)
  CALL hdf5_write(rst_info,filename,TRIM(path)//'_'//TRIM(self%field_tags(i)))
  CALL hdf5_rst_destroy(rst_info)
  !---
  CALL outfield%delete
  DEALLOCATE(outfield,valtmp)
END DO
DEBUG_STACK_POP
end subroutine fem_vec_save
!---------------------------------------------------------------------------
! SUBROUTINE: fem_vec_load
!---------------------------------------------------------------------------
!> Load a Lagrange scalar field from a HDF5 restart file.
!!
!! @param[in,out] source Destination field
!! @param[in] filen Name of source file
!! @param[in] path Field path in file
!---------------------------------------------------------------------------
subroutine fem_vec_load(self,source,filename,path)
class(oft_fem_comp_type), intent(inout) :: self
class(oft_vector), target, intent(inout) :: source
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
INTEGER(i4) :: i,version
integer(i8), pointer, dimension(:) :: lge
real(r8), pointer, dimension(:) :: valtmp
class(oft_vector), pointer :: infield
class(oft_native_vector), pointer :: invec
type(hdf5_rst) :: rst_info
logical :: legacy,success
DEBUG_STACK_PUSH
legacy=.FALSE.
CALL hdf5_read(version,filename,fem_idx_path,success)
legacy=(version<1)
!---
NULLIFY(valtmp)
IF(oft_debug_print(1))WRITE(*,'(2X,5A)')'Reading "',TRIM(path), &
  '" from file "',TRIM(filename),'"'
DO i=1,self%nfields
  IF(oft_debug_print(2))WRITE(*,'(4X,3A)')'  Field -> "', &
    TRIM(path)//'_'//TRIM(self%field_tags(i)),'"'
  CALL self%fields(i)%fe%vec_create(infield,native=.TRUE.)
  IF(native_vector_cast(invec,infield)<0)CALL oft_abort('Failed to create "infield".', &
    'fem_vec_load',__FILE__)
  lge=>self%fields(i)%fe%global%le
  IF(legacy)THEN
    CALL fem_gen_legacy(self%fields(i)%fe)
    lge=>self%fields(i)%fe%legacy_lge
  END IF
  !---
  CALL native_vector_slice_push(invec,lge,rst_info,alloc_only=.TRUE.)
  CALL hdf5_read(rst_info,filename,TRIM(path)//'_'//TRIM(self%field_tags(i)))
  CALL native_vector_slice_pop(invec,lge,rst_info)
  CALL infield%get_local(valtmp)
  CALL source%restore_local(valtmp,i)
  !---
  CALL hdf5_rst_destroy(rst_info)
  CALL infield%delete
  DEALLOCATE(infield,valtmp)
END DO
NULLIFY(lge)
DEBUG_STACK_POP
end subroutine fem_vec_load
!---------------------------------------------------------------------------
! SUBROUTINE: fem_mat_create
!---------------------------------------------------------------------------
!> Create weight vector for FE representation
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine fem_mat_create(self,new,mask)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
CLASS(oft_matrix), POINTER, INTENT(out) :: new
INTEGER(i4), OPTIONAL, INTENT(in) :: mask(:,:)
INTEGER(i4) :: i,j,k,nknown_graphs
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: mat_mask,graph_ids
CLASS(oft_vector), POINTER :: tmp_vec
TYPE(oft_graph_ptr), ALLOCATABLE :: graphs(:,:),known_graphs(:)
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Building composite FE matrix'
ALLOCATE(mat_mask(self%nfields,self%nfields))
mat_mask=1
IF(PRESENT(mask))mat_mask=mask
ALLOCATE(graphs(self%nfields,self%nfields))
!---Populate known graphs
ALLOCATE(known_graphs(self%nfields*self%nfields))
ALLOCATE(graph_ids(2,self%nfields*self%nfields))
graph_ids=0
nknown_graphs=0
DO i=1,self%nfields
  DO k=1,nknown_graphs
    IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type,self%fields(i)%fe%type/)))EXIT
  END DO
  IF(k<=nknown_graphs)CYCLE
  !---
  IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Building graph',i,i
  nknown_graphs=nknown_graphs+1
  graph_ids(:,nknown_graphs)=(/self%fields(i)%fe%type,self%fields(i)%fe%type/)
  ALLOCATE(known_graphs(nknown_graphs)%g)
  known_graphs(nknown_graphs)%g%nr=self%fields(i)%fe%ne
  known_graphs(nknown_graphs)%g%nrg=self%fields(i)%fe%global%ne
  known_graphs(nknown_graphs)%g%nc=self%fields(i)%fe%ne
  known_graphs(nknown_graphs)%g%ncg=self%fields(i)%fe%global%ne
  known_graphs(nknown_graphs)%g%nnz=self%fields(i)%fe%nee
  known_graphs(nknown_graphs)%g%kr=>self%fields(i)%fe%kee
  known_graphs(nknown_graphs)%g%lc=>self%fields(i)%fe%lee
END DO
!---Set graphs
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(mat_mask(i,j)==0)CYCLE
    IF(mat_mask(i,j)==2)THEN
      IF(i/=j)CALL oft_abort('Identity only valid on diagonal.', &
      'fem_mat_create',__FILE__)
      !---Setup identity graph
      CALL self%fields(i)%fe%vec_create(tmp_vec)
      CALL create_identity_graph(graphs(i,j)%g,tmp_vec)
      CALL tmp_vec%delete
      DEALLOCATE(tmp_vec)
      CYCLE
    END IF
    DO k=1,nknown_graphs
      IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type, &
      self%fields(j)%fe%type/)))EXIT
    END DO
    IF(k<=nknown_graphs)THEN
      IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Using known graph ',i,j
      graphs(i,j)%g=>known_graphs(k)%g
    ELSE
      IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Building graph ',i,j
      nknown_graphs=nknown_graphs+1
      graph_ids(:,nknown_graphs)=(/self%fields(i)%fe%type, &
      self%fields(j)%fe%type/)
      ALLOCATE(known_graphs(nknown_graphs)%g)
      known_graphs(nknown_graphs)%g%nr=self%fields(i)%fe%ne
      known_graphs(nknown_graphs)%g%nrg=self%fields(i)%fe%global%ne
      known_graphs(nknown_graphs)%g%nc=self%fields(j)%fe%ne
      known_graphs(nknown_graphs)%g%ncg=self%fields(j)%fe%global%ne
      CALL fem_common_linkage(self%fields(i)%fe,self%fields(j)%fe, &
        known_graphs(nknown_graphs)%g%nnz,known_graphs(nknown_graphs)%g%kr, &
        known_graphs(nknown_graphs)%g%lc)
      graphs(i,j)%g=>known_graphs(nknown_graphs)%g
    END IF
  END DO
END DO
!---
CALL self%vec_create(tmp_vec)
CALL create_matrix(new,graphs,tmp_vec,tmp_vec)
CALL tmp_vec%delete
DO i=1,nknown_graphs
  DEALLOCATE(known_graphs(i)%g)
END DO
DEALLOCATE(graphs,known_graphs,mat_mask,graph_ids,tmp_vec)
DEBUG_STACK_POP
end subroutine fem_mat_create
!---------------------------------------------------------------------------
! SUBROUTINE: fem_mat_setup_local
!---------------------------------------------------------------------------
!> Create weight vector for FE representation
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine fem_mat_setup_local(self,mloc,mask)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
TYPE(oft_local_mat), INTENT(inout) :: mloc(:,:)
INTEGER(i4), OPTIONAL, INTENT(in) :: mask(:,:)
INTEGER(i4) :: i,j,k,nknown_graphs
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: mat_mask,graph_ids,known_graphs
DEBUG_STACK_PUSH
ALLOCATE(mat_mask(self%nfields,self%nfields))
mat_mask=1
IF(PRESENT(mask))mat_mask=mask
!---Populate known graphs
ALLOCATE(known_graphs(2,self%nfields*self%nfields))
ALLOCATE(graph_ids(2,self%nfields*self%nfields))
graph_ids=0
nknown_graphs=0
DO i=1,self%nfields
  DO j=1,self%nfields
    DO k=1,nknown_graphs
      IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type,self%fields(j)%fe%type/)))EXIT
    END DO
    IF(k<=nknown_graphs)CYCLE
    nknown_graphs=nknown_graphs+1
    graph_ids(:,nknown_graphs)=(/self%fields(i)%fe%type,self%fields(j)%fe%type/)
    known_graphs(:,nknown_graphs)=(/i,j/)
  END DO
END DO
!---
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(mat_mask(i,j)==1)THEN
      ALLOCATE(mloc(i,j)%m(self%fields(i)%fe%nce,self%fields(j)%fe%nce))
      DO k=1,nknown_graphs
        IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type,self%fields(j)%fe%type/)))EXIT
      END DO
      IF(k<=nknown_graphs)THEN
        mloc(i,j)%ind=>mloc(known_graphs(1,k),known_graphs(2,k))%ind
      ELSE
        ALLOCATE(mloc(i,j)%ind(self%fields(i)%fe%nce,self%fields(j)%fe%nce))
      END IF
    END IF
  END DO
END DO
DEALLOCATE(mat_mask,known_graphs,graph_ids)
DEBUG_STACK_POP
end subroutine fem_mat_setup_local
!---------------------------------------------------------------------------
! SUBROUTINE: fem_mat_destroy_local
!---------------------------------------------------------------------------
!> Create weight vector for FE representation
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine fem_mat_destroy_local(self,mloc)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
TYPE(oft_local_mat), INTENT(inout) :: mloc(:,:)
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
!---
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(ASSOCIATED(mloc(i,j)%m))deallocate(mloc(i,j)%m)
    IF(ASSOCIATED(mloc(i,j)%ind))deallocate(mloc(i,j)%ind)
  END DO
END DO
DEBUG_STACK_POP
end subroutine fem_mat_destroy_local
!---------------------------------------------------------------------------
! SUBROUTINE: fem_mat_zero_local
!---------------------------------------------------------------------------
!> Create weight vector for FE representation
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine fem_mat_zero_local(self,mloc)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
TYPE(oft_local_mat), INTENT(inout) :: mloc(:,:)
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
!---
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(ASSOCIATED(mloc(i,j)%m))mloc(i,j)%m=0.d0
    IF(ASSOCIATED(mloc(i,j)%ind))mloc(i,j)%ind(1,1)=0
  END DO
END DO
DEBUG_STACK_POP
end subroutine fem_mat_zero_local
!---------------------------------------------------------------------------
! SUBROUTINE: fem_mat_zero_local_rows
!---------------------------------------------------------------------------
!> Zero local contributions
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine fem_mat_zero_local_rows(self,mloc,flag,irow)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
TYPE(oft_local_mat), INTENT(inout) :: mloc(:,:)
LOGICAL, DIMENSION(:), INTENT(IN) :: flag
INTEGER(i4), INTENT(IN) :: irow
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
!---
DO i=1,self%fields(irow)%fe%nce
  IF(flag(i))THEN
    DO j=1,self%nfields
      IF(ASSOCIATED(mloc(irow,j)%m))mloc(irow,j)%m(i,:)=0.d0
    END DO
  END IF
END DO
DEBUG_STACK_POP
end subroutine fem_mat_zero_local_rows
!---------------------------------------------------------------------------
! SUBROUTINE: fem_mat_add_local
!---------------------------------------------------------------------------
!> Add local contributions to full matrix
!!
!! @param[in,out] mat Full matrix
!! @param[in] mloc Local matrix
!! @param[in] iloc Local FE entries
!! @param[in,out] tlocks OpenMP row thread locks
!---------------------------------------------------------------------------
subroutine fem_mat_add_local(self,mat,mloc,iloc,tlocks)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
CLASS(oft_matrix), INTENT(inout) :: mat
TYPE(oft_local_mat), INTENT(in) :: mloc(:,:)
TYPE(oft_1d_int), INTENT(IN) :: iloc(:)
INTEGER(KIND=omp_lock_kind), INTENT(INOUT) :: tlocks(:)
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
!---Add matrix components
DO i=1,self%nfields
  CALL omp_set_lock(tlocks(i))
  DO j=1,self%nfields
    IF(ASSOCIATED(mloc(i,j)%m))THEN
      IF(ASSOCIATED(mloc(i,j)%ind))THEN
        CALL mat%add_values(iloc(i)%v,iloc(j)%v, &
        mloc(i,j)%m,self%fields(i)%fe%nce,self%fields(j)%fe%nce,i,j,mloc(i,j)%ind)
      ELSE
        CALL mat%add_values(iloc(i)%v,iloc(j)%v, &
        mloc(i,j)%m,self%fields(i)%fe%nce,self%fields(j)%fe%nce,i,j)
      END IF
    END IF
  END DO
  CALL omp_unset_lock(tlocks(i))
END DO
DEBUG_STACK_POP
end subroutine fem_mat_add_local
!---------------------------------------------------------------------------
! SUBROUTINE: ml_fem_setup
!---------------------------------------------------------------------------
!> Create weight vector for FE representation
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine ml_fem_setup(self)
class(oft_ml_fem_comp_type), intent(inout) :: self
type(oft_fem_comp_type), pointer :: fe_tmp
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
DO j=1,self%nlevels
  ALLOCATE(self%levels(j)%fe)
  fe_tmp=>self%levels(j)%fe
  fe_tmp%nfields=self%nfields
  ALLOCATE(fe_tmp%fields(fe_tmp%nfields))
  ALLOCATE(fe_tmp%field_tags(fe_tmp%nfields))
  DO i=1,self%nfields
    fe_tmp%fields(i)%fe=>self%ml_fields(i)%ml%levels(j)%fe
    fe_tmp%field_tags(i)=self%field_tags(i)
  END DO
END DO
!---Set to highest level
self%level=self%nlevels
self%current_level=>self%levels(self%level)%fe
DEBUG_STACK_POP
end subroutine ml_fem_setup
!---------------------------------------------------------------------------
! SUBROUTINE: ml_fem_vec_create
!---------------------------------------------------------------------------
!> Create weight vector for FE representation
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine ml_fem_vec_create(self,new,level,cache,native)
class(oft_ml_fem_comp_type), intent(inout) :: self
class(oft_vector), pointer, intent(out) :: new
integer(i4), optional, intent(in) :: level
logical, optional, intent(in) :: cache
logical, optional, intent(in) :: native
logical :: do_cache,force_native
DEBUG_STACK_PUSH
do_cache=.TRUE.
force_native=.FALSE.
IF(PRESENT(cache))do_cache=cache
IF(PRESENT(native))force_native=native
!---Create vector on current or new level
IF(PRESENT(level))THEN
  IF(level>self%nlevels.OR.level<=0)CALL oft_abort('Invalid FE level change requested', &
                                                   'ml_fem_vec_create',__FILE__)
  CALL self%levels(level)%fe%vec_create(new,cache=cache,native=native)
ELSE
  CALL self%current_level%vec_create(new,cache=cache,native=native)
END IF
DEBUG_STACK_POP
end subroutine ml_fem_vec_create
!---------------------------------------------------------------------------
! SUBROUTINE: ml_fem_set_level
!---------------------------------------------------------------------------
!> Set the current level for a ML Compsite-FE structure
!!
!! @param[in] level Desired level
!---------------------------------------------------------------------------
subroutine ml_fem_set_level(self,level)
class(oft_ml_fem_comp_type), intent(inout) :: self
integer(i4), intent(in) :: level
DEBUG_STACK_PUSH
IF(level>self%nlevels.OR.level<=0)CALL oft_abort('Invalid FE level change requested', &
                                                 'ml_fem_set_level',__FILE__)
!---Update level
self%level=level
self%current_level=>self%levels(self%level)%fe
self%abs_level=self%level
IF(self%level>self%blevel.AND.self%blevel>0)self%abs_level=self%level-1
DEBUG_STACK_POP
end subroutine ml_fem_set_level
!---------------------------------------------------------------------------
! SUBROUTINE: ml_fem_build_interp
!---------------------------------------------------------------------------
!> Set the current level for lagrange finite elements
!!
!! @param[in] level Desired level
!---------------------------------------------------------------------------
subroutine ml_fem_build_interp(self,minlev)
class(oft_ml_fem_comp_type), intent(inout) :: self
integer(i4), optional, intent(in) :: minlev
class(oft_vector), pointer :: fvec,cvec
type(oft_graph_ptr), POINTER :: graphs(:,:)
type(oft_matrix_ptr), POINTER :: mats(:,:)
INTEGER(i4) :: i,j,levmin
DEBUG_STACK_PUSH
levmin=2
IF(PRESENT(minlev))levmin=minlev+1
OUTER: DO j=levmin,self%nlevels
!---------------------------------------------------------------------------
! Create composite matrix
!---------------------------------------------------------------------------
  CALL self%set_level(j)
  !---Specify child graphs
  ALLOCATE(graphs(self%nfields,self%nfields))
  DO i=1,self%nfields
    IF(.NOT.ASSOCIATED(self%ml_fields(i)%ml%interp_graphs(j)%g))THEN
      DEALLOCATE(graphs)
      CYCLE OUTER
    END IF
    graphs(i,i)%g=>self%ml_fields(i)%ml%interp_graphs(j)%g
  END DO
  !---Get coarse and fine vectors
  CALL self%vec_create(cvec,level=self%level-1)
  CALL self%vec_create(fvec)
  !---Construct matrix
  CALL create_matrix(self%interp_matrices(j)%m,graphs,fvec,cvec)
  DEALLOCATE(graphs)
!---------------------------------------------------------------------------
! Combine child matrices into composite matrix
!---------------------------------------------------------------------------
  !---Specify child matrices
  ALLOCATE(mats(self%nfields,self%nfields))
  DO i=1,self%nfields
    IF(.NOT.ASSOCIATED(self%ml_fields(i)%ml%interp_matrices(j)%m))THEN
      CALL oft_abort('Sub-field matrix not allocated.','ml_fem_build_interp',__FILE__)
    END IF
    mats(i,i)%m=>self%ml_fields(i)%ml%interp_matrices(j)%m
  END DO
  !---Combine matrices
  CALL combine_matrices(mats,self%nfields,self%nfields,self%interp_matrices(j)%m)
  DEALLOCATE(mats)
  CALL self%interp_matrices(j)%m%assemble
  !---Delete temporaries
  CALL cvec%delete
  CALL fvec%delete
  DEALLOCATE(cvec,fvec)
END DO OUTER
DEBUG_STACK_POP
end subroutine ml_fem_build_interp
END MODULE fem_composite
