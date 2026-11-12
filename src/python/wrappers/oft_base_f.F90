!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_base_f.F90
!
!> @defgroup doxy_oft_python Python
!! Python interface for the Open FUSION Toolkit
!
!> Fortran part of Python wrappers for OFT
!!
!! @authors Chris Hansen
!! @date May 2023
!! @ingroup doxy_oft_python
!---------------------------------------------------------------------------------
MODULE oft_base_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_null_char, c_ptr, &
    c_f_pointer, c_bool, c_null_ptr, c_funptr, c_associated, c_f_procpointer
USE oft_base
USE oft_mesh_type, ONLY: oft_bmesh, oft_mesh, bmesh_findcell, mesh_findcell
USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct, multigrid_construct_surf
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr
USE fem_base, ONLY: oft_fem_type, oft_bfem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_interp, bfem_interp
USE oft_lag_basis, ONLY: oft_scalar_bfem
USE oft_blag_operators, ONLY: oft_lag_brinterp, oft_lag_bginterp
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_rinterp
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> Field evaluation object for FE fields exposed through the Python interface
!!
!! Holds a volume or surface interpolator, along with any vectors owned by the
!! interpolator, so that all FE fields can be evaluated and destroyed through
!! \ref oft_apply_field_eval and \ref oft_destroy_field_eval.
!---------------------------------------------------------------------------------
TYPE :: oftpy_field_eval
  CLASS(fem_interp), POINTER :: vinterp => NULL() !< Volume (3D) interpolator
  CLASS(bfem_interp), POINTER :: binterp => NULL() !< Surface (2D) interpolator
  TYPE(oft_vector_ptr), POINTER, DIMENSION(:) :: vecs => NULL() !< Vectors owned by the interpolator
END TYPE oftpy_field_eval
CONTAINS
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE copy_string(f_string,c_string)
CHARACTER(LEN=*), INTENT(in) :: f_string !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: c_string(*) !< Needs docs
INTEGER(4) :: i
DO i=1,LEN(f_string)
    c_string(i)=f_string(i:i)
END DO
c_string(LEN(f_string)+1)=c_null_char
END SUBROUTINE copy_string
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE copy_string_rev(c_string,f_string)
CHARACTER(KIND=c_char), INTENT(in) :: c_string(*) !< Needs docs
CHARACTER(LEN=*), INTENT(inout) :: f_string !< Needs docs
INTEGER(4) :: i
f_string=''
DO i=1,LEN(f_string)
    IF(c_string(i)==c_null_char)EXIT
    f_string(i:i)=c_string(i)
END DO
IF(i>LEN(f_string))CALL oft_warn("No termination character found when copying C string")
END SUBROUTINE copy_string_rev
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE oftpy_init(nthreads,quiet,input_file,slens,abort_fun) BIND(C,NAME="oftpy_init")
INTEGER(c_int), VALUE, INTENT(in) :: nthreads !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: quiet !< If `True`, do not print OFT environment information on initialization
CHARACTER(KIND=c_char), INTENT(in) :: input_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: slens !< String lengths
TYPE(c_funptr), VALUE, INTENT(in) :: abort_fun !< Abort callback for Python
INTEGER(4), POINTER, DIMENSION(:) :: slens_tmp
LOGICAL :: quiet_f
IF(oft_env%ifile/='none')RETURN
CALL copy_string_rev(input_file,oft_env%ifile)
quiet_f=quiet
CALL oft_init(nthreads=nthreads,quiet=quiet_f)
CALL c_f_pointer(slens, slens_tmp, [4])
IF(c_associated(abort_fun))CALL c_f_procpointer(abort_fun,oft_abort_cb)
slens_tmp=[OFT_MPI_PLEN,OFT_SLEN,OFT_PATH_SLEN,OFT_ERROR_SLEN]
END SUBROUTINE oftpy_init
!---------------------------------------------------------------------------------
!> Get string lengths used by the compiled library (does not require initialization)
!---------------------------------------------------------------------------------
SUBROUTINE oftpy_get_slens(slens) BIND(C,NAME="oftpy_get_slens")
TYPE(c_ptr), VALUE, INTENT(in) :: slens !< String lengths [OFT_MPI_PLEN,OFT_SLEN,OFT_PATH_SLEN,OFT_ERROR_SLEN] (4)
INTEGER(i4), POINTER, DIMENSION(:) :: slens_tmp
CALL c_f_pointer(slens, slens_tmp, [4])
slens_tmp=[OFT_MPI_PLEN,OFT_SLEN,OFT_PATH_SLEN,OFT_ERROR_SLEN]
END SUBROUTINE oftpy_get_slens
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE oftpy_load_xml(xml_file,oft_node_ptr) BIND(C,NAME="oftpy_load_xml")
CHARACTER(KIND=c_char), INTENT(in) :: xml_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), INTENT(out) :: oft_node_ptr !< Needs docs
INTEGER(i4) :: ierr
LOGICAL :: rst
CHARACTER(LEN=OFT_PATH_SLEN) :: xml_filename = 'none'
TYPE(xml_doc) :: doc
TYPE(xml_node), POINTER :: oft_node
!---Test for existence of XML file
CALL copy_string_rev(xml_file,xml_filename)
INQUIRE(FILE=TRIM(xml_filename),exist=rst)
IF(.NOT.rst)RETURN
CALL xml_parsefile(TRIM(xml_filename),doc,ierr)
IF(ierr/=0)RETURN
oft_node_ptr=C_LOC(doc%root)
END SUBROUTINE oftpy_load_xml
!---------------------------------------------------------------------------------
!> Set debug verbosity level
!---------------------------------------------------------------------------------
SUBROUTINE oftpy_set_debug(debug_level) BIND(C,NAME="oftpy_set_debug")
INTEGER(c_int), VALUE, INTENT(in) :: debug_level !< New value for debug level (must be in range [0,3])
oft_env%debug=debug_level
END SUBROUTINE oftpy_set_debug
!---------------------------------------------------------------------------------
!> Set the number of OpenMP threads to use
!---------------------------------------------------------------------------------
SUBROUTINE oftpy_set_nthreads(nthreads) BIND(C,NAME="oftpy_set_nthreads")
INTEGER(c_int), VALUE, INTENT(in) :: nthreads !< Number of threads to use for subsequent OpenMP parallel regions
CHARACTER(LEN=4) :: thrd_str,proc_str
IF(nthreads>omp_get_num_procs())THEN
    WRITE(thrd_str,'(I4.4)')nthreads
    WRITE(proc_str,'(I4.4)')omp_get_num_procs()
    CALL oft_warn("Number of requested threads ("//thrd_str//") exceeds number of available processors ("//proc_str//")")
END IF
CALL omp_set_num_threads(nthreads)
END SUBROUTINE oftpy_set_nthreads
!---------------------------------------------------------------------------------
!> Setup surface mesh from specified coordinate, connectivity, and region arrays
!---------------------------------------------------------------------------------
SUBROUTINE oft_setup_smesh(ndim,np,r_loc,npc,nc,lc_loc,reg_loc,nregs,mesh_ptr) BIND(C,NAME="oft_setup_smesh")
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc !< Pointer to point locations array `(3,np)`
TYPE(c_ptr), VALUE, INTENT(in) :: lc_loc !< Pointer to cell connectivity array `(npc,nc)`
TYPE(c_ptr), VALUE, INTENT(in) :: reg_loc !< Pointer to cell region array `(nc)`
INTEGER(c_int), VALUE, INTENT(in) :: np !< Number of points in surface mesh
INTEGER(c_int), VALUE, INTENT(in) :: ndim !< Spatial dimension of the mesh coordinates (Note: `r` array is always `(3,np)`)
INTEGER(c_int), VALUE, INTENT(in) :: npc !< Number of points per cell (3 for triangles, 4 for quads)
INTEGER(c_int), VALUE, INTENT(in) :: nc !< Number of cells in surface mesh
INTEGER(c_int), INTENT(out) :: nregs !< Number of regions in surface mesh
TYPE(c_ptr), INTENT(out) :: mesh_ptr !< Pointer to `multigrid_mesh` object containing the constructed surface mesh
integer(i4), POINTER :: lc_tmp(:,:),reg_tmp(:)
real(r8), POINTER :: r_tmp(:,:)
TYPE(multigrid_mesh), POINTER :: mg_mesh
IF(ndim>0)THEN
    ALLOCATE(r_mem(ndim,np))
    CALL c_f_pointer(r_loc, r_tmp, [ndim,np])
    r_mem=r_tmp
    ALLOCATE(lc_mem(npc,nc))
    CALL c_f_pointer(lc_loc, lc_tmp, [npc,nc])
    lc_mem=lc_tmp
    ALLOCATE(reg_mem(nc))
    CALL c_f_pointer(reg_loc, reg_tmp, [nc])
    reg_mem=reg_tmp
END IF
!---Setup Mesh
ALLOCATE(mg_mesh)
CALL multigrid_construct_surf(mg_mesh)
nregs=mg_mesh%smesh%nreg
mesh_ptr=C_LOC(mg_mesh)
END SUBROUTINE oft_setup_smesh
!---------------------------------------------------------------------------------
!> Return references to arrays defining an existing surface mesh.
!!
!! This routine exposes pointers to the coordinate, connectivity, and region
!! arrays stored inside the Fortran `multigrid_mesh` object referenced by
!! `mesh_ptr`. No data are copied and callers must not free or reallocate the objects.
!---------------------------------------------------------------------------------
SUBROUTINE oft_smesh_get(mesh_ptr,ndim,np,r_loc,npc,nc,lc_loc,reg_loc,nreg,error_str) BIND(C,NAME="oft_smesh_get")
TYPE(c_ptr), VALUE, INTENT(in) :: mesh_ptr !< Pointer to `multigrid_mesh` object containing the desired surface mesh
INTEGER(c_int), INTENT(out) :: ndim !< Spatial dimension of the mesh coordinates (Note: `r` array is always `(3,np)`)
INTEGER(c_int), INTENT(out) :: np !< Number of points in surface mesh
TYPE(c_ptr), INTENT(out) :: r_loc !< Pointer to point locations array `(3,np)`
INTEGER(c_int), INTENT(out) :: npc !< Number of points per cell (3 for triangles, 4 for quads)
INTEGER(c_int), INTENT(out) :: nc !< Number of cells in surface mesh
TYPE(c_ptr), INTENT(out) :: lc_loc !< Pointer to cell connectivity array `(npc,nc)`
TYPE(c_ptr), INTENT(out) :: reg_loc !< Pointer to cell region array `(nc)`
INTEGER(c_int), INTENT(out) :: nreg !< Number of regions in surface mesh
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error information
TYPE(multigrid_mesh), POINTER :: mg_mesh
IF(.NOT.c_associated(mesh_ptr))THEN
  CALL copy_string('Mesh object not associated',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
CALL c_f_pointer(mesh_ptr,mg_mesh)
!---Read values
ndim=mg_mesh%smesh%dim
np=mg_mesh%smesh%np
npc=mg_mesh%smesh%cell_np
nc=mg_mesh%smesh%nc
nreg=mg_mesh%smesh%nreg
!---Get pointers
r_loc=C_LOC(mg_mesh%smesh%r)
lc_loc=C_LOC(mg_mesh%smesh%lc)
reg_loc=C_LOC(mg_mesh%smesh%reg)
END SUBROUTINE oft_smesh_get
!---------------------------------------------------------------------------------
!> Setup volume mesh from specified coordinate, connectivity, and region arrays
!---------------------------------------------------------------------------------
SUBROUTINE oft_setup_vmesh(np,r_loc,npc,nc,lc_loc,reg_loc,nregs,mesh_ptr) BIND(C,NAME="oft_setup_vmesh")
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc !< Pointer to point locations array `(3,np)`
TYPE(c_ptr), VALUE, INTENT(in) :: lc_loc !< Pointer to cell connectivity array `(npc,nc)`
TYPE(c_ptr), VALUE, INTENT(in) :: reg_loc !< Pointer to cell region array `(nc)`
INTEGER(c_int), VALUE, INTENT(in) :: np !< Number of points in volume mesh
INTEGER(c_int), VALUE, INTENT(in) :: npc !< Number of points per cell (4 for tetrahedra, 8 for hexahedra)
INTEGER(c_int), VALUE, INTENT(in) :: nc !< Number of cells in volume mesh
INTEGER(c_int), INTENT(out) :: nregs !< Number of regions in volume mesh
TYPE(c_ptr), INTENT(out) :: mesh_ptr !< Pointer to `multigrid_mesh` object containing the constructed volume mesh
integer(i4), POINTER :: lc_tmp(:,:),reg_tmp(:)
real(r8), POINTER :: r_tmp(:,:)
TYPE(multigrid_mesh), POINTER :: mg_mesh
IF(np>0)THEN
    ALLOCATE(r_mem(3,np))
    CALL c_f_pointer(r_loc, r_tmp, [3,np])
    r_mem=r_tmp
    ALLOCATE(lc_mem(npc,nc))
    CALL c_f_pointer(lc_loc, lc_tmp, [npc,nc])
    lc_mem=lc_tmp
    ALLOCATE(reg_mem(nc))
    CALL c_f_pointer(reg_loc, reg_tmp, [nc])
    reg_mem=reg_tmp
END IF
!---Setup Mesh
ALLOCATE(mg_mesh)
CALL multigrid_construct(mg_mesh)
nregs=mg_mesh%mesh%nreg
mesh_ptr=C_LOC(mg_mesh)
END SUBROUTINE oft_setup_vmesh
!---------------------------------------------------------------------------------
!> Return references to arrays defining an existing volume mesh.
!!
!! This routine exposes pointers to the coordinate, connectivity, and region
!! arrays stored inside the Fortran `multigrid_mesh` object referenced by
!! `mesh_ptr`. No data are copied and callers must not free or reallocate the objects.
!---------------------------------------------------------------------------------
SUBROUTINE oft_vmesh_get(mesh_ptr,np,r_loc,npc,nc,lc_loc,reg_loc,nreg,error_str) BIND(C,NAME="oft_vmesh_get")
TYPE(c_ptr), VALUE, INTENT(in) :: mesh_ptr !< Pointer to `multigrid_mesh` object containing the desired volume mesh
INTEGER(c_int), INTENT(out) :: np !< Number of points in volume mesh
TYPE(c_ptr), INTENT(out) :: r_loc !< Pointer to point locations array `(3,np)`
INTEGER(c_int), INTENT(out) :: npc !< Number of points per cell (4 for tetrahedra, 8 for hexahedra)
INTEGER(c_int), INTENT(out) :: nc !< Number of cells in volume mesh
TYPE(c_ptr), INTENT(out) :: lc_loc !< Pointer to cell connectivity array `(npc,nc)`
TYPE(c_ptr), INTENT(out) :: reg_loc !< Pointer to cell region array `(nc)`
INTEGER(c_int), INTENT(out) :: nreg !< Number of regions in volume mesh
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error information
TYPE(multigrid_mesh), POINTER :: mg_mesh
IF(.NOT.c_associated(mesh_ptr))THEN
  CALL copy_string('Mesh object not associated',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
CALL c_f_pointer(mesh_ptr,mg_mesh)
!---Read values
np=mg_mesh%mesh%np
npc=mg_mesh%mesh%cell_np
nc=mg_mesh%mesh%nc
nreg=mg_mesh%mesh%nreg
!---Get pointers
r_loc=C_LOC(mg_mesh%mesh%r)
lc_loc=C_LOC(mg_mesh%mesh%lc)
reg_loc=C_LOC(mg_mesh%mesh%reg)
END SUBROUTINE oft_vmesh_get
!---------------------------------------------------------------------------------
!> Wrap an FE interpolator in a field evaluation object that can be passed to
!! \ref oft_apply_field_eval and \ref oft_destroy_field_eval
!!
!! This allows tool-specific interpolators (eg. Marklin's composite magnetic field)
!! to share the generic evaluation and destruction interface. Ownership of the
!! interpolator and of any vectors in `vecs` is transferred to the evaluation object,
!! so they are deleted and deallocated by \ref oft_destroy_field_eval.
!---------------------------------------------------------------------------------
SUBROUTINE oft_field_eval_wrap(int_obj,vinterp,binterp,vecs)
TYPE(c_ptr), INTENT(out) :: int_obj !< Pointer to field evaluation object
CLASS(fem_interp), OPTIONAL, TARGET, INTENT(inout) :: vinterp !< Volume (3D) interpolator (must be allocated)
CLASS(bfem_interp), OPTIONAL, TARGET, INTENT(inout) :: binterp !< Surface (2D) interpolator (must be allocated)
TYPE(oft_vector_ptr), OPTIONAL, INTENT(in) :: vecs(:) !< Vectors owned by interpolator (must be allocated)
INTEGER(i4) :: i
TYPE(oftpy_field_eval), POINTER :: eval_obj
IF(PRESENT(vinterp).EQV.PRESENT(binterp))CALL oft_abort("Exactly one of vinterp or binterp is required", &
  "oft_field_eval_wrap",__FILE__)
ALLOCATE(eval_obj)
IF(PRESENT(vinterp))eval_obj%vinterp=>vinterp
IF(PRESENT(binterp))eval_obj%binterp=>binterp
IF(PRESENT(vecs))THEN
  ALLOCATE(eval_obj%vecs(SIZE(vecs)))
  DO i=1,SIZE(vecs)
    eval_obj%vecs(i)%f=>vecs(i)%f
  END DO
END IF
int_obj=C_LOC(eval_obj)
END SUBROUTINE oft_field_eval_wrap
!---------------------------------------------------------------------------------
!> Convert C pointer to a field evaluation object, checking that it is associated
!---------------------------------------------------------------------------------
FUNCTION oft_field_eval_ccast(int_obj,eval_obj,error_str) RESULT(success)
TYPE(c_ptr), INTENT(in) :: int_obj !< C pointer to field evaluation object
TYPE(oftpy_field_eval), POINTER, INTENT(out) :: eval_obj !< Fortran pointer to field evaluation object
CHARACTER(KIND=c_char), OPTIONAL, INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
LOGICAL :: success
IF(PRESENT(error_str))CALL copy_string('',error_str)
IF(.NOT.c_associated(int_obj))THEN
  IF(PRESENT(error_str))CALL copy_string('Field evaluation object not associated',error_str)
  NULLIFY(eval_obj)
  success=.FALSE.
  RETURN
END IF
CALL c_f_pointer(int_obj,eval_obj)
success=.TRUE.
END FUNCTION oft_field_eval_ccast
!---------------------------------------------------------------------------------
!> Create a field evaluation object for a general FE field from a weight vector
!!
!! Supported field types (`imode`):
!! - `211`: Surface Lagrange scalar field (requires `lag_rep`)
!! - `212`: Gradient of surface Lagrange scalar field (requires `lag_rep`)
!! - `331`: Volume H(Curl) + Grad(H^1) vector field (requires `hcurl_grad_rep`)
!!
!! Ownership of `u` is transferred to the evaluation object (`u` is nullified on output).
!---------------------------------------------------------------------------------
SUBROUTINE oft_field_eval_setup(imode,u,int_obj,error_str,lag_rep,hcurl_grad_rep)
INTEGER(i4), INTENT(in) :: imode !< Field type
CLASS(oft_vector), POINTER, INTENT(inout) :: u !< FE weight vector for field
TYPE(c_ptr), INTENT(out) :: int_obj !< Pointer to field evaluation object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
CLASS(oft_scalar_bfem), OPTIONAL, TARGET, INTENT(inout) :: lag_rep !< Surface Lagrange FE representation
TYPE(oft_fem_comp_type), OPTIONAL, TARGET, INTENT(inout) :: hcurl_grad_rep !< Composite H(Curl) + Grad(H^1) FE representation
TYPE(oft_vector_ptr) :: vecs(1)
TYPE(oft_lag_brinterp), POINTER :: br_obj
TYPE(oft_lag_bginterp), POINTER :: bg_obj
TYPE(oft_hcurl_grad_rinterp), POINTER :: hcg_obj
CALL copy_string('',error_str)
int_obj=c_null_ptr
vecs(1)%f=>u
SELECT CASE(imode)
CASE(211)
  IF(.NOT.PRESENT(lag_rep))THEN
    CALL copy_string('Lagrange FE representation required for field type',error_str)
    RETURN
  END IF
  ALLOCATE(br_obj)
  br_obj%u=>u
  CALL br_obj%setup(lag_rep)
  CALL oft_field_eval_wrap(int_obj,binterp=br_obj,vecs=vecs)
CASE(212)
  IF(.NOT.PRESENT(lag_rep))THEN
    CALL copy_string('Lagrange FE representation required for field type',error_str)
    RETURN
  END IF
  ALLOCATE(bg_obj)
  bg_obj%u=>u
  CALL bg_obj%setup(lag_rep)
  CALL oft_field_eval_wrap(int_obj,binterp=bg_obj,vecs=vecs)
CASE(331)
  IF(.NOT.PRESENT(hcurl_grad_rep))THEN
    CALL copy_string('H(Curl) + Grad(H^1) FE representation required for field type',error_str)
    RETURN
  END IF
  ALLOCATE(hcg_obj)
  hcg_obj%u=>u
  CALL hcg_obj%setup(hcurl_grad_rep)
  CALL oft_field_eval_wrap(int_obj,vinterp=hcg_obj,vecs=vecs)
CASE DEFAULT
  CALL copy_string('Invalid field type for interpolation',error_str)
  RETURN
END SELECT
NULLIFY(u)
END SUBROUTINE oft_field_eval_setup
!---------------------------------------------------------------------------------
!> Create a field evaluation object for a general FE field
!!
!! See \ref oft_field_eval_setup for supported field types and the corresponding
!! type of FE object that `fe_ptr` must reference.
!---------------------------------------------------------------------------------
SUBROUTINE oft_get_field_eval(fe_ptr,vals,imode,int_obj,error_str) BIND(C,NAME="oft_get_field_eval")
TYPE(c_ptr), VALUE, INTENT(in) :: fe_ptr !< Pointer to FE object
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< FE weight values for field [ne]
INTEGER(KIND=c_int), VALUE, INTENT(in) :: imode !< Field type
TYPE(c_ptr), INTENT(out) :: int_obj !< Pointer to field evaluation object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CLASS(oft_vector), POINTER :: u
TYPE(oft_scalar_bfem), POINTER :: lag_rep
TYPE(oft_fem_comp_type), POINTER :: hcurl_grad_rep
int_obj=c_null_ptr
IF(.NOT.c_associated(fe_ptr))THEN
  CALL copy_string('FE object not associated',error_str)
  RETURN
END IF
SELECT CASE(imode)
CASE(211,212)
  CALL c_f_pointer(fe_ptr, lag_rep)
  CALL lag_rep%vec_create(u)
  CALL c_f_pointer(vals, vals_tmp, [u%n])
  CALL u%restore_local(vals_tmp)
  CALL oft_field_eval_setup(imode,u,int_obj,error_str,lag_rep=lag_rep)
CASE(331)
  CALL c_f_pointer(fe_ptr, hcurl_grad_rep)
  CALL hcurl_grad_rep%vec_create(u)
  CALL c_f_pointer(vals, vals_tmp, [u%n])
  CALL u%restore_local(vals_tmp)
  CALL oft_field_eval_setup(imode,u,int_obj,error_str,hcurl_grad_rep=hcurl_grad_rep)
CASE DEFAULT
  CALL copy_string('Invalid field type for interpolation',error_str)
  RETURN
END SELECT
END SUBROUTINE oft_get_field_eval
!---------------------------------------------------------------------------------
!> Evaluate a field with a field evaluation object created by \ref oft_get_field_eval
!! or a tool-specific constructor using \ref oft_field_eval_wrap
!!
!! Points outside the mesh, or whose logical coordinates fall outside the cell by more
!! than `fbary_tol`, are returned as NaN.
!---------------------------------------------------------------------------------
SUBROUTINE oft_apply_field_eval(int_obj,pt_ptr,npts,fbary_tol,dim,cell,field_ptr) BIND(C,NAME="oft_apply_field_eval")
TYPE(c_ptr), VALUE, INTENT(in) :: int_obj !< Pointer to field evaluation object
TYPE(c_ptr), VALUE, INTENT(in)  :: pt_ptr !< Locations for evaluation [3,npts], as [X,Y,Z] (3D) or [X,Y,0] (2D)
INTEGER(c_int), VALUE, INTENT(in) :: npts !< Number of points to evaluate
REAL(c_double), VALUE, INTENT(in) :: fbary_tol !< Tolerance for physical to logical mapping
INTEGER(c_int), VALUE, INTENT(in) :: dim !< Dimension of field
INTEGER(c_int), INTENT(inout) :: cell !< Guess for cell containing first point on input (`<=0` for none);
!! cell containing last point on output (`0` if not found, negative if outside of `fbary_tol`)
TYPE(c_ptr), VALUE, INTENT(in) :: field_ptr !< Field at locations specified by `pt_ptr` [dim,npts]
LOGICAL :: vol_eval
INTEGER(i4) :: i,nf,icell,guess
REAL(r8) :: f(4),goptmp(3,4),vol
REAL(r8), POINTER, DIMENSION(:,:) :: pts,field
TYPE(oftpy_field_eval), POINTER :: eval_obj
IF(.NOT.oft_field_eval_ccast(int_obj,eval_obj))CALL oft_abort("Field evaluation object not associated","oft_apply_field_eval",__FILE__)
vol_eval=ASSOCIATED(eval_obj%vinterp)
IF(vol_eval)THEN
  nf=4
ELSE IF(ASSOCIATED(eval_obj%binterp))THEN
  nf=3
ELSE
  CALL oft_abort("Interpolator not associated","oft_apply_field_eval",__FILE__)
END IF
!---Sample fields
CALL c_f_pointer(pt_ptr, pts, [3,npts])
CALL c_f_pointer(field_ptr, field, [dim,npts])
icell=cell
!$omp parallel do firstprivate(icell) lastprivate(icell) private(guess,f,goptmp,vol) if(npts>1000)
DO i=1,npts
  !---Use previous cell as initial guess, falling back to full search if not found
  guess=ABS(icell)
  icell=guess
  CALL locate_point(pts(:,i),icell,f)
  IF((guess>0).AND.(.NOT.in_cell(icell,f)))THEN
    icell=0
    CALL locate_point(pts(:,i),icell,f)
  END IF
  IF(icell==0)THEN
    field(:,i)=ieee_value(1.d0, ieee_quiet_nan)
    CYCLE
  END IF
  IF(.NOT.in_cell(icell,f))THEN
    icell=-icell
    field(:,i)=ieee_value(1.d0, ieee_quiet_nan)
    CYCLE
  END IF
  IF(vol_eval)THEN
    CALL eval_obj%vinterp%mesh%jacobian(icell,f,goptmp,vol)
    CALL eval_obj%vinterp%interp(icell,f,goptmp,field(:,i))
  ELSE
    CALL eval_obj%binterp%mesh%jacobian(icell,f(1:3),goptmp,vol)
    CALL eval_obj%binterp%interp(icell,f(1:3),goptmp(:,1:3),field(:,i))
  END IF
END DO
cell=icell
CONTAINS
!> Locate point in mesh (thread-local quantities are passed explicitly for OpenMP safety)
SUBROUTINE locate_point(pt,jcell,fl)
REAL(r8), INTENT(in) :: pt(3) !< Point to locate
INTEGER(i4), INTENT(inout) :: jcell !< Cell guess on input, containing cell on output (`0` if not found)
REAL(r8), INTENT(out) :: fl(4) !< Logical position in cell
fl=0.d0
IF(vol_eval)THEN
  CALL mesh_findcell(eval_obj%vinterp%mesh,jcell,pt,fl)
ELSE
  CALL bmesh_findcell(eval_obj%binterp%mesh,jcell,pt(1:2),fl(1:3))
END IF
END SUBROUTINE locate_point
!> Check if logical position lies in cell within tolerance
LOGICAL FUNCTION in_cell(jcell,fl)
INTEGER(i4), INTENT(in) :: jcell !< Cell index (`0` if not found)
REAL(r8), INTENT(in) :: fl(4) !< Logical position in cell
in_cell=(jcell/=0).AND.(MAXVAL(fl(1:nf))<=1.d0+fbary_tol).AND.(MINVAL(fl(1:nf))>=-fbary_tol)
END FUNCTION in_cell
END SUBROUTINE oft_apply_field_eval
!---------------------------------------------------------------------------------
!> Destroy a field evaluation object, including its interpolator and any owned vectors
!---------------------------------------------------------------------------------
SUBROUTINE oft_destroy_field_eval(int_obj) BIND(C,NAME="oft_destroy_field_eval")
TYPE(c_ptr), VALUE, INTENT(in) :: int_obj !< Pointer to field evaluation object
INTEGER(i4) :: i
TYPE(oftpy_field_eval), POINTER :: eval_obj
IF(.NOT.oft_field_eval_ccast(int_obj,eval_obj))RETURN
IF(ASSOCIATED(eval_obj%vinterp))THEN
  CALL eval_obj%vinterp%delete()
  DEALLOCATE(eval_obj%vinterp)
END IF
IF(ASSOCIATED(eval_obj%binterp))THEN
  CALL eval_obj%binterp%delete()
  DEALLOCATE(eval_obj%binterp)
END IF
IF(ASSOCIATED(eval_obj%vecs))THEN
  DO i=1,SIZE(eval_obj%vecs)
    IF(ASSOCIATED(eval_obj%vecs(i)%f))THEN
      CALL eval_obj%vecs(i)%f%delete()
      DEALLOCATE(eval_obj%vecs(i)%f)
    END IF
  END DO
  DEALLOCATE(eval_obj%vecs)
END IF
DEALLOCATE(eval_obj)
END SUBROUTINE oft_destroy_field_eval
END MODULE oft_base_f
