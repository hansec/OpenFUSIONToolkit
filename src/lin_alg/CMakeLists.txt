set( LA_OBJS
  lin_alg_base.F90
  native_la.F90
  petsc_la.F90
  deriv_matrices.F90
)
set( SOLVER_OBJS
  solver_base.F90
  native_solvers.F90
  petsc_solvers.F90
  superlu_bridge.c
  superlu_dist_bridge.c
  umfpack_bridge.c
  strumpack_bridge.cxx
  mkl_bridge.c
  oft_arpack.F90
  oft_lu.F90
)
set( MGMT_OBJS
  lin_alg_utils.F90
  solver_utils.F90
)

# Accelerate sparse interface is compiled with Apple clang, and pre-linked with the
# clang runtime (availability checks), as the Accelerate headers require clang extensions
if( OFT_ACCELERATE_SPARSE_FOUND )
  set( ACC_BRIDGE_OBJ ${CMAKE_CURRENT_BINARY_DIR}/accelerate_bridge.o )
  add_custom_command(
    OUTPUT ${ACC_BRIDGE_OBJ}
    COMMAND ${OFT_APPLE_CLANG} ${OFT_APPLE_CLANG_FLAGS} -O2 -Wall -DHAVE_ACCELERATE_SPARSE
      -c ${CMAKE_CURRENT_SOURCE_DIR}/accelerate_bridge.c -o accelerate_bridge_c.o
    COMMAND ${OFT_APPLE_LD} -r accelerate_bridge_c.o ${OFT_APPLE_CLANG_RT} -o ${ACC_BRIDGE_OBJ}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/accelerate_bridge.c
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Building Accelerate sparse interface with Apple clang"
    VERBATIM
  )
  set_source_files_properties( ${ACC_BRIDGE_OBJ} PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE )
  set( SOLVER_OBJS ${SOLVER_OBJS} ${ACC_BRIDGE_OBJ} )
endif()

add_library( oftla
  ${LA_OBJS}
  ${SOLVER_OBJS}
  ${MGMT_OBJS}
)

# STRUMPACK C++ interface requires C++17
if( OFT_STRUMPACK_FOUND )
  target_compile_features( oftla PRIVATE cxx_std_17 )
  # C++ sources are not compiled with the MPI compiler wrappers
  if( OFT_STRUMPACK_MPI )
    target_include_directories( oftla PRIVATE ${MPI_C_INCLUDE_DIRS} )
  endif()
endif()
