set( TMAKER_EQ_DIRS
  fixed_boundary
  ARC
  CUTE
  DIIID
  Dipole
  HBT
  ITER
  JET
  LTX
  MANTA
  SPARC
  AdvancedWorkflows
)

foreach( path ${TMAKER_EQ_DIRS} )
  install(DIRECTORY ${path} DESTINATION examples/TokaMaker COMPONENT examples)
  if( OFT_BUILD_TESTS )
    file(COPY ${path} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endforeach()

if( OFT_BUILD_TESTS )
  # Create test driver
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/test_examples.sh.in ${CMAKE_CURRENT_BINARY_DIR}/test_examples.sh @ONLY )
  file(
    CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_examples.sh
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ
  )
endif()