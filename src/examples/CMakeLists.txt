# Add example executables
if( OFT_BUILD_EXAMPLES )
    add_subdirectory( TokaMaker )
    add_subdirectory( ThinCurr )
    add_subdirectory( MUG )
    add_subdirectory( Marklin )
    #
    oft_add_exe( example1.F90 )
    oft_add_exe( example2.F90 )
    oft_add_exe( example3.F90 )
    oft_add_exe( example4.F90 )
endif()