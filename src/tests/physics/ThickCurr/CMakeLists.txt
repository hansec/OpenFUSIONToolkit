# ThickCurr tests
oft_add_test( test_thickcurr.F90 )
configure_file( thickcurr_testing.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
configure_file( test_ThickCurr_box.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
configure_file( test_ThickCurr_slab.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
