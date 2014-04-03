set( GMOCK_GTEST_SOURCES
"gmock-gtest-all.cc"
"gmock_main.cc"
)

include_directories ( "." ) 
include_directories ( "./gtest" ) 

# generatring the library called gmock-gtest 

add_library( gmock-gtest ${GMOCK_GTEST_SOURCES} )


