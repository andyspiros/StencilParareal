# Make headers available to all compilations
include_directories( "shared_definitions" "shared_infrastructure" "stencil_framework" )


# Add sub directories
add_subdirectory( shared_definitions )
add_subdirectory( shared_infrastructure )
add_subdirectory( stencil_framework )

# Add testing
if( STELLA_ENABLE_TESTING )
    add_subdirectory( shared_infrastructure_unittest )
    add_subdirectory( stencil_framework_unittest )
endif( STELLA_ENABLE_TESTING )

