cmake_minimum_required (VERSION 2.8)

project (HeatEquation)
enable_language (CXX)

# Boost
find_package (Boost 1.43.0)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

# OpenMP
find_package(OpenMP)
add_definitions( -D__OPENMP_BACKEND__ )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")

# Stencil library
set (STENCIL_LIBRARY_PATH "" CACHE PATH 
    "Directory where the file README of the stencil library is stored"
    )
include_directories (
    "${STENCIL_LIBRARY_PATH}/src/shared_definitions"
    "${STENCIL_LIBRARY_PATH}/src/shared_infrastructure"
    "${STENCIL_LIBRARY_PATH}/src/stencil_framework"
    )
set (SHARED_INFRASTRUCTURE_LIB
    "${STENCIL_LIBRARY_PATH}/build/lib/libSharedInfrastructure.a"
    )


# Create executable
set (SOURCES
    "main.cpp"
    "Heat.cpp"
    "MatFile.cpp"
)

add_executable (heat ${SOURCES})

# Link against shared infrastructure library for debug
target_link_libraries (heat ${SHARED_INFRASTRUCTURE_LIB})
