cmake_minimum_required (VERSION 2.8)

project (HeatEquation)
enable_language (CXX)

# Boost
find_package (Boost 1.43.0)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

# CUDA
set (CUDA_BACKEND "OFF" CACHE BOOL "Whether to use CUDA or not")

if ( "${CUDA_BACKEND}" STREQUAL "ON" )
    add_definitions("-DCUDA_BACKEND")
    find_package( CUDA 5.0 )
    add_definitions( -D__CUDA_BACKEND__ )
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_21) # -Xptxas -abi=no
else()
    # OpenMP
    find_package(OpenMP)
    add_definitions( -D__OPENMP_BACKEND__ )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
endif()

# Stencil library
set (STENCIL_LIBRARY_SRC_PATH "" CACHE PATH 
    "Source directory of the stencil library"
    )
set (STENCIL_LIBRARY_BUILD_PATH "" CACHE PATH 
    "Build directory of the stencil library"
    )
include_directories (
    "${STENCIL_LIBRARY_SRC_PATH}/shared_definitions"
    "${STENCIL_LIBRARY_SRC_PATH}/shared_infrastructure"
    "${STENCIL_LIBRARY_SRC_PATH}/stencil_framework"
    )
set (SHARED_INFRASTRUCTURE_LIB
    "${STENCIL_LIBRARY_BUILD_PATH}/src/shared_infrastructure/libSharedInfrastructure.a"
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
