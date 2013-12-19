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

# GCL
set (GCL "OFF" CACHE BOOL "Whether to use GCL or not")
if ( "${GCL}" STREQUAL "ON" )
    set (GCL_PATH "${STENCIL_LIBRARY_SRC_PATH}/../libs/gcl")
    set(USE_MPI_COMPILER "OFF" CACHE BOOL "On rosa turn this flag off since compiler takes care of mpi already" )
    find_package(MPI REQUIRED)
    add_definitions ("-D_GCL_MPI_" "-D__GCL__")
    include_directories (
        "${GCL_PATH}"
        "${GCL_PATH}/utils"
        "${GCL_PATH}/L2/include"
        "${GCL_PATH}/L3/include"
        )
    if( "${USE_MPI_COMPILER}" STREQUAL "ON" )
        INCLUDE(CMakeForceCompiler)
        CMAKE_FORCE_CXX_COMPILER(mpicxx "MPI C++ Compiler")
            include_directories(${MPI_CXX_INCLUDE_PATH})
    endif()
    if( "${CUDA_BACKEND}" STREQUAL "ON" )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GCL_GPU_")
        include_directories(${CUDA_INCLUDE_DIRS})
    endif()

    add_executable (halo
        "HaloExchange3Dtest.cpp"
        "${GCL_PATH}/L3/src/GCL.cpp"
        "MatFile.cpp"
    )
    target_link_libraries (halo ${SHARED_INFRASTRUCTURE_LIB})
endif()



# Create executable
set (SOURCES
    "main.cpp"
    "Heat.cpp"
    "MatFile.cpp"
    "${GCL_PATH}/L3/src/GCL.cpp"
)

add_executable (heat ${SOURCES})

# Link against shared infrastructure library for debug
target_link_libraries (heat ${SHARED_INFRASTRUCTURE_LIB})
