# Top Level CMakeLists.txt File
cmake_minimum_required( VERSION 2.8 )

# set the installation prefix
# Note: setting this *before* the project imposes a default that the user is able to change.
SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install/" CACHE PATH "installation prefix")

project( STELLA )
enable_language( CXX )

set ( STELLA_ENABLE_TESTING "ON" CACHE BOOL "Enables the compilation of tets" )
if ( STELLA_ENABLE_TESTING )
    add_subdirectory( libs )
    include_directories( libs/gmock-gtest )
    enable_testing()
endif ( STELLA_ENABLE_TESTING )

# BOOST
find_package(Boost 1.43.0)
if( Boost_FOUND )
    include_directories(${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )

# CUDA
set( CUDA_BACKEND "OFF" CACHE BOOL "Compile CUDA backend" )
if( CUDA_BACKEND )
    find_package( CUDA 5.0 )
    
    set( CUDA_COMPUTE_CAPABILITY "sm_35" CACHE STRING "Compute capability for CUDA" )
    set( CUDA_SEPARABLE_COMPILATION "ON" )
    
    # Do to a well-known bug in cmake, CUDA_NVCC_FLAGS are not passed correctly when the separate
    # compilation is enabled.  We fill instead a variable called STELLA_NVCC_FLAGS,
    # which has to be passed manually to all compilations
    # set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=${CUDA_COMPUTE_CAPABILITY}") # -Xptxas -abi=no
    set (STELLA_CUDA_FLAGS "-D__CUDA_BACKEND__" "-arch=${CUDA_COMPUTE_CAPABILITY}" )

    # CUDA streams
    set ( ENABLE_CUDA_STREAMS "ON" CACHE BOOL "Enable usage of CUDA streams" )
    if ( ENABLE_CUDA_STREAMS )
        set (STELLA_CUDA_FLAGS "${STELLA_CUDA_FLAGS}" "-DENABLE_CUDA_STREAMS" )
    endif( ENABLE_CUDA_STREAMS )
endif( CUDA_BACKEND )

# OpenMP
set( OMP_BACKEND "ON" CACHE BOOL "Compile OpenMP backend" )
if( OMP_BACKEND )
    find_package(OpenMP)
    set (STELLA_OMP_FLAGS "-D__OPENMP_BACKEND__  ${OpenMP_CXX_FLAGS}")
else( OMP_BACKEND )
    set (STELLA_OMP_FLAGS "")
endif( OMP_BACKEND )

# Single precision
set( SINGLEPRECISION "OFF" CACHE BOOL "Single precision" )
if( SINGLEPRECISION )
  add_definitions( -DSINGLEPRECISION )
endif( SINGLEPRECISION )

# Logging
set ( LOGGING "OFF" CACHE BOOL "Enable output of log informations" )
if ( LOGGING )
    add_definitions ( -DLOGGING )
endif( LOGGING )

# Caching
set ( ENABLE_CACHING "ON" CACHE BOOL "Enable software-managed caching" )
if ( ENABLE_CACHING )
    add_definitions ( -DENABLE_CACHING )
endif( ENABLE_CACHING )

# CUDA streams
if( CUDA_BACKEND )
endif( CUDA_BACKEND )

# Produce position independent code on Linux
if( "${CMAKE_HOST_UNIX}" AND "${SHARED_LIBRARY}" ) 
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" )
endif()

# Build sources
add_subdirectory( src )

