cmake_minimum_required (VERSION 2.8)

project (StencilParareal)
enable_language (CXX)

# Include our directories
include_directories("src")
include_directories("src/Convection")
include_directories("src/Parareal")

# Boost
find_package (Boost 1.48.0 REQUIRED COMPONENTS program_options)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

# CUDA
set (CUDA_BACKEND "OFF" CACHE BOOL "Whether to use CUDA or not")

if ("${CUDA_BACKEND}")
    find_package(CUDA 5.0 REQUIRED)
    add_definitions( -D__CUDA_BACKEND__ )
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_21) # -Xptxas -abi=no
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GCL_GPU_")
else()
    # OpenMP
    find_package(OpenMP REQUIRED)
    add_definitions( -D__OPENMP_BACKEND__ )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
endif()

# STELLA
set (STELLA_PATH "" CACHE PATH
    "Directory where STELLA is stored"
    )
set (STELLA_SRC_PATH "${STELLA_PATH}/src" )
set (STELLA_BUILD_PATH "${STELLA_PATH}/build" )
include_directories (
    "${STELLA_SRC_PATH}/shared_definitions"
    "${STELLA_SRC_PATH}/shared_infrastructure"
    "${STELLA_SRC_PATH}/stencil_framework"
    )
set (SHARED_INFRASTRUCTURE_LIB
    "${STELLA_BUILD_PATH}/src/shared_infrastructure/libSharedInfrastructure.a"
    )

# GCL
set (GCL_PATH "${STELLA_PATH}/libs/gcl")
set(USE_MPI_COMPILER "OFF" CACHE BOOL "On CSCS machines turn this flag off since compiler takes care of mpi already" )
find_package(MPI REQUIRED)
add_definitions ("-D_GCL_MPI_" "-D__GCL__")
include_directories (
    "${GCL_PATH}"
    "${GCL_PATH}/utils"
    "${GCL_PATH}/L2/include"
    "${GCL_PATH}/L3/include"
    )
if ("${USE_MPI_COMPILER}" STREQUAL "ON" )
    INCLUDE(CMakeForceCompiler)
    CMAKE_FORCE_CXX_COMPILER(mpicxx "MPI C++ Compiler")
    include_directories(${MPI_CXX_INCLUDE_PATH})
endif()
if( "${CUDA_BACKEND}" STREQUAL "ON" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GCL_GPU_")
    include_directories(${CUDA_INCLUDE_DIRS})
endif()

# Timing
set (ENABLE_TIMING "ON" CACHE BOOL "Set this on to enable timing measurements, which slow down the solver")
if (${ENABLE_TIMING})
    add_definitions ("-DENABLE_TIMING")
endif()

# Error
set (ENABLE_ERROR "ON" CACHE BOOL "Set this on to enable computation of the error after every iteration")
if (${ENABLE_ERROR})
    add_definitions ("-DENABLE_ERROR")
endif()

# Create GCL library
add_library(gcl "${GCL_PATH}/L3/src/GCL.cpp")

# Create Communication library
include_directories("src/Communication")
if (${CUDA_BACKEND})
    cuda_add_library(communication "src/Communication/Periodicity.cu")
else()
    add_library(communication "src/Communication/Periodicity.cpp")
endif()

# Create matfile library
add_library(matfile "src/MatFile.cpp")

# Upwind stencil for advection
set(UPWINDADVECTION "OFF" CACHE BOOL "Set this on to enable second-order upwind advection")
if ("${UPWINDADVECTION}")
    add_definitions("-DUPWINDADVECTION")
endif("${UPWINDADVECTION}")


###############
# Executables #
###############

## Create executable heat
#set (HEATSOURCES
#    "src/Heat/main.cpp"
#    "src/Heat/MatFile.cpp"
#)
#
#if (${CUDA_BACKEND})
#    cuda_add_executable (heat ${HEATSOURCES} "src/Heat/Heat.cu")
#else()
#    add_executable (heat ${HEATSOURCES} "src/Heat/Heat.cpp")
#endif()
#
#target_link_libraries (heat gcl "${SHARED_INFRASTRUCTURE_LIB}" "${Boost_LIBRARIES}")


# Create executable convection
set (CONVECTIONSOURCES
    "src/Convection/main.cpp"
)

if (${CUDA_BACKEND})
    cuda_add_executable (convection ${CONVECTIONSOURCES} "src/Convection/Convection.cu")
    cuda_add_executable (sendrecv "src/Communication/SendRecvtest.cu")
else()
    add_executable (convection ${CONVECTIONSOURCES} "src/Convection/Convection.cpp")
    add_executable (sendrecv "src/Communication/SendRecvtest.cpp")
endif()

target_link_libraries (convection gcl communication matfile "${SHARED_INFRASTRUCTURE_LIB}" "${Boost_LIBRARIES}")


# Periodicity test
if (${CUDA_BACKEND})
    cuda_add_executable(testperiodicity "src/Communication/testPeriodicity.cu")
    target_link_libraries(testperiodicity communication matfile "${SHARED_INFRASTRUCTURE_LIB}" )
endif()

# Parareal
set (PARAREAL_SOURCES
    "src/Parareal/RuntimeConfiguration.cpp"
)

if (${CUDA_BACKEND})
    cuda_add_executable(parareal "src/Parareal/main.cu" "src/Convection/Convection.cu" "${PARAREAL_SOURCES}")
else()
    add_executable(parareal "src/Parareal/main.cpp" "src/Convection/Convection.cpp" "${PARAREAL_SOURCES}")
endif()
target_link_libraries (parareal matfile communication "${SHARED_INFRASTRUCTURE_LIB}" "${Boost_LIBRARIES}")

