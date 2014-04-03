cmake_minimum_required (VERSION 2.8)

project (StencilParareal)
enable_language (CXX)

# Include our directories
include_directories("src")
include_directories("src/Convection")
include_directories("src/Parareal")

# Boost
set(Boost_USE_STATIC_LIBS ON)
find_package (Boost 1.48.0 REQUIRED COMPONENTS program_options)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

# CUDA
set (CUDA_BACKEND "OFF" CACHE BOOL "Whether to use CUDA or not")

if ( CUDA_BACKEND )
    find_package(CUDA 5.0 REQUIRED)
    add_definitions( -D__CUDA_BACKEND__ )
    include_directories(${CUDA_INCLUDE_DIRS})
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_35) # -Xptxas -abi=no
else( CUDA_BACKEND )
    # OpenMP
    find_package(OpenMP REQUIRED)
    add_definitions( -D__OPENMP_BACKEND__ )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
endif( CUDA_BACKEND )

# STELLA
set (STELLA_PATH "" CACHE PATH
    "Directory where STELLA is installed"
    )
set ( STELLA_PATH "STELLA-NOTFOUND" CACHE STRING "STELLA install path" )
if ( STELLA_PATH )
    get_filename_component(STELLA_PATH_ABSOLUTE "${CMAKE_BINARY_DIR}/${STELLA_PATH}"
        ABSOLUTE
    )
endif ( STELLA_PATH )
include_directories (
    "${STELLA_PATH_ABSOLUTE}/include/STELLA"
)

# STELLA libraries
if ( CUDA_BACKEND )
    set ( STELLA_SharedInfrastructure_NAME "StellaSharedInfrastructureCUDA" )
else ( CUDA_BACKEND )
    set ( STELLA_SharedInfrastructure_NAME "StellaSharedInfrastructure" )
endif ( CUDA_BACKEND )

set (STELLA_SharedInfrastructure_LIBRARY "STELLA_SharedInfrastructure_LIBRARY-NOTFOUND")
find_library ( STELLA_SharedInfrastructure_LIBRARY
    NAMES "${STELLA_SharedInfrastructure_NAME}"
    PATHS "${STELLA_PATH_ABSOLUTE}/lib"
)

if ( STELLA_SharedInfrastructure_LIBRARY )
    add_library( SharedInfrastructure STATIC IMPORTED )
    set_property( TARGET SharedInfrastructure PROPERTY IMPORTED_LOCATION "${STELLA_SharedInfrastructure_LIBRARY}" )
endif ( STELLA_SharedInfrastructure_LIBRARY )

# MPI
set(USE_MPI_COMPILER "OFF" CACHE BOOL "On Cray machines turn this flag off since compiler takes care of mpi already" )
find_package(MPI REQUIRED)

if (USE_MPI_COMPILER)
    INCLUDE(CMakeForceCompiler)
    CMAKE_FORCE_CXX_COMPILER(mpicxx "MPI C++ Compiler")
    include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

# Timing
set (ENABLE_INTERNAL_TIMING "OFF" CACHE BOOL "Set this on to enable timing measurements, which slow down the solver")
if ( ENABLE_INTERNAL_TIMING )
    add_definitions ("-DENABLE_INTERNAL_TIMING")
endif( ENABLE_INTERNAL_TIMING )

# Error
set (ENABLE_ERROR "OFF" CACHE BOOL "Set this on to enable computation of the error after every iteration")
if (ENABLE_ERROR)
    add_definitions ("-DENABLE_ERROR")
endif()

# Create Communication library
include_directories("src/Communication")
if (CUDA_BACKEND)
    cuda_add_library(communication "src/Communication/Periodicity.cu")
else()
    add_library(communication "src/Communication/Periodicity.cpp")
endif()

# Create matfile library
add_library(matfile "src/MatFile.cpp")

# Upwind stencil for advection
set(UPWINDADVECTION "ON" CACHE BOOL "Set this on to enable second-order upwind advection")
if ( UPWINDADVECTION )
    add_definitions("-DUPWINDADVECTION")
endif(UPWINDADVECTION)

#########
# Tests #
#########

#if (CUDA_BACKEND)
#    cuda_add_executable(startConvectionGPU "tests/startConvectionGPU.cu" "src/Convection/Convection.cu")
#    target_link_libraries (startConvectionGPU communication SharedInfrastructure )
#endif()

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
#target_link_libraries (heat SharedInfrastructure "${Boost_LIBRARIES}")


# Convection executable
if (CUDA_BACKEND)
    cuda_add_executable (convection "src/Convection/main.cpp" "src/Convection/Convection.cu")
else()
    add_executable (convection "src/Convection/main.cpp" "src/Convection/Convection.cpp")
endif()

target_link_libraries (convection communication matfile SharedInfrastructure "${Boost_LIBRARIES}")


# Periodicity test
#if (${CUDA_BACKEND})
#    cuda_add_executable(testperiodicity "src/Communication/testPeriodicity.cu")
#    target_link_libraries(testperiodicity communication matfile SharedInfrastructure )
#endif()

# Parareal
set (PARAREAL_SOURCES
    "src/Parareal/RuntimeConfiguration.cpp"
)

if (CUDA_BACKEND)
    cuda_add_executable(parareal "src/Parareal/main.cu" "src/Convection/Convection.cu" "${PARAREAL_SOURCES}")
else()
    add_executable(parareal "src/Parareal/main.cpp" "src/Convection/Convection.cpp" "${PARAREAL_SOURCES}")
endif()
target_link_libraries (parareal matfile communication SharedInfrastructure "${Boost_LIBRARIES}")

