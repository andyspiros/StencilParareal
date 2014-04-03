#pragma once

// define the floating point type
#ifdef SINGLEPRECISION
typedef float Real;
#else
typedef double Real;
#endif
typedef const Real ConstReal;

// define upper limits for some library constructs
#define MAX_TUPLE_SIZE 35
#define MAX_ARRAY_SIZE 10
#define MAX_GROUP_COUNT 30 
#define MAX_LOOP_COUNT 20
#define MAX_SWITCH_COUNT 10
#define MAX_CASE_COUNT 10
#define MAX_MERGE_COUNT 10
#define MAX_PARAM_COUNT 20
#define MAX_TEMP_COUNT 20
#define MAX_CACHE_COUNT 20
#define MAX_FIELD_DIMENSION 3
#define MAX_BOUNDARY_LEVELS 2


// NOTE: settings for the following preprocessor macros
// are handled by cmake in top-level CMakeLists.txt file:
//
//      Macro                       Default
//  * ENABLE_PERFORMANCE_METERS      OFF
//  * ENABLE_CACHING                 ON
//  * ENABLE_CUDA_STREAMS            ON
//  * LOGGING                        OFF

// macro defining empty copy constructors and assignment operators
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);               \
    TypeName& operator=(const TypeName&)

// define integer constants
const int cFlatLimit = 11; // constant defining the first k where the calculation is done with terrain following coordinates
const int cDefaultKSize = 60; // default k size to compute / guess the difference between two k positions at compile time 
const int cNumBoundaryLines = 3; // maximal halo size of the data fields
const int cWarpSize = 32; // warp size used for CUDA back end
const int cCacheLineSize = 128; // number of bytes in a cache line

#ifdef __CUDACC__

// define CUDA method qualifiers
#define __ACC_CPU__ __host__ __device__
#define __ACC__ __device__
#define __CPU__ __host__
#define __KERNEL__ __global__

// define an assert
#define ACC_ASSERT(x) 

#else

// define CUDA method qualifiers
#define __ACC_CPU__
#define __ACC__
#define __CPU__ 
#define __KERNEL__ 

// define an assert
#include <cassert>
#define ACC_ASSERT(x) (assert(x))

#endif
