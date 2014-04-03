set (MAIN_HEADER
   "SharedInfrastructure.h"
)

set(GENERIC_HEADERS
   "BlockSize.h"
   "IJKSize.h"
   "IJKIndex.h"
   "IJBoundary.h"
   "IJKBoundary.h"
   "KBoundary.h"
   "DataFieldStorageOrder.h"
   "DataFieldAlignment.h"
   "DataField.h"
   "DataFieldStorage.h"
   "DataFieldIJBoundary.h"
   "DataFieldStorageFormat.h"
   "DataFieldStorageStrides.h"
   "DataFieldOpenMPStorage.h"
   "DataFieldOpenMPStorageIterator.h"
   "JokerDataField.h"
   "SwapDataField.h"
   "ApplyToAll.h"
   "ApplyToAllDefinition.h"
   "TupleElements.h"
   "Tuple.h"
   "TupleDefinition.h"
   "TupleAlgorithms.h"
   "TupleAlgorithmsDefinition.h"
   "Array.h"
   "ArrayDefinition.h"
   "ArrayAlgorithms.h"
   "ArrayAlgorithmsDefinition.h"
   "TimerDummy.h"
   "Timer.h"
   "MathFunctions.h"
   "ParameterType.h"
   "PositionChecker.h"
   "ParameterTraits.h"
   "ScalarStorage.h"
   "DummyStorage.h"
   "BlockStorage.h"
   "UsageMeter.h"
   "ExternalStorage.h"
)

set(OMP_HEADERS
   "DataFieldOpenMP.h"
   "TimerOpenMP.h"
)

set(CUDA_HEADERS
   "DataFieldCUDA.h"
   "DataFieldCUDAStorage.h"
   "DataFieldCUDAStoragePointer.h"
   "DataFieldCUDAStorageIndex.h"
   "DataFieldCUDAStorageStrides.h"
   "TimerCUDA.h"
   "CUDAStreams.h"
)

set(GENERIC_SOURCES
   "PositionChecker.cpp"
)

set(CUDA_SOURCES
   "CUDAStreams.cpp"
)

# Create library for OpenMP
add_library(StellaSharedInfrastructure STATIC ${MAIN_HEADER} ${GENERIC_HEADERS} ${GENERIC_SOURCES} ${OMP_HEADERS} )
set_target_properties(StellaSharedInfrastructure PROPERTIES
    COMPILE_FLAGS "${STELLA_OMP_FLAGS}"
)
install(TARGETS StellaSharedInfrastructure DESTINATION "lib")
install (FILES ${GENERIC_HEADERS} ${OMP_HEADERS} DESTINATION "include/STELLA")
install (FILES ${MAIN_HEADER} DESTINATION "include/STELLA")

# Create library for CUDA if requested
if ( CUDA_BACKEND )
    cuda_add_library(StellaSharedInfrastructureCUDA STATIC
        ${GENERIC_HEADERS} ${MAIN_HEADER} ${GENERIC_SOURCES} ${CUDA_HEADERS} ${CUDA_SOURCES}
        OPTIONS "${STELLA_CUDA_FLAGS}"
    )
    install(TARGETS StellaSharedInfrastructureCUDA DESTINATION "lib")
    install (FILES ${CUDA_HEADERS} DESTINATION "include/STELLA")
endif( CUDA_BACKEND )

