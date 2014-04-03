set(GENERIC_HEADERS
    "Enums.h" 
    "Definitions.h"
)

set(CUDA_HEADERS
    "${STELLA_CUDA_HEADERS}"
    "CUDADefinitions.h"
)

# Install headers
install (FILES ${GENERIC_HEADERS} DESTINATION "include/STELLA")

if( CUDA_BACKEND )
    install (FILES ${CUDA_HEADERS} DESTINATION "include/STELLA")
endif( CUDA_BACKEND )
