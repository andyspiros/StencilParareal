include_directories ( "../../libs/gmock-gtest" )
include_directories ( "../shared_definitions" ) 
include_directories ( "../shared_infrastructure" ) 

set( 
	OMP_SHARED_INFRASTRUCTURES_UNITTEST_SOURCES 
	"DataFieldOpenMPUnittest.cpp"
	"DataFieldOpenMPStorageUnittest.cpp"
)

set( 
	CUDA_SHARED_INFRASTRUCTURES_UNITTEST_KERNELS
	"DataFieldCUDAUnittest.cu"
	"DataFieldCUDAStorageUnittest.cu"
)

set( 
	SHARED_INFRASTRUCTURES_UNITTEST_SOURCES 
	"DataFieldUnittest.cpp"
	"DataFieldStorageUnittest.cpp"
	"JokerDataFieldUnittest.cpp"
	"UsageMeterUnittest.cpp"
	"SwapDataFieldUnittest.cpp"
	"TupleUnittest.cpp"
	"ArrayUnittest.cpp"
	"ApplyToAllUnittest.cpp"
)

add_executable(
    SharedInfrastructureUnittest 
    ${SHARED_INFRASTRUCTURES_UNITTEST_SOURCES}
    ${OMP_SHARED_INFRASTRUCTURES_UNITTEST_SOURCES} 
)
target_link_libraries( SharedInfrastructureUnittest gmock-gtest StellaSharedInfrastructure )

if( CUDA_BACKEND )
    cuda_add_executable(
        SharedInfrastructureUnittestCUDA
        ${SHARED_INFRASTRUCTURES_UNITTEST_SOURCES}
        ${CUDA_SHARED_INFRASTRUCTURES_UNITTEST_KERNELS}
        OPTIONS "${STELLA_CUDA_FLAGS}"
    )
target_link_libraries( SharedInfrastructureUnittestCUDA gmock-gtest StellaSharedInfrastructure )
endif( CUDA_BACKEND ) 


add_test( SharedInfrastructureUnittest SharedInfrastructureUnittest )

install( TARGETS SharedInfrastructureUnittest DESTINATION bin/ )

