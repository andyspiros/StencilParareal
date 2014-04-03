set(
    GENERIC_SOURCES 
	"LaunchConfigurationUnittest.cpp"
)

set (OMP_SOURCES
	"KLoopComputationUnittest.cpp"
	"KRangeAndDomainUnittest.cpp"
	"StencilApplyUnittest.cpp" 
	"StencilParameterUnittest.cpp"
	"StencilBufferUnittest.cpp"
	"StencilCacheUnittest.cpp"
	"StencilSimpleLoopUnittest.cpp"
	"StencilComplexLoopUnittest.cpp"
	"StencilFunctionUnittest.cpp"
	"StencilStageUnittest.cpp"
	"IndirectIndexingUnittest.cpp"
	"StencilSwitchCaseUnittest.cpp"
#"StageVariableUnittest.cpp"
)

set( 
    CUDA_SOURCES 
	"KLoopComputationUnittest.cu"
	"KRangeAndDomainUnittest.cu"
	"StencilApplyUnittest.cu" 
	"StencilParameterUnittest.cu"
	"StencilBufferUnittest.cu"
	"StencilCacheUnittest.cu"
	"StencilSimpleLoopUnittest.cu"
	"StencilComplexLoopUnittest.cu"
	"StencilFunctionUnittest.cu"
	"StencilStageUnittest.cu"
	"IndirectIndexingUnittest.cu"
	"StencilSwitchCaseUnittest.cu"
#"StageVariableUnittest.cu"
)



# Add OpenMP test
add_executable(
    StencilFrameworkUnittest
    ${GENERIC_SOURCES}
    ${OMP_SOURCES}
)
target_link_libraries( StencilFrameworkUnittest gmock-gtest StellaSharedInfrastructure StellaStencilFramework )
add_test( StencilFrameworkUnittest StencilFrameworkUnittest )

# Add CUDA test if requested
if( CUDA_BACKEND )
    cuda_add_executable(
        StencilFrameworkUnittestCUDA
        ${GENERIC_SOURCES}
        ${CUDA_SOURCES}
        OPTIONS "${STELLA_CUDA_FLAGS}"
    )
    target_link_libraries( StencilFrameworkUnittestCUDA gmock-gtest StellaSharedInfrastructureCUDA StellaStencilFrameworkCUDA )
    add_test( StencilFrameworkUnittestCUDA StencilFrameworkUnittestCUDA )
endif( CUDA_BACKEND ) 

