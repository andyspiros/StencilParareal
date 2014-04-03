set (MAIN_HEADER 
    "StencilFramework.h" 
)

set(GENERIC_HEADERS
   "Temporaries.h"
   "KWindow.h"
   "Stencil.h"
   "StencilCompiler.h"
   "LaunchConfiguration.h"
   "LaunchConfigurationManager.h"
   "DummyContext.h"
   "Context.h"
   "ContextCache.h"
   "ContextCacheFunctors.h"
   "Offset.h"
   "DefineStages.h"
   "StencilCompilerFunctors.h"
   "StencilConfiguration.h"
   "StencilStage.h"
   "StencilStageEnvironment.h"
   "StencilSweepFunctor.h"
   "StencilSweepFunctorOperation.h"
   "StencilSweepFunctorOperationFunctor.h"
   "StencilSweepFunctorOperationGroup.h"
   "StencilSweepFunctorOperationGroupFunctor.h"
   "StencilSweepDescriptor.h"
   "StencilSweepGroupDescriptor.h"
   "Domains.h"
   "DataParameter.h"
   "DataParameterDynamic.h"
   "FunctionParameter.h"
   "WithWrapper.h"
   "StencilStageMacros.h"
   "Call.h"
   "StencilFunction.h"
   "StencilFunctionEnvironment.h"
   "StencilFunctionMacros.h"
   "DataParameterEvaluator.h"
   "FunctionParameterEvaluator.h"
   "HasDo.h"
   "HasParameter.h"
   "KPosition.h"
   "KRange.h"
   "IJRange.h"  
   "KBoundarySize.h"
   "KLoopCascade.h"
   "KLoopRange.h"
   "KLoopLeg.h"
   "KLoopLegFunctor.h"
   "KLoopBodyFunctors.h"
   "KLoopBody.h"
   "KLoopAdvance.h"
   "KLoopRestoreAndAdvance.h"
   "DefineTemporaries.h"
   "DefineLoops.h"
   "DefineSwitch.h"
   "DefineCase.h"
   "DefineCaches.h"
   "DefineSweep.h"
   "PackParameters.h"
   "ParameterWrapper.h"
   "YesNo.h"
   "Cache.h" 
   "CacheStorage.h"
   "CacheStorageFunctors.h"
   "ColumnBuffer.h"
   "ComparisonValue.h"
)

set(OMP_HEADERS
   "StencilSweepGroupFunctorOpenMP.h"
   "BlockLoopFunctorOpenMP.h"
   "StencilOpenMP.h"
   "StencilOpenMPFunctors.h"
   "StencilSweepFunctorOpenMP.h"
   "LaunchConfigurationOpenMP.h"
   "ContextOpenMP.h"
   "ContextOpenMPFunctors.h"
   "ColumnBufferOpenMP.h"
)

set(CUDA_HEADERS
   "StencilSweepGroupFunctorCUDA.h"
   "BlockLoopFunctorCUDA.h"
   "StencilCUDA.h"
   "StencilCUDAFunctors.h"
   "StencilSweepFunctorCUDA.h"
   "LaunchConfigurationCUDA.h"
   "ContextCUDA.h"
   "ContextCUDAFunctors.h"
   "ColumnBufferCUDA.h"
)

set(GENERIC_SOURCES
   "Dummy.cpp"
)

# Create library for OpenMP
add_library(StellaStencilFramework STATIC ${MAIN_HEADER} ${GENERIC_HEADERS} ${GENERIC_SOURCES} ${OMP_HEADERS} )
set_target_properties(StellaStencilFramework PROPERTIES
    COMPILE_FLAGS "${STELLA_OMP_FLAGS}"
)
install(TARGETS StellaStencilFramework DESTINATION "lib")
install (FILES ${GENERIC_HEADERS} ${OMP_HEADERS} DESTINATION "include/STELLA")
install (FILES ${MAIN_HEADER} DESTINATION "include/STELLA")

# Create library for CUDA if requested
if ("${CUDA_BACKEND}")
    cuda_add_library(StellaStencilFrameworkCUDA STATIC
        ${MAIN_HEADER} ${GENERIC_HEADERS} ${GENERIC_SOURCES} ${CUDA_HEADERS}
        OPTIONS "${STELLA_CUDA_FLAGS}"
    )
    install(TARGETS StellaStencilFrameworkCUDA DESTINATION "lib")
    install (FILES ${CUDA_HEADERS} DESTINATION "include/STELLA")
endif("${CUDA_BACKEND}")

