#pragma once

#include <cassert>
#include <vector>
#include <string>
#include <boost/static_assert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/fold.hpp>
#include "SharedInfrastructure.h"
#include "StencilSweepGroupFunctorCUDA.h"
#include "CUDADefinitions.h"
#include "CUDAStreams.h"
#include "ContextCUDA.h"
#include "ContextCache.h"
#include "ParameterWrapper.h"
#include "ColumnBufferCUDA.h"
#include "LaunchConfigurationCUDA.h"
#include "LaunchConfigurationManager.h"
#include "Stencil.h"
#include "StencilCUDAFunctors.h"

/**
* Kernel applying the stencil
* @param pSharedDataDevice device pointer to the parameter tuple
* @param pBlockConfigurationsDevice device pointer to the block configuration array
*/
template<
    typename TContextDescriptor,
    typename TStencilSweepGroupDescriptors>
__KERNEL__ 
void stencil_apply(float* pSharedDataDevice, float* pBlockConfigurationsDevice)
{
    // compute the maximum ij cache tuple type used in order to allocate the shared memory
    typedef typename stencil_sweep_group_descriptor_maximum_ij_cache_tuple<
        TContextDescriptor, 
        TStencilSweepGroupDescriptors
    >::type MaximumIJCacheTuple;

    // define the context shared data type (stored in shared memory)
    typedef ContextCUDASharedData<
        TContextDescriptor, 
        acc::Tuple, 
        BlockConfiguration, 
        sizeof(MaximumIJCacheTuple)
    > SharedData;
 
    // define the context shared data type without block configuration and cache tuple
    // (note this type is used in order to compute the block configuration offset
    // as the c++ offset off function only supports POD types)
    typedef ContextCUDASharedData<
        TContextDescriptor, 
        acc::Tuple, 
        boost::mpl::void_, 
        0
    > SharedDataWithoutBlockConfigAndIJCaches;

    // make sure the type sizes and offsets are a multiple of float
    BOOST_STATIC_ASSERT(sizeof(SharedData) % sizeof(float) == 0);
    BOOST_STATIC_ASSERT(sizeof(BlockConfiguration) % sizeof(float) == 0);  
    BOOST_STATIC_ASSERT(sizeof(SharedDataWithoutBlockConfigAndIJCaches) % sizeof(float) == 0);
    
    // compute the block configuration offset
    typedef boost::mpl::integral_c<int, 
        sizeof(SharedDataWithoutBlockConfigAndIJCaches) - sizeof(float)
    > OffsetOfBlockConfiguration;

    // define a shared memory variable holding the context data 
    // (allocate double buffer in order to make sure 8 byte alignment)
    __shared__ double pSharedMemory[sizeof(SharedData) / sizeof(double) + 1];
    float* pSharedData = reinterpret_cast<float*>(pSharedMemory);
   
    // setup the copyIndex which different for every thread of the block
    int copyIndex = threadIdx.x + threadIdx.y * cWarpSize;
    
    // compute the blockConfigurationOffset inside the pSharedDataArray array
    const int blockConfigurationOffset = OffsetOfBlockConfiguration::value / sizeof(float);

    // copy the shared data in parallel to the shared memory
    while(copyIndex < blockConfigurationOffset)
    {
        pSharedData[copyIndex] = pSharedDataDevice[copyIndex];
        copyIndex += warps_per_stencil<typename TContextDescriptor::BlockSize, TStencilSweepGroupDescriptors>::value * cWarpSize;
    }

    // setup the blockIndex used to extract the matching block configuration from the block configuration array
    const int blockIndex = blockIdx.x * (sizeof(BlockConfiguration) / sizeof(float));

    // copy the block configuration in parallel to the shared memory
    // (note that the copyIndex is further increment after the first while loop starting at blockConfigurationOffset)
    while(copyIndex < (OffsetOfBlockConfiguration::value + sizeof(BlockConfiguration)) / sizeof(float))
    {
        // subtract the blockConfigurationOffset from the copyIndex in order to get a copy index starting at zero
        pSharedData[copyIndex] = pBlockConfigurationsDevice[copyIndex - blockConfigurationOffset + blockIndex]; 
        copyIndex += warps_per_stencil<typename TContextDescriptor::BlockSize, TStencilSweepGroupDescriptors>::value * cWarpSize;
    }
    __syncthreads();    
      
    // apply all stencil stages to the block using the apply all functionality
    acc::apply_to_all<
        StencilSweepGroupFunctorCUDA,
        TStencilSweepGroupDescriptors, 
        SharedData*
    >(reinterpret_cast<SharedData*>(pSharedData));
}

/**
* @class StencilCUDA
* CUDA implementation of the stencil interface
*/
template<
    typename TStencilConfiguration,
    typename TParameterIJBoundary,
    typename TParameterTupleElements,
    typename TTemporaryFields,
    typename TStencilSweepGroupDescriptors>
class StencilCUDA : public StencilInterface // implements the stencil interface
{
    DISALLOW_COPY_AND_ASSIGN(StencilCUDA);
public:
    // define the block size
    typedef typename TStencilConfiguration::BlockSize BlockSize;
    typedef typename TStencilConfiguration::ValueType ValueType;
    typedef LaunchConfigurationCUDA<BlockSize, TParameterIJBoundary> LaunchConfiguration;

    // define the parameter tuple information
    typedef typename TParameterTupleElements::ElementIndexes ParameterElementIndexes;
    typedef typename TParameterTupleElements::ElementTypes ParameterElementTypes;
    typedef TParameterTupleElements ParameterTupleElements;
    typedef Tuple<ParameterTupleElements> ParameterTupleType;

    // calculate the indexes of all data field parameters
    typedef typename data_field_parameters<ParameterElementTypes>::type DataFieldParameterIndexes;

    // define the temporary field tuple information
    typedef typename boost::mpl::transform<
        TTemporaryFields,
        temporary_field_index<boost::mpl::_>
    >::type TemporaryFieldElementIndexes;
    typedef typename boost::mpl::transform<
        TTemporaryFields,
        create_temporary_field<boost::mpl::_, BlockSize, TParameterIJBoundary, TStencilSweepGroupDescriptors>
    >::type TemporaryFieldElementTypes;
    typedef TupleElements<TemporaryFieldElementIndexes, TemporaryFieldElementTypes> TemporaryFieldTupleElements;
    typedef Tuple<TemporaryFieldTupleElements> TemporaryFieldTupleType;

    // define the apply tuple
    typedef typename merge_tuple_elements<
        ParameterTupleElements,
        TemporaryFieldTupleElements,
        parameter_to_apply<boost::mpl::_, BlockSize>,
        temporary_field_to_apply<boost::mpl::_>
    >::type ApplyTupleElements;

    // define the context descriptor type
    typedef ContextDescriptor<
        TStencilConfiguration,
        ApplyTupleElements
    > ContextDescriptorType;

    // define the shared data type 
    // (set block configuration to float in order to assure spare data is available for copy in kernel)
    typedef ContextCUDASharedData<ContextDescriptorType, Tuple, boost::mpl::void_, 0> SharedData;

    __CPU__
    StencilCUDA() 
    {
        // default initializations
        threadsPerBlock_.x = 0;
        threadsPerBlock_.y = 0;
        threadsPerBlock_.z = 0;
        numberOfBlocks_.x = 0;
        numberOfBlocks_.y = 0;
        numberOfBlocks_.z = 0;

        // allocate page locked memory for the host copy of the shared data
        cudaMallocHost(reinterpret_cast<void**>(&pSharedDataHost_), sizeof(SharedData)); 

        // allocate 2 device copies of the shared data
        // (necessary to copy the data to the device while a kernel is working using the other copy)
        sharedDataDevice_.resize(2);
        cudaMalloc(reinterpret_cast<void**>(&sharedDataDevice_[0]), sizeof(SharedData)); 
        cudaMalloc(reinterpret_cast<void**>(&sharedDataDevice_[1]), sizeof(SharedData)); 
        
        // initialize the shared data index (alternating between 0 and 1)
        // it points to the currently active device copy of the shared data
        sharedDataDeviceIndex_ = 0;
     
        // initialize host copy with default values
        SharedData defaultSharedData;
        *pSharedDataHost_ = defaultSharedData;
        
        // initialize events
        cudaEventCreate(&copyDoneEvent_);
        cudaEventCreate(&oldSharedDataDeviceFreeEvent_);
        
        // increment stream user counter
        CUDAStreams::IncrementReferenceCounter();

        // verify the CUDA calls didn't fail
        assertNoCUDAError("StencilWrapperCUDA");
    }
    __CPU__
    virtual ~StencilCUDA() 
    {
        // free shared data memory
        cudaFreeHost(pSharedDataHost_);
        assert(static_cast<int>(sharedDataDevice_.size()) == 2);
        cudaFree(sharedDataDevice_[0]);
        cudaFree(sharedDataDevice_[1]);

        // free CUDA events
        cudaEventDestroy(copyDoneEvent_);
        cudaEventDestroy(oldSharedDataDeviceFreeEvent_);

        // decrement stream user counter
        CUDAStreams::DecrementReferenceCounter();

        // verify the CUDA calls didn't fail
        assertNoCUDAError("StencilWrapperCUDA");
    }

    /**
    * Initialize the wrapper
    * @param calculationDomain size of the calculation domain
    * @param parameterTuple tuple holding the stencil parameters
    */
    __CPU__
    void Init(const IJKSize& calculationDomain, Tuple<ParameterTupleElements>& parameterTuple)
    {
        // setup the launch configuration manager
        IJBoundary defaultBoundary;
        defaultBoundary.Init(0, 0, 0, 0);
        launchConfigurationManager_.Init(calculationDomain, defaultBoundary);

        // initialize the parameter tuple
        parameterTuple_ = parameterTuple;

        // initialize the temporary field tuple
        modify_tuple<
            TemporaryFieldInitFunctor,
            TemporaryFieldElementIndexes,
            const IJKSize
        >(temporaryFieldTuple_, calculationDomain);

        // initialize the apply tuple structure
        modify_2_tuples<
            ParameterToApplyFunctor<BlockSize>,
            ParameterElementIndexes
        >(parameterTuple_, pSharedDataHost_->dataTuple);

        modify_2_tuples<
            TemporaryFieldToApplyFunctor,
            TemporaryFieldElementIndexes
        >(temporaryFieldTuple_, pSharedDataHost_->dataTuple);

        // initialize the storage strides
        typedef Tuple<typename storage_strides_elements<ContextDescriptorType>::type> StorageStridesTupleType;
        modify_tuple<
            ParameterStridesInitFunctor<BlockSize, StorageStridesTupleType>,
            DataFieldParameterIndexes,
            StorageStridesTupleType
        >(parameterTuple_, pSharedDataHost_->storageStridesTuple);

        modify_tuple<
            TemporaryFieldStridesInitFunctor<StorageStridesTupleType>,
            TemporaryFieldElementIndexes,
            StorageStridesTupleType
        >(temporaryFieldTuple_, pSharedDataHost_->storageStridesTuple);

        // setup the k size member
        pSharedDataHost_->kSize = calculationDomain.kSize();

        // copy shared data to the currently active device copy
        cudaMemcpy(sharedDataDevice_[sharedDataDeviceIndex_], pSharedDataHost_, sizeof(SharedData), cudaMemcpyHostToDevice);
        
        // init the size variables
        threadsPerBlock_.x = cWarpSize;
        threadsPerBlock_.y = warps_per_stencil<BlockSize, TStencilSweepGroupDescriptors>::value;
        threadsPerBlock_.z = 1;
        numberOfBlocks_.x = 0; // set at kernel launch time
        numberOfBlocks_.y = 1;
        numberOfBlocks_.z = 1;

        // determine the ij cache tuple size
        typedef typename stencil_sweep_group_descriptor_maximum_ij_cache_tuple<
            ContextDescriptorType, 
            TStencilSweepGroupDescriptors
        >::type MaximumIJCacheTuple;

        // favor l1 cache or shared memory depending on ij cache tuple size
        cudaFuncSetCacheConfig(
            stencil_apply<ContextDescriptorType, TStencilSweepGroupDescriptors>, 
            (sizeof(MaximumIJCacheTuple) > 4096 ? cudaFuncCachePreferShared : cudaFuncCachePreferL1)
        );

        // set the memory bank configuration
        cudaDeviceSetSharedMemConfig(
            (sizeof(ValueType) == 8 ? cudaSharedMemBankSizeEightByte : cudaSharedMemBankSizeFourByte)
        );
        
        // verify the CUDA calls didn't fail
        assertNoCUDAError("StencilWrapperCUDA");
    }

    // stencil interface implementation

    /**
    * Method applying the stencil using the default boundary
    */
    __CPU__
    virtual void ApplyImpl()
    {
        // run the stencil using the default configuration
        applyImpl(launchConfigurationManager_.defaultConfiguration());
    }

    /**
    * Method applying the stencil using a variable boundary
    * @param boundary apply the stencil with a boundary around the calculation domain
    */
    __CPU__
    virtual void ApplyImpl(const IJBoundary& boundary)
    {
        // find a matching launch configuration
        const LaunchConfiguration* pLaunchConfiguration = launchConfigurationManager_.FindConfiguration(boundary);

        // if a config was found apply it 
        if(pLaunchConfiguration) 
        { 
            applyImpl(*pLaunchConfiguration);
            
        }
        // if no config was found add a new one and apply it
        else
        {
            applyImpl(launchConfigurationManager_.AddConfiguration(boundary));
        }
    }

    /**
    * Method printing the bounding box for all stencil parameters
    */
    __CPU__
    virtual std::string PrintBoundingBoxInfoImpl()
    {
        return "PrintBoundingBoxInfo not supported";
    }

private:
    // method implementing apply given a launch configuration
    __CPU__
    void applyImpl(const LaunchConfiguration& launchConfiguration)
    {
        // update the joker field pointers
        bool applyDataModified = false;
        modify_2_tuples<
            UpdateApplyFunctor<BlockSize>,
            ParameterElementIndexes,
            bool&
        >(parameterTuple_, pSharedDataHost_->dataTuple, applyDataModified);

        // if necessary copy the apply data to the device
        if(applyDataModified)
        {
            // set the index to the next available shared data device copy
            assert(static_cast<int>(sharedDataDevice_.size()) == 2);
            sharedDataDeviceIndex_ = (sharedDataDeviceIndex_ + 1) % 2;

            // pause copy until the shared data copy is not used by a kernel still in the queue anymore
            // (note the host thread is not synchronized!)
            cudaStreamWaitEvent(CUDAStreams::copyStream(), oldSharedDataDeviceFreeEvent_, 0);

            // start asynchronous copy of shared data to the device
            cudaMemcpyAsync(
                sharedDataDevice_[sharedDataDeviceIndex_], 
                pSharedDataHost_, 
                sizeof(SharedData), 
                cudaMemcpyHostToDevice, 
                CUDAStreams::copyStream()
            );
            
            // once the shared data is copied to the device raise an event
            cudaEventRecord(copyDoneEvent_, CUDAStreams::copyStream());

            // once a potentially running kernel using the old shared data device copy is done raise an event
            cudaEventRecord(oldSharedDataDeviceFreeEvent_, CUDAStreams::kernelStream());

            // pause the next kernel call until the shared data is copied to the device
            // (note the host thread is not synchronized!)
            cudaStreamWaitEvent(CUDAStreams::kernelStream(), copyDoneEvent_, 0);
        }

        // set the number of blocks
        numberOfBlocks_.x = static_cast<unsigned int>(launchConfiguration.blockConfigurations().size());
        
        // launch the stencil apply kernel
        stencil_apply<
            ContextDescriptorType, 
            TStencilSweepGroupDescriptors
        > <<<numberOfBlocks_, threadsPerBlock_, 0, CUDAStreams::kernelStream()>>> (
            sharedDataDevice_[sharedDataDeviceIndex_], 
            launchConfiguration.pBlockConfigurationsDevice()
        ); 

        // if a new shared data copy was moved to the device, stop the host thread until the copy is done
        if(applyDataModified)
        { 
            cudaStreamSynchronize(CUDAStreams::copyStream()); 
        }

        // verify the CUDA calls didn't fail
        assertNoCUDAError("StencilWrapperCUDA");
    }

    // launch configuration manager
    LaunchConfigurationManager<LaunchConfiguration> launchConfigurationManager_;

    // kernel launch configuration
    dim3 threadsPerBlock_;
    dim3 numberOfBlocks_;
    
    // CUDA events
    cudaEvent_t copyDoneEvent_;
    cudaEvent_t oldSharedDataDeviceFreeEvent_;

    // define the tuples
    ParameterTupleType parameterTuple_;
    TemporaryFieldTupleType temporaryFieldTuple_;

    // context host copy
    SharedData* pSharedDataHost_;
    
    // 2 context device copies, accessed alternatingly
    std::vector<float*> sharedDataDevice_;
    int sharedDataDeviceIndex_; 
};
  
// create a CUDA stencil
template<
    typename TStencilConfiguration, 
    typename TParameterIJBoundary,
    typename TParameterTupleElements,
    typename TTemporaryFields,
    typename TStencilSweepGroupDescriptors>
struct create_stencil
{
    typedef StencilCUDA<TStencilConfiguration, TParameterIJBoundary, TParameterTupleElements, TTemporaryFields, TStencilSweepGroupDescriptors> type;
};
