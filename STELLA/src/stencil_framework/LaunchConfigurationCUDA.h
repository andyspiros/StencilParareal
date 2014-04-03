#pragma once

#include <cassert>
#include <vector>
#include <boost/static_assert.hpp>
#include "CUDADefinitions.h"
#include "LaunchConfiguration.h"

/**
* @class LaunchConfigurationCUDA
* CUDA specific implementation of the launch configuration base class
*/
template<
    typename TBlockSize,
    typename TParameterIJBoundary>
class LaunchConfigurationCUDA : public LaunchConfiguration<TBlockSize, TParameterIJBoundary, LaunchConfigurationCUDA<TBlockSize, TParameterIJBoundary> > // CRTP
{
public:    
    using LaunchConfiguration<TBlockSize, TParameterIJBoundary, LaunchConfigurationCUDA>::Assign;
   
    __CPU__
    LaunchConfigurationCUDA() 
    { 
        pBlockConfigurationsDevice_ = NULL; 
    }
    __CPU__
    ~LaunchConfigurationCUDA() 
    {
        clearBlockConfigurations();
    }

    __CPU__
    LaunchConfigurationCUDA(const LaunchConfigurationCUDA& other) 
    { 
        // make sure pBlockConfigurationsDevice is initialized to NULL 
        // otherwise the assignment operator will try the free the memory
        pBlockConfigurationsDevice_ = NULL; 
        *this = other; 
    }
    __CPU__
    LaunchConfigurationCUDA& operator= (const LaunchConfigurationCUDA& other)
    {
        Assign(other);
        // copy data from other if necessary
        if(other.pBlockConfigurationsDevice_)
        {
            initBlockConfigurations(other.blockConfigurations());
        }
        // otherwise clear device memory
        else
        {
            clearBlockConfigurations();
        }
        // by convention
        return *this;
    }    

    /**
    * Implement CUDA launch configuration initialization
    */
    __CPU__
    void InitImpl(const std::vector<BlockConfiguration>& blockConfigurations) 
    {   
        // init device copy of block configurations
        initBlockConfigurations(blockConfigurations);
    }

    // CUDA specific interface

    /**
    * @return float pointer to the block configurations stored in device memory
    */
    __CPU__
    float* pBlockConfigurationsDevice() const 
    { 
        assert(pBlockConfigurationsDevice_); 
        return pBlockConfigurationsDevice_; 
    }

private:
    // method clearing pBlockConfigurationsDevice_ if necessary
    __CPU__
    void clearBlockConfigurations()
    {
        if(pBlockConfigurationsDevice_)
        {
            cudaFree(pBlockConfigurationsDevice_);
            pBlockConfigurationsDevice_ = NULL;
        }

        // verify the CUDA calls didn't fail
        assertNoCUDAError("LaunchConfigurationCUDA");
    }
    // method allocating device memory for block configurations
    __CPU__
    void initBlockConfigurations(const std::vector<BlockConfiguration>& blockConfigurations)
    {
        // make sure block configuration size is a multiple of float size
        BOOST_STATIC_ASSERT(sizeof(BlockConfiguration) % sizeof(float) == 0);

        // clear block configuration
        clearBlockConfigurations();

        // allocate device memory
        cudaMalloc(reinterpret_cast<void**>(&pBlockConfigurationsDevice_), sizeof(BlockConfiguration) * blockConfigurations.size());

        // copy the data to the device
        cudaMemcpy(pBlockConfigurationsDevice_, &blockConfigurations[0], sizeof(BlockConfiguration) * blockConfigurations.size(), cudaMemcpyHostToDevice);

        // verify the CUDA calls didn't fail
        assertNoCUDAError("LaunchConfigurationCUDA");
    }
    
    float* pBlockConfigurationsDevice_;
};



