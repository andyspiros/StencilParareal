#pragma once

#include <vector>
#include "LaunchConfiguration.h"

/**
* @class LaunchConfigurationOpenMP
* OpenMP specific implementation of the launch configuration base class
*/
template<
    typename TBlockSize,
    typename TParameterIJBoundary>
class LaunchConfigurationOpenMP : public LaunchConfiguration<TBlockSize, TParameterIJBoundary, LaunchConfigurationOpenMP<TBlockSize, TParameterIJBoundary> > // CRTP
{
public:    
    using LaunchConfiguration<TBlockSize, TParameterIJBoundary, LaunchConfigurationOpenMP>::Assign;
   
    LaunchConfigurationOpenMP() {}
    ~LaunchConfigurationOpenMP() {}

    LaunchConfigurationOpenMP(const LaunchConfigurationOpenMP& other) { *this = other; }
    LaunchConfigurationOpenMP& operator= (const LaunchConfigurationOpenMP& other)
    {
        Assign(other);
        // by convention
        return *this;
    }    

    /**
    * Implement OpenMP launch configuration initialization
    */
    void InitImpl(const std::vector<BlockConfiguration>& blockConfigurations) {}

    // OpenMP specific interface
};



