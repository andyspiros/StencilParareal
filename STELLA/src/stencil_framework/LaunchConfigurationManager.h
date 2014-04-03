#pragma once

#include <vector>
#include "Definitions.h"

/**
* @class LaunchConfigurationManager
* Class managing the different launch configurations of a stencil (different boundary sizes)
*/
template<typename TLaunchConfiguration>
class LaunchConfigurationManager 
{
    DISALLOW_COPY_AND_ASSIGN(LaunchConfigurationManager);
public:
    __CPU__
    LaunchConfigurationManager() {}
    __CPU__
    ~LaunchConfigurationManager() {}

    /**
    * Initialize the launch configuration manager
    * @param calculationDomain calculation domain size
    * @param defaultBoundary boundary size of default configuration
    */
    __CPU__
    void Init(const IJKSize& calculationDomain, const IJBoundary& defaultBoundary)
    {
        // store the calculation domain size
        calculationDomain_ = calculationDomain;
        
        // setup the default configuration
        defaultConfiguration_.Init(calculationDomain_, defaultBoundary);

        // clear the configuration vector
        configurations_.clear();
    }   

    /**
    * Method searching the launch configuration for a certain boundary size
    * @param boundary boundary size of the requested configuration
    * @return pointer to launch configuration or NULL if now matching configuration was found
    */
    __CPU__
    const TLaunchConfiguration* FindConfiguration(const IJBoundary& boundary) const
    { 
        // by default return NULL
        const TLaunchConfiguration* result = NULL;

        // iterate over configurations and return the pointer if a matching configuration was found
        for(typename std::vector<TLaunchConfiguration>::const_iterator iter = configurations_.begin(); iter != configurations_.end(); ++iter)
        {
            if(iter->boundary() == boundary)
            {
                result = &(*iter);
                break;
            }
        }
        
        return result;
    }
    
    /**
    * Method adding a new launch configuration 
    * @param boundary boundary of the new configuration
    * @return added configuration
    */
    __CPU__
    const TLaunchConfiguration& AddConfiguration(const IJBoundary& boundary)
    {
        // make sure there is no config with the given boundary
        assert(FindConfiguration(boundary) == NULL);

        // setup a new configuration and add it to the configurations collection
        TLaunchConfiguration configuration;
        configuration.Init(calculationDomain_, boundary);
        configurations_.push_back(configuration);

        // return pointer to the added configuration
        return configurations_.back();
    }

    /**
    * @return the default launch configuration
    */
    __CPU__
    const TLaunchConfiguration& defaultConfiguration() const { return defaultConfiguration_; }

private:
    // define the calculation domain size 
    IJKSize calculationDomain_;

    // launch configurations
    TLaunchConfiguration defaultConfiguration_;
    std::vector<TLaunchConfiguration> configurations_;
};

