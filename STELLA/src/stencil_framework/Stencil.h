#pragma once

#include <cassert>
#include <string>
#include <sstream>
#include <iostream>
#include "SharedInfrastructure.h"

/**
* @class StencilInterface
* Stencil interface
*/
class StencilInterface
{
public:
    __CPU__
    StencilInterface() {};
    __CPU__
    virtual ~StencilInterface() {};   
    
    /**
    * Method applying the stencil using the default boundary
    */
    __CPU__
    virtual void ApplyImpl() = 0;

    /**
    * Method applying the stencil using a variable boundary
    * @param boundary apply the stencil with a boundary around the calculation domain
    */
    __CPU__
    virtual void ApplyImpl(const IJBoundary& boundary) = 0;

    /**
    * Method printing the bounding box information
    */
    __CPU__
    virtual std::string PrintBoundingBoxInfoImpl() = 0;
};

/**
* @class Stencil
* Class providing a template parameter free interface for the stencil base functionality.
* The class serves as a base class for the stencil apply wrappers of the different back ends. 
* These apply wrappers provide the initialization functionality necessary to create a stencil
*/
class Stencil
{
    DISALLOW_COPY_AND_ASSIGN(Stencil);
public:
    __CPU__
    Stencil() 
    {
        pStencilWrapper_ = NULL;
    }
    __CPU__
    ~Stencil() 
    {
        if(pStencilWrapper_)
        {
            delete pStencilWrapper_;
        }
    }

    /**
    * Method used to setup the stencil
    * @param name name of the stencil
    * @param pStencilWrapper pointer to the stencil wrapper
    */
    __CPU__
    void Init(std::string name, StencilInterface* pStencilWrapper) 
    {
        // delete old stencil wrapper
        if(pStencilWrapper_)
        {
            delete pStencilWrapper_;
        }

        // set parameter fields
        pStencilWrapper_ = pStencilWrapper;
        name_ = name;

        // meters
        performanceMeter_.Init(name_ + ".PerformanceMeter");
        usageMeter_.Init(name_ + ".UsageMeter", "#");
    }

    /**
    * The method applies the stencil to the parameter data fields using the default boundary
    */
    __CPU__
    void Apply()
    {
#ifdef ENABLE_LOGGING
        std::cout << "Stencil " << name_ << " Apply() started..." << std::endl;
#endif

        assert(pStencilWrapper_);

        // call the virtual method of the interface
        // -> the call is forwarded to the actual stencil wrapper
        performanceMeter_.Start();
        pStencilWrapper_->ApplyImpl();
        performanceMeter_.Pause();
        usageMeter_.Add(1.0);

#ifdef ENABLE_LOGGING
        std::cout << "Bounding Box Info:" << std::endl << pStencilWrapper_->PrintBoundingBoxInfoImpl();
        std::cout << "Stencil " << name_ << " Apply() end" << std::endl;
#endif
    }

    /**
    * The method applies the stencil to the parameter data fields
    * @param boundary apply the stencil with a boundary around the calculation domain
    */
    __CPU__
    void Apply(const IJBoundary& boundary)
    {
#ifdef ENABLE_LOGGING
        std::cout << "Stencil " << name_ << " Apply(const IJBoundary&) started..." << std::endl;
#endif

        assert(pStencilWrapper_);

        // call the virtual method of the interface
        // -> the call is forwarded to the actual stencil wrapper
        performanceMeter_.Start();
        pStencilWrapper_->ApplyImpl(boundary);
        performanceMeter_.Pause();
        usageMeter_.Add(1.0);

#ifdef ENABLE_LOGGING
        std::cout << "Bounding Box Info:" << std::endl << pStencilWrapper_->PrintBoundingBoxInfoImpl();
        std::cout << "Stencil " << name_ << " Apply(const IJBoundary&) end" << std::endl;
#endif
    }

    /**
    * Method printing the bounding box for every stencil parameter
    */
    __CPU__
    void PrintBoundingBoxInfo() const 
    {
        assert(pStencilWrapper_);
        std::cout << pStencilWrapper_->PrintBoundingBoxInfoImpl() << std::endl;
    }

    /**
    * Method reseting all stencil meters
    */
    __CPU__
    void ResetMeters() 
    {
        performanceMeter_.Reset();
        usageMeter_.Reset();
    }

    /**
    * Method converting all stencil meters to a string
    */
    __CPU__
    std::string MetersToString() const
    {
        std::ostringstream out; 
        out << name_ 
            << "\ttime [s]: " << performanceMeter_.totalTime()
            << "\tusage [#]: " << usageMeter_.counter();
        return out.str();
    }

    /**
    * @return performance meter
    */
    __CPU__
    const PerformanceMeter& performanceMeter() const
    {
        return performanceMeter_;
    }

    /**
    * @return the name of the stencil
    */
    __CPU__
    std::string name() const
    {
        return name_;
    }

private:
    std::string name_;
   
    // meters
    PerformanceMeter performanceMeter_;
    UsageMeter usageMeter_;
    StencilInterface* pStencilWrapper_;
};
  
/**
* @struct create_stencil
* Meta function creating a stencil
*/
template<
    typename TStencilConfiguration,
    typename TParameterIJBoundary,
    typename TParameterTupleElements,
    typename TTemporaryFields,
    typename TStencilSweepGroupDescriptors>
struct create_stencil;
