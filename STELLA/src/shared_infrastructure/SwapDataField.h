#pragma once

#include <cassert>
#include <boost/static_assert.hpp>
#include "DataField.h"

/**
* @class SwapDataField
* Container for two data fields used to implement double buffering
*/
template<typename TDataField>
class SwapDataField
{
    BOOST_STATIC_ASSERT(is_data_field<TDataField>::value);
public:
    typedef typename TDataField::StorageType StorageType;

    __CPU__
    SwapDataField() { pMemorizedSwapPositionStorage_ = NULL; }
    __CPU__
    ~SwapDataField() {}

    __CPU__
    SwapDataField(const SwapDataField& other) { *this = other; }
    __CPU__
    SwapDataField& operator= (const SwapDataField& other)
    {
        in_ = other.in_;
        out_ = other.out_;
        pMemorizedSwapPositionStorage_ = other.pMemorizedSwapPositionStorage_;
        // by convention
        return *this;
    }

    /**
    * Init the swap data field and allocate the in and out data fields
    * @param name name of the data fields
    * @param calculationDomain size of the fields without any halo or boundary
    * @param kBoundary boundary around the calculation domain in k direction
    */
    __CPU__
    void Init(std::string name, const IJKSize& calculationDomain, const KBoundary& kBoundary)
    {
        name_ = name;

        // if necessary init the member data fields
        in_.Init(name + "_in", calculationDomain, kBoundary);
        out_.Init(name + "_out", calculationDomain, kBoundary);

        pMemorizedSwapPositionStorage_ = NULL;
    }

    /**
    * Init the swap data field without allocating the in and out data fields
    * @param name name of the data fields
    */
    __CPU__
    void InitWithoutDataFieldAllocation(std::string name)
    {
        name_ = name;
        pMemorizedSwapPositionStorage_ = NULL;
    }

    /**
    * Swaps input and output data field
    */
    __CPU__
    void Swap() { in_.SwapWith(out_); }

    /**
    * Method memorizing the current in data field storage
    */
    __CPU__
    void MemorizeSwapPosition() { pMemorizedSwapPositionStorage_ = &in_.storage(); }

    /**
    * @return pointer to the in data field storage memorized, or NULL if no position was memorized
    */
    __CPU__
    const StorageType* pMemorizedSwapPositionStorage() const { return pMemorizedSwapPositionStorage_; }

    /**
    * @return actual input data field
    */
    __CPU__
    TDataField& in() { return in_; } 
    
    /**
    * @return actual output data field
    */
    __CPU__
    TDataField& out() { return out_; } 

    /**
    * @return name 
    */
    __CPU__
    std::string name() { return name_; }

private:
    std::string name_;
    TDataField in_;
    TDataField out_;
    const StorageType* pMemorizedSwapPositionStorage_;
};


  
