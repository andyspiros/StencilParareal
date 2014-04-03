#pragma once

#include <boost/static_assert.hpp>
#include "DataFieldStorageFormat.h"
#include "DataFieldCUDAStorageStrides.h"

/**
* @class DataFieldCUDAStorageIndex
* Storage class which abstracts the device storage of a DataField from the host perspective.
*/
template<typename TStorageStrides>
class DataFieldCUDAStorageIndex
{
public:        
    __ACC_CPU__
    DataFieldCUDAStorageIndex() { index_ = 0; }
    __ACC_CPU__
    ~DataFieldCUDAStorageIndex() {}
    
    __ACC_CPU__
    DataFieldCUDAStorageIndex(const DataFieldCUDAStorageIndex& other) { *this = other; }
    __ACC_CPU__
    DataFieldCUDAStorageIndex& operator= (const DataFieldCUDAStorageIndex& other) 
    { 
        index_ = other.index_;
        // by convention
        return *this;
    }

    /**
    * Init the index to the thread origin and advance it in all three dimensions
    * @param storageStrides stride class used for offset calculation
    * @param iBlockIndex i index of the block
    * @param jBlockIndex j index of the block
    * @param i offset in i direction
    * @param j offset in j direction
    * @param k offset in k direction
    */
    __ACC__
    void Init(const TStorageStrides& storageStrides, const int iBlockIndex, const int jBlockIndex, const int i, const int j, const int k) 
    {
        index_ = storageStrides.ComputeBlockStride(iBlockIndex, jBlockIndex, i, j, k);
    }

    /**
    * Advance the index in all three dimensions
    * @param storageStrides stride class used for offset calculation
    * @param i offset in i direction
    * @param j offset in j direction
    * @param k offset in k direction
    */
    __ACC__
    void Advance(const TStorageStrides& storageStrides, const int i, const int j, const int k)
    {
         index_ += storageStrides.ComputeStride(i, j, k);
    }
   
    /**
    * @return the current index
    */
    __ACC__
    int index() const { return index_; }
       
private:
    int index_;
};
  
