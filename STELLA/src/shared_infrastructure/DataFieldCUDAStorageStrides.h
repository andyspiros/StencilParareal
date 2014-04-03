#pragma once

#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>
#include "Definitions.h"
#include "Enums.h"
#include "DataFieldStorageFormat.h"
#include "DataFieldStorageStrides.h"

/**
* @struct ScalarStorageStrides
* Used by DataFieldCUDAStorageStrides class to implement the strides for scalar element types
*/
template<
    typename TStorageFormat,
    typename TBlockSize>
struct ScalarStorageStrides
{
    // define the strides type
    typedef DataFieldStorageStrides<typename TStorageFormat::StorageOrder> StorageStridesType;

    /**
    * @return block internal stride
    */
    __ACC__
    static int ComputeStride(const StorageStridesType& strides, const int i, const int j, const int k)
    {
        return strides.ComputeStride(i, j, k);
    }

    /**
    * @return block stride 
    */
    __ACC__
    static int ComputeBlockStride(const StorageStridesType& strides, const int iBlockIndex, const int jBlockIndex, const int i, const int j, const int k)
    {
        // compute stride if there is no block private boundary
        return strides.ComputeStride(iBlockIndex * TBlockSize::ISize::value + i, jBlockIndex * TBlockSize::JSize::value + j, k);
    }
};

/**
* @struct BlockStorageStrides
* Used by DataFieldCUDAStorageStrides class to implement the strides for block storage element types
*/
template<
    typename TStorageFormat,
    typename TBlockStorageFormat>
struct BlockStorageStrides
{
    // define the strides type
    typedef DataFieldStorageStrides<typename TStorageFormat::StorageOrder> StorageStridesType;

    /**
    * @return block internal stride using data field strides for k dimension and block storage strides in ij dimension
    */
    __ACC__
    static int ComputeStride(const StorageStridesType& strides, const int i, const int j, const int k)
    {
        return strides.ComputeStride(0, 0, k) * TBlockStorageFormat::Size::value + TBlockStorageFormat::ComputeStride(i, j);
    }

    /**
    * @return block stride using data field strides for the block offset and block storage strides for the block internal offset
    */
    __ACC__
    static int ComputeBlockStride(const StorageStridesType& strides, const int iBlockIndex, const int jBlockIndex, const int i, const int j, const int k)
    {
        return strides.ComputeStride(iBlockIndex, jBlockIndex, k) * TBlockStorageFormat::Size::value + TBlockStorageFormat::ComputeStride(i, j);
    }
};

/**
* @class DataFieldCUDAStorageStrides
* CUDA back end storage stride class. It warps the normal strides class. There are two storage strides specializations available.
* One specialization works with scalar data field elements the other one supports block storage element types.
*/
template<typename TStorageStrides>
class DataFieldCUDAStorageStrides 
{
public:
    typedef typename TStorageStrides::StorageStridesType StorageStridesType;

    __ACC_CPU__
    DataFieldCUDAStorageStrides() {}
    __ACC_CPU__
    ~DataFieldCUDAStorageStrides() {}

    __ACC_CPU__
    DataFieldCUDAStorageStrides(const DataFieldCUDAStorageStrides& other) { *this = other; }
    __ACC_CPU__
    DataFieldCUDAStorageStrides& operator= (const DataFieldCUDAStorageStrides& other) 
    {
        strides_ = other.strides_;
        // by convention
        return *this;
    } 

    /**
    * Method initializing the strides member
    * @param strides stride instance used for initialization
    */
    __CPU__
    void Init(const StorageStridesType& strides)
    {
        strides_ = strides;
    }
    
    /**
    * Compute the stride for a given index
    * @param i index in i direction
    * @param j index in j direction
    * @param k index in k direction
    */
    __ACC__
    int ComputeStride(const int i, const int j, const int k) const
    {
        return TStorageStrides::ComputeStride(strides_, i, j, k);
    }

    /**
    * Compute the stride by explicitly specifying the offset of the block
    * @param iBlockIndex block index in i direction
    * @param jBlockIndex block index in j direction
    * @param i index inside the block in i direction
    * @param j index inside the block in j direction
    * @param k index inside the block in k direction
    */
    __ACC__
    int ComputeBlockStride(const int iBlockIndex, const int jBlockIndex, const int i, const int j, const int k) const
    {
        return TStorageStrides::ComputeBlockStride(strides_, iBlockIndex, jBlockIndex, i, j, k);
    }

    /**
    * @return storage strides
    */
    __CPU__
    const StorageStridesType& strides() const { return strides_; }

private:
    StorageStridesType strides_;
};

