#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include "Definitions.h"
#include "Enums.h"
#include "BlockSize.h"
#include "DataFieldIJBoundary.h"

/**
* @struct BlockStorageFormat
* Structure storing the format information as size and boundary if a block storage
* (assume that i minus and plus offset together add up to a cache line size
* --> add them at the beginning and end to get an aligned origin and a padded size)
*/
template<
    typename TBlockSize,
    typename TIJBoundary> 
struct BlockStorageFormat
{
    BOOST_STATIC_ASSERT(is_data_field_ij_boundary<TIJBoundary>::value);

    // define size information
    typedef boost::mpl::integral_c<int, TBlockSize::ISize::value - TIJBoundary::IMinusOffset::value + TIJBoundary::IPlusOffset::value> ISize;
    typedef boost::mpl::integral_c<int, TBlockSize::JSize::value - TIJBoundary::JMinusOffset::value + TIJBoundary::JPlusOffset::value> JSize;
    // add i minus and i plus offset at end and beginning of block in order to keep origin aligned
    typedef boost::mpl::integral_c<int, ISize::value * JSize::value + TIJBoundary::IPlusOffset::value - TIJBoundary::IMinusOffset::value> Size;
    
    // origin offset
    typedef boost::mpl::integral_c<int, TIJBoundary::IPlusOffset::value - TIJBoundary::IMinusOffset::value - TIJBoundary::JMinusOffset::value * ISize::value> OriginOffset;

    // method computing block strides 
    __ACC_CPU__
    static int ComputeStride(const int i, const int j)
    {
        ACC_ASSERT(i > TIJBoundary::IMinusOffset::value && i < TBlockSize::ISize::value + TIJBoundary::IPlusOffset::value);
        ACC_ASSERT(j > TIJBoundary::JMinusOffset::value && i < TBlockSize::JSize::value + TIJBoundary::JPlusOffset::value);
        return i + j * ISize::value;
    }
};

/**
* @class BlockStorage
* Class storing a fixed size block in i and j dimension
*/
template<
    typename TValue,
    typename TStorageFormat> 
class BlockStorage
{
public:
    typedef TValue ValueType;
    typedef ValueType* PointerType;
    typedef ValueType& ReferenceType;
    // define the return type
    typedef ReferenceType ReturnType;

    // define the storage format
    typedef TStorageFormat StorageFormat;

    __ACC_CPU__
    BlockStorage() {}
    __ACC_CPU__
    ~BlockStorage() {}
    
    // copy constructor and assignment
    __ACC_CPU__
    BlockStorage(const BlockStorage& other) { *this = other; }
    __ACC_CPU__
    BlockStorage& operator= (const BlockStorage& other)
    {
        // copy the values
        for(int i = 0; i < StorageFormat::Size::value; ++i)
        {
            values_[i] = other.values_[i];
        }
        // by convention, always return *this
        return *this;
    }
   
    /**
    * Access operator 
    * @param i position in i dimension
    * @param j position in j dimension
    * @return a reference to the value at the position
    */
    __ACC__
    ReturnType At(const int i, const int j) 
    { 
        return values_[StorageFormat::OriginOffset::value + StorageFormat::ComputeStride(i,j)];
    }

private:
    ValueType values_[StorageFormat::Size::value];
};

  
