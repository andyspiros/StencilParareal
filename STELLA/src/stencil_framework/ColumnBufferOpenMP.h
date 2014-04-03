#pragma once

#include <cassert>
#include <vector>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/if.hpp>
#include "SharedInfrastructure.h"
#include "Domains.h"
#include "ColumnBuffer.h"

/**
* @class ColumnBufferOpenMP
* OpenMP implementation of a column buffer providing memory for one block
* (note the stencil allocates one column buffer per OpenMP thread)
*/
template<
    typename TValue,
    typename TBufferStorageFormat,
    typename TIJRange,
    typename TBlockSize>
class ColumnBufferOpenMP : public ColumnBuffer<TBufferStorageFormat, ColumnBufferOpenMP<TValue, TBufferStorageFormat, TIJRange, TBlockSize> > // CRTP
{
public:
    // define storage and iterator types
    typedef DataFieldOpenMPStorage<
        TValue,
        DataFieldStorageFormat<
            DataFieldIJBoundary<
                TIJRange::IMinusOffset::value, TIJRange::IPlusOffset::value,
                TIJRange::JMinusOffset::value, TIJRange::JPlusOffset::value
            >,
            typename boost::mpl::if_<
                // note that the buffer storage format is nothing more than the k range
                is_2d_domain<typename TBufferStorageFormat::Domain>,
                StorageOrder::JI,
                StorageOrder::JIK
            >::type,
            OpenMPAlignment
        >
    > StorageType;
    typedef typename StorageType::StorageIteratorType StorageIteratorType;
 
    ColumnBufferOpenMP() {};
    ~ColumnBufferOpenMP() {};

    // copy constructor
    ColumnBufferOpenMP(const ColumnBufferOpenMP& other) { *this = other; }
    ColumnBufferOpenMP& operator= (const ColumnBufferOpenMP& other)
    {
        buffer_ = other.buffer_;
        return *this;
    }

    /**
    * Allocate and initialize the buffer data fields
    * @param calculationDomain size of the calculation domain
    * @param kMinusOffset boundary in k minus direction
    * @param kPlusOffset boundary in k plus direction
    */
    void InitImpl(const IJKSize& calculationDomain, const int kMinusOffset, const int kPlusOffset)
    {
        // calculate size of the buffer
        IJKSize bufferDomain;
        bufferDomain.Init(TBlockSize::ISize::value, TBlockSize::JSize::value, calculationDomain.kSize());

        // prepare the k boundary
        KBoundary kBoundary;
        kBoundary.Init(kMinusOffset, kPlusOffset);
        
        buffer_.Init(bufferDomain, kBoundary);
    }
    
    /**
    * @return iterator pointing to the origin of the buffer
    */
    const StorageIteratorType& originIterator() const
    {
        // return the pointer associated with the buffer
        return buffer_.originIterator();
    }

private:
    StorageType buffer_;
};

// OpenMP specialization of create column buffer
template<
    typename TValue,
    typename TBufferStorageFormat,
    typename TBlockSize,
    typename TParameterIJBoundary,
    typename TStencilSweepGroupDescriptors,
    int VBufferIndex>
struct create_column_buffer
{
    typedef ColumnBufferOpenMP<
        TValue,
        TBufferStorageFormat,
        typename maximum_parameter_ij_range<
            TStencilSweepGroupDescriptors, 
            boost::mpl::integral_c<int, VBufferIndex> 
        >::type,
        TBlockSize
    > type;
};
