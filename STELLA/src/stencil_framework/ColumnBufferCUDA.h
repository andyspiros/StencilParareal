#pragma once

#include <cassert>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/if.hpp>
#include "SharedInfrastructure.h"
#include "Domains.h"
#include "StencilSweepGroupDescriptor.h"
#include "ColumnBuffer.h"

/**
* @class ColumnBufferCUDA
* CUDA implementation of a ColumnBuffer storing a field as an array of BlockStorage elements, where each element contains an IJ block with block private boundary. 
* (note that the block private boundary avoids the need for synchronization between CUDA blocks)
*/
template<
    typename TValue,
    typename TBufferStorageFormat,
    typename TIJRange,
    typename TBlockSize,
    typename TParameterIJBoundary>
class ColumnBufferCUDA : public ColumnBuffer<TBufferStorageFormat, ColumnBufferCUDA<TValue, TBufferStorageFormat, TIJRange, TBlockSize, TParameterIJBoundary> > // CRTP
{
    DISALLOW_COPY_AND_ASSIGN(ColumnBufferCUDA);
public:    
    // make sure the boundaries are smaller than the parameter boundaries of the stencil
    BOOST_STATIC_ASSERT(
        TIJRange::IMinusOffset::value >= TParameterIJBoundary::IMinusOffset::value && 
        TIJRange::IPlusOffset::value <= TParameterIJBoundary::IPlusOffset::value 
    );
    BOOST_STATIC_ASSERT(
        TIJRange::JMinusOffset::value >= TParameterIJBoundary::JMinusOffset::value && 
        TIJRange::JPlusOffset::value <= TParameterIJBoundary::JPlusOffset::value 
    );
    
    // the i size of the block storage is padded to a multiple of alignment elements
    typedef boost::mpl::integral_c<int, cCacheLineSize / sizeof(TValue)> AlignmentElements;

    // define a boolean type which is true if no ij boundary is necessary
    typedef boost::mpl::bool_<
        TIJRange::IMinusOffset::value == 0 && TIJRange::IPlusOffset::value == 0 &&
        TIJRange::JMinusOffset::value == 0 && TIJRange::JPlusOffset::value == 0
    > HasNoBoundary;
    
    // compute the block storage ij boundary including a padding at the i plus boundary
    // (note that the smaller ij range is used and not the parameter ij boundary)
    typedef typename boost::mpl::if_<
        HasNoBoundary,
        DataFieldIJBoundary<0, 0, 0, 0>,
        DataFieldIJBoundary<
            TIJRange::IMinusOffset::value, 
            (
                (
                    (TIJRange::IPlusOffset::value - TIJRange::IMinusOffset::value + AlignmentElements::value - 1) / 
                    AlignmentElements::value
                ) *
                AlignmentElements::value + 
                TIJRange::IMinusOffset::value
            ),           
            TIJRange::JMinusOffset::value, 
            TIJRange::JPlusOffset::value
        >
    >::type BlockStorageIJBoundary;

    // define the data field ij boundary
    // (note that the parameter ij boundary needs to be covered as it is possible to execute a stencil with an apply boundary)
    typedef DataFieldIJBoundary<
        (TParameterIJBoundary::IMinusOffset::value - TBlockSize::ISize::value + 1) / TBlockSize::ISize::value, 
        (TParameterIJBoundary::IPlusOffset::value + TBlockSize::ISize::value - 1) / TBlockSize::ISize::value, 
        (TParameterIJBoundary::JMinusOffset::value - TBlockSize::JSize::value + 1) / TBlockSize::JSize::value,  
        (TParameterIJBoundary::JPlusOffset::value + TBlockSize::JSize::value - 1) / TBlockSize::JSize::value
    > StorageIJBoundary;

    // define the block storage format
    typedef BlockStorageFormat<
        TBlockSize,
        BlockStorageIJBoundary
    > BlockStorageFormatType;
    
    // check the alignment
    BOOST_STATIC_ASSERT(BlockStorageFormatType::OriginOffset::value % AlignmentElements::value == 0);
    BOOST_STATIC_ASSERT(BlockStorageFormatType::Size::value % AlignmentElements::value == 0);

    // define the buffer data field type
    typedef DataFieldStorageFormat<
        StorageIJBoundary,
        typename boost::mpl::if_<
            is_2d_domain<typename TBufferStorageFormat::Domain>, // note that the buffer storage format is nothing more than the k range
            StorageOrder::JI,
            StorageOrder::KJI
        >::type,
        DataFieldAlignment<cDimI,1>
    > StorageFormatType;

    // define storage types
    typedef BlockStorage<TValue, BlockStorageFormatType> BlockStorageType;
    typedef DataFieldCUDAStorage<BlockStorageType, StorageFormatType> StorageType;
    typedef DataFieldCUDAStorageStrides<BlockStorageStrides<StorageFormatType, BlockStorageFormatType> > StorageStridesType;
    typedef DataFieldCUDAStoragePointer<BlockStoragePointer<TValue, StorageFormatType, BlockStorageFormatType, cReadWrite> > PointerType;

    __CPU__
    ColumnBufferCUDA() {}
    __CPU__
    ~ColumnBufferCUDA() {}

    /**
    * Allocate and initialize the buffer data fields
    * @param calculationDomain size of the calculation domain
    * @param kMinusOffset boundary in k minus direction
    * @param kPlusOffset boundary in k plus direction
    */
    __CPU__
    void InitImpl(const IJKSize& calculationDomain, const int kMinusOffset, const int kPlusOffset)
    {
        // setup the origin offset
        KBoundary kBoundary;
        kBoundary.Init(kMinusOffset, kPlusOffset);

        // compute number of blocks needed to cover calculation domain size
        // (note that every IJ position holds a complete block)
        IJKSize bufferDomain;
        bufferDomain.Init(
            (calculationDomain.iSize() + TBlockSize::ISize::value - 1) / TBlockSize::ISize::value,
            (calculationDomain.jSize() + TBlockSize::JSize::value - 1) / TBlockSize::JSize::value,
            calculationDomain.kSize()
        );

        // finally initialize the storage
        storage_.Init(bufferDomain, kBoundary); 
        storage_.InitializeStoragePointerToOrigin(originStoragePointer_);
    }
    
    /**
    * @return storage pointer of the buffer
    */
    __CPU__
    const PointerType& originStoragePointer() const 
    {
        return originStoragePointer_;
    }
        
    /**
    * @return storage of the buffer
    */
    __CPU__
    const StorageType& storage() const 
    {
        return storage_;
    }

private:
    StorageType storage_;
    PointerType originStoragePointer_;
};

// CUDA specialization of create column buffer
template<
    typename TValue,
    typename TBufferStorageFormat,
    typename TBlockSize,
    typename TParameterIJBoundary,
    typename TStencilSweepGroupDescriptors,
    int VBufferIndex>
struct create_column_buffer
{
    // verify if the temporary field is accessed or if it is cached locally
    // in case the field is accessed define a column buffer type otherwise a dummy storage
    typedef typename boost::mpl::if_<
        stencil_sweep_group_descriptors_is_parameter_accessed<
            TStencilSweepGroupDescriptors,
            FullDomain,
            boost::mpl::integral_c<int, VBufferIndex>
        >,
        ColumnBufferCUDA<
            TValue,
            TBufferStorageFormat,
            typename maximum_parameter_ij_range<
                TStencilSweepGroupDescriptors, 
                boost::mpl::integral_c<int, VBufferIndex> 
            >::type,
            TBlockSize,
            TParameterIJBoundary
        >,
        DummyStorage<TValue>
    >::type type;
};


