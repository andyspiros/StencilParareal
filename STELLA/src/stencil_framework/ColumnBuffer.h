#pragma once

#include <boost/mpl/fold.hpp>
#include "SharedInfrastructure.h"
#include "IJKSize.h"
#include "ScalarStorage.h"
#include "StencilSweepGroupDescriptor.h"

/**
* @class ColumnBuffer
* Buffering class which guarantees private storage for a block of the calculation domain
*/
template<
    typename TBufferStorageFormat, 
    typename TColumnBufferImpl>
class ColumnBuffer
{
    DISALLOW_COPY_AND_ASSIGN(ColumnBuffer);
protected:
    __CPU__
    ColumnBuffer() {}
    __CPU__
    ~ColumnBuffer() {}

public:
    /**
    * Allocate and initialize the buffer
    * @param calculationDomain size of the calculation domain
    */
    __CPU__
    void Init(const IJKSize& calculationDomain)
    {
        // the storage format is nothing more than the k range of the buffer
        typedef TBufferStorageFormat BufferKRange;

        // compute the k boundary by converting the k range into a runtime k boundary
        int kMinusOffset = 
            BufferKRange::Begin::KOffset::value + 
            (BufferKRange::Begin::KReference::value == cKMaximumFlat ? cFlatLimit - 1 : 0) +
            (BufferKRange::Begin::KReference::value == cKMinimumTerrain ? cFlatLimit : 0) +
            (BufferKRange::Begin::KReference::value == cKMaximumTerrain ? calculationDomain.kSize() - 1 : 0);

        int kPlusOffset = 
            BufferKRange::End::KOffset::value - 1 + // convert half open interval in plus and minus offset
            (BufferKRange::End::KReference::value == cKMaximumFlat ? -(calculationDomain.kSize() - 1) + cFlatLimit - 1 : 0) +
            (BufferKRange::End::KReference::value == cKMinimumTerrain ? -(calculationDomain.kSize() - 1) + cFlatLimit : 0) +
            (BufferKRange::End::KReference::value == cKMinimumFlat ? -(calculationDomain.kSize() - 1) : 0);
       
        // call the implementation init
        static_cast<TColumnBufferImpl*>(this)->InitImpl(calculationDomain, kMinusOffset, kPlusOffset);
    }
};

/**
* @struct TemporaryFieldInitFunctor
* Functor used to initialize the temporary fields of a tuple
*/
struct TemporaryFieldInitFunctor
{
    // do nothing for scalar buffers
    template<typename TValue>
    static void Do(ScalarStorage<TValue, cReadWrite>&, parameter_type<const IJKSize>::type) {}
    
    // do nothing for dummy storages
    template<typename TValue>
    static void Do(DummyStorage<TValue>&, parameter_type<const IJKSize>::type) {}

    // convert column buffers to data field iterators 
    template<
        typename TBufferStorageFormat,
        typename TColumnBufferImpl>
    __CPU__
    static void Do(ColumnBuffer<TBufferStorageFormat, TColumnBufferImpl>& columnBuffer, parameter_type<const IJKSize>::type calculationDomain)
    {
        columnBuffer.Init(calculationDomain);
    }
};

/**
* @struct maximum_parameter_ij_range
* Meta function returning the maximal ij range of all stencil stages using a given parameter
*/
template<
    typename TStencilSweepGroupDescriptors,
    typename TParameterIndex>
struct maximum_parameter_ij_range :
    boost::mpl::fold<
        TStencilSweepGroupDescriptors,
        IJRange<cComplete, 0, 0, 0, 0>,
        maximum_ij_range<
            boost::mpl::_1, 
            stencil_sweep_group_descriptor_parameter_ij_range<boost::mpl::_2, TParameterIndex> 
        >
    >
{};

/**
* @struct create_column_buffer
* Meta function creating a column buffer
*/
template<
    typename TValue,
    typename TBufferStorageFormat,
    typename TBlockSize,
    typename TParameterIJBoundary,
    typename TStencilSweepGroupDescriptors,
    int VBufferIndex>
struct create_column_buffer;

/**
* @struct create_temporary_field
* Meta function computing the type of the buffer instance given a buffer type and the stencil stages using the buffer.
*/
template<
    typename TTemporaryField,
    typename TBlockSize,
    typename TParameterIJBoundary,
    typename TStencilSweepGroupDescriptors>
struct create_temporary_field;

// specialization for scalar stage variables
template<
    int VBufferIndex,
    typename TValue,
    typename TBlockSize,
    typename TParameterIJBoundary,
    typename TStencilSweepGroupDescriptors>
struct create_temporary_field<
    StageVariable<VBufferIndex, TValue>,
    TBlockSize,
    TParameterIJBoundary,
    TStencilSweepGroupDescriptors>
{
    typedef ScalarStorage<TValue, cReadWrite> type;
};

// specialization for stencil buffers
template<
    int VBufferIndex,
    typename TValue,
    typename TKRange,
    typename TBlockSize,
    typename TParameterIJBoundary,
    typename TStencilSweepGroupDescriptors>
struct create_temporary_field< 
    StencilBuffer<VBufferIndex, TValue, TKRange>,
    TBlockSize,
    TParameterIJBoundary,
    TStencilSweepGroupDescriptors>
{
    typedef typename create_column_buffer<
        TValue,
        TKRange,
        TBlockSize,
        TParameterIJBoundary,
        TStencilSweepGroupDescriptors,
        VBufferIndex
    >::type type;
};

/*
// specialization for sweep buffers
// Note: the KWindow is currently ignored, so that this is equivalent to declaring a stencil buffer
template<
    int VBufferIndex,
    typename TValue,
    typename TKRange,
    typename TKWindow,
    typename TBlockSize,
    typename TParameterIJBoundary,
    typename TStencilSweepGroupDescriptors>
struct create_temporary_field< 
    SweepBuffer<VBufferIndex, TValue, TKRange, TKWindow>,
    TBlockSize,
    TParameterIJBoundary,
    TStencilSweepGroupDescriptors>
{
    typedef typename create_column_buffer<
        TValue,
        TKRange,
        TBlockSize,
        TParameterIJBoundary,
        TStencilSweepGroupDescriptors,
        VBufferIndex
    >::type type;
};
*/

