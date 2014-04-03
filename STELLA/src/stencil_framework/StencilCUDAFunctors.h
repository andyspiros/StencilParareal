#pragma once

#include <cassert>
#include <boost/mpl/void.hpp>
#include "SharedInfrastructure.h"
#include "ParameterWrapper.h"
#include "ColumnBufferCUDA.h"

/**
* @struct parameter_to_apply
* Meta function defining the type conversion of parameter tuple elements to apply tuple elements
*/
template<
    typename T,
    typename TBlockSize>
struct parameter_to_apply;

// convert wrapped data fields to storage pointers
template<
    typename TValue,
    typename TStorageFormat, 
    ParameterIntend VParameterIntend, 
    int VParameterIndex,
    typename TBlockSize>
struct parameter_to_apply<ParameterWrapper<DataFieldCUDA<TValue, TStorageFormat>, VParameterIntend, VParameterIndex>, TBlockSize>
{
    typedef typename create_storage_pointer<DataFieldCUDA<TValue, TStorageFormat>, TBlockSize, VParameterIntend, VParameterIndex>::type type;
};

// convert wrapped data fields to storage pointers
template<
    typename TDataField, 
    ParameterIntend VParameterIntend, 
    int VParameterIndex,
    typename TBlockSize>
struct parameter_to_apply<ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex>, TBlockSize>
{
    typedef typename create_storage_pointer<TDataField, TBlockSize, VParameterIntend, VParameterIndex>::type type;
};

// convert wrapped scalars to scalar storages
template<
    typename TValue,
    int VParameterIndex,
    typename TBlockSize>
struct parameter_to_apply<ParameterWrapper<TValue, cScalar, VParameterIndex>, TBlockSize>
{
    typedef ScalarStorage<TValue, cReadOnly> type;
};

/**
* @struct ParameterToApplyFunctor
* Functor converting the parameter tuple to the apply tuple e.g. by converting fields into iterators
*/
template<typename TBlockSize>
struct ParameterToApplyFunctor
{
    // convert data fields
    template<
        typename TValue,
        typename TStorageFormat,
        ParameterIntend VParameterIntend, 
        int VParameterIndex>
    __CPU__
    static void Do(
        ParameterWrapper<DataFieldCUDA<TValue, TStorageFormat>, VParameterIntend, VParameterIndex>& in,
        typename create_storage_pointer<DataFieldCUDA<TValue, TStorageFormat>, TBlockSize, VParameterIntend, VParameterIndex>::type& out)
    {
        in.Unwrap().deviceStorage().InitializeStoragePointerToOrigin(out);
    }

    // convert joker fields
    template<
        typename TDataField,
        ParameterIntend VParameterIntend, 
        int VParameterIndex>
    __CPU__
    static void Do(
        ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex>& in, 
        typename create_storage_pointer<TDataField, TBlockSize, VParameterIntend, VParameterIndex>::type& out)
    {
        in.Unwrap().dataField().deviceStorage().InitializeStoragePointerToOrigin(out);
    }

    // unwrap scalar parameters
    template<
        typename TValue, 
        int VParameterIndex>
    __CPU__
    static void Do(
        ParameterWrapper<TValue, cScalar, VParameterIndex>& in,
        ScalarStorage<TValue, cReadOnly>& out)
    {
        out.set_value(in.Unwrap());
    }
};

/**
* @struct temporary_field_to_apply
* Meta function defining the type conversion of temporary field tuple elements to apply tuple elements
*/
template<typename TValue>
struct temporary_field_to_apply;

// convert scalar wrappers
template<typename TValue>
struct temporary_field_to_apply<ScalarStorage<TValue, cReadWrite> >
{
    typedef ScalarStorage<TValue, cReadWrite> type;
};

// convert dummy storages
template<typename TValue>
struct temporary_field_to_apply<DummyStorage<TValue> >
{
    typedef DummyStorage<TValue> type;
};

// convert column buffers to iterators
template< 
    typename TValue,
    typename TBufferStorageFormat,
    typename TIJRange, 
    typename TBlockSize,
    typename TParameterIJBoundary>
struct temporary_field_to_apply<ColumnBufferCUDA<TValue, TBufferStorageFormat, TIJRange, TBlockSize, TParameterIJBoundary> >
{
    typedef typename ColumnBufferCUDA<TValue, TBufferStorageFormat, TIJRange, TBlockSize, TParameterIJBoundary>::PointerType type;
};

/**
* @struct TemporaryFieldToApplyFunctor
* Functor converting the buffer tuple the apply tuple
*/
struct TemporaryFieldToApplyFunctor
{
    // do nothing for scalar buffers
    // (note that scalar storages are stored in the register tuple therefore the apply tuple will return void)
    template<typename TValue>
    __CPU__
    static void Do(ScalarStorage<TValue, cReadWrite>& in, boost::mpl::void_) {} 

    // do nothing for dummy storages
    template<typename TValue>
    __CPU__
    static void Do(DummyStorage<TValue>& in, DummyStorage<TValue>& out) {}

    // convert column buffers to data field iterators 
    template<
        typename TValue,
        typename TBufferStorageFormat,
        typename TIJRange, 
        typename TBlockSize,
        typename TParameterIJBoundary>
    __CPU__
    static void Do(
        ColumnBufferCUDA<TValue, TBufferStorageFormat, TIJRange, TBlockSize, TParameterIJBoundary>& in, 
        typename ColumnBufferCUDA<TValue, TBufferStorageFormat, TIJRange, TBlockSize, TParameterIJBoundary>::PointerType& out)
    {
        out = in.originStoragePointer(); 
    }
};

/**
* @struct UpdateApplyFunctor
* Functor making sure the device copies of the data fields are valid. Additionally
* the joker fields are checked for modifications. If pointers are modified 
* the boolean parameter is set to false signaling a host to device copy of the
* apply tuple is necessary.
*/
template<typename TBlockSize>
struct UpdateApplyFunctor
{
    // update data fields
    template<
        typename TValue,
        typename TStorageFormat,
        ParameterIntend VParameterIntend,
        int VParameterIndex>
    __CPU__
    static void Do(
        ParameterWrapper<DataFieldCUDA<TValue, TStorageFormat>, VParameterIntend, VParameterIndex>& in,
        typename create_storage_pointer<DataFieldCUDA<TValue, TStorageFormat>, TBlockSize, VParameterIntend, VParameterIndex>::type& out, 
        parameter_type<bool&>::type applyDataModified)
    {       
        if(!out.PointsToOrigin(in.Unwrap().deviceStorage().pStorageOrigin()))
        {
            // update the storage pointer
            in.Unwrap().deviceStorage().InitializeStoragePointerToOrigin(out);
            applyDataModified = true;
        }
        // make sure the device copy is valid
        in.Unwrap().SynchronizeDeviceStorage();
    }

    // update joker fields
    template<
        typename TDataField,
        ParameterIntend VParameterIntend,
        int VParameterIndex>
    __CPU__
    static void Do(
        ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex>& in, 
        typename create_storage_pointer<TDataField, TBlockSize, VParameterIntend, VParameterIndex>::type& out, 
        parameter_type<bool&>::type applyDataModified)
    {
        if(!out.PointsToOrigin(in.Unwrap().dataField().deviceStorage().pStorageOrigin()))
        {
            // update the storage pointer
            in.Unwrap().dataField().deviceStorage().InitializeStoragePointerToOrigin(out);
            applyDataModified = true;
        }
        // make sure the device copy is valid
        in.Unwrap().dataField().SynchronizeDeviceStorage();
    }

    // Update scalar parameters
    template<
        typename TValue,
        int VParameterIndex>
    __CPU__
    static void Do(
        ParameterWrapper<TValue, cScalar, VParameterIndex>& in, 
        ScalarStorage<TValue, cReadOnly>& out, 
        parameter_type<bool&>::type applyDataModified)
    {
        if(in.Unwrap() != out.value())
        {
            // update the storage pointer
            out.set_value(in.Unwrap());
            applyDataModified = true;
        }
    }
};

/**
* @struct ParameterStridesInitFunctor
* Functor initializing the strides used by the parameter fields of the CUDA context. 
* In case the strides for the corresponding data field type are uninitialized they are copied.
* Otherwise an assert checks that all fields of the same data field type have identical strides.
*/
template<
    typename TBlockSize,
    typename TStorageStridesTuple>
struct ParameterStridesInitFunctor
{
    // init data fields
    template<
        typename TValue,
        typename TStorageFormat,
        ParameterIntend VParameterIntend,
        int VParameterIndex>
    __CPU__
    static void Do(
        ParameterWrapper<DataFieldCUDA<TValue, TStorageFormat>, VParameterIntend, VParameterIndex>& dataField, 
        typename parameter_type<TStorageStridesTuple>::type storageStridesTuple)
    {
        typedef typename create_storage_pointer<DataFieldCUDA<TValue, TStorageFormat>, TBlockSize, VParameterIntend, VParameterIndex>::type StoragePointerType;
        typedef typename create_storage_strides<StoragePointerType>::type StorageStridesType;
        dataField.Unwrap().deviceStorage().TryToInitializeStorageStrides(storageStridesTuple(static_cast<StorageStridesType*>(0)));
    }

    // init joker fields
    template<
        typename TDataField,
        ParameterIntend VParameterIntend,
        int VParameterIndex>
    __CPU__
    static void Do(
        ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex>& dataField, 
        typename parameter_type<TStorageStridesTuple>::type storageStridesTuple)
    {
        typedef typename create_storage_pointer<TDataField, TBlockSize, VParameterIntend, VParameterIndex>::type StoragePointerType;
        typedef typename create_storage_strides<StoragePointerType>::type StorageStridesType;
        dataField.Unwrap().dataField().deviceStorage().TryToInitializeStorageStrides(storageStridesTuple(static_cast<StorageStridesType*>(0)));
    }
};

/**
* @struct TemporaryFieldStridesInitFunctor
* Functor initializing the strides used by the buffers of the CUDA context. 
* In case the strides for the corresponding data field type are uninitialized they are copied.
* Otherwise an assert checks that all fields of the same data field type have identical strides.
*/
template<typename TStorageStridesTuple>
struct TemporaryFieldStridesInitFunctor
{
    // do nothing for scalar buffers
    template<typename TValue>
    __CPU__
    static void Do(ScalarStorage<TValue, cReadWrite>& in, typename parameter_type<TStorageStridesTuple>::type) {}

    // do nothing for dummy storages
    template<typename TValue>
    __CPU__
    static void Do(DummyStorage<TValue>& in, typename parameter_type<TStorageStridesTuple>::type) {}

    // init column buffers
    template<
        typename TValue,
        typename TBufferStorageFormat,
        typename TIJRange, 
        typename TBlockSize,
        typename TParameterIJBoundary>
    __CPU__
    static void Do(
        ColumnBufferCUDA<TValue, TBufferStorageFormat, TIJRange, TBlockSize, TParameterIJBoundary>& buffer, 
        typename parameter_type<TStorageStridesTuple>::type storageStridesTuple)
    { 
        typedef typename ColumnBufferCUDA<TValue, TBufferStorageFormat, TIJRange, TBlockSize, TParameterIJBoundary>::StorageStridesType StorageStridesType;
        buffer.storage().TryToInitializeStorageStrides(storageStridesTuple(static_cast<StorageStridesType*>(0)));
    }
};





