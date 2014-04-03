#pragma once

#include <boost/mpl/void.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_back.hpp>
#include "SharedInfrastructure.h"
#include "ParameterWrapper.h"
#include "ColumnBufferOpenMP.h"
#include "StencilStage.h"
#include "StencilSweepDescriptor.h"

#include <boost/mpl/vector/vector10.hpp>

#ifndef NDEBUG
#include <string>
#include <sstream>
#endif

/**
* @struct parameter_to_apply
* Meta function defining the type conversion of parameter tuple elements to apply tuple elements
*/
template<typename T>
struct parameter_to_apply;

// convert wrapped data fields to iterators
template<
    typename TValue,
    typename TStorageFormat, 
    ParameterIntend VParameterIntend, 
    int VParameterIndex>
struct parameter_to_apply<ParameterWrapper<DataFieldOpenMP<TValue, TStorageFormat>, VParameterIntend, VParameterIndex> >
{
    typedef typename create_storage_iterator<DataFieldOpenMP<TValue, TStorageFormat>, VParameterIntend>::type type;
};

// convert wrapped data fields to iterators
template<
    typename TDataField, 
    ParameterIntend VParameterIntend, 
    int VParameterIndex>
struct parameter_to_apply<ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex> >
{
    typedef typename create_storage_iterator<TDataField, VParameterIntend>::type type;
};

// convert wrapped scalars to scalar iterators
template<
    typename TValue, 
    int VParameterIndex>
struct parameter_to_apply<ParameterWrapper<TValue, cScalar, VParameterIndex> >
{
    typedef ScalarStorage<TValue, cReadOnly> type;
};

/**
* @struct ParameterToApplyFunctor
* Functor converting the parameter tuple to the apply tuple e.g. by converting fields into iterators
*/
struct ParameterToApplyFunctor
{
    // Convert data fields to iterators
    template<
        typename TValue,
        typename TStorageFormat,
        ParameterIntend VParameterIntend, 
        int VParameterIndex>
    static void Do(
        ParameterWrapper<DataFieldOpenMP<TValue, TStorageFormat>, VParameterIntend, VParameterIndex>& in, 
        typename create_storage_iterator<DataFieldOpenMP<TValue, TStorageFormat>, VParameterIntend>::type& out)
    {
        in.Unwrap().storage().InitializeStorageIteratorToOrigin(out);
    }
    
    // Convert joker data fields to iterators
    template<
        typename TDataField,
        ParameterIntend VParameterIntend, 
        int VParameterIndex>
    static void Do(
        ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex>& in, 
        typename create_storage_iterator<TDataField, VParameterIntend>::type& out)
    {
        in.Unwrap().dataField().storage().InitializeStorageIteratorToOrigin(out);
    }

    // Unwrap scalar parameters
    template<
        typename TValue, 
        int VParameterIndex>
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

// convert scalar wrappers for stage variables
template<typename TValue>
struct temporary_field_to_apply<ScalarStorage<TValue,cReadWrite> >
{
    typedef ScalarStorage<TValue,cReadWrite> type;
};

// convert column buffers to iterators
template< 
    typename TBufferStorageFormat,
    typename TKRange,
    typename TIJRange, 
    typename TBlockSize>
struct temporary_field_to_apply<ColumnBufferOpenMP<TBufferStorageFormat, TKRange, TIJRange, TBlockSize> >
{
    typedef typename ColumnBufferOpenMP<TBufferStorageFormat, TKRange, TIJRange, TBlockSize>::StorageIteratorType type;
};

/**
* @struct TemporaryFieldToApplyFunctor
* Functor converting the temporary field tuple to the apply tuple
*/
struct TemporaryFieldToApplyFunctor
{
    // do nothing for scalar buffers
    template<typename TValue>
    static void Do(ScalarStorage<TValue, cReadWrite>& in, ScalarStorage<TValue, cReadWrite>& out) {}

    // convert column buffers to data field iterators 
    template<
        typename TBufferStorageFormat,
        typename TKRange,
        typename TIJRange, 
        typename TBlockSize>
    static void Do(
        ColumnBufferOpenMP<TBufferStorageFormat, TKRange, TIJRange, TBlockSize>& in, 
        typename ColumnBufferOpenMP<TBufferStorageFormat, TKRange, TIJRange, TBlockSize>::StorageIteratorType& out)
    {
        out = in.originIterator(); 
    }
};

/**
* @struct UpdateApplyFunctor
* Functor updating the apply tuple reading in updated positions and memorizing the start position
*/
struct UpdateApplyFunctor
{    
    // Set the data field pointers the the origin of the calculation domain 
    template<
        typename TValue,
        typename TStorageFormat,
        ParameterIntend VParameterIntend,
        int VParameterIndex>
    static void Do(
        ParameterWrapper<DataFieldOpenMP<TValue, TStorageFormat>, VParameterIntend, VParameterIndex>& in, 
        typename create_storage_iterator<DataFieldOpenMP<TValue, TStorageFormat>, VParameterIntend>::type& out)
    {
        in.Unwrap().storage().MemorizePositionOfOrigin(out);
    }
  
    // Set the data field pointers the the origin of the calculation domain
    template<
        typename TDataField,
        ParameterIntend VParameterIntend,
        int VParameterIndex>
    static void Do(
        ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex>& in, 
        typename create_storage_iterator<TDataField, VParameterIntend>::type& out)
    {
        in.Unwrap().dataField().storage().MemorizePositionOfOrigin(out);
    }

    // Update scalar parameters
    template<
        typename TValue,
        int VParameterIndex>
    static void Do(
        ParameterWrapper<TValue, cScalar, VParameterIndex>& in, 
        ScalarStorage<TValue, cReadOnly>& out)
    {
        out.set_value(in.Unwrap());
    }
};

#ifndef NDEBUG

// functors used for bounding box calculation

/**
* @struct ResetBoundingBoxFunctor
* Functor resetting the parameter access bounding boxes
*/
struct ResetBoundingBoxFunctor
{
    // Ignore all scalar types
    template<typename T>
    static void Do(T&) {}

    // Reset bounding box for all iterators
    template<
        typename TValue,
        typename TStorageStrides,
        ParameterIntend VParameterIntend, 
        AccessPolicy VAccessPolicy>
    static void Do(DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& iterator)
    {        
        iterator.positionChecker().ResetBoundingBox();
    }
};

/**
* @struct MergeBoundingBoxFunctor
* Functor merging the parameter access bounding boxes
*/
struct MergeBoundingBoxFunctor
{
    // Ignore all scalar types
    template<
        typename T1, 
        typename T2>
    static void Do(T1&, T2&) {}
    
    // Merge bounding box for all iterators
    template<
        typename TValue,
        typename TStorageStrides,
        AccessPolicy VAccessPolicy>
    static void Do(
        DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& in, 
        DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& out)
    {        
        out.positionChecker().MergeBoundingBox(in.positionChecker());
    }
};

/**
* @struct PrintBoundingBoxFunctor
* Functor printing the parameter bounding boxes
*/
struct PrintBoundingBoxFunctor
{
    // Ignore all scalar types
    template<
        typename T1, 
        typename T2>
    static void Do(T1&, T2&, parameter_type<std::ostream>::type) {}

    // Print bounding box for all iterators
    template<
        typename TDataField,
        typename TValue,
        typename TStorageStrides,
        AccessPolicy VAccessPolicy>
    static void Do(
        TDataField& dataField,
        DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& iterator,
        parameter_type<std::ostringstream>::type boundingBoxInfo)
    {    
        // get the positions and subtract origin
        IJKIndex floor = iterator.positionChecker().boundingBoxFloor();
        floor.AddOffset(
            -iterator.positionChecker().currentPosition().iIndex(),
            -iterator.positionChecker().currentPosition().jIndex(),
            -iterator.positionChecker().currentPosition().kIndex()
        );
        IJKIndex ceiling = iterator.positionChecker().boundingBoxCeiling();
        ceiling.AddOffset(
            -iterator.positionChecker().currentPosition().iIndex(),
            -iterator.positionChecker().currentPosition().jIndex(),
            -iterator.positionChecker().currentPosition().kIndex()
        );

        // format the string
        std::string whiteSpaces;
        whiteSpaces.insert(0, std::max(0,(int)(30-dataField.Unwrap().name().length())), ' ');
        boundingBoxInfo 
            << dataField.Unwrap().name() << whiteSpaces;
        if (VAccessPolicy == cReadOnly) {
          boundingBoxInfo << " R  ";
        } else {
          boundingBoxInfo << " RW ";
        }
        boundingBoxInfo 
            << " from\t" << floor.ToString() << "\t" 
            << "to\t" << ceiling.ToString() << std::endl;
    }
};

#endif




  
