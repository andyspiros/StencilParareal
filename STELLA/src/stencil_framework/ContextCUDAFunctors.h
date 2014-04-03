#pragma once

#include <boost/mpl/void.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/and.hpp>
#include "SharedInfrastructure.h"

/**
* @struct is_register_data
* Meta function returning true for all elements that shall be stored in a register
*/
template<typename T>
struct is_register_data : boost::mpl::false_ {};

template<typename T>
struct is_register_data<ScalarStorage<T, cReadWrite> > : boost::mpl::true_ {};

/**
* @struct storage_strides_elements
* Meta function calculating the storage strides set needed for a certain parameter set
*/
template<typename TContextDescriptor>
struct storage_strides_elements 
{
    // compute a unique set of storage strides
    typedef typename boost::mpl::fold<
        typename boost::mpl::fold<
            typename TContextDescriptor::ElementTypes,
            boost::mpl::set0<>, // first insert strides in a set (a set only contains one element per type)
            boost::mpl::if_<
                is_iterable<boost::mpl::_2>,
                boost::mpl::insert<boost::mpl::_1, create_storage_strides<boost::mpl::_2> >,
                boost::mpl::_1
            >
        >::type,
        boost::mpl::vector0<>,
        boost::mpl::push_back<boost::mpl::_1, boost::mpl::_2> // fill the set into a vector as a random access sequence is expected by the tuple
    >::type StorageStridesSet;

    // define the tuple elements
    typedef TupleElements<StorageStridesSet, StorageStridesSet> type;
};

/**
* @struct shared_data_elements
* Meta function calculating the shared data elements needed to store a certain parameter set
*/
template<typename TContextDescriptor>
struct shared_data_elements 
{
    // partition the element indexes into register and const tuple indexes
    typedef typename boost::mpl::fold<
        typename TContextDescriptor::ElementIndexes,
        boost::mpl::vector0_c<int>,
        boost::mpl::if_<
            is_register_data<boost::mpl::at<typename TContextDescriptor::IndexToElementMap, boost::mpl::_2> >, 
            boost::mpl::_1,
            boost::mpl::push_back<boost::mpl::_1, boost::mpl::_2>
        >
    >::type SharedDataIndexes;

    // define the data type stored in shared memory
    typedef TupleElements<
        SharedDataIndexes,
        typename boost::mpl::transform<
            SharedDataIndexes,
            boost::mpl::at<typename TContextDescriptor::IndexToElementMap, boost::mpl::_>
        >::type            
    > type;
};

/**
* @struct register_data_elements
* Meta function computing the register data tuple element type
*/
template<
    typename TContextDescriptor,
    typename TStencilSweepDescriptor>
struct register_data_elements 
{
    // define the register tuple indexes
    typedef typename boost::mpl::fold<
        typename TContextDescriptor::ElementIndexes,
        boost::mpl::vector0_c<int>,
        boost::mpl::if_<
            boost::mpl::and_<
                is_register_data<boost::mpl::at<typename TContextDescriptor::IndexToElementMap, boost::mpl::_2> >,    
                stencil_sweep_descriptor_uses_parameter<TStencilSweepDescriptor, FullDomain, boost::mpl::_2>
            >,
            boost::mpl::push_back<boost::mpl::_1, boost::mpl::_2>,
            boost::mpl::_1
        >
    >::type RegisterIndexes;

    // compute the register tuple elements
    typedef TupleElements<
        RegisterIndexes,
        typename boost::mpl::transform<
            RegisterIndexes,
            boost::mpl::at<typename TContextDescriptor::IndexToElementMap, boost::mpl::_>
        >::type
    > type;
};

/**
* @struct IJOffset
* Structure storing an ij offset
*/
struct IJOffset
{
    int iBlockIndex;
    int jBlockIndex;
    int i;
    int j;
};

/**
* @struct IJKOffset
* Structure storing an ijk offset
*/
struct IJKOffset
{
    int iBlockIndex;
    int jBlockIndex;
    int i;
    int j;
    int k;
};

/**
* @struct AdvanceFunctor
* Functor advancing the index position by a constant offset
*/
template<int VI, int VJ, int VK>
struct AdvanceFunctor
{
    // as we work on strides no empty method needed
    template<typename TStorageStrides>
    __ACC__
    static void Do(
        const DataFieldCUDAStorageStrides<TStorageStrides>& storageStrides, 
        DataFieldCUDAStorageIndex<DataFieldCUDAStorageStrides<TStorageStrides> >& storageIndex)
    {
        storageIndex.Advance(storageStrides, VI, VJ, VK);
    }
};

/**
* @struct AdvanceFunctor
* Functor advancing index positions by a k offset
*/
template<int VI, int VJ>
struct AdvanceInKFunctor
{
    // as we work on strides no empty method needed
    template<typename TStorageStrides>
    __ACC__
    static void Do(
        const DataFieldCUDAStorageStrides<TStorageStrides>& storageStrides, 
        DataFieldCUDAStorageIndex<DataFieldCUDAStorageStrides<TStorageStrides> >& storageIndex,
        parameter_type<const int>::type kOffset)
    {
        storageIndex.Advance(storageStrides, VI, VJ, kOffset);
    }
};

/**
* @struct RestoreAndAdvanceInIJFunctor
* Functor restoring the storage positions and advancing the iterators by an ij offset 
*/
template<int VK>
struct RestoreAndAdvanceInIJFunctor
{
    // as we work on strides no empty method needed
    template<typename TStorageStrides>
    __ACC__
    static void Do(
        const DataFieldCUDAStorageStrides<TStorageStrides>& storageStrides, 
        DataFieldCUDAStorageIndex<DataFieldCUDAStorageStrides<TStorageStrides> >& storageIndex, 
        parameter_type<const IJOffset>::type ijOffset)
    {
        storageIndex.Init(storageStrides, ijOffset.iBlockIndex, ijOffset.jBlockIndex, ijOffset.i, ijOffset.j, VK);
    }
};

/**
* @struct RestoreAndAdvanceInIJKFunctor
* Functor restoring the storage positions and advancing the iterators by an ijk offset 
*/
struct RestoreAndAdvanceInIJKFunctor
{
    // as we work on strides no empty method needed
    template<typename TStorageStrides>
    __ACC__
    static void Do(
        const DataFieldCUDAStorageStrides<TStorageStrides>& storageStrides, 
        DataFieldCUDAStorageIndex<DataFieldCUDAStorageStrides<TStorageStrides> >& storageIndex, 
        parameter_type<const IJKOffset>::type ijkOffset)
    {
        storageIndex.Init(storageStrides, ijkOffset.iBlockIndex, ijkOffset.jBlockIndex, ijkOffset.i, ijkOffset.j, ijkOffset.k);
    }
};




