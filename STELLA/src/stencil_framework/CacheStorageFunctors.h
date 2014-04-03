#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/next.hpp>
#include "SharedInfrastructure.h"
#include "Offset.h"

/**
* @struct compute_k_cache_storage_array_index
* Meta function computing a k cache storage array index given a data parameter k offset
*/
template<
    typename TKOffset, 
    typename TKWindow,
    KLoopDirection VKLoopDirection>
struct compute_k_cache_storage_array_index
{
    // check that the offset is in a valid range
    BOOST_STATIC_ASSERT(
        TKOffset::value >= TKWindow::KMinusOffset::value && 
        TKOffset::value <= TKWindow::KPlusOffset::value
    );

    // compute the cache array index
    typedef typename boost::mpl::if_c<
        VKLoopDirection == cKIncrement,
        boost::mpl::integral_c<int, TKWindow::KPlusOffset::value - TKOffset::value>,
        boost::mpl::integral_c<int, TKOffset::value - TKWindow::KMinusOffset::value>
    >::type type;
};

/**
* @struct KCacheStorageFillFunctor
* Functor filling an array element given a data parameter
*/
template<
    typename TDataParameter,
    typename TKWindow,
    KLoopDirection VKLoopDirection>
struct KCacheStorageFillFunctor
{
    template<
        typename TCacheArray,
        typename TContext,
        typename TKOffset>
    __ACC__
    static void Do(
        typename parameter_type<TCacheArray>::type cacheArray,
        typename parameter_type<TContext>::type context)
    {
        // compute the data parameter offset
        typedef Offset<0, 0, TKOffset::value> DataParameterOffset;

        // compute the array index
        typedef typename compute_k_cache_storage_array_index<
            TKOffset,
            TKWindow,
            VKLoopDirection
        >::type ArrayIndex;
        
        // fill the cache array at the given index
        cacheArray.At(static_cast<ArrayIndex*>(0)) = 
            context.BypassCache(TDataParameter::At(DataParameterOffset()));
    }
};

/**
* @struct KCacheStorageFillFunctor
* Functor filling an array element given a data parameter
*/
template<
    typename TDataParameter,
    typename TKWindow,
    KLoopDirection VKLoopDirection>
struct KCacheStorageFlushFunctor
{
    template<
        typename TCacheArray,
        typename TContext,
        typename TKOffset>
    __ACC__
    static void Do(
        typename parameter_type<TCacheArray>::type cacheArray,
        typename parameter_type<TContext>::type context)
    {
        // compute the data parameter offset
        typedef Offset<0, 0, TKOffset::value> DataParameterOffset;

        // compute the array index
        typedef typename compute_k_cache_storage_array_index<
            TKOffset,
            TKWindow,
            VKLoopDirection
        >::type ArrayIndex;

        // flush the cache array at the given index
        context.BypassCache(TDataParameter::At(DataParameterOffset())) =
            cacheArray.At(static_cast<ArrayIndex*>(0));
    }
};

/**
* @struct KCacheStorageSlideFunctor
* Functor copying array element to the next position in the array
*/
struct KCacheStorageSlideFunctor
{
    template<
        typename TCacheArray, 
        typename TArrayIndex>
    __ACC__
    static void Do(typename parameter_type<TCacheArray>::type cacheArray)
    {
        // compute the array index of the next element
        typedef typename boost::mpl::next<TArrayIndex>::type NextArrayIndex;

        // slide the element to the next array index
        cacheArray.At(static_cast<NextArrayIndex*>(0)) = 
            cacheArray.At(static_cast<TArrayIndex*>(0));
    }
};

