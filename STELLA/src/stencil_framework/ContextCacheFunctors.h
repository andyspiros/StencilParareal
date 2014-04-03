#pragma once

#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/pair.hpp>
#include "SharedInfrastructure.h"
#include "StencilStage.h"
#include "CacheStorage.h"

#include <boost/mpl/vector/vector20.hpp> // depends on MAX_CACHE_COUNT

/**
* @struct k_cache_elements
* Meta function computing the k cache tuple element type, used to setup the k cache tuple in the CUDA context
*/
template<
    typename TCaches,
    typename TIndexToElementMap,
    KLoopDirection VKLoopDirection>
struct k_cache_elements
{
    // compute index and element types
    typedef typename boost::mpl::transform<
        TCaches,
        cache_index<boost::mpl::_>
    >::type KCacheIndexes;
    typedef typename boost::mpl::transform<
        TCaches,
        create_k_cache_strorage<
            boost::mpl::_, 
            boost::mpl::at<TIndexToElementMap, cache_index<boost::mpl::_> >,
            boost::mpl::integral_c<KLoopDirection, VKLoopDirection>
        >
    >::type KCacheElementTypes;

    // define the tuple element type
    typedef TupleElements<
        KCacheIndexes, 
        KCacheElementTypes
    > type;
};

/**
* @struct ij_cache_elements
* Meta function computing the ij cache tuple element type, used to setup the ij cache tuple in the CUDA context
*/
template<
    typename TCaches,
    typename TIndexToElementMap,
    typename TStencilStages,
    typename TBlockSize>
struct ij_cache_elements
{
    // compute index and element types
    typedef typename boost::mpl::transform<
        TCaches,
        cache_index<boost::mpl::_>
    >::type IJCacheIndexes;
    typedef typename boost::mpl::transform<
        TCaches,
        create_ij_cache_strorage<
            boost::mpl::_, 
            boost::mpl::at<TIndexToElementMap, cache_index<boost::mpl::_> >,
            stencil_stages_parameter_ij_range<TStencilStages, cache_index<boost::mpl::_> >,
            TBlockSize
        >
    >::type IJCacheElementTypes;

    // define the tuple element type
    typedef TupleElements<
        IJCacheIndexes, 
        IJCacheElementTypes
    > type;
};

/**
* @struct index_to_cache_element_map_add
* Meta function adding all elements of a cache tuple to an index to cache element map
*/
template<
    typename TIndexToCacheElementMap,
    typename TCacheElements>
struct index_to_cache_element_map_add
{
    typedef typename boost::mpl::fold<
        boost::mpl::range_c<int, 0, boost::mpl::size<typename TCacheElements::ElementIndexes>::value>,
        TIndexToCacheElementMap,
        boost::mpl::insert<
            boost::mpl::_1,
            boost::mpl::pair<
                boost::mpl::at<typename TCacheElements::ElementIndexes, boost::mpl::_2>,
                boost::mpl::at<typename TCacheElements::ElementTypes, boost::mpl::_2>
            >
        >
    >::type type;
};

/**
* @struct create_cache_element_map
* Meta function creating a map for all elements cached in the register cache
*/
template<
    typename TKCacheElements,
    typename TIJCacheElements>
struct create_index_to_cache_element_map
{
    // sequentially add the k and the ij cache elements to the result map
    typedef typename index_to_cache_element_map_add<
        typename index_to_cache_element_map_add<
            boost::mpl::map0<>,
            TKCacheElements
        >::type,
        TIJCacheElements
    >::type type;
};

/**
* @struct ContextCacheSlideFunctor
* Functor sliding an individual k cache storage
*/
struct ContextCacheSlideFunctor
{
    template<
        typename TValue,
        typename TKWindow,
        typename TKRange, 
        int VParameterIndex, 
        CacheIOPolicy VIOPolicy,
        KLoopDirection VKLoopDirection>
    __ACC__
    static void Do(KCacheStorage<TValue, TKWindow, TKRange, VParameterIndex, VIOPolicy, VKLoopDirection>& cacheStorage)
    {
        cacheStorage.Slide();
    }
};

/**
* @struct ContextCacheFillFrontFunctor
* Functor filling the front of a k cache storage
*/
template<
    typename TContext,
    typename TDomain>
struct ContextCacheFillFrontFunctor
{
    template<
        typename TValue,
        typename TKWindow,
        typename TKRange, 
        int VParameterIndex, 
        CacheIOPolicy VIOPolicy,
        KLoopDirection VKLoopDirection>
    __ACC__
    static void Do(
        KCacheStorage<TValue, TKWindow, TKRange, VParameterIndex, VIOPolicy, VKLoopDirection>& cacheStorage, 
        typename parameter_type<TContext>::type context)
    {
        cacheStorage.template FillFront<TDomain>(context);
    }
};

/**
* @struct ContextCacheFillFunctor
* Functor filling a k cache storage
*/
template<
    typename TContext,
    typename TDomain>
struct ContextCacheFillFunctor
{
    template<
        typename TValue,
        typename TKWindow,
        typename TKRange, 
        int VParameterIndex, 
        CacheIOPolicy VIOPolicy,
        KLoopDirection VKLoopDirection>
    __ACC__
    static void Do(
        KCacheStorage<TValue, TKWindow, TKRange, VParameterIndex, VIOPolicy, VKLoopDirection>& cacheStorage, 
        typename parameter_type<TContext>::type context)
    {
        cacheStorage.template Fill<TDomain>(context);
    }
};

/**
* @struct ContextCacheFlushBackFunctor
* Functor flushing the back of a k cache storage
*/
template<
    typename TContext,
    typename TDomain>
struct ContextCacheFlushBackFunctor
{
    template<
        typename TValue,
        typename TKWindow,
        typename TKRange, 
        int VParameterIndex, 
        CacheIOPolicy VIOPolicy,
        KLoopDirection VKLoopDirection>
    __ACC__
    static void Do(
        KCacheStorage<TValue, TKWindow, TKRange, VParameterIndex, VIOPolicy, VKLoopDirection>& cacheStorage, 
        typename parameter_type<TContext>::type context)
    {
        cacheStorage.template FlushBack<TDomain>(context);
    }
};

/**
* @struct ContextCacheFlushFunctor
* Functor flushing a k cache storage
*/
template<
    typename TContext,
    typename TDomain>
struct ContextCacheFlushFunctor
{
    template<
        typename TValue,
        typename TKWindow,
        typename TKRange, 
        int VParameterIndex, 
        CacheIOPolicy VIOPolicy,
        KLoopDirection VKLoopDirection>
    __ACC__
    static void Do(
        KCacheStorage<TValue, TKWindow, TKRange, VParameterIndex, VIOPolicy, VKLoopDirection>& cacheStorage, 
        typename parameter_type<TContext>::type context)
    {
        cacheStorage.template Flush<TDomain>(context);
    }
};
