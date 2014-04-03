#pragma once

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include "KRange.h"
#include "Domains.h"
#include "KWindow.h"

/**
* @struct KCache
* Structure holding information necessary to cache a data field in k direction
*/
template<
    int VParameterIndex,
    CacheIOPolicy VIOPolicy,
    typename TKWindow,
    typename TKRange>
struct KCache
{
    BOOST_STATIC_ASSERT(is_k_window<TKWindow>::value);
    BOOST_STATIC_ASSERT(is_k_range<TKRange>::value);

    // make sure no caches for a single k level are defined
    BOOST_STATIC_ASSERT(is_3d_domain<typename TKRange::Domain>::value);

    typedef boost::mpl::integral_c<int, VParameterIndex> ParameterIndex;
    typedef boost::mpl::integral_c< ::CacheIOPolicy, VIOPolicy> CacheIOPolicy;
    typedef TKRange KRange;
    typedef TKWindow KWindow;
};

/**
* @struct is_k_cache
* Meta function returning true if the parameter is a k cache
*/
template<typename T>
struct is_k_cache : boost::mpl::false_ {};

template<
    int VParameterIndex,
    CacheIOPolicy VIOPolicy,
    typename TKWindow,
    typename TKRange>
struct is_k_cache<KCache<VParameterIndex, VIOPolicy, TKWindow, TKRange> > : boost::mpl::true_ {};

/**
* @struct IJCache
* Structure holding information necessary to cache a data field in ij direction
*/
template<
    int VParameterIndex,
    typename TKRange>
struct IJCache
{
    BOOST_STATIC_ASSERT(is_k_range<TKRange>::value);

    // make sure no caches for a single k level are defined
    BOOST_STATIC_ASSERT(is_3d_domain<typename TKRange::Domain>::value);

    typedef boost::mpl::integral_c<int, VParameterIndex> ParameterIndex;
    typedef boost::mpl::integral_c< ::CacheIOPolicy, cLocal> CacheIOPolicy;
    typedef TKRange KRange;
    typedef ::KWindow<0,0> KWindow;
};

/**
* @struct is_ij_cache
* Meta function returning true if the parameter is an ij cache
*/
template<typename T>
struct is_ij_cache : boost::mpl::false_ {};

template<
    int VParameterIndex,
    typename TKRange>
struct is_ij_cache<IJCache<VParameterIndex, TKRange> > : boost::mpl::true_ {};

/**
* @struct is_cache
* Meta function returning true if the parameter is a k or ij cache
*/
template<typename T>
struct is_cache : boost::mpl::or_<is_k_cache<T>, is_ij_cache<T> > {};

/**
* @struct cache_index
* Meta function computing the index of a cache
*/
template<typename TCache>
struct cache_index;

template<
    int VParameterIndex,
    CacheIOPolicy VIOPolicy,
    typename TKWindow,
    typename TKRange>
struct cache_index<KCache<VParameterIndex, VIOPolicy, TKWindow, TKRange> > : boost::mpl::integral_c<int, VParameterIndex> {};

template<
    int VParameterIndex,
    typename TKRange>
struct cache_index<IJCache<VParameterIndex, TKRange> > : boost::mpl::integral_c<int, VParameterIndex> {};

/**
* @struct cache_k_range
* Meta function returning the cache k range
*/
template<typename TCache>
struct cache_k_range;

template<    
    int VParameterIndex,
    CacheIOPolicy VIOPolicy,
    typename TKWindow,
    typename TKRange>
struct cache_k_range<KCache<VParameterIndex, VIOPolicy, TKWindow, TKRange> > 
{
    typedef TKRange type;
};

template<    
    int VParameterIndex,
    typename TKRange>
struct cache_k_range<IJCache<VParameterIndex, TKRange> > 
{
    typedef TKRange type;
};

/**
* @struct cache_is_flush_needed
* Meta function returning true if a cache flush is needed
*/
template<typename TCache>
struct cache_is_flush_needed;

template<
    int VParameterIndex,
    CacheIOPolicy VIOPolicy,
    typename TKWindow,
    typename TKRange>
struct cache_is_flush_needed<KCache<VParameterIndex, VIOPolicy, TKWindow, TKRange> > :
    boost::mpl::bool_<VIOPolicy == cFlush || VIOPolicy == cFillAndFlush>
{};

// currently flush is not supported for ij caches
template<
    int VParameterIndex,
    typename TKRange>
struct cache_is_flush_needed<IJCache<VParameterIndex, TKRange> > : boost::mpl::false_ {};

/**
* @struct cache_is_fill_needed
* Meta function returning true if a cache fill is needed
*/
template<typename TCache>
struct cache_is_fill_needed;

template<
    int VParameterIndex,
    CacheIOPolicy VIOPolicy,
    typename TKWindow,
    typename TKRange>
struct cache_is_fill_needed<KCache<VParameterIndex, VIOPolicy, TKWindow, TKRange> > :
    boost::mpl::bool_<VIOPolicy == cFill || VIOPolicy == cFillAndFlush>
{};

// currently fill is not supported for ij caches
template<
    int VParameterIndex,
    typename TKRange>
struct cache_is_fill_needed<IJCache<VParameterIndex, TKRange> > : boost::mpl::false_ {};

/**
* @struct cache_is_slide_needed
* Meta function returning true if a cache slide is needed
*/
template<typename TCache>
struct cache_is_slide_needed;

template<
    int VParameterIndex,
    CacheIOPolicy VIOPolicy,
    int VKMinusOffset,
    int VKPlusOffset,
    typename TKRange>
struct cache_is_slide_needed<KCache<VParameterIndex, VIOPolicy, KWindow<VKMinusOffset, VKPlusOffset>, TKRange> > : 
    boost::mpl::bool_<VKMinusOffset != VKPlusOffset> 
{};

// slides are not needed for single k level caches
template<
    int VParameterIndex,
    typename TKRange>
struct cache_is_slide_needed<IJCache<VParameterIndex, TKRange> > : boost::mpl::false_ {};

/**
* @struct cache_is_local
* Meta function returning true if a cache is local
* (note a cache is local if no fill or flush is performed)
*/
template<typename TCache>
struct cache_is_local :
    boost::mpl::not_<
        typename boost::mpl::or_<
            typename cache_is_fill_needed<TCache>::type,
            typename cache_is_flush_needed<TCache>::type
        >::type
    >
{};

/**
* @struct cache_full_update_k_range_minimum_offset
* Meta function computing the full update k range minimum offset 
* (above this offset all cache operations as fill and flush are performed regularly while slide is not considered)
*/
template<
    typename TCache,
    typename TKLoopDirection>
struct cache_full_update_k_range_minimum_offset
{
    // return the maximum fill or flush index
    typedef typename boost::mpl::if_c<
        (cache_is_fill_needed<TCache>::value && TKLoopDirection::value == cKIncrement) ||
        (cache_is_flush_needed<TCache>::value && TKLoopDirection::value == cKDecrement),
        boost::mpl::integral_c<int, TCache::KWindow::KPlusOffset::value>,
        boost::mpl::void_ 
    >::type MaximumFillOrFlushOffset;  
    
    // return the minimum fill or flush index
    typedef typename boost::mpl::if_c<
        (cache_is_flush_needed<TCache>::value && TKLoopDirection::value == cKIncrement) ||
        (cache_is_fill_needed<TCache>::value && TKLoopDirection::value == cKDecrement),
        boost::mpl::integral_c<int, TCache::KWindow::KMinusOffset::value>,
        MaximumFillOrFlushOffset
    >::type MinimumFillOrFlushOffset;

    // if MinimumFillOrFlushOffset is not void compute the absolute operation offset otherwise return the max boundary level offset
    // (note in case any operation is performed MinimumFillOrFlushOffset is not void)
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::eval_if<
                boost::mpl::is_void_<MinimumFillOrFlushOffset>,
                boost::mpl::integral_c<int, -MAX_BOUNDARY_LEVELS>,
                boost::mpl::minus<
                    typename TCache::KRange::KMinimumOffset,
                    MinimumFillOrFlushOffset            
                >
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct cache_full_update_k_range_maxium_offset
* Meta function computing the full update k range maximum offset 
* (below this offset all cache operations as fill and flush are performed regularly while slide is not considered)
*/
template<
    typename TCache,
    typename TKLoopDirection>
struct cache_full_update_k_range_maximum_offset
{
    // return the minimum fill or flush index
    typedef typename boost::mpl::if_c<
        (cache_is_flush_needed<TCache>::value && TKLoopDirection::value == cKIncrement) ||
        (cache_is_fill_needed<TCache>::value && TKLoopDirection::value == cKDecrement),
        boost::mpl::integral_c<int, TCache::KWindow::KMinusOffset::value>,
        boost::mpl::void_
    >::type MinimumFillOrFlushOffset;  

    // return the maximum fill or flush index
    typedef typename boost::mpl::if_c<
        (cache_is_fill_needed<TCache>::value && TKLoopDirection::value == cKIncrement) ||
        (cache_is_flush_needed<TCache>::value && TKLoopDirection::value == cKDecrement),
        boost::mpl::integral_c<int, TCache::KWindow::KPlusOffset::value>,
        MinimumFillOrFlushOffset
    >::type MaximumFillOrFlushOffset;

    // if MaximumFillOrFlushOffset is not void compute the absolute operation offset otherwise return the max boundary level offset
    // (note in case any operation is performed MaximumFillOrFlushOffset is not void)
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::eval_if<
                boost::mpl::is_void_<MaximumFillOrFlushOffset>,
                boost::mpl::integral_c<int, MAX_BOUNDARY_LEVELS>, 
                boost::mpl::minus<
                    typename TCache::KRange::KMaximumOffset,
                    MaximumFillOrFlushOffset            
                >
            >::type::value 
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct cache_full_update_k_range
* Meta function returning a k range in which all cache fill and flush operations are performed regularly
* (note outside this range it is k loop boundary level specific which cache operations need to be performed)
*/
template<
    typename TCache,
    typename TKLoopDirection>
struct cache_full_update_k_range
{
    typedef KRange<
        typename TCache::KRange::Domain,        
        cache_full_update_k_range_minimum_offset<TCache, TKLoopDirection>::value,
        cache_full_update_k_range_maximum_offset<TCache, TKLoopDirection>::value
    > type;
};