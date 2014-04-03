#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/mpl/erase_key.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include "IJRange.h"
#include "StencilStage.h"

#include <boost/mpl/map/map10.hpp> // depends on MAX_MERGE_COUNT
#include <boost/mpl/vector/vector30.hpp> // depends on MAX_MERGE_COUNT * 3 (slide, fill, and flush per stage / ij range)

/**
* @struct StencilSweepFunctorOperationSlide
* Structure defining the slide caches stencil sweep functor operation
*/
template<
    typename TSlideMask,
    typename TIJRange>
struct StencilSweepFunctorOperationSlide {};

/**
* @struct StencilSweepFunctorOperationFill
* Structure defining the fill caches stencil sweep functor operation
*/
template<
    typename TFillMask,
    typename TIJRange>
struct StencilSweepFunctorOperationFill {};

/**
* @struct StencilSweepFunctorOperationStage
* Structure defining the apply stencil stage stencil sweep functor operation
*/
template<typename TStencilStage>
struct StencilSweepFunctorOperationStage {};

/**
* @struct StencilSweepFunctorFlushCaches
* Structure defining the flush caches stencil sweep functor operation
*/
template<
    typename TFlushMask,
    typename TIJRange>
struct StencilSweepFunctorOperationFlush {};  

/**
* @struct stencil_sweep_functor_operation_ij_range
* Meta function computing the ij range of a stencil sweep functor operation
*/
template<typename TStencilSweepFunctorOperation>
struct stencil_sweep_functor_operation_ij_range;

template<
    typename TSlideMask,
    typename TIJRange>
struct stencil_sweep_functor_operation_ij_range<StencilSweepFunctorOperationSlide<TSlideMask, TIJRange> > 
{
    BOOST_STATIC_ASSERT(is_ij_range<TIJRange>::value);
    typedef TIJRange type;
};

template<
    typename TFillMask,
    typename TIJRange>
struct stencil_sweep_functor_operation_ij_range<StencilSweepFunctorOperationFill<TFillMask, TIJRange> > 
{
    BOOST_STATIC_ASSERT(is_ij_range<TIJRange>::value);
    typedef TIJRange type;
};

template<
    typename TFlushMask,
    typename TIJRange>
struct stencil_sweep_functor_operation_ij_range<StencilSweepFunctorOperationFlush<TFlushMask, TIJRange> > 
{
    BOOST_STATIC_ASSERT(is_ij_range<TIJRange>::value);
    typedef TIJRange type;
};

template<typename TStencilStage>
struct stencil_sweep_functor_operation_ij_range<StencilSweepFunctorOperationStage<TStencilStage> > 
{
    typedef typename stencil_stage_ij_range<TStencilStage>::type type;
};

/**
* @struct stencil_sweep_functor_operation_is_sync_needed
* Meta function returning true if the operation requires a synchronization after execution and false otherwise
*/
template<typename TStencilSweepFunctorOperation>
struct stencil_sweep_functor_operation_is_sync_needed : boost::mpl::false_ {};

// currently only stage apply operations require a synchronization
template<typename TStencilStage>
struct stencil_sweep_functor_operation_is_sync_needed<StencilSweepFunctorOperationStage<TStencilStage> > : boost::mpl::true_ {};

// meta functions used in order to create a map which associates an ij range with all caches updated in the given ij range

/**
* @struct cache_ij_range_map_insert
* Meta function inserting a cache into a cache ij range map
* (the cache ij range map associates an ij range with all caches updated in the given ij range)
*/
template<
    typename TCacheIJRangeMap,
    typename TCache,
    typename TStencilStages>
struct cache_ij_range_map_insert
{
    // compute the ij range of the given cache
    typedef typename stencil_stages_parameter_ij_range<
        TStencilStages, 
        typename cache_index<TCache>::type 
    >::type IJRange;

    // compute a vector of caches updated in the given ij range
    // (either extend an existing vector of caches or create a new vector containing the current cache)
    typedef typename boost::mpl::eval_if<
        boost::mpl::has_key<TCacheIJRangeMap, IJRange>,
        boost::mpl::push_back<typename boost::mpl::at<TCacheIJRangeMap, IJRange>::type, TCache>,
        boost::mpl::vector1<TCache>
    >::type Caches;

    // insert the cache into the map
    typedef typename boost::mpl::insert<
        typename boost::mpl::erase_key<TCacheIJRangeMap, IJRange>::type,
        boost::mpl::pair<IJRange, Caches>
    >::type type;
};

/**
* @struct create_cache_ij_range_map
* Meta function creating a cache ij range map given the cache usage map and all stencil stages of the sweep
* (the cache ij range map associates an ij range with all caches updated in the given ij range)
*/
template<
    typename TCacheUsageMap,   
    typename TStencilStages>
struct create_cache_ij_range_map
{
    typedef typename boost::mpl::fold<
        TCacheUsageMap,
        boost::mpl::map0<>,
        cache_ij_range_map_insert<boost::mpl::_1, boost::mpl::_2, TStencilStages>
    >::type type;
};

// meta functions used in order to compute a vector of stencil sweep functor operations for a given stencil stage

/**
* @struct stencil_sweep_functor_operations_add_slide
* Meta function adding all slide operations for a given ij range
*/
template<
    typename TStencilSweepFunctorOperations,
    typename TFirstCacheUsage,
    typename TIJRange>
struct stencil_sweep_functor_operations_add_slide
{
    // compute the cache slide mask
    typedef typename boost::mpl::fold<
        TFirstCacheUsage,
        boost::mpl::vector0<>,
        boost::mpl::if_<
            cache_is_slide_needed<boost::mpl::_2>,
            boost::mpl::push_back<boost::mpl::_1, cache_index<boost::mpl::_2> >,
            boost::mpl::_1
        >
    >::type SlideMask;

    // add the slide operations
    typedef typename boost::mpl::eval_if<
        boost::mpl::empty<SlideMask>,
        TStencilSweepFunctorOperations,
        boost::mpl::push_back<
            TStencilSweepFunctorOperations,
            StencilSweepFunctorOperationSlide<SlideMask, TIJRange>
        >
    >::type type;
};

/**
* @struct stencil_sweep_functor_operations_add_fill
* Meta function adding all fill operations for a given ij range
*/
template<
    typename TStencilSweepFunctorOperations,
    typename TFirstCacheUsage,
    typename TIJRange>
struct stencil_sweep_functor_operations_add_fill
{
    // compute the cache fill mask
    typedef typename boost::mpl::fold<
        TFirstCacheUsage,
        boost::mpl::vector0<>,
        boost::mpl::if_<
            cache_is_fill_needed<boost::mpl::_2>,
            boost::mpl::push_back<boost::mpl::_1, cache_index<boost::mpl::_2> >,
            boost::mpl::_1
        >
    >::type FillMask;

    // add the fill operations
    typedef typename boost::mpl::eval_if<
        boost::mpl::empty<FillMask>,
        TStencilSweepFunctorOperations,
        boost::mpl::push_back<
            TStencilSweepFunctorOperations,
            StencilSweepFunctorOperationFill<FillMask, TIJRange>
        >
    >::type type;
};

/**
* @struct stencil_sweep_functor_operations_add_flush
* Meta function adding all flush operations for a given ij range
*/
template<
    typename TStencilSweepFunctorOperations,
    typename TLastCacheUsage,
    typename TIJRange>
struct stencil_sweep_functor_operations_add_flush
{
    // compute the cache flush mask
    typedef typename boost::mpl::fold<
        TLastCacheUsage,
        boost::mpl::vector0<>,
        boost::mpl::if_<
            cache_is_flush_needed<boost::mpl::_2>,
            boost::mpl::push_back<boost::mpl::_1, cache_index<boost::mpl::_2> >,
            boost::mpl::_1
        >
    >::type FlushMask;

    // add the flush operations
    typedef typename boost::mpl::eval_if<
        boost::mpl::empty<FlushMask>,
        TStencilSweepFunctorOperations,
        boost::mpl::push_back<
            TStencilSweepFunctorOperations,
            StencilSweepFunctorOperationFlush<FlushMask, TIJRange>
        >
    >::type type;
};

/**
* @struct stencil_sweep_functor_operations_add_slide_and_fill
* Meta function adding all slide and fill operations for a given ij range
*/
template<
    typename TStencilSweepFunctorOperations,
    typename TFirstCacheUsage,
    typename TIJRange>
struct stencil_sweep_functor_operations_add_slide_and_fill
{
    typedef typename stencil_sweep_functor_operations_add_fill<
        typename stencil_sweep_functor_operations_add_slide<
            TStencilSweepFunctorOperations,
            TFirstCacheUsage,
            TIJRange
        >::type,
        TFirstCacheUsage,
        TIJRange
    >::type type;
};

/**
* @struct compute_stencil_stage_operations
* Meta function computing a vector containing all stencil sweep functor operations associated to stencil stage
*/
template<
    typename TStencilStage,
    typename TStencilStages,
    typename TFirstCacheUsage,
    typename TLastCacheUsage>
struct compute_stencil_sweep_functor_operations
{
    // define the stencil stage ij range
    typedef typename TStencilStage::IJRange IJRange;

    // compute maps with associate ij ranges with the caches updated in the given ij range
    typedef typename create_cache_ij_range_map<
        TFirstCacheUsage, 
        TStencilStages
    >::type FirstUsageCacheIJRangeMap;
    typedef typename create_cache_ij_range_map<
        TLastCacheUsage, 
        TStencilStages
    >::type LastUsageCacheIJRangeMap;

    // compute all slide and fill operations with an ij range different from the stencil stage ij range
    typedef typename boost::mpl::fold<
        typename boost::mpl::erase_key<FirstUsageCacheIJRangeMap, IJRange>::type,
        boost::mpl::vector0<>,
        stencil_sweep_functor_operations_add_slide_and_fill<
            boost::mpl::_1,
            boost::mpl::second<boost::mpl::_2>,
            boost::mpl::first<boost::mpl::_2>
        >
    >::type JaggedSlidesAndFills;

    // add slide an fill operations with an ij range equivalent to the stage ij range
    typedef typename boost::mpl::eval_if<
        boost::mpl::has_key<FirstUsageCacheIJRangeMap, IJRange>,
        stencil_sweep_functor_operations_add_slide_and_fill<
            JaggedSlidesAndFills,
            typename boost::mpl::at<FirstUsageCacheIJRangeMap, IJRange>::type,
            IJRange
        >,
        JaggedSlidesAndFills
    >::type AllSlidesAndFills;
  
    // add the apply stage operation
    typedef typename boost::mpl::push_back<
        AllSlidesAndFills,
        StencilSweepFunctorOperationStage<TStencilStage>
    >::type SlidesFillsAndApplyStage;

    // add all flush operations with an ij range equivalent to the stage ij range
    typedef typename boost::mpl::eval_if<
        boost::mpl::has_key<LastUsageCacheIJRangeMap, IJRange>,
        stencil_sweep_functor_operations_add_flush<
            SlidesFillsAndApplyStage,
            typename boost::mpl::at<LastUsageCacheIJRangeMap, IJRange>::type,
            IJRange
        >,
        SlidesFillsAndApplyStage
    >::type SlidesFillsApplyStageAndNonJaggedFlushes;

    // finally add all flush operations with an ij range different from the stencil stage ij range
    typedef typename boost::mpl::fold<
        typename boost::mpl::erase_key<LastUsageCacheIJRangeMap, IJRange>::type,
        SlidesFillsApplyStageAndNonJaggedFlushes,
        stencil_sweep_functor_operations_add_flush<
            boost::mpl::_1,
            boost::mpl::second<boost::mpl::_2>,
            boost::mpl::first<boost::mpl::_2>
        >
    >::type type;
};