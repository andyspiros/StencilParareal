#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/pop_back.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/mpl/erase_key.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include "IJRange.h"
#include "KRange.h"
#include "StencilStage.h"
#include "StencilSweepFunctorOperation.h"

#include <boost/mpl/map/map10.hpp> // depends on MAX_MERGE_COUNT 
#include <boost/mpl/set/set10.hpp> // depends on MAX_MERGE_COUNT
#include <boost/mpl/vector/vector30.hpp> // depends on MAX_MERGE_COUNT * 3 (slide, fill, and flush per stage / ij range)

/**
* @struct StencilSweepFunctorOperationGroup
* Structure holding a group of stencil sweep functor operations
* (note that all operations of a group are applied with the same ij range)
*/
template<
    typename TStencilSweepFunctorOperations,
    typename TIJRange>
struct StencilSweepFunctorOperationGroup
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TStencilSweepFunctorOperations>::value);
    BOOST_STATIC_ASSERT(is_ij_range<TIJRange>::value);

    typedef TStencilSweepFunctorOperations StencilSweepFunctorOperations;
    typedef TIJRange IJRange;
}; 

// meta functions used in order to acces stencil sweep functor operation group members

/**
* @struct stencil_sweep_functor_operation_group_ij_range
* Meta function returning the ij range of a stencil sweep functor operation group
*/
template<typename TStencilSweepFunctorOperationGroup>
struct stencil_sweep_functor_operation_group_ij_range;

template<
    typename TStencilSweepFunctorOperations,
    typename TIJRange>
struct stencil_sweep_functor_operation_group_ij_range<StencilSweepFunctorOperationGroup<TStencilSweepFunctorOperations, TIJRange> >
{
    BOOST_STATIC_ASSERT(is_ij_range<TIJRange>::value);
    typedef TIJRange type;
};

/**
* @struct stencil_sweep_functor_operation_group_operations
* Meta function returning the operations of a stencil sweep functor operation group
*/
template<typename TStencilSweepFunctorOperationGroup>
struct stencil_sweep_functor_operation_group_operations;

template<
    typename TStencilSweepFunctorOperations,
    typename TIJRange>
struct stencil_sweep_functor_operation_group_operations<StencilSweepFunctorOperationGroup<TStencilSweepFunctorOperations, TIJRange> >
{
    typedef TStencilSweepFunctorOperations type;
};

// meta functions used in order to create a map which associates a stencil stage with the caches it uses first or last

/**
* @struct cache_usage_map_insert
* Meta function inserting a cache into a cache usage map
* (the cache usage map associates a stencil stage with a vector of caches used first or last by the stencil stage)
*/
template<
    typename TCacheUsageMap,
    typename TStencilStage,
    typename TCache>
struct cache_usage_map_insert
{
    // define the vector of caches which shall be added to the map
    // (either extend an existing vector of caches or create a new vector containing the current cache)
    typedef typename boost::mpl::eval_if<
        boost::mpl::has_key<TCacheUsageMap, TStencilStage>,
        boost::mpl::push_back<typename boost::mpl::at<TCacheUsageMap, TStencilStage>::type, TCache>,
        boost::mpl::vector1<TCache>
    >::type Caches;
    
    // insert the new entry into the map
    typedef typename boost::mpl::insert<
        typename boost::mpl::erase_key<TCacheUsageMap, TStencilStage>::type,
        boost::mpl::pair<TStencilStage, Caches>
    >::type type;
};

/**
* @struct cache_usage_map_at
* Meta function returning all caches used by a given stencil stage
* (return an empty vector if there is no map entry for the given stencil stage)
*/
template<
    typename TCacheUsageMap,
    typename TStencilStage>
struct cache_usage_map_at
{
    typedef typename boost::mpl::eval_if<
        boost::mpl::has_key<TCacheUsageMap, TStencilStage>,
        boost::mpl::at<TCacheUsageMap, TStencilStage>,
        boost::mpl::vector0<>
    >::type type;
};

/**
* @struct create_cache_first_usage_map
* Meta function creating a cache usage map containing the first usages
* (the cache usage map associates a stencil stage with a vector of caches used first or last by the stencil stage)
*/
template<
    typename TStencilStages,
    typename TCaches>
struct create_cache_first_usage_map
{
	// iterate over all caches and add them to the map
    typedef typename boost::mpl::fold<
        TCaches,
        boost::mpl::map0<>,
        cache_usage_map_insert<
            boost::mpl::_1,
            stencil_stages_first_parameter_usage<TStencilStages, cache_index<boost::mpl::_2> >,
            boost::mpl::_2
        >
    >::type type;
};

/**
* @struct create_cache_last_usage_map
* Meta function creating a cache usage map containing the last usages
* (the cache usage map associates a stencil stage with a vector of caches used first or last by the stencil stage)
*/
template<
    typename TStencilStages,
    typename TCaches>
struct create_cache_last_usage_map
{
	// iterate over all caches and add them to the map
    typedef typename boost::mpl::fold<
        TCaches,
        boost::mpl::map0<>,
        cache_usage_map_insert<
            boost::mpl::_1,
            stencil_stages_last_parameter_usage<TStencilStages, cache_index<boost::mpl::_2> >,
            boost::mpl::_2
        >
    >::type type;
};

// meta functions used in order to compute a vector of stencil sweep functor operation groups

/**
* @struct create_stencil_sweep_functor_operation_group
* Meta function creating a new stencil sweep operation group given an operation
*/
template<typename TStencilSweepFunctorOperation>
struct create_stencil_sweep_functor_operation_group
{
    typedef StencilSweepFunctorOperationGroup<
        boost::mpl::vector1<TStencilSweepFunctorOperation>,
        typename stencil_sweep_functor_operation_ij_range<TStencilSweepFunctorOperation>::type
    > type;
};

/**
* @struct stencil_sweep_functor_operation_group_add_operation
* Meta function adding an operation to a stencil sweep functor operation group
*/
template<
    typename TStencilSweepFunctorOperationGroup,
    typename TStencilSweepFunctorOperation>
struct stencil_sweep_functor_operation_group_add_operation
{
    // add the operation at the back of the operation group
    typedef StencilSweepFunctorOperationGroup<
        typename boost::mpl::push_back<
            typename stencil_sweep_functor_operation_group_operations<TStencilSweepFunctorOperationGroup>::type,
            TStencilSweepFunctorOperation
        >::type,
        typename stencil_sweep_functor_operation_group_ij_range<TStencilSweepFunctorOperationGroup>::type
    > type;
};

/**
* @struct stencil_sweep_functor_operation_groups_add_operation
* Meta function adding an operation to a vector of stencil sweep functor operations
* (either the operation is added to the last group in case the ij ranges match or a new operation group is created)
*/
template<
    typename TStencilSweepFunctorOperationGroups,
    typename TStencilSweepFunctorOperations>
struct stencil_sweep_functor_operation_groups_add_operations
{
    // add a new operation to the stencil sweep functor operation groups vector
    typedef typename boost::mpl::fold<
        TStencilSweepFunctorOperations,
        TStencilSweepFunctorOperationGroups,
        boost::mpl::if_<
            boost::mpl::empty<boost::mpl::_1>, 
            // in case there is no group insert a new one
            boost::mpl::push_back<
                boost::mpl::_1,
                create_stencil_sweep_functor_operation_group<boost::mpl::_2>
            >,
            // check if the ij range of the last group matches the operation ij range
            boost::mpl::if_< 
                boost::is_same<
                    stencil_sweep_functor_operation_group_ij_range<boost::mpl::back<boost::mpl::_1> >,
                    stencil_sweep_functor_operation_ij_range<boost::mpl::_2>
                >, 
				// insert the operation into the last group
                boost::mpl::push_back<
                    boost::mpl::pop_back<boost::mpl::_1>,
                    stencil_sweep_functor_operation_group_add_operation<boost::mpl::back<boost::mpl::_1>, boost::mpl::_2>
                >,
				// create a new group and insert the operation as first element
                boost::mpl::push_back<
                    boost::mpl::_1,
                    create_stencil_sweep_functor_operation_group<boost::mpl::_2>
                >
            >
        >
    >::type type;
};

/**
* @struct compute_stencil_sweep_functor_operation_groups
* Meta function computing a vector of sweep functor operation groups given the stages and caches which need to be updated by the sweep
*/
template<
    typename TStencilStages,
    typename TCaches,
    typename TDomain>
struct compute_stencil_sweep_functor_operation_groups
{
    // filter out all stencil stages which are not active in the current domain
    // (we do currently not filter the caches as the expected performance gain is negligible)
    typedef typename boost::mpl::copy_if<
        TStencilStages,
        in_k_range<stencil_stage_k_range<boost::mpl::_>, TDomain>
    >::type ActiveStencilStages;
    
    // compute a set holding all stencil stages of a sweep
    typedef typename boost::mpl::fold<
        TStencilStages,
        boost::mpl::set0<>,
        boost::mpl::insert<boost::mpl::_1, boost::mpl::_2>
    >::type StencilStageSet;

    // compare the stencil stage set size with the stencil stage vector size
    // (if the sizes match the stencil stages are unique and can be used as key element)
    BOOST_MPL_ASSERT_MSG( 
        (
			boost::mpl::size<StencilStageSet>::value == 
            boost::mpl::size<TStencilStages>::value
		),
        STENCIL_SWEEP_FUNCTOR_STENCIL_STAGE_DEFINITIONS_NOT_UNIQUE,
        (TStencilStages) 
    );

    // define two maps which associate the caches with the stencil stage using the caches first or last
    typedef typename create_cache_first_usage_map<ActiveStencilStages, TCaches>::type FirstCacheUsage;
    typedef typename create_cache_last_usage_map<ActiveStencilStages, TCaches>::type LastCacheUsage;

    // iterate over all active stencil stages and add the associated operations
    typedef typename boost::mpl::fold<
        ActiveStencilStages,
        boost::mpl::vector0<>,
        // add the operations associated to a single stage into the final vector of operations
        // (merge successive operations with identical ij range)
        stencil_sweep_functor_operation_groups_add_operations<
            boost::mpl::_1,
            // compute a vector of operations associated to the given stencil stage
            compute_stencil_sweep_functor_operations<
                boost::mpl::_2, 
                TStencilStages,
                cache_usage_map_at<FirstCacheUsage, boost::mpl::_2>,
                cache_usage_map_at<LastCacheUsage, boost::mpl::_2>
            >
        >
    >::type type;
};

