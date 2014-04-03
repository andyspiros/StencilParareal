#pragma once

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/count_if.hpp>
#include <boost/mpl/less.hpp>
#include "Enums.h"
#include "Domains.h"
#include "Cache.h"
#include "StencilStage.h"

/**
* @struct StencilSweepDescriptor
* Structure defining a stencil sweep given the caches, stencil stages and the k loop direction
*/
template<
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection> 
struct StencilSweepDescriptor 
{
    // check TCaches is either void or a list of caches
    BOOST_STATIC_ASSERT( 
        (
            boost::mpl::eval_if<
                boost::mpl::is_sequence<TCaches>,
                boost::mpl::fold<
                    TCaches,
                    boost::mpl::true_,
                    boost::mpl::if_<
                        is_cache<boost::mpl::_2>,
                        boost::mpl::_1,
                        boost::mpl::false_
                    >
                >,
                boost::mpl::is_void_<TCaches>
            >::type::value
        )
    );
    
    // check TStencilStages is a sequence of stencil stages
    BOOST_STATIC_ASSERT( 
        (
            boost::mpl::eval_if<
                boost::mpl::is_sequence<TStencilStages>,
                boost::mpl::fold<
                    TStencilStages,
                    boost::mpl::true_,
                    boost::mpl::if_<
                        is_stencil_stage<boost::mpl::_2>,
                        boost::mpl::_1,
                        boost::mpl::false_
                    >
                >,
                boost::mpl::false_
            >::type::value
        )
    );
};

/**
* @struct is_stencil_sweep_descriptor
* Meta function returning true if the parameter is a stencil sweep descriptor
*/
template<typename T>
struct is_stencil_sweep_descriptor : boost::mpl::false_ {};

template<
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection>
struct is_stencil_sweep_descriptor<StencilSweepDescriptor<TCaches, TStencilStages, VKLoopDirection> > : boost::mpl::true_ {};

/**
* @struct stencil_sweep_descriptor_stages
* Meta function returning stencil stages of a stencil sweep descriptor
*/
template<typename T>
struct stencil_sweep_descriptor_stages;

template<
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection>
struct stencil_sweep_descriptor_stages<StencilSweepDescriptor<TCaches, TStencilStages, VKLoopDirection> >
{
    typedef TStencilStages type; 
};

/**
* @struct stencil_sweep_descriptor_uses_parameter
* Meta function returning true if a stencil sweep descriptor uses a given parameter
*/
template<
    typename TStencilSweepDescriptor, 
    typename TDomain,
    typename TParameterIndex>
struct stencil_sweep_descriptor_uses_parameter : 
    stencil_stages_have_parameter<
        typename stencil_sweep_descriptor_stages<TStencilSweepDescriptor>::type, 
        TDomain, 
        TParameterIndex
    > 
{};

/**
* @struct stencil_sweep_descriptor_maximum_ij_range
* Meta function returning the maximal ij range of all sweep stencil stages
*/
template<typename TStencilSweepDescriptor>
struct stencil_sweep_descriptor_maximum_ij_range : 
    stencil_stages_maximum_ij_range<
        typename stencil_sweep_descriptor_stages<TStencilSweepDescriptor>::type
    > 
{};

/**
* @struct stencil_sweep_descriptor_parameter_ij_range
* Meta function returning the maximal ij range of all sweep stencil stages using a given parameter
*/
template<
    typename TStencilSweepDescriptor, 
    typename TParameterIndex>
struct stencil_sweep_descriptor_parameter_ij_range : 
    stencil_stages_parameter_ij_range<
        typename stencil_sweep_descriptor_stages<TStencilSweepDescriptor>::type, 
        TParameterIndex
    > 
{};

/**
* @struct stencil_sweep_descriptor_is_parameter_usage_unique
* Meta function that tests if a parameter is used in more than one stage of a stencil sweep in which
*   - returns mpl::false_ if the variable is used in more than one stencil stage
*   - returns mpl::true otherwise
*/
template<
    typename TStencilSweepDescriptor,
    typename TParameterIndex>
struct stencil_sweep_descriptor_is_parameter_usage_unique
{
    // check if usage count is less than 2
    typedef typename boost::mpl::less<
        typename boost::mpl::count_if< // count how many stages use the parameter
            typename stencil_sweep_descriptor_stages<TStencilSweepDescriptor>::type, // get the sequence of stencil stages
            stencil_stage_has_parameter<boost::mpl::_, FullDomain, TParameterIndex>
        >::type,
        boost::mpl::int_<2>
    >::type type;
};

/**
* @struct stencil_sweep_descriptors_is_parameter_usage_unique
* Meta function returning false if a parameter is used by more than one stage of a stencil sweep
*/
template<
    typename TStencilSweepDescriptors,
    typename TParameterIndex>
struct stencil_sweep_descriptors_is_parameter_usage_unique
{
    typedef typename boost::mpl::fold<
        TStencilSweepDescriptors,
        boost::mpl::true_,
        boost::mpl::if_<
            stencil_sweep_descriptor_is_parameter_usage_unique<
                boost::mpl::_2,
                TParameterIndex
            >,
            boost::mpl::_1,
            boost::mpl::false_
        >
    >::type type;
};

/**
* @struct stencil_sweep_descriptor_is_parameter_cached_locally
* Meta function returning true if a given parameter is cached locally without cache fills or flushes
*/
template<
    typename TStencilSweepDescriptor,
    typename TParameterIndex>
struct stencil_sweep_descriptor_is_parameter_cached_locally;

template<
    typename TParameterIndex,
    typename TStencilStages,
    KLoopDirection VKLoopDirection>
struct stencil_sweep_descriptor_is_parameter_cached_locally<
    StencilSweepDescriptor<boost::mpl::void_, TStencilStages, VKLoopDirection>,
    TParameterIndex> : 
    boost::mpl::false_ 
{};

template<
    typename TParameterIndex,
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection>
struct stencil_sweep_descriptor_is_parameter_cached_locally<
    StencilSweepDescriptor<TCaches, TStencilStages, VKLoopDirection>,
    TParameterIndex> 
{
    // find if the given parameter is cached
    typedef typename boost::mpl::find_if<
        TCaches, 
        boost::is_same<cache_index<boost::mpl::_>, TParameterIndex> 
    >::type CacheIterator;

    // compute the cache value respectively void if no cache was found
    typedef typename boost::mpl::eval_if<
        boost::is_same<CacheIterator, typename boost::mpl::end<TCaches>::type>,
        boost::mpl::void_,
        boost::mpl::deref<CacheIterator>
    >::type Cache;
    
    // return true if the parameter is cached locally
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            boost::mpl::eval_if<
                boost::mpl::is_void_<Cache>,
                boost::mpl::false_,
                cache_is_local<Cache>
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<bool, bool(value)> type;
};

/**
* @struct stencil_sweep_descriptor_is_parameter_accessed
* Meta function returning true if a given parameter is accessed not counting local cache accesses
* (note that the meta function is used to decide if a temporary field needs to be allocated,
* therefore local cache accesses are not considered as the underlying temporary field is not accessed)
*/
template<
    typename TStencilSweepDescriptor,
    typename TDomain,
    typename TParameterIndex>
struct stencil_sweep_descriptor_is_parameter_accessed
{
    // check if the parameter is used by the sweep and if yes not cached locally
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            stencil_stages_have_parameter<
                typename stencil_sweep_descriptor_stages<TStencilSweepDescriptor>::type, 
                TDomain, 
                TParameterIndex
            >::value &&
            !stencil_sweep_descriptor_is_parameter_cached_locally<
                TStencilSweepDescriptor,
                TParameterIndex
            >::value
        )
    );
    typedef boost::mpl::integral_c<bool, bool(value)> type;
};

/**
* @struct stencil_sweep_descriptors_is_parameter_accessed
* Meta function returning true if a given parameter is accessed not counting local cache accesses
* (note that the meta function is used to decide if a temporary field needs to be allocated,
* therefore local cache accesses are not considered as the underlying temporary field is not accessed)
*/
template<
    typename TStencilSweepDescriptors,
    typename TDomain,
    typename TParameterIndex>
struct stencil_sweep_descriptors_is_parameter_accessed
{
    // iterate over all sweeps and search of a parameter access
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            boost::mpl::fold<
                TStencilSweepDescriptors,
                boost::mpl::false_,
                boost::mpl::if_<
                    stencil_sweep_descriptor_is_parameter_accessed<boost::mpl::_2, TDomain, TParameterIndex>,
                    boost::mpl::true_,
                    boost::mpl::_1
                >
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<bool, bool(value)> type;
};
