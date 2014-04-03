#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/is_sequence.hpp>
#include "Cache.h"
#include "StencilSweepDescriptor.h"

#ifdef ENABLE_CACHING

/**
* Method defining a stencil sweep given a loop direction, caches and stencil stages
*/
template<
    KLoopDirection VKLoopDirection,
    typename TCaches,
    typename TStencilStages>
__CPU__
StencilSweepDescriptor<TCaches, TStencilStages, VKLoopDirection> define_sweep(TCaches, TStencilStages)
{
    StencilSweepDescriptor<TCaches, TStencilStages, VKLoopDirection> result;
    return result;
};

#else

/**
* Method defining a stencil sweep given a loop direction, caches and stencil stages
*/
template<
    KLoopDirection VKLoopDirection,
    typename TCaches,
    typename TStencilStages>
__CPU__
StencilSweepDescriptor<boost::mpl::void_, TStencilStages, VKLoopDirection> define_sweep(TCaches, TStencilStages)
{
    // check TCaches is either void or a list of caches
    // (note TCaches is not checked by the stencil sweep descriptor as void instead of TCaches is passed)
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

    // return a sweep descriptor with TCaches set to void
    StencilSweepDescriptor<boost::mpl::void_, TStencilStages, VKLoopDirection> result;
    return result;
};

#endif

/**
* Method defining a stencil sweep given a loop direction and stencil stages
*/
template<
    KLoopDirection VKLoopDirection,
    typename TStencilStages>
__CPU__
StencilSweepDescriptor<boost::mpl::void_, TStencilStages, VKLoopDirection> define_sweep(TStencilStages)
{
    StencilSweepDescriptor<boost::mpl::void_, TStencilStages, VKLoopDirection> result;
    return result;
};
