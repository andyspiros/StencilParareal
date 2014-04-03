#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/pop_back.hpp>
#include "Definitions.h"
#include "Enums.h"
#include "KPosition.h"
#include "KRange.h"
#include "ApplyToAll.h"
#include "KLoopBodyFunctors.h"

/**
* @struct KLoop2DBody
* Structure implementing the k loop body for 2d domains
*/
template<
    typename TDomain,
    bool VFillOrFlushCaches>
struct KLoop2DBody
{
    BOOST_STATIC_ASSERT(is_2d_domain<TDomain>::value);

    template<
        typename TContext,
        typename TStencilSweepFunctor>
    __ACC__
    static void Do(TContext& context)
    {
        TStencilSweepFunctor::template Do<TDomain, VFillOrFlushCaches, VFillOrFlushCaches>(context);
    }
};

/**
* @struct KLoop3DBody
* Structure implementing the k loop body for 3d domains
*/
template<
    typename TDomain,
    typename TIterationMask,
    typename TBeginKPosition,
    int VBeginKBoundarySize,
    typename TEndKPosition,
    int VEndKBoundarySize,
    bool VFillOrFlushCaches,
    KLoopDirection VKLoopDirection>
struct KLoop3DBody;

// specialization for k loops which support caching
// (there is a specific )
template<
    typename TDomain,
    typename TIterationMask,
    typename TBeginKPosition,
    int VBeginKBoundarySize,
    typename TEndKPosition,
    int VEndKBoundarySize,
    KLoopDirection VKLoopDirection>
struct KLoop3DBody<
    TDomain, 
    TIterationMask, 
    TBeginKPosition, 
    VBeginKBoundarySize, 
    TEndKPosition, 
    VEndKBoundarySize, 
    true, // fill or flush caches
    VKLoopDirection>
{
    BOOST_STATIC_ASSERT(is_domain<TDomain>::value);
    BOOST_STATIC_ASSERT(is_k_position<TBeginKPosition>::value);
    BOOST_STATIC_ASSERT(is_k_position<TEndKPosition>::value);

    template<
        typename TContext,
        typename TStencilSweepFunctor>
    __ACC__
    static void Do(TContext& context)
    {
        // compute a type list containing all boundary domains
        // (make sure there is at least one boundary domain used for cache initialization)
        typedef typename compute_begin_boundary_domains<
            TDomain,
            VKLoopDirection,
            TBeginKPosition::KOffset::value,
            (VBeginKBoundarySize > 1 ? VBeginKBoundarySize : 1)
        >::type BeginBoundaryDomains;

        // execute the first boundary level
        TStencilSweepFunctor::template Do<typename boost::mpl::front<BeginBoundaryDomains>::type, true, false>(context);
        context.template Advance<TIterationMask, 0, 0, VKLoopDirection>();

        // execute all from levels
        acc::apply_to_all<
            KLoopBoundaryLevelFunctor<TIterationMask, TStencilSweepFunctor, VKLoopDirection>,
            typename boost::mpl::pop_front<BeginBoundaryDomains>::type,
            TContext
        >(context);

        // init the loop counter with the compile time length of the loop range
        // (subtract one as the last iteration is done outside the loop) 
        int loopCounter = 
            k_positions_compile_time_length_component<
                TBeginKPosition, 
                TEndKPosition
            >::value - 
            (VBeginKBoundarySize > 1 ? VBeginKBoundarySize : 1) - 
            (VEndKBoundarySize > 1 ? VEndKBoundarySize : 1);

        // add the run time length deviation if necessary
        if(have_k_positions_run_time_length_component<TBeginKPosition, TEndKPosition>::value)
        {
            loopCounter += (context.kSize() - cDefaultKSize); 
        }
       
        // execute the main loop
        ACC_ASSERT(loopCounter >= 0);
        for(; loopCounter > 0; --loopCounter)
        {
            TStencilSweepFunctor::template Do<TDomain, false, false>(context);
            context.template Advance<TIterationMask, 0, 0, VKLoopDirection>();
        }

        // compute a type list containing all boundary domains
        // (make sure there is at least one boundary domain used for cache initialization)
        typedef typename compute_end_boundary_domains<
            TDomain,
            VKLoopDirection,
            TEndKPosition::KOffset::value,
            (VEndKBoundarySize > 1 ? VEndKBoundarySize : 1)
        >::type EndBoundaryDomains;

        // execute all to levels except the last one
        acc::apply_to_all<
            KLoopBoundaryLevelFunctor<TIterationMask, TStencilSweepFunctor, VKLoopDirection>,
            typename boost::mpl::pop_back<EndBoundaryDomains>::type,
            TContext
        >(context);      

        // updated the top most level
        TStencilSweepFunctor::template Do<typename boost::mpl::back<EndBoundaryDomains>::type, false, true>(context);
        // do not advance the iterators after the last loop level
    }
};

// specialization for k loops which do not require caching
template<
    typename TDomain,
    typename TIterationMask,
    typename TBeginKPosition,
    int VBeginKBoundarySize,
    typename TEndKPosition,
    int VEndKBoundarySize,
    KLoopDirection VKLoopDirection>
struct KLoop3DBody<
    TDomain, 
    TIterationMask, 
    TBeginKPosition, 
    VBeginKBoundarySize, 
    TEndKPosition, 
    VEndKBoundarySize, 
    false, // no cache fill or flush necessary
    VKLoopDirection>
{
    BOOST_STATIC_ASSERT(is_domain<TDomain>::value);
    BOOST_STATIC_ASSERT(is_k_position<TBeginKPosition>::value);
    BOOST_STATIC_ASSERT(is_k_position<TEndKPosition>::value);

    template<
        typename TContext,
        typename TStencilSweepFunctor>
    __ACC__
    static void Do(TContext& context)
    {
        // compute a type list containing all boundary domains
        // (make sure there is at least one boundary domain used for cache initialization)
        typedef typename compute_begin_boundary_domains<
            TDomain,
            VKLoopDirection,
            TBeginKPosition::KOffset::value,
            VBeginKBoundarySize
        >::type BeginBoundaryDomains;

        // execute all from levels
        acc::apply_to_all<
            KLoopBoundaryLevelFunctor<TIterationMask, TStencilSweepFunctor, VKLoopDirection>,
            BeginBoundaryDomains,
            TContext
        >(context);

        // init the loop counter with the compile time length of the loop range
        // (subtract one as the last iteration is done outside the loop) 
        int loopCounter = 
            k_positions_compile_time_length_component<
                TBeginKPosition, 
                TEndKPosition
            >::value - 
            VBeginKBoundarySize - 
            (VEndKBoundarySize > 1 ? VEndKBoundarySize : 1);

        // add the run time length deviation if necessary
        if(have_k_positions_run_time_length_component<TBeginKPosition, TEndKPosition>::value)
        {
            loopCounter += (context.kSize() - cDefaultKSize); 
        }
       
        // execute the main loop
        ACC_ASSERT(loopCounter >= 0);
        for(; loopCounter > 0; --loopCounter)
        {
            TStencilSweepFunctor::template Do<TDomain, false, false>(context);
            context.template Advance<TIterationMask, 0, 0, VKLoopDirection>();
        }

        // compute a type list containing all boundary domains
        // (make sure there is at least one boundary domain used for cache initialization)
        typedef typename compute_end_boundary_domains<
            TDomain,
            VKLoopDirection,
            TEndKPosition::KOffset::value,
            (VEndKBoundarySize > 1 ? VEndKBoundarySize : 1)
        >::type EndBoundaryDomains;

        // execute all to levels except the last one
        acc::apply_to_all<
            KLoopBoundaryLevelFunctor<TIterationMask, TStencilSweepFunctor, VKLoopDirection>,
            typename boost::mpl::pop_back<EndBoundaryDomains>::type,
            TContext
        >(context);      

        // updated the top most level
        TStencilSweepFunctor::template Do<typename boost::mpl::back<EndBoundaryDomains>::type, false, false>(context);
        // do not advance the iterators after the last loop level
    }
};

// specialization for k loops neither boundary layers nor caching
// (should simplify compilation of stencils which horizontal dependencies only)
template<
    typename TDomain,
    typename TIterationMask,
    typename TBeginKPosition,
    typename TEndKPosition,
    KLoopDirection VKLoopDirection>
struct KLoop3DBody<
    TDomain, 
    TIterationMask, 
    TBeginKPosition, 0, // no begin boundary
    TEndKPosition, 0, // no end boundary
    false, // no cache fill or flush necessary
    VKLoopDirection>
{
    BOOST_STATIC_ASSERT(is_domain<TDomain>::value);
    BOOST_STATIC_ASSERT(is_k_position<TBeginKPosition>::value);
    BOOST_STATIC_ASSERT(is_k_position<TEndKPosition>::value);

    template<
        typename TContext,
        typename TStencilSweepFunctor>
    __ACC__
    static void Do(TContext& context)
    {
        // init the loop counter with the compile time length of the loop range
        // (subtract one as the last iteration is done outside the loop) 
        int loopCounter = 
            k_positions_compile_time_length_component<
                TBeginKPosition, 
                TEndKPosition
            >::value - 1; 
            
        // add the run time length deviation if necessary
        if(have_k_positions_run_time_length_component<TBeginKPosition, TEndKPosition>::value)
        {
            loopCounter += (context.kSize() - cDefaultKSize); 
        }
           
        // execute the main loop
        ACC_ASSERT(loopCounter >= 0);
        for(; loopCounter > 0; --loopCounter)
        {
            TStencilSweepFunctor::template Do<TDomain, false, false>(context);
            context.template Advance<TIterationMask, 0, 0, VKLoopDirection>();
        }
            
        // updated the top most level
        TStencilSweepFunctor::template Do<TDomain, false, false>(context);
        // do not advance the iterators after the last loop level
    }
};




