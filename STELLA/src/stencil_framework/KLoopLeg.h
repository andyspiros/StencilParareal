#pragma once

#include <boost/config.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/partition.hpp>
#include <boost/mpl/reverse_fold.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/if.hpp>
#include "Enums.h"
#include "KLoopRange.h"
#include "Cache.h"
#include "StencilStage.h"

/**
* @struct KLoopLeg
* Structure defining a k loop leg used by the k loop leg functor in order to execute a k loop cascade element
*/
template<
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages,
    typename TAdvanceToKPosition,
    KLoopDirection VKLoopDirection>
struct KLoopLeg 
{
    // check the input parameters
    BOOST_STATIC_ASSERT(is_k_position<TAdvanceToKPosition>::value || boost::mpl::is_void_<TAdvanceToKPosition>::value); 

    typedef TKLoopRange KLoopRange;
    typedef TCaches Caches;
    typedef TStencilStages StencilStages;
    typedef TAdvanceToKPosition AdvanceToKPosition;
    typedef boost::mpl::integral_c< ::KLoopDirection, VKLoopDirection> KLoopDirection;
};

/**
* @struct k_loop_leg_from_k_position
* Meta function computing the from position of the k loop leg given
*/
template<typename TKLoopLeg>
struct k_loop_leg_from_k_position;

template<
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages,
    typename TAdvanceToKPosition,
    KLoopDirection VKLoopDirection>
struct k_loop_leg_from_k_position<KLoopLeg<TKLoopRange, TCaches, TStencilStages, TAdvanceToKPosition, VKLoopDirection> >
{
    typedef typename boost::mpl::eval_if_c<
        VKLoopDirection == cKIncrement,
        k_loop_range_k_minimum<TKLoopRange>,
        k_loop_range_k_maximum<TKLoopRange>        
    >::type type;
};

/**
* @struct k_loop_leg_to_k_position
* Meta function computing the to position of the k loop leg given
*/
template<typename TKLoopLeg>
struct k_loop_leg_to_k_position;

template<
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages,
    typename TAdvanceToKPosition,
    KLoopDirection VKLoopDirection>
struct k_loop_leg_to_k_position<KLoopLeg<TKLoopRange, TCaches, TStencilStages, TAdvanceToKPosition, VKLoopDirection> >
{
    typedef typename boost::mpl::eval_if_c<
        VKLoopDirection == cKIncrement,
        k_loop_range_k_maximum<TKLoopRange>,
        k_loop_range_k_minimum<TKLoopRange>
    >::type type;
};

/**
* @struct k_loop_leg_from_k_boundary_size
* Meta function computing the from boundary size of the k loop leg given
*/
template<typename TKLoopLeg>
struct k_loop_leg_from_k_boundary_size;

template<
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages,
    typename TAdvanceToKPosition,
    KLoopDirection VKLoopDirection>  
struct k_loop_leg_from_k_boundary_size<KLoopLeg<TKLoopRange, TCaches, TStencilStages, TAdvanceToKPosition, VKLoopDirection> >
{
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::if_c<
                VKLoopDirection == cKIncrement,
                typename TKLoopRange::KBoundarySize::KMinimumSize,
                typename TKLoopRange::KBoundarySize::KMaximumSize
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct k_loop_leg_to_k_boundary_size
* Meta function computing the to boundary size of the k loop leg given
*/
template<typename TKLoopLeg>
struct k_loop_leg_to_k_boundary_size;

template<
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages,
    typename TAdvanceToKPosition,
    KLoopDirection VKLoopDirection>
struct k_loop_leg_to_k_boundary_size<KLoopLeg<TKLoopRange, TCaches, TStencilStages, TAdvanceToKPosition, VKLoopDirection> >
{
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::if_c<
                VKLoopDirection == cKIncrement,
                typename TKLoopRange::KBoundarySize::KMaximumSize,
                typename TKLoopRange::KBoundarySize::KMinimumSize                    
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct k_loop_legs_add_k_loop_leg
* Meta function creating and adding a new k loop leg at the beginning of the k loop legs sequence.
* We assume that the k loop legs are added in reverse execution order which means the front element 
* of k loop legs is executed right after the added element.
*/
template<
    typename TKLoopLegs,
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection>
struct k_loop_legs_add_k_loop_leg
{
    // method partitioning the caches enabled and disabled once
    // the k range of the disabled caches does not overlap with the k range of the leg
    typedef typename boost::mpl::eval_if<
        boost::mpl::is_void_<TCaches>,
        boost::mpl::pair<boost::mpl::vector0<>, boost::mpl::vector0<> >,
        boost::mpl::partition<
            TCaches,
            in_k_loop_range<TKLoopRange, cache_k_range<boost::mpl::_> >,
            boost::mpl::back_inserter<boost::mpl::vector0<> >,
            boost::mpl::back_inserter<boost::mpl::vector0<> >
        >
    >::type PartitionedCaches;
    typedef typename PartitionedCaches::first Caches;
    typedef typename PartitionedCaches::second DisabledCaches;
    
    // filter out all stages which are not touched by the k loop leg
    typedef typename boost::mpl::copy_if<
        TStencilStages,
        in_k_loop_range<TKLoopRange, stencil_stage_k_range<boost::mpl::_> >
    >::type StencilStages;

    // make sure that caches which are disabled in the k loop leg are not accesses
    // (note that a cached parameter shall be cached during the full sweep 
    // and no uncached accesses in a specific k loop leg shall be possible)
    BOOST_MPL_ASSERT_MSG(
        (
            boost::is_same<
                typename boost::mpl::find_if<
                    DisabledCaches,
                    stencil_stages_have_parameter<
                        StencilStages,
                        typename TKLoopRange::KRange::Domain,
                        cache_index<boost::mpl::_>
                    >
                >::type,
                typename boost::mpl::end<DisabledCaches>::type
            >::value
        ),
        CACHE_ACCESSED_OUT_OF_KRANGE,
        (DisabledCaches)
    );

    // extend the k loop range boundaries 
    // (necessary in order to handle the cache fill and flush operations properly)
    typedef typename boost::mpl::fold<
        Caches,
        TKLoopRange,
        extend_k_loop_range_boundary<
            boost::mpl::_1,
            cache_full_update_k_range<
                boost::mpl::_2,
                boost::mpl::integral_c<KLoopDirection, VKLoopDirection>
            >
        >
    >::type KLoopRange;

    // compute the following k loop leg or void if there is none
    typedef typename boost::mpl::eval_if<
        boost::mpl::empty<TKLoopLegs>,
        boost::mpl::void_,
        boost::mpl::front<TKLoopLegs>
    >::type FollowingKLoopLeg;

    // compute the advance to k position depending on the following k loop leg
    typedef typename boost::mpl::eval_if<
        boost::mpl::is_void_<FollowingKLoopLeg>,
        boost::mpl::void_, // if there is no subsequent k loop leg return void
        k_loop_leg_from_k_position<FollowingKLoopLeg> // otherwise return the from k position of the next leg
    >::type AdvanceToKPosition;
    
    // create and add a new k loop leg at beginning of the k loop legs sequence
    typedef typename boost::mpl::push_front<
        TKLoopLegs,
        KLoopLeg<
            KLoopRange,
            Caches,
            StencilStages,
            AdvanceToKPosition,
            VKLoopDirection
        >
    >::type type;
};

/**
* @struct k_loop_legs_add_k_increment_k_loop_leg
* Meta function calling k_loop_legs_add_k_loop_leg with k increment loop direction
*/
template<
    typename TKLoopLegs,
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages>
struct k_loop_legs_add_k_increment_k_loop_leg 
{
    typedef typename k_loop_legs_add_k_loop_leg<
        TKLoopLegs, 
        TKLoopRange, 
        TCaches, 
        TStencilStages, 
        cKIncrement
    >::type type;
};

/**
* @struct k_loop_legs_add_k_decrement_k_loop_leg
* Meta function calling k_loop_legs_add_k_loop_leg with k decrement loop direction
*/
template<
    typename TKLoopLegs,
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages>
struct k_loop_legs_add_k_decrement_k_loop_leg 
{
    typedef typename k_loop_legs_add_k_loop_leg<
        TKLoopLegs, 
        TKLoopRange, 
        TCaches, 
        TStencilStages, 
        cKDecrement
    >::type type;
};

/**
* @struct compute_k_loop_legs
* Meta function computing the k loop legs of k loop cascade
*/
template<
    typename TKLoopRanges,
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection>
struct compute_k_loop_legs
{
    // assemble the k loop cascade legs in reversed execution order 
    // as we need to compute the advance to k position based on the following k loop leg
    typedef typename boost::mpl::eval_if_c<
        VKLoopDirection == cKIncrement,
        // if the cascade loops in k increment direction create loops backwards
        // (note we need the start position of the following loop in the cascade)
        boost::mpl::reverse_fold<
            TKLoopRanges,
            boost::mpl::vector0<>,
            k_loop_legs_add_k_increment_k_loop_leg<
                boost::mpl::_1,
                boost::mpl::_2,
                TCaches,
                TStencilStages
            >
        >,
        // if the cascade loops in k decrement direction create loops forward
        // (note we need the start position of the following loop in the cascade)
        boost::mpl::fold<
            TKLoopRanges,
            boost::mpl::vector0<>,
            k_loop_legs_add_k_decrement_k_loop_leg<
                boost::mpl::_1,    
                boost::mpl::_2,
                TCaches,
                TStencilStages
            >
        >
    >::type type;
};