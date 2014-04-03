#pragma once

#include <boost/mpl/front.hpp>
#include <boost/mpl/back.hpp>
#include "Definitions.h"
#include "ApplyToAll.h"
#include "StencilSweepDescriptor.h"
#include "KLoopRestoreAndAdvance.h"
#include "KLoopAdvance.h"
#include "KLoopRange.h"
#include "KLoopLegFunctor.h"

/**
* @struct KLoopCascade
* Structure implementing a loop over the full k domain given a sweep descriptor
* (one or more k loops covering the k domain)
*/
template<
    typename TContext,
    typename TStencilSweepDescriptor>
struct KLoopCascade;

template<
    typename TContext,
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection> 
struct KLoopCascade<TContext, StencilSweepDescriptor<TCaches, TStencilStages, VKLoopDirection> >
{
    // define the k loops of the cascade
    typedef typename compute_k_loop_ranges<
        TStencilStages
    >::type KLoopRanges;

    // compute the domain containing all loop ranges
    typedef typename compute_k_loop_ranges_base_domain<
        KLoopRanges
    >::type Domain;

    // compute the k loop legs
    typedef typename compute_k_loop_legs<
        KLoopRanges, 
        TCaches, 
        TStencilStages, 
        VKLoopDirection
    >::type KLoopLegs;

    // compute from and to positions of the loop cascade
    // (using the first and last k loop of the cascade)
    typedef typename k_loop_leg_from_k_position<
        typename boost::mpl::front<KLoopLegs>::type
    >::type FromKPosition;
    typedef typename k_loop_leg_to_k_position<
        typename boost::mpl::back<KLoopLegs>::type
    >::type ToKPosition;

    // compute the iteration mask covering all stencil sweep parameters
    typedef typename TContext::template compute_iteration_mask<
        TStencilStages, 
        Domain
    >::type IterationMask;

    /**
    * Restore the context position and move the position by an ij offset
    * @param context context object
    * @param iOffset offset in i direction 
    * @param jOffset offset in j direction
    */
    __ACC__
    static void RestoreAndAdvance(TContext& context, const int iOffset, const int jOffset)
    {
        // move the iterators / indexes to the from-position of the first k loop
        KLoopRestoreAndAdvance<
            FromKPosition,
            IterationMask
        >::Do(context, iOffset, jOffset);
    }

    /**
    * Advance the context position to the next column given a constant ij offset
    * @param context context object
    */
    template<int VIOffset, int VJOffset>
    __ACC__
    static void ColumnAdvance(TContext& context)
    {
        // move the iterators / indexes from the to-position of the last k loop to the from-position of the first k loop
        KLoopAdvance<
            ToKPosition,
            FromKPosition,
            IterationMask,
            VIOffset,
            VJOffset
        >::Do(context);
    }

    /**
    * Method executing the k loop
    * @param context context object
    */
    __ACC__
    static void Do(TContext& context)
    {
        // apply the k loops of the cascade consecutively
        acc::apply_to_all<
            KLoopLegFunctor<IterationMask>,
            KLoopLegs,
            TContext
        >(context);  
    }
};
