#pragma once

#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/erase_key.hpp>
#include <boost/mpl/if.hpp>
#include "Definitions.h"
#include "KLoopAdvance.h"
#include "KLoopLeg.h"
#include "KLoopBody.h"

// include the matching stencil stage functor
#ifdef __CUDA_BACKEND__
#include "StencilSweepFunctorCUDA.h"
#else
#include "StencilSweepFunctorOpenMP.h"
#endif

/**
* @struct KLoopLegFunctor
* Functor executing a k loop leg of a k loop cascade
*/
template<typename TIterationMask>
struct KLoopLegFunctor
{
    template<
        typename TContext,
        typename TKLoopLeg>
    __ACC__
    static void Do(TContext& context)
    {
        // define k loop leg properties
        typedef typename TKLoopLeg::KLoopDirection KLoopDirection;
        typedef typename TKLoopLeg::StencilStages StencilStages; 
        typedef typename TKLoopLeg::Caches Caches;
        typedef typename TKLoopLeg::KLoopRange::KRange::Domain Domain;
        
        // compute the loop iteration mask covering all parameters processed by the k loop
        // and the idle iteration mask covering all parameters not touched by the k loop
        typedef typename TContext::template compute_iteration_mask<
            StencilStages, 
            Domain
        >::type LoopIterationMask;
        // erase all loop iteration mask elements from the iteration mask
        typedef typename boost::mpl::fold<
            LoopIterationMask,
            TIterationMask, 
            boost::mpl::erase_key<boost::mpl::_1, boost::mpl::_2>
        >::type IdleIterationMask;

        // define closed loop interval from and to positions are updated by the k loop
        // (in order to avoid by-one-off-error) 
        typedef typename k_loop_leg_from_k_position<TKLoopLeg>::type FromKPosition;
        typedef typename k_loop_leg_to_k_position<TKLoopLeg>::type ToKPosition;

        // define the k loop boundary sizes
        typedef typename k_loop_leg_from_k_boundary_size<TKLoopLeg>::type BeginKBoundarySize;
        typedef typename k_loop_leg_to_k_boundary_size<TKLoopLeg>::type EndKBoundarySize;

        // define k loop begin and end positions plus the corresponding k offsets
        // (note that a half open interval is defined which is extended by 1 in loop direction)
        typedef FromKPosition BeginKPosition;
        typedef typename shift_k_position_by<
            ToKPosition, 
            KLoopDirection::value
        >::type EndKPosition;
        typedef typename FromKPosition::KOffset FromKOffset;
        typedef typename ToKPosition::KOffset ToKOffset;
        
        // define target k position used to advance iterators after the k loop
        typedef typename boost::mpl::if_<
            boost::mpl::is_void_<typename TKLoopLeg::AdvanceToKPosition>,
            ToKPosition,
            typename TKLoopLeg::AdvanceToKPosition
        >::type AdvanceToKPosition; 
        
        // iterate over all caches and search for a cache which requires either a fill or a flush
        // (used in order to instantiate a loop body with specific update function at the loop boundaries)
        typedef typename boost::mpl::not_<
            typename boost::is_same<
                typename boost::mpl::find_if<
                    Caches,
                    boost::mpl::or_<
                        cache_is_fill_needed<boost::mpl::_>, 
                        cache_is_flush_needed<boost::mpl::_> 
                    >
                >::type, 
                typename boost::mpl::end<Caches>::type
            >::type
        >::type FillOrFlushCaches;

        // compute the k loop body type
        typedef typename boost::mpl::if_<
            boost::is_same<FromKPosition, ToKPosition>,
            KLoop2DBody<
                Domain, 
                FillOrFlushCaches::value
            >,
            KLoop3DBody<
                Domain, 
                LoopIterationMask, 
                FromKPosition, 
                BeginKBoundarySize::value, 
                EndKPosition, 
                EndKBoundarySize::value, 
                FillOrFlushCaches::value,
                KLoopDirection::value
            >
        >::type KLoopBody;
    
        // define stencil stage functor
        typedef typename create_stencil_sweep_functor<
            StencilStages,
            Caches
        >::type StencilSweepFunctor;

        // execute the k loop body
        KLoopBody::template Do<
            TContext,
            StencilSweepFunctor            
        >(context);

        // move the iterators to the next loop
        KLoopAdvance<
            ToKPosition,
            AdvanceToKPosition,
            LoopIterationMask,
            0, 0
        >::Do(context);
        
        // move the iterators not used by the k loop to the start position of the next loop
        KLoopAdvance<
            FromKPosition,
            AdvanceToKPosition,
            IdleIterationMask,
            0, 0
        >::Do(context);
    }
};


