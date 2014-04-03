#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/void.hpp>
#include "KRange.h"
#include "StencilStage.h"
#include "StencilSweepFunctor.h"

/**
* @struct StencilSweepFunctorOpenMP
* Functor executing all stencil stages of a stencil sweep
* (the OpenMP backend always runs just one stencil stage per sweep) 
*/
template<typename TStencilStages>
struct StencilSweepFunctorOpenMP : StencilSweepFunctor<StencilSweepFunctorOpenMP<TStencilStages> > // CRTP
{
    // implement OpenMP apply of stencil stages
    template<
        typename TDomain,
        bool VFillCaches,
        bool VFlushCaches,
        typename TContext>
    __ACC__
    static void DoImpl(TContext& context)
    {
        // assert that no stencil stage sequences are passed
        BOOST_STATIC_ASSERT(boost::mpl::size<TStencilStages>::value == 1);
        typedef typename boost::mpl::front<TStencilStages>::type StencilStage;
        
        // extract stencil stage definition functor 
        typedef typename stencil_stage_definition<
            StencilStage, 
            TContext, 
            TDomain
        >::type StencilStageDefinition;

        // find best matching do method of the stage
        typedef typename stencil_stage_lookup_do_domain<
            StencilStage, 
            TDomain
        >::type DoDomain;
        // make sure a do method was found
        BOOST_STATIC_ASSERT(boost::mpl::is_not_void_<DoDomain>::value);
        
        // make sure the domain is inside the k range of the stencil stage
        BOOST_STATIC_ASSERT((in_k_range<typename stencil_stage_k_range<StencilStage>::type, TDomain>::value));

        // finally execute the stencil stage
        StencilStageDefinition::Do(context, DoDomain());
    }
};

// create a OpenMP stencil sweep functor
template<
    typename TStencilStages,
    typename TCaches>
struct create_stencil_sweep_functor
{
    typedef StencilSweepFunctorOpenMP<TStencilStages> type;
};


