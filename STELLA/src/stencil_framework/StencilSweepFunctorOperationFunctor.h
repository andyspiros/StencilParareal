#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/void.hpp>
#include "KRange.h"
#include "StencilStage.h"
#include "StencilSweepFunctorOperation.h"

template<
    typename TDomain,
    bool VFillCaches,
    bool VFlushCaches>
struct StencilSweepFunctorOperationFunctor
{
    template<
        typename TContext, 
        typename TStencilSweepFunctorOperation>
    __ACC__
    static void Do(typename parameter_type<TContext>::type context)
    {
        // select the operation specific do implementation
        doImpl(context, static_cast<TStencilSweepFunctorOperation*>(0));
    }

private:
    // operation specific do method implementations
    template<
        typename TContext,
        typename TSlideMask,
        typename TIJRange>
    __ACC__
    static void doImpl(TContext& context, StencilSweepFunctorOperationSlide<TSlideMask, TIJRange>*)
    {
        // slide the caches except on the first k loop leg level
        if(!VFillCaches)
        {
            context.template SlideCaches<TSlideMask>();
        }
    }

    template<
        typename TContext,
        typename TFillMask,
        typename TIJRange>
    __ACC__
    static void doImpl(TContext& context, StencilSweepFunctorOperationFill<TFillMask, TIJRange>*)
    {
        // on the first k loop level fill the full cache 
        // (otherwise only fill the cache front)
        if(VFillCaches)
        {
            context.template FillCaches<TFillMask, TDomain>();
        }
        else
        {
            context.template FillCachesFront<TFillMask, TDomain>();
        }
    }

    template<
        typename TContext,
        typename TStencilStage>
    __ACC__
    static void doImpl(TContext& context, StencilSweepFunctorOperationStage<TStencilStage>*)
    {
        // define the stage stencil and apply it
        typedef typename stencil_stage_definition<
            TStencilStage, 
            TContext, 
            TDomain
        >::type StencilStageDefinition;

        // define the stencil stage do domain
        typedef typename stencil_stage_lookup_do_domain<
            TStencilStage,
            TDomain
        >::type DoDomain;
        BOOST_STATIC_ASSERT(boost::mpl::is_not_void_<DoDomain>::value);

        // make sure the domain is inside the k range of the stencil stage
        BOOST_STATIC_ASSERT((in_k_range<typename stencil_stage_k_range<TStencilStage>::type, TDomain>::value));

        // run the stencil stage
        StencilStageDefinition::Do(context, DoDomain());
    }

    template<
        typename TContext,
        typename TFlushMask,
        typename TIJRange>
    __ACC__
    static void doImpl(TContext& context, StencilSweepFunctorOperationFlush<TFlushMask, TIJRange>*)
    {
        // on the last k loop level flush the full cache 
        // (otherwise only flush the cache back)
        if(VFlushCaches)
        {
            context.template FlushCaches<TFlushMask, TDomain>();
        }
        else
        {
            context.template FlushCachesBack<TFlushMask, TDomain>();
        }
    }
};




