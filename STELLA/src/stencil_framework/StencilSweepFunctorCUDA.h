#pragma once

#include <boost/mpl/copy_if.hpp>
#include "ApplyToAll.h"
#include "KRange.h"
#include "StencilStage.h"
#include "StencilSweepFunctor.h"
#include "StencilSweepFunctorOperationGroupFunctor.h"

/**
* @struct StencilSweepFunctorCUDA
* Functor executing the stencil stages of a stencil sweep
*/
template<
    typename TStencilStages,
    typename TCaches>
struct StencilSweepFunctorCUDA : StencilSweepFunctor<StencilSweepFunctorCUDA<TStencilStages, TCaches> > // CRTP
{
    template<
        typename TDomain,
        bool VFillCaches,
        bool VFlushCaches,
        typename TContext>
    __ACC__
    static void DoImpl(TContext& context)
    {
        // compute a vector of stencil sweep functor operation groups
        // (every group of operations is applied to a common ij range)
        typedef typename compute_stencil_sweep_functor_operation_groups<
            TStencilStages,
            TCaches,
            TDomain
        >::type StencilSweepFunctorOperationGroups;

        // apply all stencil sweep functor operation groups
        acc::apply_to_all<
            StencilSweepFunctorOperationGroupFunctor<
                TDomain, 
                VFillCaches,
                VFlushCaches,
                (boost::mpl::size<StencilSweepFunctorOperationGroups>::value > 1)
            >, 
            StencilSweepFunctorOperationGroups, 
            TContext
        >(context);
    }
};

// create a CUDA stencil sweep functor
template<
    typename TStencilStages,
    typename TCaches>
struct create_stencil_sweep_functor
{
    typedef StencilSweepFunctorCUDA<
        TStencilStages,
        TCaches
    > type;
};


