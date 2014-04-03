#pragma once

#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/deref.hpp>
#include "ApplyToAll.h"
#include "StencilSweepFunctorOperationFunctor.h"
#include "StencilSweepFunctorOperationGroup.h"

/**
* @struct StencilSweepFunctorOperationGroupFunctor
* Functor applying a group of stencil sweep functor operations which are applied to a common ij range
*/
template<
    typename TDomain,
    bool VFillCaches,
    bool VFlushCaches,
    bool VSyncThreads>
struct StencilSweepFunctorOperationGroupFunctor
{
    template<
        typename TContext,
        typename TStencilSweepFunctorOperationGroup>
    __ACC__
    static void Do(typename parameter_type<TContext>::type context)
    {
        // extract the operations and their ij range from the operation group
        typedef typename TStencilSweepFunctorOperationGroup::IJRange IJRange;
        typedef typename TStencilSweepFunctorOperationGroup::StencilSweepFunctorOperations StencilSweepFunctorOperations;

        // pause all threads which are outside the operation ij range
        if(context.template IsInRange<IJRange>())
        {
            acc::apply_to_all<
                StencilSweepFunctorOperationFunctor<TDomain, VFillCaches, VFlushCaches>, 
                StencilSweepFunctorOperations, 
                TContext
            >(context);

            // note that operations with identical ij ranges cannot have dependencies in ij direction
            // therefore no synchronization is necessary between the individual operations of a group
        }

        // check if none of the stencil sweep functors requires a synchronization of the group
        // (after applying a stencil stage a sync is necessary while k cache updates require no sync)
        typedef typename boost::is_same<
            typename boost::mpl::find_if<
                StencilSweepFunctorOperations, 
                stencil_sweep_functor_operation_is_sync_needed<boost::mpl::_> 
            >::type,
            typename boost::mpl::end<StencilSweepFunctorOperations>::type
        >::type NoSyncRequired;

        // synchronize all threads of the block if necessary
        if(VSyncThreads && !NoSyncRequired::value)
        {
            __syncthreads();
        }
    }
};


