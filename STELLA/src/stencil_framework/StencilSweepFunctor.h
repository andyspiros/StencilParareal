#pragma once

/**
* @struct StencilSweepFunctor
* Stencil sweep functor interface
*/
template<typename TStencilSweepFunctorImpl>
struct StencilSweepFunctor
{
    template<
        typename TDomain,
        bool VFillCaches,
        bool VFlushCaches,
        typename TContext>
    __ACC__
    static void Do(TContext& context)
    {
        TStencilSweepFunctorImpl::template DoImpl<TDomain, VFillCaches, VFlushCaches, TContext>(context);
    }
};

/**
* @struct create_stencil_sweep_functor
* Meta function returning the stencil sweep functor type
*/
template<
    typename TStencilStages,
    typename TCaches>
struct create_stencil_sweep_functor;





