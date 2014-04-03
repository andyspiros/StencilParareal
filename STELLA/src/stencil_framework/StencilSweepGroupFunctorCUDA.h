#pragma once

#include <boost/mpl/void.hpp>
#include <boost/mpl/if.hpp>
#include "SharedInfrastructure.h"
#include "ParameterType.h"
#include "BlockLoopFunctorCUDA.h"
#include "StencilSweepDescriptor.h"
#include "StencilSweepGroupDescriptor.h"

/**
* @struct StencilSweepGroupFunctorCUDA
* CUDA functor which executes the stencil sweeps of the group in case the stencil sweep group condition is fulfilled
*/
struct StencilSweepGroupFunctorCUDA
{
    template<
        typename TSharedDataPointer,
        typename TStencilSweepGroupDescriptor>
    __ACC__
    static void Do(typename parameter_type<TSharedDataPointer>::type pSharedData)
    {
        // check if there is a valid comparison value 
        // which needs to be evaluated before the sweep group execution
        typedef typename boost::mpl::is_not_void_<
            typename stencil_sweep_group_descriptor_comparison_value<TStencilSweepGroupDescriptor>::type
        >::type PerformComparison;

        // select the do implementation with or without comparison
        doImpl<
            TSharedDataPointer, 
            TStencilSweepGroupDescriptor
        >(pSharedData, static_cast<PerformComparison*>(0));
    }

private:
    // implementation which directly executes the stencil sweeps 
    // (chosen if the comparison value is void)
    template<
        typename TSharedDataPointer, 
        typename TStencilSweepGroupDescriptor>
    __ACC__
    static void doImpl(TSharedDataPointer pSharedData, boost::mpl::false_*)
    {
        // apply the stencil sweeps
        applyStencilSweepGroup<
            TSharedDataPointer, 
            TStencilSweepGroupDescriptor
        >(pSharedData);
    }

    // implementation which executes the stencil sweeps conditionally
    // (chosen if the comparison value is not void)
    template<
        typename TSharedDataPointer, 
        typename TStencilSweepGroupDescriptor>
    __ACC__
    static void doImpl(TSharedDataPointer pSharedData, boost::mpl::true_*)
    {
        // extract the type information necessary for the comparison operation
        typedef typename TStencilSweepGroupDescriptor::ComparisonValue ComparisonValue;
        typedef typename TStencilSweepGroupDescriptor::ParameterIndex ParameterIndex;

        // make sure the shared data pointer is valid
        ACC_ASSERT(pSharedData_);

        // compare the compile time comparison value to the parameter with the given index
        if(pSharedData->dataTuple(static_cast<ParameterIndex*>(0)).value() == ComparisonValue::value)
        {
            // apply the stencil sweeps
            applyStencilSweepGroup<
                TSharedDataPointer, 
                TStencilSweepGroupDescriptor
            >(pSharedData);
        }
    }

    // method applying the stencil sweeps of the group
    template<
        typename TSharedDataPointer,
        typename TStencilSweepGroupDescriptor>
    __ACC__
    static void applyStencilSweepGroup(TSharedDataPointer pSharedData)
    {
        // extract the stencil sweep descriptors
        typedef typename stencil_sweep_group_descriptor_sweeps<
            TStencilSweepGroupDescriptor
        >::type StencilSweepDescriptors;

        // finally apply the sweeps using the block loop functor
        acc::apply_to_all<
            BlockLoopFunctorCUDA,
            StencilSweepDescriptors,
            TSharedDataPointer
        >(pSharedData);
    }
};