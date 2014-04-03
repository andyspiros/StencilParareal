#pragma once

#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_back.hpp>
#include "SharedInfrastructure.h"
#include "ParameterType.h"
#include "BlockLoopFunctorOpenMP.h"
#include "StencilSweepDescriptor.h"
#include "StencilSweepGroupDescriptor.h"

/**
* @struct set_stencil_sweep_descriptor_stages
* Meta function stripping a sweep descriptor by replacing the stencil stages by a single stage and removing the caches
*/
template<
    typename TStencilSweepDescriptor,
    typename TStencilStage>
struct strip_stencil_sweep_descriptor;

template<
    typename TStencilStage,
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection>
struct strip_stencil_sweep_descriptor<StencilSweepDescriptor<TCaches, TStencilStages, VKLoopDirection>, TStencilStage>
{
    typedef StencilSweepDescriptor<
        boost::mpl::void_, 
        boost::mpl::vector1<TStencilStage>, 
        VKLoopDirection
    > type;
};

/**
* @struct expand_stencil_sweep_descriptor
* Meta function expanding a stencil sweep with multiple stages into multiple stripped stencil sweeps (single stencil stage no cache)
*/
template<
    typename TStencilSweepDescriptors,
    typename TStencilSweepDescriptor>
struct expand_stencil_sweep_descriptor :
    boost::mpl::fold<
        typename stencil_sweep_descriptor_stages<TStencilSweepDescriptor>::type,
        TStencilSweepDescriptors,
        boost::mpl::push_back<
            boost::mpl::_1, 
            strip_stencil_sweep_descriptor<TStencilSweepDescriptor, boost::mpl::_2>
        >
    >
{};

/**
* @struct StencilSweepGroupFunctorOpenMP
* OpenMP functor which executes the stencil sweeps of the group in case the stencil sweep group condition is fulfilled
*/
struct StencilSweepGroupFunctorOpenMP
{
    template<
        typename TContext,
        typename TStencilSweepGroupDescriptor>
    static void Do(typename parameter_type<TContext>::type context)
    {
        // check if there is a valid comparison value 
        // which needs to be evaluated before the sweep group execution
        typedef typename boost::mpl::is_not_void_<
            typename stencil_sweep_group_descriptor_comparison_value<TStencilSweepGroupDescriptor>::type
        >::type PerformComparison;

        // select the do implementation with or without comparison
        doImpl<
            TContext, 
            TStencilSweepGroupDescriptor
        >(context, static_cast<PerformComparison*>(0));
    }

private:
    // implementation which directly executes the stencil sweeps 
    // (chosen if the comparison value is void)
    template<
        typename TContext, 
        typename TStencilSweepGroupDescriptor>
    static void doImpl(TContext& context, boost::mpl::false_*)
    {
        // apply the stencil sweeps
        applyStencilSweepGroup<
            TContext, 
            TStencilSweepGroupDescriptor
        >(context);
    }

    // implementation which executes the stencil sweeps conditionally
    // (chosen if the comparison value is not void)
    template<
        typename TContext, 
        typename TStencilSweepGroupDescriptor>
    static void doImpl(TContext& context, boost::mpl::true_*)
    {
        // extract the type information necessary for the comparison operation
        typedef typename TStencilSweepGroupDescriptor::ComparisonValue ComparisonValue;
        typedef DataParameter<
            TStencilSweepGroupDescriptor::ParameterIndex::value,
            Offset<0,0,0>
        > ParameterIndex;

        // compare the compile time comparison value to the parameter with the given index
        if(context[ParameterIndex()] == ComparisonValue::value)
        {
            // apply the stencil sweeps
            applyStencilSweepGroup<
                TContext, 
                TStencilSweepGroupDescriptor
            >(context);
        }
    }

    // method applying the stencil sweeps of the group
    template<
        typename TContext,
        typename TStencilSweepGroupDescriptor>
    static void applyStencilSweepGroup(TContext& context)
    {
        // extract the stencil sweep descriptors
        typedef typename stencil_sweep_group_descriptor_sweeps<
            TStencilSweepGroupDescriptor
        >::type StencilSweepDescriptors;

        // expand merged stages in the OpenMP implementation
        typedef typename boost::mpl::fold<
            StencilSweepDescriptors,
            boost::mpl::vector0<>,
            expand_stencil_sweep_descriptor<boost::mpl::_1, boost::mpl::_2>
        >::type ExpandedStencilSweepDescriptors;

        // finally apply the sweeps using the block loop functor
        apply_to_all<
            BlockLoopFunctorOpenMP,
            ExpandedStencilSweepDescriptors,
            TContext
        >(context);
    }
};

