#pragma once

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/identity.hpp>
#include "StencilSweepDescriptor.h"
#include "IJRange.h"
#include "ComparisonValue.h"

/**
* @struct StencilSweepGroupDescriptor
* Structure defining a stencil sweep group which is group of stencil sweeps which are executed conditionally.
* The condition is satisfied if either the comparison value is void and the parameter index -1 or if
* the comparison value is equivalent to the stencil parameter with the given parameter index.
*/
template<
    typename TStencilSweepDescriptors,
    typename TComparisonValue,
    int VParameterIndex>
struct StencilSweepGroupDescriptor
{
    // check TStencilSweepDescriptors is a sequence of stencil sweep descriptors 
    BOOST_STATIC_ASSERT( 
        (
            boost::mpl::eval_if<
                boost::mpl::is_sequence<TStencilSweepDescriptors>,
                boost::mpl::fold<
                    TStencilSweepDescriptors,
                    boost::mpl::true_,
                    boost::mpl::if_<
                        is_stencil_sweep_descriptor<boost::mpl::_2>,
                        boost::mpl::_1,
                        boost::mpl::false_
                    >
                >,
                boost::mpl::false_
            >::type::value
        )
    );

    // check if TComparisonValue is void or a compile time value
    BOOST_STATIC_ASSERT(
        (
            is_comparison_value<TComparisonValue>::value || 
            (boost::mpl::is_void_<TComparisonValue>::value && VParameterIndex == -1)
        )
    );

    // define comparison type information
    typedef TComparisonValue ComparisonValue;
    typedef boost::mpl::integral_c<int, VParameterIndex> ParameterIndex;
};

/**
* @struct is_stencil_sweep_group_descriptor
* Meta function returning true if the parameter is a stencil sweep group descriptor
*/
template<typename T>
struct is_stencil_sweep_group_descriptor : boost::mpl::false_ {};

template<
    typename TStencilSweepDescriptors,
    typename TComparisonValue,
    int VParameterIndex>
struct is_stencil_sweep_group_descriptor<
    StencilSweepGroupDescriptor<TStencilSweepDescriptors, TComparisonValue, VParameterIndex> > : boost::mpl::true_ {};

/**
* @struct stencil_sweep_group_descriptor_sweeps
* Meta function returning the stencil sweep descriptors of a stencil sweep group descriptor
*/
template<typename T>
struct stencil_sweep_group_descriptor_sweeps;

template<
    typename TStencilSweepDescriptors,
    typename TComparisonValue,
    int VParameterIndex>
struct stencil_sweep_group_descriptor_sweeps<
    StencilSweepGroupDescriptor<TStencilSweepDescriptors, TComparisonValue, VParameterIndex> >
{
    typedef TStencilSweepDescriptors type;
};

/**
* @struct stencil_sweep_group_descriptor_parameter_index
* Meta function returning the parameter index of a stencil sweep group descriptor
*/
template<typename T>
struct stencil_sweep_group_descriptor_parameter_index;

template<
    typename TStencilSweepDescriptors,
    typename TComparisonValue,
    int VParameterIndex>
struct stencil_sweep_group_descriptor_parameter_index<
    StencilSweepGroupDescriptor<TStencilSweepDescriptors, TComparisonValue, VParameterIndex> > : boost::mpl::integral_c<int, VParameterIndex> 
{};

/**
* @struct stencil_sweep_group_descriptor_comparison_value
* Meta function returning the comparison value of a stencil sweep group descriptor
*/
template<typename T>
struct stencil_sweep_group_descriptor_comparison_value;

template<
    typename TStencilSweepDescriptors,
    typename TComparisonValue,
    int VParameterIndex>
struct stencil_sweep_group_descriptor_comparison_value<
    StencilSweepGroupDescriptor<TStencilSweepDescriptors, TComparisonValue, VParameterIndex> >
{
    typedef TComparisonValue type;
};

/**
* @struct stencil_sweep_group_descriptor_parameter_ij_range
* Meta function returning the maximal ij range of all sweep group stencil stages using a given parameter
*/
template<
    typename TStencilSweepGroupDescriptor, 
    typename TParameterIndex>
struct stencil_sweep_group_descriptor_parameter_ij_range 
{
    typedef typename boost::mpl::fold<
        typename stencil_sweep_group_descriptor_sweeps<TStencilSweepGroupDescriptor>::type,
        IJRange<cComplete, 0, 0, 0, 0>,
        maximum_ij_range<
            boost::mpl::_1,
            stencil_sweep_descriptor_parameter_ij_range<boost::mpl::_2, TParameterIndex>
        >
    >::type type;
};

/**
* @struct stencil_sweep_group_descriptors_is_parameter_usage_unique
* Meta function returning false if a parameter is used by more than one stage of a stencil sweep
*/
template<
    typename TStencilSweepGroupDescriptors,
    typename TParameterIndex>
struct stencil_sweep_group_descriptors_is_parameter_usage_unique
{
    typedef typename boost::mpl::fold<
        TStencilSweepGroupDescriptors,
        boost::mpl::true_,
        boost::mpl::if_<
            stencil_sweep_descriptors_is_parameter_usage_unique<
                stencil_sweep_group_descriptor_sweeps<boost::mpl::_2>,
                TParameterIndex
            >,
            boost::mpl::_1,
            boost::mpl::false_
        >
    >::type type;
};

/**
* @struct stencil_sweep_group_descriptors_is_parameter_usage_unique
* Meta function returning true if a given parameter is accessed not counting local cache accesses
* (note that the meta function is used to decide if a temporary field needs to be allocated,
* therefore local cache accesses are not considered as the underlying temporary field is not accessed)
*/
template<
    typename TStencilSweepGroupDescriptors,
    typename TDomain,
    typename TParameterIndex>
struct stencil_sweep_group_descriptors_is_parameter_accessed
{
    // iterate over all sweeps and search of a parameter access
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            boost::mpl::fold<
                TStencilSweepGroupDescriptors,
                boost::mpl::false_,
                boost::mpl::if_<
                    stencil_sweep_descriptors_is_parameter_accessed<
                        stencil_sweep_group_descriptor_sweeps<boost::mpl::_2>,
                        TDomain,
                        TParameterIndex
                    >,
                    boost::mpl::true_,
                    boost::mpl::_1
                >
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<bool, bool(value)> type;
};

/**
* @struct create_stencil_sweep_group_descriptor_from_sweep
* Meta function creating a stencil sweep group descriptor with void comparison value given a stencil sweep descriptor
*/
template<typename TStencilSweepDescriptor>
struct create_stencil_sweep_group_descriptor_from_sweep
{
    BOOST_STATIC_ASSERT(is_stencil_sweep_descriptor<TStencilSweepDescriptor>::value);

    // define the stencil sweep group descriptor
    typedef StencilSweepGroupDescriptor<
        boost::mpl::vector1<TStencilSweepDescriptor>, 
        boost::mpl::void_, 
        -1
    > type;
};

/**
* @struct stencil_sweep_group_descriptors_is_merge_possible
* Meta function returning true if it is possible to merge two stencil sweep group descriptors
*/
template<
    typename TStencilSweepGroupDescriptor1,
    typename TStencilSweepGroupDescriptor2>
struct stencil_sweep_group_descriptors_is_merge_possible : boost::mpl::false_ {};

template<
    typename TStencilSweepDescriptors1,
    typename TStencilSweepDescriptors2,
    typename TComparisonValue,
    int VParameterIndex>
struct stencil_sweep_group_descriptors_is_merge_possible<
    StencilSweepGroupDescriptor<TStencilSweepDescriptors1, TComparisonValue, VParameterIndex>,
    StencilSweepGroupDescriptor<TStencilSweepDescriptors2, TComparisonValue, VParameterIndex> > : boost::mpl::true_ {};

/**
* @struct merge_stencil_sweep_group_descriptors
* Meta function merging two stencil sweep group descriptors with identical comparison value and parameter index
*/
template<
    typename TStencilSweepGroupDescriptor1,
    typename TStencilSweepGroupDescriptor2>
struct merge_stencil_sweep_group_descriptors;

template<
    typename TStencilSweepDescriptors1,
    typename TStencilSweepDescriptors2,
    typename TComparisonValue,
    int VParameterIndex>
struct merge_stencil_sweep_group_descriptors<
    StencilSweepGroupDescriptor<TStencilSweepDescriptors1, TComparisonValue, VParameterIndex>,
    StencilSweepGroupDescriptor<TStencilSweepDescriptors2, TComparisonValue, VParameterIndex> >
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TStencilSweepDescriptors1>::value);
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TStencilSweepDescriptors2>::value);

    // specialization concatenating the two sweep descriptor vectors
    // the comparison value and the parameter index have to be identical
    typedef StencilSweepGroupDescriptor<
        typename boost::mpl::fold<
            TStencilSweepDescriptors2,
            TStencilSweepDescriptors1,
            boost::mpl::push_back<boost::mpl::_1, boost::mpl::_2>
        >::type,
        TComparisonValue,
        VParameterIndex
    > type;
};

/**
* @struct stencil_sweep_group_descriptors_add_stencil_sweep_group_descriptor
* Meta function adding a stencil sweep group descriptor to a vector of stencil sweep group descriptors
* (note that successive stencil sweep group descriptors if possible are merged into a single stencil sweep group descriptor)
*/
template<
    typename TStencilSweepGroupDescriptors,
    typename TStencilSweepGroupDescriptor>
struct stencil_sweep_group_descriptors_add_stencil_sweep_group_descriptor
{
    // extract the preceding stencil sweep group descriptor
    typedef typename boost::mpl::eval_if<
        boost::mpl::empty<TStencilSweepGroupDescriptors>,
        boost::mpl::void_,
        boost::mpl::back<TStencilSweepGroupDescriptors>
    >::type LastStencilSweepGroupDescriptor;
    
    // figure out if it is possible to merge the last stencil sweep group descriptor with the new one
    typedef typename stencil_sweep_group_descriptors_is_merge_possible<
        LastStencilSweepGroupDescriptor,
        TStencilSweepGroupDescriptor
    >::type IsMergePossible;
    
    // compute a new stencil sweep group descriptors vector
    // merge the new stencil sweep group descriptor with the last from the vector if possible
    // otherwise add the new stencil sweep group descriptor at the end of the vector
    typedef typename boost::mpl::push_back<
        typename boost::mpl::eval_if<
            IsMergePossible,
            boost::mpl::pop_back<TStencilSweepGroupDescriptors>, 
            boost::mpl::identity<TStencilSweepGroupDescriptors>
        >::type,
        typename boost::mpl::eval_if<
            IsMergePossible,
            merge_stencil_sweep_group_descriptors<LastStencilSweepGroupDescriptor, TStencilSweepGroupDescriptor>,
            boost::mpl::identity<TStencilSweepGroupDescriptor>
        >::type
    >::type type;
};