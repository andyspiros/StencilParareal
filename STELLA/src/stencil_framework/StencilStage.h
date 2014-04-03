#pragma once

#include <boost/config.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/reverse.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/end.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include "IJRange.h"
#include "DummyContext.h"
#include "HasParameter.h"
#include "HasDo.h"
#include "Domains.h"
#include "KRange.h"
#include "StencilStageEnvironment.h"

#include <boost/mpl/map/map10.hpp> // depends on MAX_MERGE_COUNT

/**
* @struct StencilSubStage
* Structure defining a stencil stage given the stencil stage definition functor and the associated ranges
*/
template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    typename TKRange> 
struct StencilStage
{
    // verify the ij and k range parameters 
    BOOST_STATIC_ASSERT(is_ij_range<TIJRange>::value);
    BOOST_STATIC_ASSERT(is_k_range<TKRange>::value);

    typedef TIJRange IJRange;
    typedef TKRange KRange;
};

/**
* @struct is_stencil_stage
* Meta function returning true if the parameter is a stencil stage
*/
template<typename T>
struct is_stencil_stage : boost::mpl::false_ {};

template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    typename TKRange>
struct is_stencil_stage<StencilStage<TStencilStageDefinition, TIJRange, TKRange> > : boost::mpl::true_ {};

/**
* @struct stencil_stage_definition
* Meta function returning the stencil stage definition for a given context and domain
*/
template<
    typename TStencilStage,
    typename TContext,
    typename TDomain>
struct stencil_stage_definition;

template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    typename TKRange,
    typename TContext,
    typename TDomain>
struct stencil_stage_definition<StencilStage<TStencilStageDefinition, TIJRange, TKRange>, TContext, TDomain> 
{
    // define a stencil stage environment given a context and a domain
    typedef TStencilStageDefinition<StencilStageEnvironment<TContext, TDomain> > type;
};

/**
* @struct stencil_stage_ij_range
* Meta function returning the ij range of a stencil stage 
*/
template<typename TStencilStage>
struct stencil_stage_ij_range; 

template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    typename TKRange>
struct stencil_stage_ij_range<StencilStage<TStencilStageDefinition, TIJRange, TKRange> > 
{
    typedef TIJRange type;
};

/**
* @struct stencil_stage_k_range
* Meta function returning the k range of a stencil stage 
*/
template<typename TStencilStage>
struct stencil_stage_k_range; 

template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    typename TKRange>
struct stencil_stage_k_range<StencilStage<TStencilStageDefinition, TIJRange, TKRange> > 
{
    typedef TKRange type;
};

/**
* @struct stencil_stage_has_parameter
* Meta function returning true if a stencil stage has a parameter with the given index and domain
*/
template<
    typename TStencilStage,
    typename TDomain,
    typename TParameterIndex>
struct stencil_stage_has_parameter 
{
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            has_parameter<
                typename stencil_stage_definition<TStencilStage, DummyContext, TDomain>::type,
                TParameterIndex
            >::value
        )
    );
    typedef boost::mpl::bool_<bool(value)> type;
};

/**
* @struct stencil_stage_has_do
* Meta function returning true if a stencil stage has a do method for a given domain
*/
template<
    typename TStencilStage,
    typename TDomain>
struct stencil_stage_has_do 
{
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            has_do<
                typename stencil_stage_definition<TStencilStage, DummyContext, FullDomain>::type, 
                void, 
                DummyContext, 
                TDomain
            >::value
        )
    );
    typedef boost::mpl::bool_<bool(value)> type;
};

/**
* @struct stencil_stage_lookup_do_domain
* Meta function searching the stencil stage do method which matches best a given domain
* The function starts the search with the parameter domain and the iterates over all base domains
* (Note that void is returned in case no matching do method is found)
*/
template<
    typename TStencilStage,
    typename TDomain>
struct stencil_stage_lookup_do_domain
{
    // compute vector containing all candidate domains, e.g. for TDomain=KMaximumCenter 
    // the vector will contain KMaximumCenter, TerrainCoordinates and FullDomain
    typedef typename boost::mpl::push_front<
        typename base_domains<TDomain>::type, 
        TDomain
    >::type BaseDomains;

    // find first base domain with a matching do method
    typedef typename boost::mpl::find_if<
        BaseDomains, 
        stencil_stage_has_do<TStencilStage, boost::mpl::_> 
    >::type DoDomainIterator;

    // if a do method was found return the domain otherwise return void
    typedef typename boost::mpl::eval_if<
        boost::is_same<DoDomainIterator, typename boost::mpl::end<BaseDomains>::type>,
        boost::mpl::void_,        
        boost::mpl::deref<DoDomainIterator>
    >::type type;
};

/**
* @struct stencil_stage_parameter_ij_range
* Meta function returning the stencil stage ij range if it uses a given parameter and the empty ij range otherwise
*/
template<
    typename TStencilStage,
    typename TParameterIndex>
struct stencil_stage_parameter_ij_range
{
    // if the stage uses the parameter return its ij range
    // otherwise return the empty range
    typedef typename boost::mpl::if_<
        stencil_stage_has_parameter<TStencilStage, FullDomain, TParameterIndex>, 
        typename stencil_stage_ij_range<TStencilStage>::type,
        IJRange<cIndented, 0, 0, 0, 0> 
    >::type type;
};

/**
* @struct stencil_stage_find_first_boundary_domain
* Meta function searching for the first matching stencil stage do method given a sequence of boundary domains
* (Note that void is returned in case no matching do method is found)
*/
template<
    typename TStencilStage,
    typename TBoundaryDomains>
struct stencil_stage_find_first_boundary_domain
{
    // find the first boundary level
    typedef typename boost::mpl::find_if<
        TBoundaryDomains,
        stencil_stage_has_do<TStencilStage, boost::mpl::_>
    >::type BoundaryDomainIterator;

    // if a boundary level was found return the domain level otherwise return void
    typedef typename boost::mpl::eval_if<
        boost::is_same<BoundaryDomainIterator, typename boost::mpl::end<TBoundaryDomains>::type>,
        boost::mpl::void_,        
        boost::mpl::deref<BoundaryDomainIterator>
    >::type type;
};

// meta functions working on multiple stencil stages 

/**
* @struct stencil_stages_maximum_ij_range
* Meta function returning the maximal ij range of all stencil stage ij ranges
*/
template<typename TStencilStages>
struct stencil_stages_maximum_ij_range
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TStencilStages>::value);
    typedef typename boost::mpl::fold<
        TStencilStages,
        IJRange<cIndented, 0, 0, 0, 0>,
        maximum_ij_range<
            boost::mpl::_1,
            stencil_stage_ij_range<boost::mpl::_2>
        >
    >::type type;
};

/**
* @struct stencil_stages_parameter_ij_range
* Meta function returning the maximal ij range of all stencil stages using a given parameter
*/
template<
    typename TStencilStages,
    typename TParameterIndex>
struct stencil_stages_parameter_ij_range
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TStencilStages>::value);
    typedef typename boost::mpl::fold<
        TStencilStages,
        IJRange<cIndented, 0, 0, 0, 0>,
        maximum_ij_range<
            boost::mpl::_1,
            stencil_stage_parameter_ij_range<boost::mpl::_2, TParameterIndex>
        >
    >::type type;
};

/**
* @struct stencil_stages_have_parameter
* Meta function returning true if the stencil stages use a given parameter
*/
template<
    typename TStencilStages,
    typename TDomain,
    typename TParameterIndex>
struct stencil_stages_have_parameter 
{
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            !boost::is_same<
                typename boost::mpl::find_if<
                    TStencilStages,
                    stencil_stage_has_parameter<boost::mpl::_, TDomain, TParameterIndex>
                >::type,
                typename boost::mpl::end<TStencilStages>::type
            >::value
        )
    );
    typedef boost::mpl::bool_<bool(value)> type;
};

/**
* @struct stencil_stages_first_parameter_usage
* Meta function returning the first stencil stage using a given parameter
* (if no stage uses the given parameter the first stage is returned)
*/
template<
    typename TStencilStages,
    typename TParameterIndex>
struct stencil_stages_first_parameter_usage
{
    // make sure the stencil stages vector contains at least one stage
    BOOST_STATIC_ASSERT(
        boost::mpl::is_sequence<TStencilStages>::value && 
        !boost::mpl::empty<TStencilStages>::value
    );

    // find the first parameter usage
    typedef typename boost::mpl::find_if<
        TStencilStages,
        stencil_stage_has_parameter<boost::mpl::_, FullDomain, TParameterIndex>
    >::type FirstUsageIterator;

    // if a usage was found return the stage the iterator is pointing to
    // otherwise return the first stage of the stencil stages vector
    typedef typename boost::mpl::eval_if<
        boost::is_same<FirstUsageIterator, typename boost::mpl::end<TStencilStages>::type>,
        boost::mpl::front<TStencilStages>,        
        boost::mpl::deref<FirstUsageIterator>
    >::type type;
};

/**
* @struct stencil_stages_last_parameter_usage
* Meta function returning the last stencil stage using a given parameter
* (if no stage uses the given parameter the last stage is returned)
*/
template<
    typename TStencilStages,
    typename TParameterIndex>
struct stencil_stages_last_parameter_usage
{
    // find the first usage for a reversed stencil stage list
    typedef typename stencil_stages_first_parameter_usage<
        typename boost::mpl::reverse<TStencilStages>::type,
        TParameterIndex
    >::type type;
};

