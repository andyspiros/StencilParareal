#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include "Domains.h"

/**
* @struct StencilStageEnvironment
* Structure holding the type information used for the stencil stage execution
* (e.g. context type and domain type)
*/
template<
    typename TContext,
    typename TDomain> 
struct StencilStageEnvironment
{
    BOOST_STATIC_ASSERT(is_domain<TDomain>::value == true);
 
    typedef typename TContext::ValueType ValueType;
    typedef TContext& Context;
    typedef TDomain Domain;
};

/**
* @struct is_stencil_stage_environment
* Meta function returning true the parameter is a stencil stage environment
*/
template<typename T>
struct is_stencil_stage_environment : boost::mpl::false_ {};

template<
    typename TContext,
    typename TDomain> 
struct is_stencil_stage_environment<StencilStageEnvironment<TContext, TDomain> > : boost::mpl::true_ {};
