#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/is_sequence.hpp>

/**
* @struct Environment
* Structure holding the type information necessary for stencil function execution, e.g. context type and parameter list
*/
template<
    typename TContext,
    typename TParameterList> 
struct StencilFunctionEnvironment
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TParameterList>::value);
    
    typedef typename TContext::ValueType ValueType;
    typedef TContext& Context;
    typedef TParameterList ParameterList;
};

/**
* @struct is_stencil_function_environment
* Meta function returning true the parameter is a stencil function environment
*/
template<typename T>
struct is_stencil_function_environment : boost::mpl::false_ {};

template<
    typename TContext,
    typename TParameterList> 
struct is_stencil_function_environment<StencilFunctionEnvironment<TContext, TParameterList> > : boost::mpl::true_ {};
  
