#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>
#include "YesNo.h"
#include "DataParameter.h"
#include "WithWrapper.h"
#include "StencilFunction.h"

// macro used to define a stencil function
#define STENCIL_FUNCTION(environment) \
    BOOST_STATIC_ASSERT(is_stencil_function_environment<environment>::value); \
    typedef environment StencilFunctionEnvironment; \
    typedef typename environment::ValueType T; \
    typedef typename environment::Context Context; 

// macro used for function parameter definition
// it redefines the parameter in case it is in the environment parameter list
// otherwise the parameter is -1 (used by function stencil wrapper which instantiates
// the function stencil with the dummy environment for memory usage reasons)
#define FUNCTION_PARAMETER(number, param) \
    typedef typename boost::mpl::eval_if_c< \
        (number < boost::mpl::size<typename StencilFunctionEnvironment::ParameterList>::value), \
        boost::mpl::at_c<typename StencilFunctionEnvironment::ParameterList, number>, \
        boost::mpl::identity<DataParameter<-1, Offset<0, 0, 0> > > \
    >::type param; \
    static yes Parameter(boost::mpl::integral_c<int, number>*); 



