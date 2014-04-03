#pragma once

#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/vector/vector0.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/next.hpp>
#include <boost/mpl/void.hpp>
#include "DummyContext.h"
#include "StencilFunctionEnvironment.h"
#include "Definitions.h"
#include "HasParameter.h"

/**
* @struct count_stencil_function_parameters
* Meta function returning the number of parameters of a stencil function
* (note that void is returned in case the parameters are not ok)
*/
template<template <typename> class TStencilFunction>
struct count_stencil_function_parameters
{
    // define stencil function with a standard environment
    typedef TStencilFunction<
        StencilFunctionEnvironment<
            DummyContext,
            boost::mpl::vector0<>
        >
    > StencilFunction;  

    // compute the number of parameters
    // return void in case the parameter numbers are not incremented continuously
    typedef typename boost::mpl::fold<
        boost::mpl::range_c<int, 0, MAX_PARAM_COUNT>,
        boost::mpl::integral_c<int, 0>,
        boost::mpl::if_<
            boost::is_same<boost::mpl::_1, boost::mpl::_2>, 
            boost::mpl::if_<
                has_parameter<StencilFunction, boost::mpl::_2>,
                boost::mpl::next<boost::mpl::_1>, // if there is a parameter increment the counter
                boost::mpl::_1 // otherwise set number of valid parameters to the current number
            >,
            boost::mpl::if_<
                has_parameter<StencilFunction, boost::mpl::_2>,
                boost::mpl::void_, // if there is a parameter lager than number of valid parameters return void
                boost::mpl::_1               
            >
        >
    >::type type;
};





