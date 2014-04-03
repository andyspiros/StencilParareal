#pragma once

#include <boost/mpl/void.hpp>
#include "StencilFunction.h"
#include "Domains.h"
#include "WithWrapper.h"

/**
* @struct Call
* Class used to wrap a stencil function before handing it over to the context
* (Note that the main functionality is provided by the WithWrapper)
*/
template<
    template <typename> class TStencilFunction,
    typename TDomain = boost::mpl::void_>
struct Call : WithWrapper<TStencilFunction, TDomain, typename count_stencil_function_parameters<TStencilFunction>::type> {};



