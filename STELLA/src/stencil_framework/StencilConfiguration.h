#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include "BlockSize.h"

/**
* @struct StencilConfiguration
* Structure holding global configuration information of a stencil.
* The block size defines the size of the individual columns the stencil application strategy works on.
* The value type is the type used to define local variables in the stencil
*/
template<
    typename TValue,
    typename TBlockSize>
struct StencilConfiguration
{
    BOOST_STATIC_ASSERT(is_block_size<TBlockSize>::value);

    typedef TValue ValueType;
    typedef TBlockSize BlockSize; 
};

/**
* @struct is_stencil_configuration
* Meta function returning true if the parameter is a StencilConfiguration
*/
template<typename T>
struct is_stencil_configuration : boost::mpl::false_ {};

template<
    typename TValue,
    typename TBlockSize>
struct is_stencil_configuration<StencilConfiguration<TValue, TBlockSize> > : boost::mpl::true_ {};

  
