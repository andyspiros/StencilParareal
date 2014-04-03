#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>
#include "KRange.h"

/**
* @struct KWindow
* Structure defining the caching window in k direction
*/
template<
    int VKMinusOffset,
    int VKPlusOffset>
struct KWindow
{
    BOOST_STATIC_ASSERT(VKMinusOffset <= VKPlusOffset);
    
    typedef boost::mpl::integral_c<int, VKMinusOffset> KMinusOffset;    
    typedef boost::mpl::integral_c<int, VKPlusOffset> KPlusOffset;
};

/**
* @struct is_k_window
* Meta function returning true if the parameter is a cache window
*/
template<typename T>
struct is_k_window : boost::mpl::false_ {};

template<
    int VKMinusOffset,
    int VKPlusOffset>
struct is_k_window<KWindow<VKMinusOffset, VKPlusOffset> > : boost::mpl::true_ {};

