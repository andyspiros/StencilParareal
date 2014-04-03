#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>

/**
* @struct KBoundarySize
* Structure defining how many k maximum and k minimum levels are supported by a loop
*/
template<
    int VKMinimumSize,
    int VKMaximumSize>
struct KBoundarySize
{
    // sizes cannot be negative
    BOOST_STATIC_ASSERT((VKMinimumSize >= 0) && (VKMaximumSize >= 0));

    typedef boost::mpl::integral_c<int, VKMinimumSize> KMinimumSize;
    typedef boost::mpl::integral_c<int, VKMaximumSize> KMaximumSize;
};

/**
* @struct is_kboundary_size
* Meta function returning true if the parameter is a valid k boundary size
*/
template<typename T>
struct is_k_boundary_size : boost::mpl::false_ {};

template<int VKMinimumSize, int VKMaximumSize>
struct is_k_boundary_size<KBoundarySize<VKMinimumSize, VKMaximumSize> > : boost::mpl::true_ {};
