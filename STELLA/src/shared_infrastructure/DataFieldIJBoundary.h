#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>

/**
* @struct DataFieldIJBoundary
* Structure defining the ij boundary size of a data field
*/
template<
    int VIMinusOffset,
    int VIPlusOffset,
    int VJMinusOffset,
    int VJPlusOffset>
struct DataFieldIJBoundary
{
    BOOST_STATIC_ASSERT(VIMinusOffset <= 0 && VJMinusOffset <= 0);
    BOOST_STATIC_ASSERT(VIPlusOffset >= 0 && VJPlusOffset >= 0);

    // width in negative and positive i direction
    typedef boost::mpl::integral_c<int, VIMinusOffset> IMinusOffset;    
    typedef boost::mpl::integral_c<int, VIPlusOffset> IPlusOffset;
    
    // width in negative and positive j direction
    typedef boost::mpl::integral_c<int, VJMinusOffset> JMinusOffset;
    typedef boost::mpl::integral_c<int, VJPlusOffset> JPlusOffset;
};

/**
* @struct is_data_field_ij_boundary
* Meta function returning true if the parameter is a DataFieldIJBoundary
*/
template<typename T>
struct is_data_field_ij_boundary : boost::mpl::false_ {};

template<
    int VIMinusOffset,
    int VIPlusOffset,
    int VJMinusOffset,
    int VJPlusOffset>
struct is_data_field_ij_boundary<DataFieldIJBoundary<VIMinusOffset, VIPlusOffset, VJMinusOffset, VJPlusOffset> > : boost::mpl::true_ {};

