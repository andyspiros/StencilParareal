#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>

/**
* @struct Alignment
* Structure defining the alignment properties of a data field.
* The properties are the dimension of the alignment and the size of the alignment/padding in number of elements (not bytes!)
*/
template<
    Dimension VDimension,
    int VElements>
struct DataFieldAlignment 
{
    BOOST_STATIC_ASSERT(VElements > 0); // make sure the alignment is at least 1

    typedef boost::mpl::integral_c< ::Dimension, VDimension> Dimension;
    typedef boost::mpl::integral_c<int, VElements> Elements;
};

/**
* @struct is_data_field_alignment
* Meta function returning true if the parameter an alignment type
*/
template<typename T>
struct is_data_field_alignment : boost::mpl::false_ {};

template<
    Dimension VDimension,
    int VElements>
struct is_data_field_alignment<DataFieldAlignment<VDimension, VElements> > : boost::mpl::true_ {};

