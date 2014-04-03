#pragma once

#include <boost/config.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/count.hpp>
#include "Enums.h"

#include <boost/mpl/vector/vector10_c.hpp>

// define the storage orders
namespace StorageOrder {
    // 1D orders
    typedef boost::mpl::vector1_c<Dimension, cDimI> I;
    typedef boost::mpl::vector1_c<Dimension, cDimJ> J;
    typedef boost::mpl::vector1_c<Dimension, cDimK> K;

    // 2D orders
    typedef boost::mpl::vector2_c<Dimension, cDimI, cDimJ> JI;
    typedef boost::mpl::vector2_c<Dimension, cDimI, cDimK> KI;
    typedef boost::mpl::vector2_c<Dimension, cDimJ, cDimI> IJ;
    typedef boost::mpl::vector2_c<Dimension, cDimJ, cDimK> KJ;
    typedef boost::mpl::vector2_c<Dimension, cDimK, cDimI> IK;
    typedef boost::mpl::vector2_c<Dimension, cDimK, cDimJ> JK;

    // 3D orders
    typedef boost::mpl::vector3_c<Dimension, cDimI, cDimJ, cDimK> KJI;
    typedef boost::mpl::vector3_c<Dimension, cDimI, cDimK, cDimJ> JKI;
    typedef boost::mpl::vector3_c<Dimension, cDimJ, cDimI, cDimK> KIJ;
    typedef boost::mpl::vector3_c<Dimension, cDimJ, cDimK, cDimI> IKJ;
    typedef boost::mpl::vector3_c<Dimension, cDimK, cDimJ, cDimI> IJK;
    typedef boost::mpl::vector3_c<Dimension, cDimK, cDimI, cDimJ> JIK;
}

// is_data_field_storage_order implementation class
template<
    typename TStorageOrder, 
    bool VIsSequence>
struct is_data_field_storage_order_impl : boost::mpl::false_ {};

template<typename TStorageOrder>
struct is_data_field_storage_order_impl<TStorageOrder, true>
{
    typedef typename boost::mpl::size<TStorageOrder>::type Size;

    // compute dimension counts
    typedef typename boost::mpl::count<TStorageOrder, boost::mpl::integral_c<Dimension, cDimI> >::type ICount;
    typedef typename boost::mpl::count<TStorageOrder, boost::mpl::integral_c<Dimension, cDimJ> >::type JCount;
    typedef typename boost::mpl::count<TStorageOrder, boost::mpl::integral_c<Dimension, cDimK> >::type KCount;

    // check the storage order
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            (Size::value <= 3 && Size::value >= 0) &&
            (ICount::value <= 1 && KCount::value <= 1 && JCount::value <= 1) && 
            (ICount::value + JCount::value + KCount::value == Size::value)  
        )
    );
    typedef boost::mpl::bool_<bool(value)> type;
};

/**
* @struct is_data_field_storage_order
* Meta function returning true if the parameter is a valid storage order
*/
template<typename TStorageOrder>
struct is_data_field_storage_order : is_data_field_storage_order_impl<TStorageOrder, boost::mpl::is_sequence<TStorageOrder>::value> {};


  

