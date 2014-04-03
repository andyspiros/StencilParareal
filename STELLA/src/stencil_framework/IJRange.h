#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>
#include "Enums.h"

/**
* @struct IJRange
* Structure storing the ij range of a stencil stage
*/
template<
    CornerPolicy VCornerPolicy,
    int VIMinusOffset,
    int VIPlusOffset,
    int VJMinusOffset,
    int VJPlusOffset>
struct IJRange
{
    BOOST_STATIC_ASSERT(VIMinusOffset <= 0 && VJMinusOffset <= 0);
    BOOST_STATIC_ASSERT(VIPlusOffset >= 0 && VJPlusOffset >= 0);

    // define corner policy
    typedef boost::mpl::integral_c< ::CornerPolicy, VCornerPolicy> CornerPolicy;

    // width in negative and positive i direction
    typedef boost::mpl::integral_c<int, VIMinusOffset> IMinusOffset;    
    typedef boost::mpl::integral_c<int, VIPlusOffset> IPlusOffset;
    
    // width in negative and positive j direction
    typedef boost::mpl::integral_c<int, VJMinusOffset> JMinusOffset;
    typedef boost::mpl::integral_c<int, VJPlusOffset> JPlusOffset;
};

/**
* @struct is_ij_boundary
* Meta function returning true if the parameter is a IJRange
*/
template<typename T>
struct is_ij_range : boost::mpl::false_ {};

template<
    CornerPolicy VCornerPolicy,
    int VIMinusOffset,
    int VIPlusOffset,
    int VJMinusOffset,
    int VJPlusOffset>
struct is_ij_range<IJRange<VCornerPolicy, VIMinusOffset, VIPlusOffset, VJMinusOffset, VJPlusOffset> > : boost::mpl::true_ {};

/**
* @struct maximum_ij_range
* Meta function computing the minimal IJRange containing both parameter IJRanges
*/
template<
    typename TIJRange1,
    typename TIJRange2>
struct maximum_ij_range
{
    typedef IJRange<
        (TIJRange1::CornerPolicy::value == cIndented && TIJRange2::CornerPolicy::value == cIndented ? cIndented : cComplete),
        (TIJRange1::IMinusOffset::value <= TIJRange2::IMinusOffset::value ? TIJRange1::IMinusOffset::value : TIJRange2::IMinusOffset::value),
        (TIJRange1::IPlusOffset::value >= TIJRange2::IPlusOffset::value ? TIJRange1::IPlusOffset::value : TIJRange2::IPlusOffset::value),
        (TIJRange1::JMinusOffset::value <= TIJRange2::JMinusOffset::value ? TIJRange1::JMinusOffset::value : TIJRange2::JMinusOffset::value),
        (TIJRange1::JPlusOffset::value >= TIJRange2::JPlusOffset::value ? TIJRange1::JPlusOffset::value : TIJRange2::JPlusOffset::value)
    > type;
};
  
