#pragma once

#include <boost/config.hpp> 
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>

/**
* @struct ComparisonValue
* Structure storing a compile time integral value which can be used for comparison operations
*/
template<
    typename TValue,
    TValue VValue>
struct ComparisonValue 
{
    typedef TValue ValueType;

    BOOST_STATIC_CONSTANT(TValue, value = VValue);
    typedef boost::mpl::integral_c<TValue, VValue> type;
};

/**
* @struct is_comparison_value
* Meta function returning true if the parameter is a comparison value
*/
template<typename T>
struct is_comparison_value : boost::mpl::false_ {};

template<
    typename TValue,
    TValue VValue>
struct is_comparison_value<ComparisonValue<TValue, VValue> > : boost::mpl::true_ {};

/**
* @struct comparison_value_is_comparison_possible
* Meta function returning true if a comparison value can be compared to a given value type
*/
template<
    typename TComparisonValue,
    typename TValue>
struct comparison_value_is_comparison_possible : boost::mpl::false_ {};

template<
    typename TValue,
    TValue VValue>
struct comparison_value_is_comparison_possible<ComparisonValue<TValue, VValue>, TValue> : boost::mpl::true_ {};