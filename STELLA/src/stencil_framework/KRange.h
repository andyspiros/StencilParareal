#pragma once

#include <boost/config.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include "Enums.h"
#include "Definitions.h"
#include "Domains.h"
#include "KPosition.h"

/**
* @struct KRange
* Structure defining a half open k range associated to a domain plus additional k minimum and k maximum offsets
*/
template<
    typename TDomain,
    int VKMinimumOffset = 0,
    int VKMaximumOffset = 0>
struct KRange
{
    // check that the 2d domain offsets are 0
    BOOST_STATIC_ASSERT( 
        is_3d_domain<TDomain>::value ||
        (is_2d_domain<TDomain>::value && VKMinimumOffset == 0 && VKMaximumOffset == 0)
    );

    // make sure VKMinimumOffset respects MAX_BOUNDARY_LEVELS
    // additionally guarantee that terrain coordinates k ranges do not overlap with flat coordinates k ranges
    // (a k range can start between -MAX_BOUNDARY_LEVELS and MAX_BOUNDARY_LEVELS+1)
    BOOST_STATIC_ASSERT(
        (
            VKMinimumOffset >= (boost::is_same<TDomain, TerrainCoordinates>::value ? 0 : -MAX_BOUNDARY_LEVELS) && 
            VKMinimumOffset <= (MAX_BOUNDARY_LEVELS + 1) // it's fine to reduce the k range until there is no overlap with a boundary level
        )        
    );
    
    // make sure VKMaximumOffset respects MAX_BOUNDARY_LEVELS 
    // additionally guarantee that flat coordinates k ranges do not overlap with terrain coordinates k ranges
    // (a k range can start between -1-MAX_BOUNDARY_LEVELS and MAX_BOUNDARY_LEVELS)
    BOOST_STATIC_ASSERT(
        (
            VKMaximumOffset >= -(MAX_BOUNDARY_LEVELS + 1) && // it's fine to reduce the k range until there is no overlap with a boundary level
            VKMaximumOffset <= (boost::is_same<TDomain, FlatCoordinates>::value ? 0 : MAX_BOUNDARY_LEVELS) 
        )        
    );
   
    typedef TDomain Domain;
    typedef boost::mpl::integral_c<int, VKMinimumOffset> KMinimumOffset;
    typedef boost::mpl::integral_c<int, VKMaximumOffset> KMaximumOffset;

    // define domain and the half open interval around the domain
    typedef typename shift_k_position_by<typename domain_begin<TDomain>::type, VKMinimumOffset>::type Begin;
    typedef typename shift_k_position_by<typename domain_end<TDomain>::type, VKMaximumOffset>::type End;
};

/**
* @struct is_krange
* Meta function returning true if the parameter is a valid k range.
*/
template<typename T>
struct is_k_range : boost::mpl::false_ {};

template<typename TDomain, int VKMinimumOffset, int VKMaximumOffset>
struct is_k_range<KRange<TDomain, VKMinimumOffset, VKMaximumOffset> > : boost::mpl::true_ {};

/**
* @struct in_k_range
* Meta function returning true if a certain domain plus optional offset is inside a k range
*/
template<
    typename TKRange,
    typename TDomain,
    typename TDomainOffset = boost::mpl::integral_c<int, 0> >
struct in_k_range
{
    BOOST_STATIC_ASSERT(is_k_range<TKRange>::value);
    BOOST_STATIC_ASSERT(is_domain<TDomain>::value);

    // compute domain begin and end shifted by the domain offset
    typedef typename shift_k_position_by<
        typename domain_begin<TDomain>::type,
        TDomainOffset::value
    >::type DomainBegin;
    typedef typename shift_k_position_by<
        typename domain_end<TDomain>::type,
        TDomainOffset::value
    >::type DomainEnd;

    // if the domain is a 3d domain check if it is a sub domain of the k range domain
    // otherwise compare the compile time origin offsets to check if the domain is inside the k range
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            is_3d_domain<TDomain>::value ?
            is_sub_domain<TDomain, typename TKRange::Domain>::value :
            (
                (
                    k_position_compile_time_origin_offset<typename TKRange::Begin>::value <= 
                    k_position_compile_time_origin_offset<DomainBegin>::value
                ) && (
                    k_position_compile_time_origin_offset<typename TKRange::End>::value >= 
                    k_position_compile_time_origin_offset<DomainEnd>::value
                ) 
            )
        )
    );
    typedef boost::mpl::bool_<bool(value)> type;
};




