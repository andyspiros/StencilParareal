#pragma once 

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/if.hpp>
#include "KPosition.h"

#include <boost/mpl/vector/vector10.hpp>

/**
* @struct FullDomain
* Tag used to identify update functions and parameters which are valid 
* in the complete domain.
*/
struct FullDomain {};

/**
* @struct TerrainCoordinates
* Tag used to identify update functions and parameters which are valid 
* in the area of terrain following coordinates.
*/
struct TerrainCoordinates {};

/**
* @struct FlatCoordinates
* Tag used to identify update functions and parameters which are valid 
* in the area of flat coordinates.
*/
struct FlatCoordinates {};

/**
* @struct KMaximum
* Tag used to identify update functions and parameters which are valid at the k maximum boundary
*/
template<
    typename TBaseDomain, 
    int VLevel>
struct KMaximum
{
    // make sure the base domain is valid
    BOOST_STATIC_ASSERT( 
        (
            boost::is_same<TBaseDomain, TerrainCoordinates>::value || 
            boost::is_same<TBaseDomain, FlatCoordinates>::value 
        )
    );
    // make sure the level is in range between -MAX_BOUNDARY_LEVELS and MAX_BOUNDARY_LEVELS
    // (additionally make sure terrain and flat coordinates do not overlap --> no levels > 0 at the flat boundary)
    BOOST_STATIC_ASSERT(
        (
            VLevel <= (boost::is_same<TBaseDomain, FlatCoordinates>::value ? 0 : MAX_BOUNDARY_LEVELS) &&
            VLevel >= -MAX_BOUNDARY_LEVELS
        )        
    );
};

/**
* @struct KMinimum
* Tag used to identify update functions and parameters which are valid at the k minimum boundary
*/
template<
    typename TBaseDomain, 
    int VLevel>
struct KMinimum
{
    // make sure the base domain is valid
    BOOST_STATIC_ASSERT( 
        (
            boost::is_same<TBaseDomain, TerrainCoordinates>::value || 
            boost::is_same<TBaseDomain, FlatCoordinates>::value 
        )
    );
    // make sure the level is in range between -MAX_BOUNDARY_LEVELS and MAX_BOUNDARY_LEVELS
    // (additionally make sure terrain and flat coordinates do not overlap --> no levels > 0 at the flat boundary)
    BOOST_STATIC_ASSERT(
        (
            VLevel <= MAX_BOUNDARY_LEVELS &&
            VLevel >= (boost::is_same<TBaseDomain, TerrainCoordinates>::value ? 0 : -MAX_BOUNDARY_LEVELS)
        )        
    );
};

// predefined levels for k maximum and minimum
// (depends on MAX_BOUNDARY_LEVELS)
typedef KMinimum<FlatCoordinates,-2> KMinimumMinus2;
typedef KMinimum<FlatCoordinates,-1> KMinimumMinus1;
typedef KMinimum<FlatCoordinates, 0> KMinimumCenter;
typedef KMinimum<FlatCoordinates, 1> KMinimumPlus1;
typedef KMinimum<FlatCoordinates, 2> KMinimumPlus2;
// flat coordinates
typedef KMaximum<FlatCoordinates, -2> KMaximumFlatMinus2;
typedef KMaximum<FlatCoordinates, -1> KMaximumFlatMinus1;
typedef KMaximum<FlatCoordinates, 0> KMaximumFlatCenter;
typedef KMinimum<TerrainCoordinates, 0> KMinimumTerrainCenter;
typedef KMinimum<TerrainCoordinates, 1> KMinimumTerrainPlus1;
typedef KMinimum<TerrainCoordinates, 2> KMinimumTerrainPlus2;
// terrain coordinates
typedef KMaximum<TerrainCoordinates,-2> KMaximumMinus2;
typedef KMaximum<TerrainCoordinates,-1> KMaximumMinus1;
typedef KMaximum<TerrainCoordinates, 0> KMaximumCenter;
typedef KMaximum<TerrainCoordinates, 1> KMaximumPlus1;
typedef KMaximum<TerrainCoordinates, 2> KMaximumPlus2;

// define vectors containing all minimum and maximum levels
// (depends on MAX_BOUNDARY_LEVELS)
typedef boost::mpl::vector5<KMinimumMinus2, KMinimumMinus1, KMinimumCenter, KMinimumPlus1, KMinimumPlus2> FlatCoordinatesKMinimumLevels;
typedef boost::mpl::vector3<KMaximumFlatMinus2, KMaximumFlatMinus1, KMaximumFlatCenter> FlatCoordinatesKMaximumLevels;

typedef boost::mpl::vector3<KMinimumTerrainCenter, KMinimumTerrainPlus1, KMinimumTerrainPlus2> TerrainCoordinatesKMinimumLevels;
typedef boost::mpl::vector5<KMaximumMinus2, KMaximumMinus1, KMaximumCenter, KMaximumPlus1, KMaximumPlus2> TerrainCoordinatesKMaximumLevels;

/**
* @struct domain_k_minimum_levels
* Meta function returning a vector containing all domain minimum levels
*/
template<typename TDomain>
struct domain_k_minimum_levels;

template<>
struct domain_k_minimum_levels<FlatCoordinates> : FlatCoordinatesKMinimumLevels {};

template<>
struct domain_k_minimum_levels<TerrainCoordinates> : TerrainCoordinatesKMinimumLevels {};

template<>
struct domain_k_minimum_levels<FullDomain> : FlatCoordinatesKMinimumLevels {};

/**
* @struct domain_k_maximum_levels
* Meta function returning a vector containing all domain maximum levels
*/
template<typename TDomain>
struct domain_k_maximum_levels;

template<>
struct domain_k_maximum_levels<FlatCoordinates> : FlatCoordinatesKMaximumLevels {};

template<>
struct domain_k_maximum_levels<TerrainCoordinates> : TerrainCoordinatesKMaximumLevels {};

template<>
struct domain_k_maximum_levels<FullDomain> : TerrainCoordinatesKMaximumLevels {};

/**
* @struct is_2d_domain
* Meta function returning true if the parameter is a 2d domain object
*/
template<typename T>
struct is_2d_domain : boost::mpl::false_ {};

template<typename TBaseDomain, int VLevel>
struct is_2d_domain<KMaximum<TBaseDomain, VLevel> > : boost::mpl::true_ {};

template<typename TBaseDomain, int VLevel>
struct is_2d_domain<KMinimum<TBaseDomain, VLevel> > : boost::mpl::true_ {};

/**
* @struct is_3d_domain
* Meta function returning true if the parameter is a 3d domain object.
*/
template<typename T>
struct is_3d_domain : boost::mpl::false_ {};

template<>
struct is_3d_domain<FullDomain> : boost::mpl::true_ {};

template<>
struct is_3d_domain<FlatCoordinates> : boost::mpl::true_ {};

template<>
struct is_3d_domain<TerrainCoordinates> : boost::mpl::true_ {};

/**
* @struct is_domain
* Meta function returning true if the parameter is a valid domain object
*/
template<typename T>
struct is_domain 
{
    BOOST_STATIC_CONSTANT(bool, value = (is_3d_domain<T>::value || is_2d_domain<T>::value));
    typedef boost::mpl::bool_<bool(value)> type;
};

/**
* @struct base_domains
* Meta function computing a vector containing all the base domains of domain
*/
template<typename T>
struct base_domains;

template<>
struct base_domains<FullDomain> : boost::mpl::vector0<> {};

template<>
struct base_domains<TerrainCoordinates> : boost::mpl::vector1<FullDomain> {};

template<>
struct base_domains<FlatCoordinates> : boost::mpl::vector1<FullDomain> {};

template<typename TBaseDomain, int VLevel>
struct base_domains<KMaximum<TBaseDomain, VLevel> > : boost::mpl::push_front<typename base_domains<TBaseDomain>::type, TBaseDomain> {};

template<typename TBaseDomain, int VLevel>
struct base_domains<KMinimum<TBaseDomain, VLevel> > : boost::mpl::push_front<typename base_domains<TBaseDomain>::type, TBaseDomain> {};

/**
* @struct is_sub_domain
* Meta function returning true if the sub domain is contained inside the domain
*/
template<
    typename TSubDomain,
    typename TDomain>
struct is_sub_domain : boost::mpl::contains<typename base_domains<TSubDomain>::type, TDomain> {};

// return true if sub domain and domain are equivalent
template<typename TDomain>
struct is_sub_domain<TDomain, TDomain> : boost::mpl::true_ {};

/**
* @struct domain_level
* Meta function returning the level of a domain for 3d domains 0 is returned
*/
template<typename TDomain>
struct domain_level : boost::mpl::integral_c<int, 0> 
{ 
    BOOST_STATIC_ASSERT(is_3d_domain<TDomain>::value);
};

template<
    typename TBaseDomain, 
    int VLevel>
struct domain_level<KMinimum<TBaseDomain, VLevel> > : boost::mpl::integral_c<int, VLevel> {};

template<
    typename TBaseDomain, 
    int VLevel>
struct domain_level<KMaximum<TBaseDomain, VLevel> > : boost::mpl::integral_c<int, VLevel> {};

/**
* @struct domain_begin
* Meta function returning the k position at the beginning of the domain
* (domain_begin and domain_end define the half open interval around the domain)
*/
template<typename TDomain>
struct domain_begin;

/**
* @struct domain_end
* Meta function returning the k position after the end of the domain
* (domain_begin and domain_end define the half open interval around the domain)
*/
template<typename TDomain>
struct domain_end;

// define 3d domain begin values
template<>
struct domain_begin<FullDomain>
{
    typedef KPosition<cKMinimumFlat, 0> type;
};
template<>
struct domain_begin<TerrainCoordinates>
{
    typedef KPosition<cKMinimumTerrain, 0> type;
};
template<>
struct domain_begin<FlatCoordinates>
{
    typedef KPosition<cKMinimumFlat, 0> type;
};

// define 3d domain end values
template<>
struct domain_end<FullDomain>
{
    typedef KPosition<cKMaximumTerrain, 1> type;
};
template<>
struct domain_end<TerrainCoordinates>
{
    typedef KPosition<cKMaximumTerrain, 1> type;
};
template<>
struct domain_end<FlatCoordinates>
{
    typedef KPosition<cKMaximumFlat, 1> type;
};

// define 2d domain begin values
template<int VLevel>
struct domain_begin<KMaximum<TerrainCoordinates, VLevel> >
{
    typedef KPosition<cKMaximumTerrain, VLevel> type;
};
template<int VLevel>
struct domain_begin<KMaximum<FlatCoordinates, VLevel> >
{
    typedef KPosition<cKMaximumFlat, VLevel> type;
};
template<int VLevel>
struct domain_begin<KMinimum<TerrainCoordinates, VLevel> >
{
    typedef KPosition<cKMinimumTerrain, VLevel> type;
};
template<int VLevel>
struct domain_begin<KMinimum<FlatCoordinates, VLevel> >
{
    typedef KPosition<cKMinimumFlat, VLevel> type;
};

// define 2d domain end values
template<int VLevel>
struct domain_end<KMaximum<TerrainCoordinates, VLevel> >
{
    typedef KPosition<cKMaximumTerrain, VLevel + 1> type;
};
template<int VLevel>
struct domain_end<KMaximum<FlatCoordinates, VLevel> >
{
    typedef KPosition<cKMaximumFlat, VLevel + 1> type;
};
template<int VLevel>
struct domain_end<KMinimum<TerrainCoordinates, VLevel> >
{
    typedef KPosition<cKMinimumTerrain, VLevel + 1> type;
};
template<int VLevel>
struct domain_end<KMinimum<FlatCoordinates, VLevel> >
{
    typedef KPosition<cKMinimumFlat, VLevel + 1> type;
};

/**
* @struct compute_base_domain
* Meta function computing the smallest base domain covering two given input domains
*/
template<
    typename TDomain1,
    typename TDomain2>
struct compute_base_domain
{
    // find the first base domain of TDomain1 which is a sub domain of TDomain2
    // note that TDomain1 is added to the front of the base domain vector as it is a candidate as well
    typedef typename boost::mpl::deref<
        typename boost::mpl::find_if<
            typename boost::mpl::push_front<
                typename base_domains<TDomain1>::type,
                TDomain1                
            >::type,
            is_sub_domain<TDomain2, boost::mpl::_>
        >::type
    >::type type;
};

/**
* @struct compute_k_minimum_boundary_domain
* Meta function computing a k minimum boundary domain given a base domain and a level
*/
template<
    typename TBaseDomain,
    typename TLevel>
struct compute_k_minimum_boundary_domain
{
    BOOST_STATIC_ASSERT(is_3d_domain<TBaseDomain>::value);

    typedef typename boost::mpl::if_<
        boost::is_same<FullDomain, TBaseDomain>,
        KMinimum<FlatCoordinates, TLevel::value>,
        KMinimum<TBaseDomain, TLevel::value>
    >::type type;
};

/**
* @struct compute_k_maximum_boundary_domain
* Meta function computing a k maximum boundary domain given a base domain and a level
*/
template<
    typename TBaseDomain,
    typename TLevel>
struct compute_k_maximum_boundary_domain
{
    BOOST_STATIC_ASSERT(is_3d_domain<TBaseDomain>::value);

    typedef typename boost::mpl::if_<
        boost::is_same<FullDomain, TBaseDomain>,
        KMaximum<TerrainCoordinates, TLevel::value>,
        KMaximum<TBaseDomain, TLevel::value>
    >::type type;
};