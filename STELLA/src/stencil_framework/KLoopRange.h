#pragma once

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/mpl/erase_key.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/prior.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/minus.hpp>
#include <boost/mpl/min.hpp>
#include <boost/mpl/max.hpp>
#include "Enums.h"
#include "Domains.h"
#include "KPosition.h"
#include "KRange.h"
#include "Cache.h"
#include "KBoundarySize.h"
#include "StencilStage.h"

// there are at most 6 * MAX_BOUNDARY_LEVELS 2d domains plus 3 3d domains
#include <boost/mpl/vector/vector20.hpp> // depends on MAX_BOUNDARY_LEVELS 
#include <boost/mpl/map/map20.hpp> // depends on MAX_BOUNDARY_LEVELS

/**
* @struct KLoopRange
* Structure storing the range and boundary size of a k loop
*/
template<
    typename TKRange,
    typename TKBoundarySize>
struct KLoopRange
{
    // check the parameter types
    BOOST_STATIC_ASSERT(is_k_range<TKRange>::value);
    BOOST_STATIC_ASSERT(is_k_boundary_size<TKBoundarySize>::value);

    // check that boundary is 0,0 for non 3d domain
    BOOST_STATIC_ASSERT( 
        (TKBoundarySize::KMinimumSize::value == 0 && TKBoundarySize::KMaximumSize::value == 0) ||
        (is_3d_domain<typename TKRange::Domain>::value)
    );

    // define range and boundary size
    typedef TKRange KRange;
    typedef TKBoundarySize KBoundarySize;
};

// meta function which compute basic k loop range features / traits

/**
* @struct is_k_loop_range
* Meta function returning true if the parameter is a k loop range
*/
template<typename TKLoopRange>
struct is_k_loop_range : boost::mpl::false_ {};

template<
    typename TKRange,
    typename TKBoundarySize>
struct is_k_loop_range<KLoopRange<TKRange, TKBoundarySize> > : boost::mpl::true_ {};

/**
* @struct k_loop_range_k_maximum
* Meta function computing k loop range maximum k position
*/
template<typename TKLoopRange>
struct k_loop_range_k_maximum;

template<    
    typename TKRange,
    typename TKBoundarySize>
struct k_loop_range_k_maximum<KLoopRange<TKRange, TKBoundarySize> >
{
    typedef typename shift_k_position_by<typename TKRange::End, -1>::type type;
};

/**
* @struct k_loop_range_k_minimum
* Meta function computing k loop range minimum k position
*/
template<typename TKLoopRange>
struct k_loop_range_k_minimum;

template<    
    typename TKRange,
    typename TKBoundarySize>
struct k_loop_range_k_minimum<KLoopRange<TKRange, TKBoundarySize> >
{
    typedef typename TKRange::Begin type;
};

/**
* @struct k_loop_range_domain
* Meta function computing the k loop range domain
*/
template<typename TKLoopRange>
struct k_loop_range_domain;

template<
    typename TKRange,
    typename TKBoundarySize>
struct k_loop_range_domain<KLoopRange<TKRange, TKBoundarySize> >
{
    typedef typename TKRange::Domain type;
};

/**
* @struct in_k_loop_range
* Meta function returning true if a given k range overlaps withs the k loop range
*/
template<
    typename TKLoopRange, 
    typename TKRange>
struct in_k_loop_range
{
    BOOST_STATIC_ASSERT(is_k_loop_range<TKLoopRange>::value);
    BOOST_STATIC_ASSERT(is_k_range<TKRange>::value);

    // extract ranges and domains
    typedef typename TKLoopRange::KRange KLoopRangeKRange;
    typedef typename KLoopRangeKRange::Domain KLoopRangeDomain;
    typedef typename TKRange::Domain KRangeDomain;
    
    // make sure the k range overlaps with the k loop range
    // 1. make sure one domain is a sub domain of the other one otherwise return false
    // 2. if the sub domain is a 3d domain we know that it overlaps for sure with the base domain return true
    // 3. if the sub domain is a 2d domain make sure that it is in the k range of the base domain using the in_k_range method
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            is_sub_domain<KRangeDomain, KLoopRangeDomain>::value ?
            (
                is_2d_domain<KRangeDomain>::value ? in_k_range<KLoopRangeKRange, KRangeDomain>::value : true
            ) :
            (
                is_sub_domain<KLoopRangeDomain, KRangeDomain>::value ?
                (
                    is_2d_domain<KLoopRangeDomain>::value ? in_k_range<TKRange, KLoopRangeDomain>::value : true
                ) :
                (
                    false // the domains do not overlap no sub domain found
                )
            ) 
        )
    );
    typedef boost::mpl::integral_c<bool, bool(value)> type;
};

// meta functions computing the k loop ranges needed to process a stencil stage

/**
* @struct create_k_loop_range
* Meta function creating a k loop range for a given stencil stage and k range
*/
template<
    typename TStencilStage,
    typename TKRange> 
struct create_k_loop_range
{
    // assert that a 3d domain is passed and that the stage provides a matching do method
    BOOST_STATIC_ASSERT(is_3d_domain<typename TKRange::Domain>::value);
        
    // define an open interval with a minimal and maximal boundary limit
    typedef typename boost::mpl::prior<typename TKRange::Begin::KOffset>::type KMinimumLimitLevel;
    typedef typename TKRange::End::KOffset KMaximumLimitLevel;

    // compute the last k minimum boundary domain (or void if there is none)
    typedef typename stencil_stage_find_first_boundary_domain<
        TStencilStage,
        typename boost::mpl::reverse_copy_if< // reverse the domain levels in order to find the last boundary domain
            typename domain_k_minimum_levels<typename TKRange::Domain>::type,
            boost::mpl::greater<
                domain_level<boost::mpl::_>, 
                KMinimumLimitLevel
            >
        >::type
    >::type KMinimumLastBoundaryDomain;
    
    // compute the first k maximum boundary domain (or void if there is none)
    typedef typename stencil_stage_find_first_boundary_domain<
        TStencilStage,
        typename boost::mpl::copy_if<
            typename domain_k_maximum_levels<typename TKRange::Domain>::type,
            boost::mpl::less<
                domain_level<boost::mpl::_>, 
                KMaximumLimitLevel
            >
        >::type
    >::type KMaximumFirstBoundaryDomain;

    // compute the k minimum boundary size given the last boundary level
    typedef typename boost::mpl::minus<
        typename boost::mpl::eval_if<
            boost::mpl::is_void_<KMinimumLastBoundaryDomain>,
            KMinimumLimitLevel, // if no boundary domain was found compute size 0
            domain_level<KMinimumLastBoundaryDomain>
        >::type,
        KMinimumLimitLevel
    >::type KMinimumBoundarySize;
      
    // compute the k minimum boundary size given the last boundary level
    typedef typename boost::mpl::minus<
        KMaximumLimitLevel,
        typename boost::mpl::eval_if<
            boost::mpl::is_void_<KMaximumFirstBoundaryDomain>,
            KMaximumLimitLevel, // if no boundary domain was found compute size 0
            domain_level<KMaximumFirstBoundaryDomain>
        >::type
    >::type KMaximumBoundarySize;

    // define the k loop range
    typedef KLoopRange<
        TKRange,
        KBoundarySize<
            KMinimumBoundarySize::value,
            KMaximumBoundarySize::value
        >
    > type;
};

/**
* @struct create_k_loop_ranges
* Meta function creating a vector containing all k loop ranges for a given stencil stage
*/
template<typename TStencilStage>
struct create_k_loop_ranges;

// specialization for flat or terrain coordinate loops
template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    typename TKRange> 
struct create_k_loop_ranges<StencilStage<TStencilStageDefinition, TIJRange, TKRange> >
{
    // define the stencil stage
    typedef StencilStage<
        TStencilStageDefinition, 
        TIJRange, 
        TKRange
    > StencilStageType;

    // make sure the k range spans the flat or the terrain coordinates
    BOOST_STATIC_ASSERT(
        (
            boost::is_same<FlatCoordinates, typename TKRange::Domain>::value ||
            boost::is_same<TerrainCoordinates, typename TKRange::Domain>::value
        )
    );
    // check that the stencil stage provides a do method matching the k range
    BOOST_MPL_ASSERT_MSG(
        (        
            boost::mpl::is_not_void_<
                typename stencil_stage_lookup_do_domain<
                    StencilStageType,
                    typename TKRange::Domain
                >::type 
            >::value
        ),
        STENCIL_STAGE_DEFINTION_NO_DO_FOR_THE_REQUESTED_KRANGE,
        (StencilStageType, TKRange)
    );

    // define the k loop range vector
    typedef boost::mpl::vector1<
        typename create_k_loop_range< 
            StencilStageType, 
            TKRange
        >::type
    > type;
};

// specialization for full domain loops handling the split into flat and terrain loops if necessary
template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    int VKMinimumOffset,
    int VKMaximumOffset>
struct create_k_loop_ranges<StencilStage<TStencilStageDefinition, TIJRange, KRange<FullDomain, VKMinimumOffset, VKMaximumOffset> > >
{
    // define the stencil stage
    typedef StencilStage<
        TStencilStageDefinition, 
        TIJRange, 
        KRange<FullDomain, VKMinimumOffset, VKMaximumOffset>
    > StencilStageType;

    // check that the stencil stage provides a do method matching the k range
    BOOST_MPL_ASSERT_MSG(
        (
            stencil_stage_has_do<StencilStageType, FullDomain>::value ^
            (
                stencil_stage_has_do<StencilStageType, FlatCoordinates>::value && 
                stencil_stage_has_do<StencilStageType, TerrainCoordinates>::value
            )            
        ),
        STENCIL_STAGE_DEFINTION_NO_DO_FOR_THE_REQUESTED_KRANGE,
        (StencilStageType, KRange<FullDomain, VKMinimumOffset, VKMaximumOffset>)
    );
  
    // distinguish between stages with a single full domain do method 
    // and stages which provide 2 do methods for terrain and flat domains
    // -> in case there are two do methods return two k loop ranges for the two domains
    typedef typename boost::mpl::transform<
        typename boost::mpl::if_<
            stencil_stage_has_do<StencilStageType, FullDomain>,
            boost::mpl::vector1<
                KRange<FullDomain, VKMinimumOffset, VKMaximumOffset> 
            >,
            boost::mpl::vector2<
                KRange<FlatCoordinates, VKMinimumOffset, 0>,
                KRange<TerrainCoordinates, 0, VKMaximumOffset>
            >
        >::type,
        create_k_loop_range<StencilStageType, boost::mpl::_>
    >::type type;
};

// specialization for k minimum levels
template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    typename TBaseDomain, 
    int VLevel>
struct create_k_loop_ranges<StencilStage<TStencilStageDefinition, TIJRange, KRange<KMinimum<TBaseDomain, VLevel>, 0, 0> > > 
{
    // define the stencil stage
    typedef StencilStage<
        TStencilStageDefinition, 
        TIJRange, 
        KRange<KMinimum<TBaseDomain, VLevel>, 0, 0>
    > StencilStageType;

    // check that the stencil stage provides a do method matching the k range
    BOOST_MPL_ASSERT_MSG(
        (        
            boost::mpl::is_not_void_<
                typename stencil_stage_lookup_do_domain<
                    StencilStageType, 
                    KMinimum<TBaseDomain, VLevel> 
                >::type 
            >::value
        ),
        STENCIL_STAGE_DEFINTION_NO_DO_FOR_THE_REQUESTED_KRANGE,
        (StencilStageType, KRange<KMinimum<TBaseDomain, VLevel>, 0, 0>)
    );

    // define the k loop range vector
    typedef boost::mpl::vector1<
        KLoopRange<
            KRange<KMinimum<TBaseDomain, VLevel>, 0, 0>, 
            KBoundarySize<0, 0> 
        > 
    > type;
};

// specialization for k maximum levels
template<
    template<typename> class TStencilStageDefinition,
    typename TIJRange,
    typename TBaseDomain, 
    int VLevel>
struct create_k_loop_ranges<StencilStage<TStencilStageDefinition, TIJRange, KRange<KMaximum<TBaseDomain, VLevel>, 0, 0> > > 
{
    // define the stencil stage
    typedef StencilStage<
        TStencilStageDefinition, 
        TIJRange, 
        KRange<KMaximum<TBaseDomain, VLevel>, 0, 0>
    > StencilStageType;

    // check that the stencil stage provides a do method matching the k range
    BOOST_MPL_ASSERT_MSG(
        (        
            boost::mpl::is_not_void_<
                typename stencil_stage_lookup_do_domain<
                    StencilStageType, 
                    KMaximum<TBaseDomain, VLevel> 
                >::type 
            >::value
        ),
        STENCIL_STAGE_DEFINTION_NO_DO_FOR_THE_REQUESTED_KRANGE,
        (StencilStageType, KRange<KMaximum<TBaseDomain, VLevel>, 0, 0>)
    );

    // define the k loop range vector
    typedef boost::mpl::vector1<
        KLoopRange<
            KRange<KMaximum<TBaseDomain, VLevel>, 0, 0>, 
            KBoundarySize<0, 0> 
        > 
    > type;
};

// meta functions useful for merging and extending k loop ranges

/**
* @struct compute_minimum_offset
* Meta function computing the minimum k loop range offset given two offsets
*/
template<
    int VKMinimumOffset1,
    int VKMinimumOffset2>
struct compute_minimum_offset : 
    boost::mpl::integral_c<int, (VKMinimumOffset1 < VKMinimumOffset2 ? VKMinimumOffset1 : VKMinimumOffset2)> 
{};

/**
* @struct compute_maximum_offset
* Meta function computing the maximum k loop range offset given two offsets
*/
template<
    int VKMaximumOffset1,
    int VKMaximumOffset2>
struct compute_maximum_offset : 
    boost::mpl::integral_c<int, (VKMaximumOffset1 > VKMaximumOffset2 ? VKMaximumOffset1 : VKMaximumOffset2)>
{};

/**
* @struct compute_minimum_boundary_size
* Meta function computing the minimum k loop range boundary size given two boundary sizes and minimum offsets
*/
template<
    int VKMinimumOffset1,
    int VKMinimumSize1,
    int VKMinimumOffset2,
    int VKMinimumSize2>
struct compute_minimum_boundary_size 
{
    // compute the absolute minimum boundary maximums
    typedef boost::mpl::integral_c<int, VKMinimumOffset1 + VKMinimumSize1> KMinimumBoundaryMaximum1;
    typedef boost::mpl::integral_c<int, VKMinimumOffset2 + VKMinimumSize2> KMinimumBoundaryMaximum2;

    // compute the boundary size maximal boundary level minus minimal offset
    // (note this is true for boundary size 0 as well as we need the boundary levels before the start of the smaller k range 
    // -> otherwise the framework would call the stencil stage associated to the smaller range too early)
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::minus<
                typename boost::mpl::max<KMinimumBoundaryMaximum1, KMinimumBoundaryMaximum2>::type,
                typename compute_minimum_offset<VKMinimumOffset1, VKMinimumOffset2>::type
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};
    
/**
* @struct compute_maximum_boundary_size
* Meta function computing the maximum k loop range boundary size given two boundary sizes and maximum offsets
*/
template<
    int VKMaximumOffset1,
    int VKMaximumSize1,
    int VKMaximumOffset2,
    int VKMaximumSize2>
struct compute_maximum_boundary_size 
{
    // compute the absolute maximum boundary minimums
    typedef boost::mpl::integral_c<int, VKMaximumOffset1 - VKMaximumSize1> KMaximumBoundaryMinimum1;
    typedef boost::mpl::integral_c<int, VKMaximumOffset2 - VKMaximumSize2> KMaximumBoundaryMinimum2;

    // compute the boundary size as maximal offset minus minimal boundary level
    // (note this is true for boundary size 0 as well as we need the boundary levels after the end of the smaller k range 
    // -> otherwise the framework would call the stencil stage associated to the smaller range too late)
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::minus<
                typename compute_maximum_offset<VKMaximumOffset1, VKMaximumOffset2>::type,
                typename boost::mpl::min<KMaximumBoundaryMinimum1, KMaximumBoundaryMinimum2>::type                
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct extend_k_loop_range
* Meta function extending a k loop range by another k loop range
* (there are specializations extending a 3d domain by a k minimum or maximum range)
*/
template<
    typename TKLoopRange1,
    typename TKLoopRange2>
struct extend_k_loop_range;

// specialization extending the minimum of a 3d domain
template<
    typename TDomain,
    typename TBaseDomainKMminimum,    
    int VKMinimumOffset,
    int VKMaximumOffset,
    int VKMinimumSize,
    int VKMaximumSize,
    int VLevelKMinimum>
struct extend_k_loop_range<
    KLoopRange<KRange<TDomain, VKMinimumOffset, VKMaximumOffset>, KBoundarySize<VKMinimumSize, VKMaximumSize> >,
    KLoopRange<KRange<KMinimum<TBaseDomainKMminimum, VLevelKMinimum>, 0, 0>, KBoundarySize<0, 0> > >
{
    // make sure the base domain is a sub domain of the 3d domain
    BOOST_STATIC_ASSERT((is_sub_domain<TBaseDomainKMminimum, TDomain>::value));

    // compute the extended k loop range adding the 2d domain
    typedef KLoopRange<
        KRange<
            TDomain, 
            compute_minimum_offset<VKMinimumOffset, VLevelKMinimum>::value, // include VLevel in the k loop range
            VKMaximumOffset
        >,
        KBoundarySize< 
            compute_minimum_boundary_size<VKMinimumOffset, VKMinimumSize, VLevelKMinimum, 1>::value, // make sure a boundary size of 1 is supported
            VKMaximumSize
        >
    > type;
};

// specialization extending the maximum of a 3d domain
template<
    typename TDomain,
    typename TBaseDomainKMaximum,    
    int VKMinimumOffset,
    int VKMaximumOffset,
    int VKMinimumSize,
    int VKMaximumSize,
    int VLevelKMaximum>
struct extend_k_loop_range<
    KLoopRange<KRange<TDomain, VKMinimumOffset, VKMaximumOffset>, KBoundarySize<VKMinimumSize, VKMaximumSize> >,
    KLoopRange<KRange<KMaximum<TBaseDomainKMaximum, VLevelKMaximum>, 0, 0>, KBoundarySize<0, 0> > >
{
    // make sure the base domain is a sub domain of the 3d domain
    BOOST_STATIC_ASSERT((is_sub_domain<TBaseDomainKMaximum, TDomain>::value));

    // compute the extended k loop range adding the 2d domain
    typedef KLoopRange<
        KRange<
            TDomain, 
            VKMinimumOffset,
            compute_maximum_offset<VKMaximumOffset, VLevelKMaximum>::value // include VLevel in the k loop range
        >,
        KBoundarySize< 
            VKMinimumSize,
            compute_maximum_boundary_size<VKMaximumOffset, VKMaximumSize, VLevelKMaximum, 1>::value // make sure a boundary size of 1 is supported
        >
    > type;
};

/**
* @struct split_k_loop_range
* Meta function splitting of a sub domain of a k loop range
* (there are specializations splitting the full domain)
*/
template<
    typename TKLoopRange,
    typename TSubDomain>
struct split_k_loop_range;

// split of the flat coordinates of a full domain k loop range
template<
    int VKMinimumOffset,
    int VKMaximumOffset,
    int VKMinimumSize,
    int VKMaximumSize>
struct split_k_loop_range<
    KLoopRange<KRange<FullDomain, VKMinimumOffset, VKMaximumOffset>, KBoundarySize<VKMinimumSize, VKMaximumSize> >,
    FlatCoordinates>
{
    typedef KLoopRange<
        KRange<FlatCoordinates, VKMinimumOffset, 0>,
        KBoundarySize<VKMinimumSize, 0>
    > type;
};

// split of the terrain coordinates of a full domain k loop range
template<
    int VKMinimumOffset,
    int VKMaximumOffset,
    int VKMinimumSize,
    int VKMaximumSize>
struct split_k_loop_range<
    KLoopRange<KRange<FullDomain, VKMinimumOffset, VKMaximumOffset>, KBoundarySize<VKMinimumSize, VKMaximumSize> >,
    TerrainCoordinates>
{
    typedef KLoopRange<
        KRange<TerrainCoordinates, 0, VKMaximumOffset>,
        KBoundarySize<0, VKMaximumSize>
    > type;
};

/**
* @struct merge_k_loop_ranges
* Meta function merging two k loop ranges which cover the same base domain
* (note that the resulting k loop range covers is the maximum of both domains providing the full boundary needed to process both k loop ranges)
*/
template<
    typename TKLoopRange1,
    typename TKLoopRange2>
struct merge_k_loop_ranges;
    
// specialization merging k loop ranges with different offsets and boundaries
template<
    typename TDomain,
    int VKMinimumOffset1,
    int VKMaximumOffset1,
    int VKMinimumSize1,
    int VKMaximumSize1,
    int VKMinimumOffset2,
    int VKMaximumOffset2,
    int VKMinimumSize2,
    int VKMaximumSize2>
struct merge_k_loop_ranges<
    KLoopRange<KRange<TDomain, VKMinimumOffset1, VKMaximumOffset1>, KBoundarySize<VKMinimumSize1, VKMaximumSize1> >,
    KLoopRange<KRange<TDomain, VKMinimumOffset2, VKMaximumOffset2>, KBoundarySize<VKMinimumSize2, VKMaximumSize2> > >
{
    // make sure this specialization is only used for 3d domains -> 2d domains have no boundaries and offsets
    BOOST_STATIC_ASSERT(is_3d_domain<TDomain>::value);
    
    // compute the merged range
    typedef KLoopRange<
        KRange<
            TDomain, 
            compute_minimum_offset<VKMinimumOffset1, VKMinimumOffset2>::value,
            compute_maximum_offset<VKMaximumOffset1, VKMaximumOffset2>::value
        >,
        KBoundarySize< 
            compute_minimum_boundary_size<VKMinimumOffset1, VKMinimumSize1, VKMinimumOffset2, VKMinimumSize2>::value,
            compute_maximum_boundary_size<VKMaximumOffset1, VKMaximumSize1, VKMaximumOffset2, VKMaximumSize2>::value
        >
    > type;
};

// specialization merging two identical k loop ranges 
// (note covers the case where two 2d k loop ranges are merged)
template<
    typename TDomain,
    int VKMinimumOffset,
    int VKMaximumOffset,
    int VKMinimumSize,
    int VKMaximumSize>
struct merge_k_loop_ranges<
    KLoopRange<KRange<TDomain, VKMinimumOffset, VKMaximumOffset>, KBoundarySize<VKMinimumSize, VKMaximumSize> >,
    KLoopRange<KRange<TDomain, VKMinimumOffset, VKMaximumOffset>, KBoundarySize<VKMinimumSize, VKMaximumSize> > >
{
    // return the input k loop range as we are merging two identical k loop ranges
    typedef KLoopRange<
        KRange<TDomain, VKMinimumOffset, VKMaximumOffset>, 
        KBoundarySize<VKMinimumSize, VKMaximumSize> 
    > type;
};

// meta functions computing the k loops needed to process a set of stencil stages

/**
* @struct k_loop_range_map_add_k_loop_range
* Add a k loop range to a k loop range map which associates a domain to a k loop range
* (if there is already a k loop range in the map covering the same domain the k loop ranges are merged)
*/
template<
    typename TKLoopRangeMap,
    typename TKLoopRange>
struct k_loop_range_map_add_k_loop_range
{
    // define the k loop range domain
    typedef typename k_loop_range_domain<TKLoopRange>::type Domain;

    // define the k loop range which shall be added to the map
    // (merge the k loop range with an existing k loop range if necessary)
    typedef typename boost::mpl::eval_if<
        boost::mpl::has_key<TKLoopRangeMap, Domain>,
        merge_k_loop_ranges<TKLoopRange, typename boost::mpl::at<TKLoopRangeMap, Domain>::type>,
        boost::mpl::identity<TKLoopRange>
    >::type KLoopRange;

    // insert the updated k loop range into the map
    typedef typename boost::mpl::insert<
        typename boost::mpl::erase_key<TKLoopRangeMap, Domain>::type,
        boost::mpl::pair<Domain, KLoopRange>
    >::type type;
};

/**
* @struct k_loop_range_map_add_stencil_stage
* Meta function adding all k loop ranges of a certain stencil stage to a k loop range map which associates a domain to a k loop range
*/
template<
    typename TKLoopRangeMap,
    typename TStencilStage>
struct k_loop_range_map_add_stencil_stage
{
    // insert k loops into the t k loop map
    typedef typename boost::mpl::fold<
        typename create_k_loop_ranges<TStencilStage>::type,
        TKLoopRangeMap,
        k_loop_range_map_add_k_loop_range<boost::mpl::_1, boost::mpl::_2>
    >::type type;
};

/**
* @struct k_loop_range_map_lookup_domains
* Meta function searching all k loop ranges matching a list of domains
*/
template<
    typename TKLoopRangeMap,
    typename TDomains>
struct k_loop_range_map_lookup_domains
{
    // for all domains present in the k loop range map return the associated k loop range
    typedef typename boost::mpl::fold<
        TDomains, 
        boost::mpl::vector0<>,
        boost::mpl::if_<
            boost::mpl::has_key<TKLoopRangeMap, boost::mpl::_2>, 
            boost::mpl::push_back<
                boost::mpl::_1, 
                boost::mpl::at<TKLoopRangeMap, boost::mpl::_2> 
            >,
            boost::mpl::_1
        >
    >::type type;
};

/**
* @struct k_loop_ranges_add_domain
* Meta function adding all k loop ranges associated to a certain domain to a k loop range vector
*/
template<
    typename TKLoopRanges,
    typename TDomain,
    typename TKLoopRangeMap>
struct k_loop_ranges_add_domain
{
    // intended for 3d domains
    BOOST_STATIC_ASSERT(is_3d_domain<TDomain>::value);

    // in case we compute the k loop ranges of terrain or flat coordinates 
    // check if there is a full domain k loop range which has to be split up
    typedef typename boost::mpl::eval_if_c<
        boost::mpl::has_key<TKLoopRangeMap, FullDomain>::value &&
        (
            boost::is_same<FlatCoordinates, TDomain>::value ||
            boost::is_same<TerrainCoordinates, TDomain>::value
        ),
        split_k_loop_range<
            typename boost::mpl::at<TKLoopRangeMap, FullDomain>::type, 
            TDomain
        >,
        boost::mpl::void_
    >::type SplitFullDomainKLoopRange;

    // if there is a split of add it to the k loop range map
    typedef typename boost::mpl::eval_if<
        boost::mpl::is_void_<SplitFullDomainKLoopRange>,
        boost::mpl::identity<TKLoopRangeMap>,
        k_loop_range_map_add_k_loop_range<TKLoopRangeMap, SplitFullDomainKLoopRange>
    >::type KLoopRangeMap;
    
    // define the main domain k loop range 
    // (set to void if there is no k loop range matching the domain)
    typedef typename boost::mpl::at<KLoopRangeMap, TDomain>::type MainKLoopRange;

    // define k loop ranges at the domain boundaries
    typedef typename boost::mpl::fold<
        typename domain_k_maximum_levels<TDomain>::type, 
        typename domain_k_minimum_levels<TDomain>::type,
        boost::mpl::push_back<boost::mpl::_1, boost::mpl::_2>
    >::type BoundaryDomains;
    typedef typename k_loop_range_map_lookup_domains<
        KLoopRangeMap, 
        BoundaryDomains
    >::type BoundaryKLoopRanges;

    // if there is a main k loop range merge in all boundary ranges
    typedef typename boost::mpl::eval_if<
        boost::mpl::is_void_<MainKLoopRange>,
        boost::mpl::void_,
        boost::mpl::fold<
            BoundaryKLoopRanges,
            MainKLoopRange,
            extend_k_loop_range<boost::mpl::_1, boost::mpl::_2>
        >
    >::type ExtendedMainKLoopRange;
    
    // if the main domain k loop range is void return the boundary k loop ranges
    // otherwise extend the domain k loop range by all the boundary k loop ranges
    typedef typename boost::mpl::eval_if<
        boost::mpl::is_void_<MainKLoopRange>,
        boost::mpl::fold<
            BoundaryKLoopRanges,
            TKLoopRanges,
            boost::mpl::push_back<boost::mpl::_1, boost::mpl::_2>
        >,            
        boost::mpl::push_back<
            TKLoopRanges, 
            ExtendedMainKLoopRange
        >
    >::type type;
};

/**
* @struct compute_k_loop_ranges
* Meta function computing the k loop ranges given the stencil stages of the sweep
*/
template<typename TStencilStages> 
struct compute_k_loop_ranges
{
    // define a map containing all loop domains used by the stencil stages
    typedef typename boost::mpl::fold<
        TStencilStages,
        boost::mpl::map0<>,
        k_loop_range_map_add_stencil_stage<boost::mpl::_1, boost::mpl::_2>
    >::type KLoopRangeMap;

    // define terrain and flat coordinate specific boundary domains
    typedef typename boost::mpl::fold<
        typename domain_k_maximum_levels<FlatCoordinates>::type,
        typename domain_k_minimum_levels<TerrainCoordinates>::type,
        boost::mpl::push_back<boost::mpl::_1, boost::mpl::_2>
    >::type BoundaryDomains;
    typedef typename k_loop_range_map_lookup_domains<
        KLoopRangeMap,
        BoundaryDomains
    >::type BoundaryKLoopRanges;

    // compute the list of base domains needed to run all stencil stages
    // if there are k loop ranges for flat or terrain coordinates work with flat and terrain coordinates
    // by default work with the full domain only
    typedef typename boost::mpl::if_c<
        (
            boost::mpl::has_key<KLoopRangeMap, TerrainCoordinates>::value ||
            boost::mpl::has_key<KLoopRangeMap, FlatCoordinates>::value ||
            !boost::mpl::empty<BoundaryKLoopRanges>::value 
        ),
        boost::mpl::vector2<FlatCoordinates, TerrainCoordinates>,
        boost::mpl::vector1<FullDomain>
    >::type Domains;
    
    // compute the ordered list of k loop ranges
    typedef typename boost::mpl::fold<
        Domains,
        boost::mpl::vector0<>,
        k_loop_ranges_add_domain<boost::mpl::_1, boost::mpl::_2, KLoopRangeMap>
    >::type type;
};

/**
* @struct compute_k_loop_ranges_base_domain
* Meta function computing the domain covering all k loops
*/
template<typename TKLoopRanges>
struct compute_k_loop_ranges_base_domain
{
    // make sure a sequence is passed
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TKLoopRanges>::value);

    // iterate over all ranges and compute the base domain
    // (first transform the k loop ranges vector into a domain vector)
    typedef typename boost::mpl::fold<
        typename boost::mpl::transform<
            TKLoopRanges,
            k_loop_range_domain<boost::mpl::_>
        >::type,
        boost::mpl::void_,
        boost::mpl::if_<
            boost::mpl::is_void_<boost::mpl::_1>,
            boost::mpl::_2,
            compute_base_domain<boost::mpl::_1, boost::mpl::_2>
        >
    >::type type;
};

/**
* @struct compute_extended_minimum_boundary_size
* Meta function extending the minimum boundary to the inside given the minimum offset of an inner k range
*/
template<
    int VKMinimumOffset1,
    int VKMinimumSize1,
    int VKMinimumOffset2>
struct compute_extended_minimum_boundary_size : 
    boost::mpl::integral_c<int, (VKMinimumOffset1 + VKMinimumSize1 < VKMinimumOffset2 ? VKMinimumOffset2 - VKMinimumOffset1 : VKMinimumSize1)> 
{};

/**
* @struct compute_extended_maximum_boundary_size
* Meta function extending the maximum boundary to the inside given the maximum offset of an inner k range
*/
template<
    int VKMaximumOffset1,
    int VKMaximumSize1,
    int VKMaximumOffset2>
struct compute_extended_maximum_boundary_size : 
    boost::mpl::integral_c<int, (VKMaximumOffset1 - VKMaximumSize1 > VKMaximumOffset2 ? VKMaximumOffset1 - VKMaximumOffset2 : VKMaximumSize1)> 
{};

/**
* @struct extend_k_loop_range_boundary
* Meta function extending the k loop range boundary given an inner k range
* (note that the function is used in order to adapt the k loop range boundary size to the caches)
*/
template<
    typename TKLoopRange,
    typename TKRange>
struct extend_k_loop_range_boundary
{
    // error while computing the k loop boundary levels necessary in order to update the caches
    // this can happen in case a cache covers only a sub domain of the k loop domain 
    // (e.g. if you loop over the full domain and specify a cache for the terrain or flat domain)
    BOOST_MPL_ASSERT_MSG( 
        false,
        K_LOOP_COMPUTATION_FAILED_WITH_INVALID_CACHE_K_RANGE,
        (TKLoopRange, TKRange) 
    );
};

// specialization updating minimum and maximum boundary sizes for an arbitrary 3d domain
template<
    typename TDomain,
    int KMinimumOffset1,
    int KMaximumOffset1,
    int KMinimumSize1,
    int KMaximumSize1,
    int KMinimumOffset2,
    int KMaximumOffset2>
struct extend_k_loop_range_boundary<
    KLoopRange<KRange<TDomain, KMinimumOffset1, KMaximumOffset1>, KBoundarySize<KMinimumSize1, KMaximumSize1> >, 
    KRange<TDomain, KMinimumOffset2, KMaximumOffset2> >
{
    BOOST_STATIC_ASSERT(is_3d_domain<TDomain>::value);

    typedef KLoopRange<
        KRange<TDomain, KMinimumOffset1, KMaximumOffset1>,
        KBoundarySize<
            compute_extended_minimum_boundary_size<KMinimumOffset1, KMinimumSize1, KMinimumOffset2>::value,
            compute_extended_maximum_boundary_size<KMaximumOffset1, KMaximumSize1, KMaximumOffset2>::value
        >
    > type;
};

// specialization updating the full domain minimum boundary
template<
    int KMinimumOffset1,
    int KMaximumOffset1,
    int KMinimumSize1,
    int KMaximumSize1,
    int KMinimumOffset2,
    int KMaximumOffset2>
struct extend_k_loop_range_boundary<
    KLoopRange<KRange<FlatCoordinates, KMinimumOffset1, KMaximumOffset1>, KBoundarySize<KMinimumSize1, KMaximumSize1> >, 
    KRange<FullDomain, KMinimumOffset2, KMaximumOffset2> >
{
    typedef KLoopRange<
        KRange<FlatCoordinates, KMinimumOffset1, KMaximumOffset1>,
        KBoundarySize<
            compute_extended_minimum_boundary_size<KMinimumOffset1, KMinimumSize1, KMinimumOffset2>::value,
            KMaximumSize1
        >
    > type;
};

// specialization updating the full domain maximum boundary
template<
    int KMinimumOffset1,
    int KMaximumOffset1,
    int KMinimumSize1,
    int KMaximumSize1,
    int KMinimumOffset2,
    int KMaximumOffset2>
struct extend_k_loop_range_boundary<
    KLoopRange<KRange<TerrainCoordinates, KMinimumOffset1, KMaximumOffset1>, KBoundarySize<KMinimumSize1, KMaximumSize1> >, 
    KRange<FullDomain, KMinimumOffset2, KMaximumOffset2> >
{
    typedef KLoopRange<
        KRange<TerrainCoordinates, KMinimumOffset1, KMaximumOffset1>,
        KBoundarySize<
            KMinimumSize1,
            compute_extended_maximum_boundary_size<KMaximumOffset1, KMaximumSize1, KMaximumOffset2>::value
        >
    > type;
};

// specialization for a k minimum domain k loop range
template<
    typename TKRange,
    typename TBaseDomain,
    int VLevel>
struct extend_k_loop_range_boundary<KLoopRange<KRange<KMinimum<TBaseDomain, VLevel>, 0, 0>, KBoundarySize<0, 0> >, TKRange>
{
    typedef KLoopRange<
        KRange<KMinimum<TBaseDomain, VLevel>, 0, 0>, 
        KBoundarySize<0, 0> 
    > type;
};

// specialization for a k maximum domain k loop range
template<
    typename TKRange,
    typename TBaseDomain,
    int VLevel>
struct extend_k_loop_range_boundary<KLoopRange<KRange<KMaximum<TBaseDomain, VLevel>, 0, 0>, KBoundarySize<0, 0> >, TKRange>
{
    typedef KLoopRange<
        KRange<KMaximum<TBaseDomain, VLevel>, 0, 0>, 
        KBoundarySize<0, 0> 
    > type;
};
