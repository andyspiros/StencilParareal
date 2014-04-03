#include <boost/type_traits/is_same.hpp>
#include "gtest/gtest.h"
#include "StencilFramework.h"

// test if base domain computation works
TEST(KRangeAndDomainUnittest, ComputeBaseDomain)
{   
    // test a number of domain combinations and check the right base domain is computed
    ASSERT_TRUE((
        boost::is_same<
            compute_base_domain<KMinimumCenter, FullDomain>::type,
            FullDomain
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            compute_base_domain<FlatCoordinates, TerrainCoordinates>::type,
            FullDomain
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            compute_base_domain<FlatCoordinates, KMaximumCenter>::type,
            FullDomain
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            compute_base_domain<KMaximumMinus2, KMaximumCenter>::type,
            TerrainCoordinates
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            compute_base_domain<KMinimumCenter, KMaximumFlatMinus2>::type,
            FlatCoordinates
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            compute_base_domain<KMaximumFlatMinus2, KMaximumFlatMinus2>::type,
            KMaximumFlatMinus2
        >::value
    ));
}

// test if in k range computation works
TEST(KRangeAndDomainUnittest, InKRange)
{
    // test 3d domains
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, TerrainCoordinates>::value));
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, FlatCoordinates>::value));
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, FullDomain>::value));
    ASSERT_TRUE((in_k_range<KRange<FlatCoordinates,0,0>, FlatCoordinates>::value));
    ASSERT_FALSE((in_k_range<KRange<FlatCoordinates,0,0>, TerrainCoordinates>::value));
    ASSERT_FALSE((in_k_range<KRange<FlatCoordinates,0,0>, FullDomain>::value));

    // test 2d domains
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, KMinimumCenter>::value));
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, KMaximumCenter>::value));
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, KMinimumPlus1>::value));
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, KMaximumMinus1>::value));
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, KMaximum<FlatCoordinates,-1> >::value));

    // should be outside
    ASSERT_FALSE((in_k_range<KRange<FullDomain,0,0>, KMinimumMinus1>::value));
    ASSERT_FALSE((in_k_range<KRange<FullDomain,0,0>, KMaximumPlus1>::value));

    // test offset
    ASSERT_TRUE((in_k_range<KRange<FullDomain,0,0>, KMaximumPlus1, boost::mpl::integral_c<int,-1> >::value));

    // test used defined rang
    ASSERT_TRUE((in_k_range<KRange<FlatCoordinates,1,-1>, KMinimumPlus1>::value));
    ASSERT_TRUE((in_k_range<KRange<FlatCoordinates,1,-1>, KMaximum<FlatCoordinates,-1> >::value));
    ASSERT_FALSE((in_k_range<KRange<FlatCoordinates,1,-1>, KMinimumMinus1>::value));
    ASSERT_FALSE((in_k_range<KRange<FlatCoordinates,1,-1>, KMinimumCenter>::value));
    ASSERT_FALSE((in_k_range<KRange<FlatCoordinates,1,-1>, KMaximum<FlatCoordinates,0> >::value));

    // test 2d domains with 2d ranges
    ASSERT_TRUE((in_k_range<KRange<KMinimumPlus1>, KMinimumPlus1>::value));
    ASSERT_TRUE((in_k_range<KRange<KMaximumFlatCenter>, KMaximumFlatCenter>::value));
    ASSERT_FALSE((in_k_range<KRange<KMinimumPlus1>, KMaximumFlatCenter>::value));
    ASSERT_FALSE((in_k_range<KRange<KMaximumFlatCenter>, KMinimumPlus1>::value));
}

template<typename TEnv>
struct DummyStencilStage
{
    STENCIL_STAGE(TEnv)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
};

// test if in k loop range computation works
TEST(KRangeAndDomainUnittest, InKLoopRange)
{
    // test 2d k loop ranges 
    ASSERT_FALSE((
        in_k_loop_range<
            KLoopRange<KRange<KMinimumMinus1>, KBoundarySize<0,0> >,
            KRange<FullDomain,0,0>
        >::value
    ));
    ASSERT_FALSE((
        in_k_loop_range<
            KLoopRange<KRange<KMinimumMinus1>, KBoundarySize<0,0> >,
            KRange<KMinimumMinus2>
        >::value
    ));
    ASSERT_TRUE((
        in_k_loop_range<
            KLoopRange<KRange<KMinimumMinus1>, KBoundarySize<0,0> >,
            KRange<FullDomain,-1,0>
        >::value
    ));
    ASSERT_TRUE((
        in_k_loop_range<
            KLoopRange<KRange<KMinimumMinus1>, KBoundarySize<0,0> >,
            KRange<KMinimumMinus1> 
        >::value
    ));

    // test 2d stencil stages
    ASSERT_FALSE((
        in_k_loop_range<
            KLoopRange<KRange<FullDomain,0,0>, KBoundarySize<0,0> >,
            KRange<KMinimumMinus1>
        >::value
    ));
    ASSERT_FALSE((
        in_k_loop_range<
            KLoopRange<KRange<KMinimumMinus2>, KBoundarySize<0,0> >,
            KRange<KMinimumMinus1>
        >::value
    ));
    ASSERT_TRUE((
        in_k_loop_range<
            KLoopRange<KRange<FullDomain,-1,0>, KBoundarySize<0,0> >,
            KRange<KMinimumMinus1> 
        >::value
    ));
    ASSERT_TRUE((
        in_k_loop_range<
            KLoopRange<KRange<KMinimumMinus1>, KBoundarySize<0,0> >,
            KRange<KMinimumMinus1>
        >::value
    ));

    // test 3d k loop ranges and stencil stages
    ASSERT_TRUE((
        in_k_loop_range<
            KLoopRange<KRange<FullDomain,0,0>, KBoundarySize<0,0> >,
            KRange<TerrainCoordinates,0,0> 
        >::value
    ));
    ASSERT_FALSE((
        in_k_loop_range<
            KLoopRange<KRange<FlatCoordinates,0,0>, KBoundarySize<0,0> >,
            KRange<TerrainCoordinates,0,0>
        >::value
    ));
}

// test if merging k loop ranges works
TEST(KRangeAndDomainUnittest, MergeKLoopRanges)
{
    // test two identical ranges merge
    ASSERT_TRUE((
        boost::is_same<
            merge_k_loop_ranges<
                KLoopRange<KRange<KMaximumMinus1>, KBoundarySize<0,0> >,
                KLoopRange<KRange<KMaximumMinus1>, KBoundarySize<0,0> >
            >::type,
            KLoopRange<KRange<KMaximumMinus1>, KBoundarySize<0,0> > // -> return the same k loop range
        >::value
    ));

    // test merging of 3d domains with various k ranges and boundary sizes
    ASSERT_TRUE((
        boost::is_same<
            merge_k_loop_ranges<
                KLoopRange<KRange<TerrainCoordinates,-1,-2>, KBoundarySize<0,1> >,
                KLoopRange<KRange<TerrainCoordinates,1,0>, KBoundarySize<1,2> >
            >::type,
            KLoopRange<KRange<TerrainCoordinates,-1,0>, KBoundarySize<3,3> > // --> maximize range and adapt boundaries so that all boundary levels are contained
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            merge_k_loop_ranges<
                KLoopRange<KRange<FullDomain,0,1>, KBoundarySize<0,1> >,
                KLoopRange<KRange<FullDomain,-1,-1>, KBoundarySize<2,0> >
            >::type,
            KLoopRange<KRange<FullDomain,-1,1>, KBoundarySize<2,2> > // --> increase the boundary size up to the start of the smaller k loop range
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            merge_k_loop_ranges<
                KLoopRange<KRange<FlatCoordinates,1,0>, KBoundarySize<0,1> >,
                KLoopRange<KRange<FlatCoordinates,-2,0>, KBoundarySize<0,2> >
            >::type,
            KLoopRange<KRange<FlatCoordinates,-2,0>, KBoundarySize<3,2> > // note we need the minimum boundary levels in order to know where to start with the first loop 
        >::value
    ));
}

// test if k loop range extension works
TEST(KRangeAndDomainUnittest, ExtendKLoopRanges)
{
    // try to extend full domains by a k minimum range either inside or outside the k loop range
    ASSERT_TRUE((
        boost::is_same<
            extend_k_loop_range<
                KLoopRange<KRange<FullDomain,1,0>, KBoundarySize<0,0> >,
                KLoopRange<KRange<KMinimumCenter>, KBoundarySize<0,0> > // minimum level right before the full domain
            >::type,
            KLoopRange<KRange<FullDomain,0,0>, KBoundarySize<1,0> > // extend range and create a boundary of size one
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            extend_k_loop_range<
                KLoopRange<KRange<FlatCoordinates,1,0>, KBoundarySize<0,0> >,
                KLoopRange<KRange<KMinimumMinus1,0,0> , KBoundarySize<0,0> > // minimum level 2 levels before the flat coordinates
            >::type,
            KLoopRange<KRange<FlatCoordinates,-1,0>, KBoundarySize<2,0> > // minimum levels up the the start of the other range needed
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            extend_k_loop_range<
                KLoopRange<KRange<TerrainCoordinates,0,0>, KBoundarySize<1,0> >,
                KLoopRange<KRange<KMinimumTerrainPlus2,0,0> , KBoundarySize<0,0> > // minimum level inside the terrain coordinates
            >::type,
            KLoopRange<KRange<TerrainCoordinates,0,0>, KBoundarySize<3,0> > // extend the boundary size 1 -> 3
        >::value
    ));

    // try to extend full domains by a k maximum range either inside or outside the k loop range
    ASSERT_TRUE((
        boost::is_same<
            extend_k_loop_range<
                KLoopRange<KRange<FullDomain,1,0>, KBoundarySize<0,0> >,
                KLoopRange<KRange<KMaximumCenter>, KBoundarySize<0,0> > // maximum level right insight full domain
            >::type,
            KLoopRange<KRange<FullDomain,1,0>, KBoundarySize<0,1> > // --> extend boundary size to 1
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            extend_k_loop_range<
                KLoopRange<KRange<FlatCoordinates,0,-1>, KBoundarySize<1,1> >,
                KLoopRange<KRange<KMaximumFlatCenter,0,0> , KBoundarySize<0,0> > // maximum level right after the flat coordinates range
            >::type,
            KLoopRange<KRange<FlatCoordinates,0,0>, KBoundarySize<1,2> > // extend range and boundary
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            extend_k_loop_range<
                KLoopRange<KRange<TerrainCoordinates,0,1>, KBoundarySize<1,0> >,
                KLoopRange<KRange<KMaximumCenter,0,0> , KBoundarySize<0,0> > // minimum level inside the terrain coordinates
            >::type,
            KLoopRange<KRange<TerrainCoordinates,0,1>, KBoundarySize<1,2> > // extend the boundary size in order to cover k maximum center
        >::value
    ));
}
