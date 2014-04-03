#include "gtest/gtest.h"
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>
#include "StencilFramework.h"

template<typename TEnv>
struct KLoopDescriptorStage1
{
    STENCIL_STAGE(TEnv)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
    __ACC__
    static void Do(Context ctx, KMinimumPlus2) {}
    __ACC__
    static void Do(Context ctx, KMinimumPlus1) {}   
    __ACC__
    static void Do(Context ctx, KMinimumCenter) {}
};

template<typename TEnv>
struct KLoopDescriptorStage2
{
    STENCIL_STAGE(TEnv)

    __ACC__
    static void Do(Context ctx, KMaximumPlus2) {}
    __ACC__
    static void Do(Context ctx, KMaximumCenter) {}
    __ACC__
    static void Do(Context ctx, KMaximumMinus1) {}
    __ACC__
    static void Do(Context ctx, TerrainCoordinates) {}
    __ACC__
    static void Do(Context ctx, KMinimumTerrainCenter) {}
    __ACC__
    static void Do(Context ctx, KMaximumFlatCenter) {}
    __ACC__
    static void Do(Context ctx, KMaximumFlatMinus2) {}
    __ACC__
    static void Do(Context ctx, FlatCoordinates) {}
    __ACC__
    static void Do(Context ctx, KMinimumMinus1) {}
};

template<typename TEnv>
struct KLoopDescriptorStage3
{
    STENCIL_STAGE(TEnv)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
};

// meta function verifying the stencil stages of a k loop leg
template<
    typename TKLoopLeg,
    typename TExpectedStencilStages>
struct check_k_loop_leg_stencil_stages;

template<
    typename TExpectedStencilStages,
    typename TKLoopRange,
    typename TCaches,
    typename TStencilStages,
    typename TAdvanceToKPosition,
    KLoopDirection VKLoopDirection>
struct check_k_loop_leg_stencil_stages<KLoopLeg<TKLoopRange, TCaches, TStencilStages, TAdvanceToKPosition, VKLoopDirection>, TExpectedStencilStages> : 
    boost::mpl::equal<TStencilStages, TExpectedStencilStages> 
{};

// make sure the loop computation works 
TEST(KLoopComputationUnittest, ComputeKLoops)
{
    // simple test case without custom boundary no split in flat and terrain
    typedef boost::mpl::vector1<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >
    > Stages1FullDomain; 
   
    typedef compute_k_loop_ranges<
        Stages1FullDomain
    >::type Stage1FullDomainKLoopRanges;
    typedef compute_k_loop_legs<
        Stage1FullDomainKLoopRanges,
        boost::mpl::void_,
        Stages1FullDomain,
        cKIncrement
    >::type Stage1FullDomainKLoopLegs;
       
    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<Stage1FullDomainKLoopRanges>::value == 1);

    // check the loop definition
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<Stage1FullDomainKLoopRanges, 0>::type,
                KLoopRange<KRange<FullDomain,0,0>, KBoundarySize<3,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<Stage1FullDomainKLoopLegs, 0>::type,
                Stages1FullDomain
            >::value
        )
    );

    // test case without custom boundary but with split 
    typedef boost::mpl::vector1<
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> > 
    > Stages2FullDomain; 
       
    typedef compute_k_loop_ranges<Stages2FullDomain>::type Stage2FullDomainKLoopRanges;
    typedef compute_k_loop_legs<
        Stage2FullDomainKLoopRanges,
        boost::mpl::void_,
        Stages2FullDomain,
        cKIncrement
    >::type Stage2FullDomainKLoopLegs;

    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<Stage2FullDomainKLoopRanges>::value == 2);

    // flat loop first
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<Stage2FullDomainKLoopRanges, 0>::type,
                KLoopRange<KRange<FlatCoordinates,0,0>, KBoundarySize<0,3> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<Stage2FullDomainKLoopLegs, 0>::type,
                Stages2FullDomain
            >::value
        )
    );    

    // terrain loop second
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<Stage2FullDomainKLoopRanges, 1>::type,
                KLoopRange<KRange<TerrainCoordinates,0,0>, KBoundarySize<1,2> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<Stage2FullDomainKLoopLegs, 1>::type,
                Stages2FullDomain
            >::value
        )
    );   
   
    // setup flat coordinates test with custom boundary
    typedef boost::mpl::vector1<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,-1,-2> >
    > Stages1FlatCoordinates; 
    
    typedef compute_k_loop_ranges<Stages1FlatCoordinates>::type Stages1FlatCoordinatesKLoopRanges;
    typedef compute_k_loop_legs<
        Stages1FlatCoordinatesKLoopRanges,
        boost::mpl::void_,
        Stages1FlatCoordinates,
        cKIncrement
    >::type Stages1FlatCoordinatesKLoopLegs;

    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<Stages1FlatCoordinatesKLoopRanges>::value == 1);

    // check the loop definition
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<Stages1FlatCoordinatesKLoopRanges, 0>::type,
                KLoopRange<KRange<FlatCoordinates,-1,-2>, KBoundarySize<4,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<Stages1FlatCoordinatesKLoopLegs, 0>::type,
                Stages1FlatCoordinates
            >::value
        )
    );  
        
    // setup terrain coordinates test with custom boundary
    typedef boost::mpl::vector1<
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,1,1> > 
    > Stages2TerrainCoordinates; 
        
    typedef compute_k_loop_ranges<Stages2TerrainCoordinates>::type Stages2TerrainCoordinatesKLoopRanges;
    typedef compute_k_loop_legs<
        Stages2TerrainCoordinatesKLoopRanges,
        boost::mpl::void_,
        Stages2TerrainCoordinates,
        cKIncrement
    >::type Stages2TerrainCoordinatesKLoopLegs;

    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<Stages2TerrainCoordinatesKLoopRanges>::value == 1);

    // check the forward loop definition
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<Stages2TerrainCoordinatesKLoopRanges, 0>::type,
                KLoopRange<KRange<TerrainCoordinates,1,1>, KBoundarySize<0,3> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<Stages2TerrainCoordinatesKLoopLegs, 0>::type,
                Stages2TerrainCoordinates
            >::value
        )
    );  
}

// make sure full domain loops are extended at the boundary
TEST(KLoopComputationUnittest, ComputeKLoopsFullDomain)
{
    // define 2 stencil stages one full domain and extend it with two boundary domains
    typedef boost::mpl::vector3<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,-1> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus2> >
    > FullStencilStages;

    typedef compute_k_loop_ranges<FullStencilStages>::type FullKLoopRanges;
    typedef compute_k_loop_legs<
        FullKLoopRanges,
        boost::mpl::void_,
        FullStencilStages,
        cKIncrement
    >::type FullKLoopLegs;

    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<FullKLoopRanges>::value == 1);

    // full domain loop
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullKLoopRanges, 0>::type,
                KLoopRange<KRange<FullDomain,-1,0>, KBoundarySize<4,1> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullKLoopLegs, 0>::type,
                FullStencilStages
            >::value
        )
    );  
}

// make sure flat domain loops are extended 
TEST(KLoopComputationUnittest, ComputeKLoopsFlatDomain)
{
    // define 2 stencil stages one full domain and extend it with two boundary domains
    typedef boost::mpl::vector3<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,2,-1> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus1> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >
    > FlatStencilStages;

    typedef compute_k_loop_ranges<FlatStencilStages>::type FlatKLoopRanges;
    typedef compute_k_loop_legs<
        FlatKLoopRanges,
        boost::mpl::void_,
        FlatStencilStages,
        cKIncrement
    >::type FlatKLoopLegs;
    
    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<FlatKLoopRanges>::value == 2);

    // flat domain loop
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FlatKLoopRanges, 0>::type,
                KLoopRange<KRange<FlatCoordinates,1,-1>, KBoundarySize<2,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FlatKLoopLegs, 0>::type,
                boost::mpl::vector2<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,2,-1> >,
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus1> >
                >
            >::value
        )
    );  

    // k maximum center loop
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FlatKLoopRanges, 1>::type,
                KLoopRange<KRange<KMaximumCenter>, KBoundarySize<0,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FlatKLoopLegs, 1>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >
                >
            >::value
        )
    );  
}

// make sure terrain domain loops are extended 
TEST(KLoopComputationUnittest, ComputeKLoopsTerrainDomain)
{
    // define 2 stencil stages one full domain and extend it with two boundary domains
    typedef boost::mpl::vector3<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,2,0> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainPlus1> >
    > TerrainStencilStages;

    typedef compute_k_loop_ranges<TerrainStencilStages>::type TerrainKLoopRanges;
    typedef compute_k_loop_legs<
        TerrainKLoopRanges,
        boost::mpl::void_,
        TerrainStencilStages,
        cKIncrement
    >::type TerrainKLoopLegs;
    
    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<TerrainKLoopRanges>::value == 2);

    // k flat maximum domain loop
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<TerrainKLoopRanges, 0>::type,
                KLoopRange<KRange<KMaximumFlatCenter>, KBoundarySize<0,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<TerrainKLoopLegs, 0>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> > 
                >
            >::value
        )
    );

    // k terrain loop
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<TerrainKLoopRanges, 1>::type,
                KLoopRange<KRange<TerrainCoordinates,1,0>, KBoundarySize<1,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<TerrainKLoopLegs, 1>::type,
                boost::mpl::vector2<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,2,0> >,
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainPlus1> >
                >
            >::value
        )
    );
}

// make sure boundary only loops work
TEST(KLoopComputationUnittest, ComputeKLoopsBoundaries)
{
    // define 2 stencil stages with multiple boundary levels
    typedef boost::mpl::vector3<
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumMinus1> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainPlus1> >
    > BoundaryStencilStages;

    typedef compute_k_loop_ranges<BoundaryStencilStages>::type BoundaryKLoopRanges;
    typedef compute_k_loop_legs<
        BoundaryKLoopRanges,
        boost::mpl::void_,
        BoundaryStencilStages,
        cKIncrement
    >::type BoundaryKLoopLegs;
        
    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<BoundaryKLoopRanges>::value == 3);

    // test all boundary levels
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<BoundaryKLoopRanges, 0>::type,
                KLoopRange<KRange<KMinimumMinus1>, KBoundarySize<0,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<BoundaryKLoopLegs, 0>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumMinus1> >
                >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<BoundaryKLoopRanges, 1>::type,
                KLoopRange<KRange<KMaximumFlatCenter>, KBoundarySize<0,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<BoundaryKLoopLegs, 1>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >
                >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<BoundaryKLoopRanges, 2>::type,
                KLoopRange<KRange<KMinimumTerrainPlus1>, KBoundarySize<0,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<BoundaryKLoopLegs, 2>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainPlus1> >
                >
            >::value
        )
    );
}

// make sure full domain loops are split up
TEST(KLoopComputationUnittest, ComputeKLoopsFullDomainSplit)
{
    // define 2 stencil stages one full domain and one flat coordinates
    typedef boost::mpl::vector2<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,0> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,1,-2> >
    > FullAndFlatStencilStages;

    typedef compute_k_loop_ranges<FullAndFlatStencilStages>::type FullAndFlatKLoopRanges;
    typedef compute_k_loop_legs<
        FullAndFlatKLoopRanges,
        boost::mpl::void_,
        FullAndFlatStencilStages,
        cKIncrement
    >::type FullAndFlatKLoopLegs;

    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<FullAndFlatKLoopRanges>::value == 2);

    // flat loop first
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullAndFlatKLoopRanges, 0>::type,
                KLoopRange<KRange<FlatCoordinates,-1,0>, KBoundarySize<4,3> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullAndFlatKLoopLegs, 0>::type,
                boost::mpl::vector2<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,0> >,
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,1,-2> >
                >
            >::value
        )
    );
    
    // terrain loop second
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullAndFlatKLoopRanges, 1>::type,
                KLoopRange<KRange<TerrainCoordinates,0,0>, KBoundarySize<0,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullAndFlatKLoopLegs, 1>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,0> >
                >
            >::value
        )
    );
    
    // define 2 stencil stages one full domain and one terrain coordinates
    typedef boost::mpl::vector2<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,0> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,2,1> >
    > FullAndTerrainStencilStages;

    typedef compute_k_loop_ranges<FullAndTerrainStencilStages>::type FullAndTerrainKLoopRanges;
    typedef compute_k_loop_legs<
        FullAndTerrainKLoopRanges,
        boost::mpl::void_,
        FullAndTerrainStencilStages,
        cKIncrement
    >::type FullAndTerrainKLoopLegs;

    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<FullAndTerrainKLoopRanges>::value == 2);

    // flat loop first
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullAndTerrainKLoopRanges, 0>::type,
                KLoopRange<KRange<FlatCoordinates,-1,0>, KBoundarySize<4,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullAndTerrainKLoopLegs, 0>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,0> >
                >
            >::value
        )
    );
    
    // terrain loop second
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullAndTerrainKLoopRanges, 1>::type,
                KLoopRange<KRange<TerrainCoordinates,0,1>, KBoundarySize<2,3> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullAndTerrainKLoopLegs, 1>::type,
                boost::mpl::vector2<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,0> >,
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,2,1> >
                >
            >::value
        )
    );

    // define 2 stencil stages one full domain the other one at k maximum flat minus 1 boundary
    typedef boost::mpl::vector2<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,1> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatMinus1> >
    > FullAndFlatBoundaryStencilStages;

    typedef compute_k_loop_ranges<FullAndFlatBoundaryStencilStages>::type FullAndFlatBoundaryKLoopRanges;
    typedef compute_k_loop_legs<
        FullAndFlatBoundaryKLoopRanges,
        boost::mpl::void_,
        FullAndFlatBoundaryStencilStages,
        cKIncrement
    >::type FullAndFlatBoundaryKLoopLegs;
    
    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<FullAndFlatBoundaryKLoopRanges>::value == 2);

    // flat loop first
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullAndFlatBoundaryKLoopRanges, 0>::type,
                KLoopRange<KRange<FlatCoordinates,-1,0>, KBoundarySize<4,2> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullAndFlatBoundaryKLoopLegs, 0>::type,
                boost::mpl::vector2<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,1> >,
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatMinus1> >
                >
            >::value
        )
    );

    // terrain loop second
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullAndFlatBoundaryKLoopRanges, 1>::type,
                KLoopRange<KRange<TerrainCoordinates,0,1>, KBoundarySize<0,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullAndFlatBoundaryKLoopLegs, 1>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,1> >
                >
            >::value
        )
    );
    
    // define 2 stencil stages one full domain the other one at k minimum terrain center boundary
    typedef boost::mpl::vector2<
        StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,-1> >,
        StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainCenter> >
    > FullAndTerrainBoundaryStencilStages;

    typedef compute_k_loop_ranges<FullAndTerrainBoundaryStencilStages>::type FullAndTerrainBoundaryKLoopRanges;
    typedef compute_k_loop_legs<
        FullAndTerrainBoundaryKLoopRanges,
        boost::mpl::void_,
        FullAndTerrainBoundaryStencilStages,
        cKIncrement
    >::type FullAndTerrainBoundaryKLoopLegs;
     
    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<FullAndTerrainBoundaryKLoopRanges>::value == 2);

    // flat loop first
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullAndTerrainBoundaryKLoopRanges, 0>::type,
                KLoopRange<KRange<FlatCoordinates,-1,0>, KBoundarySize<4,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullAndTerrainBoundaryKLoopLegs, 0>::type,
                boost::mpl::vector1<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,-1> >
                >
            >::value
        )
    );

    // terrain loop second
    ASSERT_TRUE(
        (
            boost::is_same<
                boost::mpl::at_c<FullAndTerrainBoundaryKLoopRanges, 1>::type,
                KLoopRange<KRange<TerrainCoordinates,0,-1>, KBoundarySize<1,0> >
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_k_loop_leg_stencil_stages<
                boost::mpl::at_c<FullAndTerrainBoundaryKLoopLegs, 1>::type,
                boost::mpl::vector2<
                    StencilStage<KLoopDescriptorStage1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,-1,-1> >,
                    StencilStage<KLoopDescriptorStage2, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainCenter> >
                >
            >::value
        )
    );
}

// test that the cache k range computation works properly
TEST(KLoopComputationUnittest, ComputeCacheInnerKRange)
{
    // test slide is not considered for the update range calculation
    // (return the max boundary levels range!)
    ASSERT_TRUE(
        (
            boost::is_same<
                cache_full_update_k_range<
                    KCache<0, cLocal, KWindow<-1,1>, KRange<FullDomain,1,-1> >,
                    boost::mpl::integral_c<KLoopDirection, cKIncrement>
                >::type,
                KRange<FullDomain,-MAX_BOUNDARY_LEVELS,MAX_BOUNDARY_LEVELS>
            >::value 
        )
    );

    // test fill and flush cache
    ASSERT_TRUE(
        (
            boost::is_same<
                cache_full_update_k_range<
                    KCache<0, cFillAndFlush, KWindow<-1,2>, KRange<FlatCoordinates,-1,0> >,
                    boost::mpl::integral_c<KLoopDirection, cKIncrement>
                >::type,
                KRange<FlatCoordinates,0,-2>
            >::value &&
            boost::is_same<
                cache_full_update_k_range<
                    KCache<0, cFillAndFlush, KWindow<-1,2>, KRange<FlatCoordinates,-1,0> >,
                    boost::mpl::integral_c<KLoopDirection, cKDecrement>
                >::type,
                KRange<FlatCoordinates,0,-2>
            >::value
        )
    );

    // test fill only cache
    ASSERT_TRUE(
        (
            boost::is_same<
                cache_full_update_k_range<
                    KCache<0, cFill, KWindow<-2,-1>, KRange<FlatCoordinates,0,-1> >,
                    boost::mpl::integral_c<KLoopDirection, cKIncrement>
                >::type,
                KRange<FlatCoordinates,1,0>
            >::value &&
            boost::is_same<
                cache_full_update_k_range<
                    KCache<0, cFill, KWindow<-2,-1>, KRange<FlatCoordinates,0,-1> >,
                    boost::mpl::integral_c<KLoopDirection, cKDecrement>
                >::type,
                KRange<FlatCoordinates,2,1>
            >::value
        )
    );

    // test flush only cache
    ASSERT_TRUE(
        (
            boost::is_same<
                cache_full_update_k_range<
                    KCache<0, cFlush, KWindow<-1,2>, KRange<FlatCoordinates,0,-1> >,
                    boost::mpl::integral_c<KLoopDirection, cKIncrement>
                >::type,
                KRange<FlatCoordinates,1,0>
            >::value &&
            boost::is_same<
                cache_full_update_k_range<
                    KCache<0, cFlush, KWindow<-1,2>, KRange<FlatCoordinates,0,-1> >,
                    boost::mpl::integral_c<KLoopDirection, cKDecrement>
                >::type,
                KRange<FlatCoordinates,-2,-3>
            >::value
        )
    );
}

// meta function verifying two k loop legs are identical
template<
    typename TExpectedKLoopLeg,
    typename TActualKLoopLeg>
struct check_k_loop_leg : boost::mpl::false_ {};

template<
    typename TKLoopRange,
    typename TExpectedCaches,
    typename TExpectedStencilStages,
    typename TActualCaches,
    typename TActualStencilStages,
    typename TAdvanceToKPosition,
    KLoopDirection VKLoopDirection>
struct check_k_loop_leg<
    KLoopLeg<TKLoopRange, TExpectedCaches, TExpectedStencilStages, TAdvanceToKPosition, VKLoopDirection>,
    KLoopLeg<TKLoopRange, TActualCaches, TActualStencilStages, TAdvanceToKPosition, VKLoopDirection> > : 
    boost::mpl::and_<
        typename boost::mpl::equal<TExpectedCaches, TActualCaches>::type, 
        typename boost::mpl::equal<TExpectedStencilStages, TActualStencilStages>::type 
    > 
{};
    
// test that caches extend the k loop boundaries correctly
TEST(KLoopComputationUnittest, ExtendKLoopBoundariesForCaches)
{
    // test caches with a single large full domain
    typedef boost::mpl::vector1<
        StencilStage<KLoopDescriptorStage3, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >
    > Stages1FullDomain; 
    typedef boost::mpl::vector1<
        KCache<0, cFillAndFlush, KWindow<-1,2>, KRange<FullDomain,0,0> >
    > Caches1;

    typedef compute_k_loop_ranges<
        Stages1FullDomain
    >::type Stage1FullDomainKLoopRanges;
    typedef compute_k_loop_legs<
        Stage1FullDomainKLoopRanges,
        Caches1,
        Stages1FullDomain,
        cKIncrement
    >::type Stage1FullDomainKLoopLegs;
       
    // assert there is one loop
    ASSERT_TRUE(boost::mpl::size<Stage1FullDomainKLoopRanges>::value == 1);

    // check the loop definition
    ASSERT_TRUE(
        (
            check_k_loop_leg<
                KLoopLeg<
                    KLoopRange<KRange<FullDomain,0,0>, KBoundarySize<1,2> >,
                    Caches1,
                    Stages1FullDomain,
                    boost::mpl::void_,
                    cKIncrement
                >,
                boost::mpl::at_c<Stage1FullDomainKLoopLegs, 0>::type
            >::value
        )
    );

    // test that full domain caches work well with flat coordinates domains
    typedef boost::mpl::vector1<
        StencilStage<KLoopDescriptorStage3, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,-1,0> >
    > Stages2FlatCoordinates; 
    typedef boost::mpl::vector2<
        KCache<0, cFill, KWindow<-1,2>, KRange<FullDomain,0,0> >,
        KCache<1, cFlush, KWindow<-2,2>, KRange<FullDomain,0,0> >
    > Caches2;

    typedef compute_k_loop_ranges<
        Stages2FlatCoordinates
    >::type Stage2FlatCoordinatesKLoopRanges;
    typedef compute_k_loop_legs<
        Stage2FlatCoordinatesKLoopRanges,
        Caches2,
        Stages2FlatCoordinates,
        cKDecrement
    >::type Stage2FlatCoordinatesKLoopLegs;
       
    // assert there are 2 loops
    ASSERT_TRUE(boost::mpl::size<Stage2FlatCoordinatesKLoopLegs>::value == 1);

    // check the loop definition
    ASSERT_TRUE(
        (
            check_k_loop_leg<
                KLoopLeg<
                    KLoopRange<KRange<FlatCoordinates,-1,0>, KBoundarySize<2,0> >,
                    Caches2,
                    Stages2FlatCoordinates,
                    boost::mpl::void_,
                    cKDecrement
                >,
                boost::mpl::at_c<Stage2FlatCoordinatesKLoopLegs, 0>::type
            >::value
        )
    );

    // test that full domain caches work well with terrain coordinates
    typedef boost::mpl::vector1<
        StencilStage<KLoopDescriptorStage3, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,-1> >
    > Stages3TerrainCoordinates; 
    typedef boost::mpl::vector2<
        KCache<0, cFill, KWindow<-0,2>, KRange<FullDomain,0,0> >,
        KCache<1, cFlush, KWindow<-0,1>, KRange<FullDomain,0,0> >
    > Caches3;

    typedef compute_k_loop_ranges<
        Stages3TerrainCoordinates
    >::type Stage3TerrainCoordinatesKLoopRanges;
    typedef compute_k_loop_legs<
        Stage3TerrainCoordinatesKLoopRanges,
        Caches3,
        Stages3TerrainCoordinates,
        cKIncrement
    >::type Stage3TerrainCoordinatesKLoopLegs;
       
    // assert there are 2 loops
    ASSERT_TRUE(boost::mpl::size<Stage3TerrainCoordinatesKLoopLegs>::value == 1);

    // check the loop definition
    ASSERT_TRUE(
        (
            check_k_loop_leg<
                KLoopLeg<
                    KLoopRange<KRange<TerrainCoordinates,0,-1>, KBoundarySize<0,1> >,
                    Caches3,
                    Stages3TerrainCoordinates,
                    boost::mpl::void_,
                    cKIncrement
                >,
                boost::mpl::at_c<Stage3TerrainCoordinatesKLoopLegs, 0>::type
            >::value
        )
    );

    // test that caches overlapping a 2d domain are fine
    typedef boost::mpl::vector1<
        StencilStage<KLoopDescriptorStage3, IJRange<cComplete,0,0,0,0>, KRange<KMaximum<TerrainCoordinates,-1>,0,0> >
    > Stages4KMaximum; 
    typedef boost::mpl::vector2<
        KCache<0, cFill, KWindow<-0,2>, KRange<FullDomain,0,0> >,
        KCache<1, cFlush, KWindow<-0,1>, KRange<FullDomain,0,0> >
    > Caches4;

    typedef compute_k_loop_ranges<
        Stages4KMaximum
    >::type Stage4KMaximumKLoopRanges;
    typedef compute_k_loop_legs<
        Stage4KMaximumKLoopRanges,
        Caches4,
        Stages4KMaximum,
        cKIncrement
    >::type Stage4KMaximumKLoopLegs;
       
    // assert there are 2 loops
    ASSERT_TRUE(boost::mpl::size<Stage4KMaximumKLoopLegs>::value == 1);

    // check the loop definition
    ASSERT_TRUE(
        (
            check_k_loop_leg<
                KLoopLeg<
                    KLoopRange<KRange<KMaximum<TerrainCoordinates,-1>,0,0>, KBoundarySize<0,0> >,
                    Caches4,
                    Stages4KMaximum,
                    boost::mpl::void_,
                    cKIncrement
                >,
                boost::mpl::at_c<Stage4KMaximumKLoopLegs, 0>::type
            >::value
        )
    );
}

