#include "gtest/gtest.h"
#include <boost/type_traits/is_same.hpp>
#include "StencilFramework.h"

enum { par0, par1, par2, par3, par4, par5 };

template<typename TEnv>
struct Stage1
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, par1)
    STAGE_PARAMETER(FullDomain, par2)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
};

template<typename TEnv>
struct Stage2
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, par0)
    STAGE_PARAMETER(FullDomain, par1)
    STAGE_PARAMETER(FullDomain, par2)
    STAGE_PARAMETER(FullDomain, par3)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
};

template<typename TEnv>
struct Stage3
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, par0)
    STAGE_PARAMETER(FullDomain, par1)
    STAGE_PARAMETER(FullDomain, par4)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
};

// define stencil stage list and single stage
typedef boost::mpl::vector3<
    StencilStage<Stage1, IJRange<cComplete, 0, 0, 0, 0>, KRange<FullDomain,0,0> >,
    StencilStage<Stage2, IJRange<cComplete, -1, 2, -3, 4>, KRange<FullDomain,0,0> >,
    StencilStage<Stage3, IJRange<cComplete, -5, 6, -7, 8>, KRange<FullDomain,0,0> >
> AllStages;
typedef boost::mpl::vector1<
    StencilStage<Stage2, IJRange<cComplete, -1, 2, -3, 4>, KRange<FullDomain,0,0> >
> Stage2Only;

// compute the maximal ij range
TEST(StencilStageUnittest, MaximumIJRange)
{
    typedef stencil_stages_maximum_ij_range<AllStages>::type AllStagesMaxRange;
    typedef stencil_stages_maximum_ij_range<Stage2Only>::type Stage2MaxRange;

    ASSERT_EQ(-5, AllStagesMaxRange::IMinusOffset::value);
    ASSERT_EQ(6, AllStagesMaxRange::IPlusOffset::value);
    ASSERT_EQ(-7, AllStagesMaxRange::JMinusOffset::value);
    ASSERT_EQ(8, AllStagesMaxRange::JPlusOffset::value);

    ASSERT_EQ(-1, Stage2MaxRange::IMinusOffset::value);
    ASSERT_EQ(2, Stage2MaxRange::IPlusOffset::value);
    ASSERT_EQ(-3, Stage2MaxRange::JMinusOffset::value);
    ASSERT_EQ(4, Stage2MaxRange::JPlusOffset::value);
}

// compute parameter ij range
TEST(StencilStageUnittest, ParameterIJRange)
{
    typedef stencil_stages_parameter_ij_range<AllStages, boost::mpl::integral_c<int, par0> >::type AllStagesPar0;
    typedef stencil_stages_parameter_ij_range<AllStages, boost::mpl::integral_c<int, par1> >::type AllStagesPar1;
    typedef stencil_stages_parameter_ij_range<AllStages, boost::mpl::integral_c<int, par4> >::type AllStagesPar4;
    
    ASSERT_TRUE((boost::is_same<AllStagesPar0, AllStagesPar1>::value));
    ASSERT_TRUE((boost::is_same<AllStagesPar0, AllStagesPar4>::value));
    ASSERT_EQ(-5, AllStagesPar0::IMinusOffset::value);
    ASSERT_EQ(6, AllStagesPar0::IPlusOffset::value);
    ASSERT_EQ(-7, AllStagesPar0::JMinusOffset::value);
    ASSERT_EQ(8, AllStagesPar0::JPlusOffset::value);

    typedef stencil_stages_parameter_ij_range<AllStages, boost::mpl::integral_c<int, par2> >::type AllStagesPar2;
    typedef stencil_stages_parameter_ij_range<AllStages, boost::mpl::integral_c<int, par3> >::type AllStagesPar3;

    ASSERT_TRUE((boost::is_same<AllStagesPar2, AllStagesPar3>::value));
    ASSERT_EQ(-1, AllStagesPar2::IMinusOffset::value);
    ASSERT_EQ(2, AllStagesPar2::IPlusOffset::value);
    ASSERT_EQ(-3, AllStagesPar2::JMinusOffset::value);
    ASSERT_EQ(4, AllStagesPar2::JPlusOffset::value);
}

// compute stages first and last usage
TEST(StencilStageUnittest, StageFirstLastUsage)
{
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_first_parameter_usage<AllStages, boost::mpl::integral_c<int, par0> >::type,
            StencilStage<Stage2, IJRange<cComplete, -1, 2, -3, 4>, KRange<FullDomain,0,0> >
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_last_parameter_usage<AllStages, boost::mpl::integral_c<int, par0> >::type,
            StencilStage<Stage3, IJRange<cComplete, -5, 6, -7, 8>, KRange<FullDomain,0,0> >
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_first_parameter_usage<AllStages, boost::mpl::integral_c<int, par1> >::type,
            StencilStage<Stage1, IJRange<cComplete, 0, 0, 0, 0>, KRange<FullDomain,0,0> >
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_last_parameter_usage<AllStages, boost::mpl::integral_c<int, par1> >::type,
            StencilStage<Stage3, IJRange<cComplete, -5, 6, -7, 8>, KRange<FullDomain,0,0> >
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_first_parameter_usage<AllStages, boost::mpl::integral_c<int, par2> >::type,
            StencilStage<Stage1, IJRange<cComplete, 0, 0, 0, 0>, KRange<FullDomain,0,0> >
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_last_parameter_usage<AllStages, boost::mpl::integral_c<int, par2> >::type,
            StencilStage<Stage2, IJRange<cComplete, -1, 2, -3, 4>, KRange<FullDomain,0,0> >
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_first_parameter_usage<AllStages, boost::mpl::integral_c<int, par3> >::type,
            StencilStage<Stage2, IJRange<cComplete, -1, 2, -3, 4>, KRange<FullDomain,0,0> >
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_last_parameter_usage<AllStages, boost::mpl::integral_c<int, par3> >::type,
            StencilStage<Stage2, IJRange<cComplete, -1, 2, -3, 4>, KRange<FullDomain,0,0> >
        >::value
    ));

    // check the first / last stage is returned if the parameter is not used
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_first_parameter_usage<AllStages, boost::mpl::integral_c<int, par5> >::type,
            StencilStage<Stage1, IJRange<cComplete, 0, 0, 0, 0>, KRange<FullDomain,0,0> >
        >::value
    ));
    ASSERT_TRUE((
        boost::is_same<
            stencil_stages_last_parameter_usage<AllStages, boost::mpl::integral_c<int, par5> >::type,
            StencilStage<Stage3, IJRange<cComplete, -5, 6, -7, 8>, KRange<FullDomain,0,0> >
        >::value
    ));
}

// check if stages use parameter
TEST(StencilStageUnittest, UseParameter)
{
    ASSERT_TRUE((stencil_stages_have_parameter<AllStages, FullDomain, boost::mpl::integral_c<int, par0> >::value));
    ASSERT_TRUE((stencil_stages_have_parameter<AllStages, FullDomain, boost::mpl::integral_c<int, par1> >::value));
    ASSERT_TRUE((stencil_stages_have_parameter<AllStages, FullDomain, boost::mpl::integral_c<int, par2> >::value));
    ASSERT_TRUE((stencil_stages_have_parameter<AllStages, FullDomain, boost::mpl::integral_c<int, par3> >::value));
    ASSERT_TRUE((stencil_stages_have_parameter<AllStages, FullDomain, boost::mpl::integral_c<int, par4> >::value));
    ASSERT_FALSE((stencil_stages_have_parameter<AllStages, FullDomain, boost::mpl::integral_c<int, par5> >::value));

    ASSERT_TRUE((stencil_stages_have_parameter<Stage2Only, FullDomain, boost::mpl::integral_c<int, par0> >::value));
    ASSERT_TRUE((stencil_stages_have_parameter<Stage2Only, FullDomain, boost::mpl::integral_c<int, par1> >::value));
    ASSERT_TRUE((stencil_stages_have_parameter<Stage2Only, FullDomain, boost::mpl::integral_c<int, par2> >::value));
    ASSERT_TRUE((stencil_stages_have_parameter<Stage2Only, FullDomain, boost::mpl::integral_c<int, par3> >::value));
    ASSERT_FALSE((stencil_stages_have_parameter<Stage2Only, FullDomain, boost::mpl::integral_c<int, par4> >::value));
    ASSERT_FALSE((stencil_stages_have_parameter<Stage2Only, FullDomain, boost::mpl::integral_c<int, par5> >::value));
}

template<typename TEnv>
struct DomainsStage
{
    STENCIL_STAGE(TEnv)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
    __ACC__
    static void Do(Context ctx, TerrainCoordinates) {}
    __ACC__
    static void Do(Context ctx, KMinimumCenter) {}
};

// check do domain computation
TEST(StencilStageUnittest, LookupDoDomain)
{
    // define the test stage
    typedef StencilStage<DomainsStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> > StencilStage;

    ASSERT_TRUE((boost::is_same<stencil_stage_lookup_do_domain<StencilStage, KMaximumCenter>::type, TerrainCoordinates>::value));
    ASSERT_TRUE((boost::is_same<stencil_stage_lookup_do_domain<StencilStage, KMaximumMinus1>::type, TerrainCoordinates>::value));
    ASSERT_TRUE((boost::is_same<stencil_stage_lookup_do_domain<StencilStage, FlatCoordinates>::type, FullDomain>::value));
    ASSERT_TRUE((boost::is_same<stencil_stage_lookup_do_domain<StencilStage, KMinimumMinus1>::type, FullDomain>::value));
    ASSERT_TRUE((boost::is_same<stencil_stage_lookup_do_domain<StencilStage, KMinimumCenter>::type, KMinimumCenter>::value));
}

