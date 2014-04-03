#include <boost/mpl/integral_c.hpp>
#include "gtest/gtest.h"
#include "StencilFramework.h"

// define the parameter enum
enum { a, b, c, fullbuffer, terrainbuffer, flatbuffer,  blockbuffer, scalarbuffer, result };

template<typename TEnv>
struct WriteStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, fullbuffer)
    STAGE_PARAMETER(FullDomain, terrainbuffer)
    STAGE_PARAMETER(FullDomain, flatbuffer)
    STAGE_PARAMETER(FullDomain, blockbuffer)

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        ctx[terrainbuffer::At(kplus1)] = 1.0;
        ctx[terrainbuffer::Center()] = 0.0;
        ctx[fullbuffer::At(kplus1)] = 1.0;
        ctx[fullbuffer::Center()] = 0.0;
    }

    __ACC__
    static void Do(Context ctx, TerrainCoordinates)
    {
        ctx[fullbuffer::Center()] = 22.0;
        ctx[terrainbuffer::Center()] = 44.0;
    }

    __ACC__
    static void Do(Context ctx, FlatCoordinates)
    {
        ctx[fullbuffer::Center()] = 22.0;
        ctx[flatbuffer::Center()] = 33.0;
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        ctx[blockbuffer::Center()] = 7.0;
        ctx[fullbuffer::Center()] = 0.0;
        ctx[fullbuffer::At(kminus1)] = -1.0;
        ctx[flatbuffer::Center()] = 0.0;
        ctx[flatbuffer::At(kminus1)] = -1.0;
    }
};

template<typename TEnv>
struct ReadStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, result)
    STAGE_PARAMETER(FullDomain, fullbuffer)
    STAGE_PARAMETER(FullDomain, terrainbuffer)
    STAGE_PARAMETER(FullDomain, flatbuffer)
    STAGE_PARAMETER(FullDomain, blockbuffer)
    STAGE_PARAMETER(FullDomain, scalarbuffer)

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        ctx[result::At(kplus2)] = ctx[scalarbuffer::Center()];
        ctx[result::At(kplus1)] = ctx[fullbuffer::At(kplus1)] + ctx[terrainbuffer::At(kplus1)];
        ctx[result::Center()] = ctx[fullbuffer::Center()] + ctx[terrainbuffer::Center()] + ctx[blockbuffer::Center()];
    }

    __ACC__
    static void Do(Context ctx, TerrainCoordinates)
    {
        ctx[scalarbuffer::Center()] += 1.0;
        ctx[result::Center()] = ctx[fullbuffer::Center()] + ctx[terrainbuffer::Center()];
    }

    __ACC__
    static void Do(Context ctx, FlatCoordinates)
    {
        ctx[scalarbuffer::Center()] += 1.0;
        ctx[result::Center()] = ctx[fullbuffer::Center()] + ctx[flatbuffer::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        ctx[scalarbuffer::Center()] = 0;
        ctx[result::Center()] = ctx[fullbuffer::Center()] + ctx[flatbuffer::Center()];
        ctx[result::At(kminus1)] = ctx[fullbuffer::At(kminus1)] + ctx[flatbuffer::At(kminus1)];
    }
};

// check that stencil buffers work
TEST(StencilBufferUnittest, Buffers)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);
    
    // define a result field
    IJKRealField result;
    result.Init("result", size, boundary);

    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::result, cInOut>(result)
        ),
        define_temporaries(
            StencilBuffer<fullbuffer, Real, KRange<FullDomain,-1,1> >(),
            StencilBuffer<terrainbuffer, Real, KRange<TerrainCoordinates,0,1> >(),
            StencilBuffer<flatbuffer, Real, KRange<FlatCoordinates,-1,0> >(),
            StencilBuffer<blockbuffer, Real,  KRange<KMinimumCenter> >(),
            StageVariable<scalarbuffer, Real>()
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<WriteStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<ReadStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // check the sum was done right
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            ASSERT_DOUBLE_EQ(38.0, result(i,j,41));
            ASSERT_DOUBLE_EQ(1.0 + 1.0, result(i,j,40));
            ASSERT_DOUBLE_EQ(7.0, result(i,j,39));
            for(int k=cFlatLimit; k<38; k++)
            {
                ASSERT_DOUBLE_EQ(22.0 + 44.0, result(i,j,k));
            }
            for(int k=1; k<cFlatLimit; k++)
            {
                ASSERT_DOUBLE_EQ(22.0 + 33.0, result(i,j,k));
            }
            ASSERT_DOUBLE_EQ(0.0, result(i,j,0));
            ASSERT_DOUBLE_EQ(-1.0 + -1.0, result(i,j,-1));
        }
    }
}

template<typename TEnv>
struct InitStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, fullbuffer)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[fullbuffer::Center()] = 0.0;
    }
};

template<typename TEnv>
struct AddStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, fullbuffer)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[fullbuffer::Center()] += 2.0;
    }
};

template<typename TEnv>
struct AverageStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, result)
    STAGE_PARAMETER(FullDomain, fullbuffer)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[result::Center()] = (T)0.25 *
            (
                ctx[fullbuffer::At(Offset<-3,-3,0>())] +
                ctx[fullbuffer::At(Offset<3,-3,0>())] +
                ctx[fullbuffer::At(Offset<-3,3,0>())] +
                ctx[fullbuffer::At(Offset<3,3,0>())]
            );
    }
};

// check that stencil buffers boundaries are thread/block private
TEST(StencilBufferUnittest, BlockPrivateAccess)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);

    // define a result field
    IJKRealField result;
    result.Init("result", size, boundary);

    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::result, cInOut>(result)
        ),
        define_temporaries(
            StencilBuffer<fullbuffer, Real, KRange<FullDomain,0,0> >()
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<InitStage, IJRange<cComplete,-3,3,-3,3>, KRange<FullDomain,0,0> >(),
                    StencilStage<AddStage, IJRange<cComplete,-3,3,-3,3>, KRange<FullDomain,0,0> >(),
                    StencilStage<AverageStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // check that none of the buffer boundary elements was updated by two threads?
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(2.0, result(i,j,k));
            }
        }
    }
}

// check the ij boundary size calculation 
template<typename TEnv>
struct StageA
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, a)
    STAGE_PARAMETER(FullDomain, b)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
};

template<typename TEnv>
struct StageB
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, b)
    STAGE_PARAMETER(FullDomain, c)

    __ACC__
    static void Do(Context ctx, FullDomain) {}
};

template<typename TStencilSweepDescriptors>
void check_maximum_range(TStencilSweepDescriptors)
{
    // compute stencil sweep group descriptors list
    typedef typename boost::mpl::transform<
        TStencilSweepDescriptors,
        create_stencil_sweep_group_descriptor_from_sweep<boost::mpl::_>
    >::type StencilSweepGroupDescriptors;

    typedef typename maximum_parameter_ij_range<
        StencilSweepGroupDescriptors,
        boost::mpl::integral_c<int,a>
    >::type Param1Range;

    ASSERT_EQ(-1, Param1Range::IMinusOffset::value);
    ASSERT_EQ(2, Param1Range::IPlusOffset::value);
    ASSERT_EQ(-3, Param1Range::JMinusOffset::value);
    ASSERT_EQ(4, Param1Range::JPlusOffset::value);

    typedef typename maximum_parameter_ij_range<
        StencilSweepGroupDescriptors,
        boost::mpl::integral_c<int,b>
    >::type Param2Range;

    ASSERT_EQ(-4, Param2Range::IMinusOffset::value);
    ASSERT_EQ(3, Param2Range::IPlusOffset::value);
    ASSERT_EQ(-3, Param2Range::JMinusOffset::value);
    ASSERT_EQ(4, Param2Range::JPlusOffset::value);

    typedef typename maximum_parameter_ij_range<
        StencilSweepGroupDescriptors,
        boost::mpl::integral_c<int,c>
    >::type Param3Range;

    ASSERT_EQ(-4, Param3Range::IMinusOffset::value);
    ASSERT_EQ(3, Param3Range::IPlusOffset::value);
    ASSERT_EQ(-2, Param3Range::JMinusOffset::value);
    ASSERT_EQ(1, Param3Range::JPlusOffset::value);
}

TEST(StencilBufferUnittest, MaxParameterBoundary)
{
    // check if the boundaries for the individual parameters are correct
    check_maximum_range(
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StageA, IJRange<cComplete,-1,2,-3,4>, KRange<FullDomain,0,0> >()
                )
            ),
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StageB, IJRange<cComplete,-4,3,-2,1>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
}
  
TEST(StencilBufferUnittest, MaximumIJRange)
{
    typedef IJRange<cComplete,-1,0,-2,3> Boundary1;
    typedef IJRange<cComplete,-3,1,-1,0> Boundary2;

    // test if the boundary is maximized
    typedef maximum_ij_range<
        Boundary1,
        Boundary2
    >::type MaximumIJRange;

    ASSERT_EQ(-3, MaximumIJRange::IMinusOffset::value);
    ASSERT_EQ(1, MaximumIJRange::IPlusOffset::value);
    ASSERT_EQ(-2, MaximumIJRange::JMinusOffset::value);
    ASSERT_EQ(3, MaximumIJRange::JPlusOffset::value);
}
