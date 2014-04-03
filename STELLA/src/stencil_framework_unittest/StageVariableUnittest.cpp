#include <boost/mpl/integral_c.hpp>
#include "gtest/gtest.h"
#include "StencilFramework.h"

// define the parameter enum
enum { buffer,input, output, stagevar };

// stage that copies input to buffer, using stagevar as intermediate storage
template<typename TEnv>
struct StageOne
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, input)
    STAGE_PARAMETER(FullDomain, buffer)
    STAGE_PARAMETER(FullDomain, stagevar)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[stagevar::Center()]  = ctx[input::Center()];
        ctx[buffer::Center()] = ctx[stagevar::Center()];
    }

};

// stage that writes buffer to output, using stagevar as intermediate storage
template<typename TEnv>
struct StageTwo
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, output)
    STAGE_PARAMETER(FullDomain, buffer)
    STAGE_PARAMETER(FullDomain, stagevar)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[stagevar::Center()]  = ctx[buffer::Center()];
        ctx[output::Center()]    = ctx[stagevar::Center()];
    }

};

// check that stage buffers work
TEST(StageVariableUnittest, Buffers)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(0, 0);

    // define input and output fields
    IJKRealField input, output;
    input.Init("input",   size, boundary);
    output.Init("output", size, boundary);

    // this should fail at compile time
    // because the StageVariable stagevar is used in more than one stage of the sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "StageVariableTestStencil",
        size,
        StencilConfiguration<double, BlockSize<4,4> >(),
        pack_parameters(
            Param< ::input,  cInOut>(input),
            Param< ::output, cInOut>(output)
        ),
        define_temporaries(
            StencilBuffer<buffer, Real, KRange<FullDomain,0,0> >(),
            StageVariable<stagevar, Real>()
        ),
        concatenate_sweeps(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StageOne, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<StageTwo, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
}

