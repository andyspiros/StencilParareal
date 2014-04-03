#include "gtest/gtest.h"
#include "StencilFramework.h"

#include <boost/mpl/vector/vector10.hpp>
#include <boost/mpl/print.hpp>

// define a parameter enum
enum { data, variant, sign };

template<typename TEnv>
struct InitStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data::Center()] = 0.0;
    }
};

template<typename TEnv>
struct AddTwoStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data::Center()] += 2.0;
    }
};

template<typename TEnv>
struct AddFourStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data::Center()] += 4.0;
    }
};

template<typename TEnv>
struct TimesTenStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data::Center()] *= 10.0;
    }
};

template<typename TEnv>
struct SignStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data::Center()] *= -1.0;
    }
};

// define a variant enum
enum Variant 
{
    cAddNothing,
    cAddTwo,
    cAddFour,
    cAddFourAndTwo
};

TEST(StencilSwitchCaseUnittest, SwitchCase)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);

    // define a result field
    IJKRealField data;
    data.Init("data", size, boundary);

    // define a variant
    Variant variant = cAddNothing;
    bool sign = false;

    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::data, cInOut>(data),
            Param< ::variant, cScalar>(variant),
            Param< ::sign, cScalar>(sign)
        ),
        define_loops(
            // always initialize to 0
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<InitStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            ),
            // first switch case using an enum, either add 2, 4 or 6
            define_switch< ::variant>(
                define_case<Variant, cAddTwo>(
                    define_sweep<cKIncrement>(
                        define_stages(
                            StencilStage<AddTwoStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                        )
                    )
                ),
                define_case<Variant, cAddFour>(
                    define_sweep<cKIncrement>(
                        define_stages(
                            StencilStage<AddFourStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                        )
                    )
                ),
                define_case<Variant, cAddFourAndTwo>(
                    define_sweep<cKIncrement>(
                        define_stages(
                            StencilStage<AddTwoStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >(),
                            StencilStage<AddFourStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                        )
                    )
                )
            ),
            // always multiply by 10
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<TimesTenStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >() 
                )
            ),
            // second switch case changing the sign if sign is set to true
            define_switch< ::sign>(
                define_case<bool, true>(
                    define_sweep<cKIncrement>(
                        define_stages(
                            StencilStage<SignStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                        )
                    )
                )
            )
        )
    );

    // check data is set to 4 * 10 * -1
    variant = cAddFour;
    sign = true;
    testStencil.Apply();
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(-40.0, data(i,j,k));
            }
        }
    }

    // check data is set to 6 * 10
    variant = cAddFourAndTwo;
    sign = false;
    testStencil.Apply();
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(60.0, data(i,j,k));
            }
        }
    }

    // check data is set to 2 * 10
    variant = cAddTwo;
    sign = false;
    testStencil.Apply();
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(20.0, data(i,j,k));
            }
        }
    }

    // check data is set to 0 * 10 * -1
    variant = cAddNothing;
    sign = true;
    testStencil.Apply();
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(-0.0, data(i,j,k));
            }
        }
    }
}