#include <boost/mpl/vector/vector0.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>
#include "gtest/gtest.h"
#include "DummyContext.h"
#include "StencilFramework.h"
#include "HasDo.h"

enum { randomField, out1, out2, out3 };

template<typename TEnv>
struct NoDo
{  
    STENCIL_FUNCTION(TEnv)
};

template<typename TEnv>
struct Average
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, offset)
    FUNCTION_PARAMETER(1, data)

    __ACC__
    static T Do(Context ctx)
    {
        return (Real)0.5 * (ctx[data::At(offset())] + ctx[data::Center()]);
    }
};

template<typename TEnv>
struct AverageFlatCoordinates
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, offset)
    FUNCTION_PARAMETER(1, data)

    __ACC__
    static T Do(Context ctx, FlatCoordinates)
    {
        return (Real)0.5 * (ctx[data::At(offset())] + ctx[data::Center()]);
    }
};

template<typename TEnv>
struct EmptyAverage
{  
    STENCIL_FUNCTION(TEnv)

    __ACC__
    static void Do(Context ctx) {}
};

TEST(StencilFunctionUnittest, HasDo)
{
    typedef StencilFunctionEnvironment<DummyContext, boost::mpl::vector0<> > DummyStencilFunctionEnvironment;

    // check that has do method works    
    ASSERT_TRUE(has_do_member<Average<DummyStencilFunctionEnvironment> >::value);
    ASSERT_TRUE(has_do_member<EmptyAverage<DummyStencilFunctionEnvironment> >::value);
    ASSERT_FALSE(has_do_member<NoDo<DummyStencilFunctionEnvironment> >::value);

    // check no do method is detected
    ASSERT_FALSE( (has_do<NoDo<DummyStencilFunctionEnvironment>, void, DummyContext, boost::mpl::void_>::value) );
    ASSERT_FALSE( (has_do<NoDo<DummyStencilFunctionEnvironment>, double, DummyContext, boost::mpl::void_>::value) );
    ASSERT_FALSE( (has_do<NoDo<DummyStencilFunctionEnvironment>, void, DummyContext, FullDomain>::value) );
    ASSERT_FALSE( (has_do<NoDo<DummyStencilFunctionEnvironment>, double, DummyContext, FullDomain>::value) );
  
    // verify a method returning void is detected
    ASSERT_TRUE( (has_do<EmptyAverage<DummyStencilFunctionEnvironment>, void, DummyContext, boost::mpl::void_>::value) );
    ASSERT_FALSE( (has_do<EmptyAverage<DummyStencilFunctionEnvironment>, double, DummyContext, boost::mpl::void_>::value) );
    ASSERT_FALSE( (has_do<EmptyAverage<DummyStencilFunctionEnvironment>, void, DummyContext, FullDomain>::value) );
    ASSERT_FALSE( (has_do<EmptyAverage<DummyStencilFunctionEnvironment>, double, DummyContext, FullDomain>::value) );

    // verify a do method returning a value is detected
    ASSERT_TRUE( (has_do<Average<DummyStencilFunctionEnvironment>, double, DummyContext, boost::mpl::void_>::value) );
    ASSERT_FALSE( (has_do<Average<DummyStencilFunctionEnvironment>, void, DummyContext, boost::mpl::void_>::value) );
    ASSERT_FALSE( (has_do<Average<DummyStencilFunctionEnvironment>, double, DummyContext, FullDomain>::value) );
    ASSERT_FALSE( (has_do<Average<DummyStencilFunctionEnvironment>, void, DummyContext, FullDomain>::value) );

    // check do method returning a value at a given domain is detected
    ASSERT_TRUE( (has_do<AverageFlatCoordinates<DummyStencilFunctionEnvironment>, double, DummyContext, FlatCoordinates>::value) );
    ASSERT_FALSE( (has_do<AverageFlatCoordinates<DummyStencilFunctionEnvironment>, double, DummyContext, TerrainCoordinates>::value) );
    ASSERT_FALSE( (has_do<AverageFlatCoordinates<DummyStencilFunctionEnvironment>, double, DummyContext, boost::mpl::void_>::value) );
    ASSERT_FALSE( (has_do<AverageFlatCoordinates<DummyStencilFunctionEnvironment>, void, DummyContext, TerrainCoordinates>::value) );
    ASSERT_FALSE( (has_do<AverageFlatCoordinates<DummyStencilFunctionEnvironment>, void, DummyContext, boost::mpl::void_>::value) );
}

TEST(StencilFunctionUnittest, HasParameter)
{
    typedef StencilFunctionEnvironment<DummyContext, boost::mpl::vector0<> > DummyStencilFunctionEnvironment;

    // check has parameter returns true
    ASSERT_TRUE( (has_parameter<Average<DummyStencilFunctionEnvironment>, boost::mpl::integral_c<int,0> >::value) );
    ASSERT_TRUE( (has_parameter<Average<DummyStencilFunctionEnvironment>, boost::mpl::integral_c<int,1> >::value) );

    // check has parameter returns false
    ASSERT_FALSE( (has_parameter<Average<DummyStencilFunctionEnvironment>, boost::mpl::integral_c<int,2> >::value) );
    ASSERT_FALSE( (has_parameter<Average<DummyStencilFunctionEnvironment>, boost::mpl::integral_c<int,3> >::value) );
}

template<typename TEnv>
struct AverageWrongParameters1
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, offset)
    FUNCTION_PARAMETER(2, data)
    FUNCTION_PARAMETER(3, nix)
};

template<typename TEnv>
struct AverageWrongParameters2
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(1, offset)
    FUNCTION_PARAMETER(2, data)
};

TEST(StencilFunctionUnittest, CountParameters)
{
    // test parameter list length
    ASSERT_EQ(2, count_stencil_function_parameters<Average>::type::value);

    ASSERT_TRUE(boost::mpl::is_void_<count_stencil_function_parameters<AverageWrongParameters1>::type>::type::value);
    ASSERT_TRUE(boost::mpl::is_void_<count_stencil_function_parameters<AverageWrongParameters2>::type>::type::value);
}

template<typename TEnv>
struct StencilStageUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, randomField)
    STAGE_PARAMETER(FullDomain, out1)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[out1::Center()] = (Real)0.5 * (
            (Real)0.5 * (ctx[randomField::At(Offset<1,1,0>())] +  ctx[randomField::At(iplus1)]) +  
            (Real)0.5 * (ctx[randomField::At(jplus1)] +  ctx[randomField::Center()])
        );
    }
};

template<typename TEnv>
struct StencilFunctionUpdate1
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, randomField)
    STAGE_PARAMETER(FullDomain, out2)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[out2::Center()] = ctx[Call<Average>::With(iplus1, Call<Average>::With(jplus1, randomField::Center()))];
    }
};

template<typename TEnv>
struct StencilFunctionUpdate2
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, randomField)
    STAGE_PARAMETER(FullDomain, out3)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[out3::Center()] = (Real)0.5 * (
            ctx[Call<Average>::With(jplus1, randomField::At(iplus1))] +  
            ctx[Call<Average>::With(jplus1, randomField::Center())] 
        );
    }
};

TEST(StencilFunctionUnittest, Execution)
{
    IJKSize size;
    size.Init(3, 7, 25);

    KBoundary boundary;
    boundary.Init(0, 0);

    // setup three data fields
    IJKRealField randomField;
    IJKRealField out1;
    IJKRealField out2;
    IJKRealField out3;
 
    randomField.Init("randomField", size, boundary);
    out1.Init("out1", size, boundary);
    out2.Init("out2", size, boundary);
    out3.Init("out3", size, boundary);
    
    for(int i = -cNumBoundaryLines; i < size.iSize() + cNumBoundaryLines; ++i)
    {
        for(int j = -cNumBoundaryLines; j < size.jSize() + cNumBoundaryLines; ++j)
        {
            for(int k = boundary.kMinusOffset(); k < size.kSize() + boundary.kPlusOffset(); ++k)
            {
                randomField(i,j,k) = i * (Real)0.34 + j * (Real)23.43 + k * (Real)1324.3;
            }
        }
    }

    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param< ::randomField, cInOut>(randomField),
            Param< ::out1, cInOut>(out1),
            Param< ::out2, cInOut>(out2),
            Param< ::out3, cInOut>(out3)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StencilStageUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            ),
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StencilFunctionUpdate1, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            ),
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StencilFunctionUpdate2, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // run the same logic by hand and verify the result
    for(int i=0; i<size.iSize(); i++)
    {
        for(int j=0; j<size.jSize(); j++)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                // calculate result
                Real result = (Real)0.5 * (
                    (Real)0.5 * (randomField(i+1,j+1,k) + randomField(i+1,j,k)) +
                    (Real)0.5 * (randomField(i,j+1,k) + randomField(i,j,k)) 
                );
                
                ASSERT_DOUBLE_EQ(result, out1(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
                ASSERT_DOUBLE_EQ(result, out2(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
                ASSERT_DOUBLE_EQ(result, out3(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
            }
        }
    }
}

template<typename TEnv>
struct DomainSpecificAverageFunction
{
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, randomField)
    FUNCTION_PARAMETER(1, out1)

    __ACC__
    static void Do(Context ctx, TerrainCoordinates)
    {
        ctx[out1::Center()] = (Real)0.1 * (
            ctx[Call<Average>::With(jplus1, randomField::At(iplus1))] +  
            ctx[Call<Average>::With(jplus1, randomField::Center())] 
        );
    }

    __ACC__
    static void Do(Context ctx, FlatCoordinates)
    {
        ctx[out1::Center()] = (Real)0.2 * (
            ctx[Call<Average>::With(jplus1, randomField::At(iplus1))] +  
            ctx[Call<Average>::With(jplus1, randomField::Center())] 
        );
    }
};

template<typename TEnv>
struct DomainSpecificAverageStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, randomField)
    STAGE_PARAMETER(FullDomain, out1)

    __ACC__
    static void Do(Context ctx, TerrainCoordinates)
    {
        ctx[Call<DomainSpecificAverageFunction, TerrainCoordinates>::With(randomField::Center(), out1::Center())];
    }

    __ACC__
    static void Do(Context ctx, FlatCoordinates)
    {
        ctx[Call<DomainSpecificAverageFunction, FlatCoordinates>::With(randomField::Center(), out1::Center())];
    }
};

TEST(StencilFunctionUnittest, Domains)
{
    IJKSize size;
    size.Init(3, 7, 25);

    KBoundary boundary;
    boundary.Init(0, 0);

    // setup three data fields
    IJKRealField randomField;
    IJKRealField out1;
 
    randomField.Init("randomField", size, boundary);
    out1.Init("out1", size, boundary);
    
    for(int i = -cNumBoundaryLines; i < size.iSize() + cNumBoundaryLines; ++i)
    {
        for(int j = -cNumBoundaryLines; j < size.jSize() + cNumBoundaryLines; ++j)
        {
            for(int k = boundary.kMinusOffset(); k < size.kSize() + boundary.kPlusOffset(); ++k)
            {
                randomField(i,j,k) = i * (Real)0.34 + j * (Real)23.43 + k * (Real)1324.3;
            }
        }
    }

    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param< ::randomField, cInOut>(randomField),
            Param< ::out1, cInOut>(out1)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<DomainSpecificAverageStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // run the same logic by hand and verify the result
    for(int i=0; i<size.iSize(); i++)
    {
        for(int j=0; j<size.jSize(); j++)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                Real result = 
                    (Real)0.5 * (randomField(i+1,j+1,k) + randomField(i+1,j,k)) +
                    (Real)0.5 * (randomField(i,j+1,k) + randomField(i,j,k)); 

                if(k < cFlatLimit)
                {
                    result *= 0.2;
                }
                else
                {
                    result *= 0.1;
                }
                
                ASSERT_DOUBLE_EQ(result, out1(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
            }
        }
    }
}


