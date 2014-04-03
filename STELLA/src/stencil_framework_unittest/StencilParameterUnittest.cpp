#include "gtest/gtest.h"
#include "StencilFramework.h"

enum { sum, result3d, result2d, result1d, par3d, par2d, par1d, parscalar, par_in, par_out, data };

// access test stages
template<typename TEnv>
struct Delta 
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, offset)
    FUNCTION_PARAMETER(1, data)

    __ACC__
    static T Do(Context ctx)
    {
        return ctx[data::At(offset())] - ctx[data::Center()];
    }
};

template<typename TEnv>
struct DataAccessStencil
{
    STENCIL_STAGE(TEnv)
 
    STAGE_PARAMETER(FullDomain, sum)
    STAGE_PARAMETER(FullDomain, result3d)
    STAGE_PARAMETER(FullDomain, result2d)
    STAGE_PARAMETER(FullDomain, result1d)
    STAGE_PARAMETER(FullDomain, par3d)
    STAGE_PARAMETER(FullDomain, par2d)
    STAGE_PARAMETER(FullDomain, par1d)
    STAGE_PARAMETER(FullDomain, parscalar)

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        // sum up a scalar field
        ctx[sum::At(kplus1)] += ctx[par3d::Center()];
        
        // read in the scalar for test purposes
        ctx[result3d::Center()] = ctx[parscalar::Center()];
        ctx[result2d::Center()] = ctx[parscalar::At(Offset<3,2,1>())];
    }

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        // sum up a scalar field
        ctx[sum::At(kminus1)] += ctx[par3d::Center()];
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        // check that the 3d field varies in all dimensions
        ctx[result3d::Center()] =
            ctx[Call<Delta>::With(kplus1,par3d::Center())] + // 100
            ctx[Call<Delta>::With(jplus1,par3d::Center())] + // 10
            ctx[Call<Delta>::With(iplus1,par3d::Center())];  // 1

        // check that the 2d field varies in all dimensions
        ctx[result2d::Center()] =
            ctx[Call<Delta>::With(kplus1,par2d::Center())] + // 0
            ctx[Call<Delta>::With(jplus1,par2d::Center())] + // 10
            ctx[Call<Delta>::With(iplus1,par2d::Center())];  // 1

        // check that the 1d field varies in all dimensions
        ctx[result1d::Center()] =
            ctx[Call<Delta>::With(kplus1,par1d::Center())] + // 0
            ctx[Call<Delta>::With(jplus1,par1d::Center())] + // 0
            ctx[Call<Delta>::With(iplus1,par1d::Center())];  // 1

        // sum up a scalar field
        ctx[sum::Center()] += ctx[par3d::Center()];
    }
};

// test that the framework is able to access data fields 1-3d and scalar parameters
TEST(StencilParameterUnittest, Access)
{
    IJKSize size;
    size.Init(8, 12, 20);
    
    KBoundary boundary; 
    boundary.Init(0,0); 

    // setup three data fields
    IJKRealField par3d;
    IJRealField par2d;
    IRealField par1d;
    IJKRealField result3d;
    IJKRealField result2d;
    IJKRealField result1d;
    IJRealField sum;
    IJRealField checksum;

    par3d.Init("par3d", size, boundary);
    par2d.Init("par2d", size, boundary);
    par1d.Init("par1d", size, boundary);
    result3d.Init("result3d", size, boundary);
    result2d.Init("result2d", size, boundary);
    result1d.Init("result1d", size, boundary);  
    sum.Init("sum", size, boundary);
    checksum.Init("checksum", size, boundary);

    Real parscalar;
    parscalar = (Real)3.14;
    for(int i=0; i<size.iSize()+1; i++)
    {
        for(int j=0; j<size.jSize()+1; j++)
        {
            checksum(i,j,0) = (Real)0.0;
            sum(i,j,0) = (Real)0.0;
            
            for(int k=0; k<size.kSize(); k++)
            {
                par3d(i,j,k) = (Real)i + (Real)j * (Real)10.0 + (Real)k * (Real)100.0;
                par2d(i,j,k) = (Real)i + (Real)j * (Real)10.0;
                par1d(i,j,k) = (Real)i;
                
                checksum(i,j,0) += par3d(i,j,k);
            }
        }
    }
    
    // setup and run the test stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::sum, cInOut>(sum),
            Param< ::result3d, cInOut>(result3d),
            Param< ::result2d, cInOut>(result2d),
            Param< ::result1d, cInOut>(result1d),
            Param< ::par3d, cInOut>(par3d),
            Param< ::par2d, cInOut>(par2d),
            Param< ::par1d, cInOut>(par1d),
            Param< ::parscalar, cScalar>(parscalar)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<DataAccessStencil, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
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
            // test the check sum
            ASSERT_DOUBLE_EQ(checksum(i,j,0), sum(i,j,0)) << "i=" << i << "; j=" << j;
            ASSERT_DOUBLE_EQ(parscalar, result3d(i,j,0)) << "i=" << i << "; j=" << j;
            ASSERT_DOUBLE_EQ(parscalar, result2d(i,j,0)) << "i=" << i << "; j=" << j;

            for(int k=1; k<size.kSize()-1; k++)
            {
                // check the deltas
                ASSERT_DOUBLE_EQ(111.0, result3d(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
                ASSERT_DOUBLE_EQ(11.0, result2d(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
                ASSERT_DOUBLE_EQ(1.0, result1d(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
            }
        }
    }
}

// double buffering test stages
template<typename TEnv>
struct ParameterInit
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, par_in)
    STAGE_PARAMETER(FullDomain, par_out)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[par_out::Center()] = (T)1.0;
    }
};

template<typename TEnv>
struct DoubleBufferFunction
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, data)

    __ACC__
    static T Do(Context ctx)
    {
        if( (ctx[data::At(iminus1)] == 0) || (ctx[data::At(iplus1)] == 0) ||
            (ctx[data::At(jminus1)] == 0) || (ctx[data::At(jplus1)] == 0) ||
            (ctx[data::At(kminus1)] == 0) || (ctx[data::At(kplus1)] == 0)   )
        {
            return 0;
        }
        else
        {
            return (
                ctx[data::Center()]
                + ctx[data::At(iminus1)] + ctx[data::At(iplus1)]
                + ctx[data::At(jminus1)] + ctx[data::At(jplus1)]
                + ctx[data::At(kminus1)] + ctx[data::At(kplus1)] 
            ) / (Real)7.0 + 1;
        }
    }
};

template<typename TEnv>
struct DoubleBufferUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, par_in)
    STAGE_PARAMETER(FullDomain, par_out)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[par_out::Center()] = ctx[Call<DoubleBufferFunction>::With(par_in::Center())];
    }
};

// test if the double buffering works
TEST(StencilParameterUnittest, DoubleBuffering)
{
    const int isize = 17;
    const int jsize = 21;
    const int ksize = 31;
    const int kboundary = 1;
    const int cNumRounds = 4;

    IJKSize size;
    size.Init(isize, jsize, ksize);
    KBoundary boundary;
    boundary.Init(-kboundary, kboundary);

    SwapDataField<IJKRealField> par;
    par.Init("par", size, boundary);

    // init the fields to 0
    for(int i = -cNumBoundaryLines; i < isize + cNumBoundaryLines; ++i)
    {
        for(int j = -cNumBoundaryLines; j < jsize + cNumBoundaryLines; ++j)
        {
            for(int k = -kboundary; k < ksize + kboundary; ++k)
            {
                par.in()(i,j,k) = 0.0;
                par.out()(i,j,k) = 0.0;
            }
        }
    }

    // fill par out with 1.0
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<par_in, cIn>(par.in()),
            Param<par_out, cInOut>(par.out())
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<ParameterInit, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )  
        )
    );
    testStencil.Apply();
    par.Swap();

    // apply filter stencil
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(), 
        pack_parameters(
            Param<par_in, cIn>(par.in()),
            Param<par_out, cInOut>(par.out())
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<DoubleBufferUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )           
        )
    );

    // apply DoubleBufferTestUpdate cNumRounds times
    // par_in:  1  1  1  1  par_out:  0  2  2  2  2  round:  0
    // par_in:  0  2  2  2  par_out:  0  0  3  3  3  round:  1
    // par_in:  0  0  3  3  par_out:  0  0  0  4  4  round:  2
    // par_in:  0  0  0  4  par_out:  0  0  0  0  5  round:  3
    // etc.
    for (int round = 0; round < cNumRounds; ++round) 
    {
        testStencil.Apply();
        par.Swap();
    }

    // checking results
    for(int i = -cNumBoundaryLines; i < isize + cNumBoundaryLines; ++i) 
    {
        for(int j = -cNumBoundaryLines; j < jsize + cNumBoundaryLines; ++j)
        {
            for(int k = -kboundary; k < ksize + kboundary; ++k)
            {
                if( (i < cNumRounds) || (i >= isize - cNumRounds) ||
                    (j < cNumRounds) || (j >= jsize - cNumRounds) ||
                    (k < cNumRounds) || (k >= ksize - cNumRounds)  )
                {
                    ASSERT_DOUBLE_EQ(0, par.in()(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
                }
                else
                {
                    ASSERT_DOUBLE_EQ(cNumRounds + 1, par.in()(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
                }
            }
        }
    }
}

// test that joker fields work
TEST(StencilParameterUnittest, JokerFields)
{
    const int isize = 17;
    const int jsize = 21;
    const int ksize = 31;

    IJKSize size;
    size.Init(isize, jsize, ksize);
    KBoundary boundary;
    boundary.Init(0, 0);

    JokerDataField<IJKRealField> joker;
    IJKRealField par0;
    IJKRealField par1;
    
    // init parameter fields
    par0.Init("par0", size, boundary);
    par1.Init("par1", size, boundary);
    joker.Init("joker", par0);

    // init the fields to 0
    for(int i = -cNumBoundaryLines; i < isize + cNumBoundaryLines; ++i)
    {
        for(int j = -cNumBoundaryLines; j < jsize + cNumBoundaryLines; ++j)
        {
            for(int k = 0; k < ksize; ++k)
            {
                par0(i,j,k) = 0.0;
                par1(i,j,k) = 0.0;
            }
        }
    }

    // fill joker out with 1.0
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<par_in, cIn>(par0),
            Param<par_out, cInOut>(joker)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<ParameterInit, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )  
        )
    );
    testStencil.Apply();

    // make sure par0 == joker is filled with ones, while par 1 is still set to 1.0
    // init the fields to 0
    for(int i = 0; i < isize; ++i)
    {
        for(int j = 0; j < jsize; ++j)
        {
            for(int k = 0; k < ksize; ++k)
            {
                ASSERT_DOUBLE_EQ(1.0, par0(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
                ASSERT_DOUBLE_EQ(0.0, par1(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
            }
        }
    }

    // apply the stencil to par1 as well
    joker.set_dataField(par1);
    testStencil.Apply();

    // init the fields to 1
    for(int i = 0; i < isize; ++i)
    {
        for(int j = 0; j < jsize; ++j)
        {
            for(int k = 0; k < ksize; ++k)
            {
                ASSERT_DOUBLE_EQ(1.0, par0(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
                ASSERT_DOUBLE_EQ(1.0, par1(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
            }
        }
    }
}

template<typename TEnv>
struct OffsetUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, par_in)
    STAGE_PARAMETER(FullDomain, par_out)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[par_out::Center()] = 
            ctx[par_in::At(center)] +
            ctx[par_in::At(iminus3)] + ctx[par_in::At(iminus2)] + ctx[par_in::At(iminus1)] +
            ctx[par_in::At(iplus3)] + ctx[par_in::At(iplus2)] + ctx[par_in::At(iplus1)] +
            ctx[par_in::At(jminus3)] + ctx[par_in::At(jminus2)] + ctx[par_in::At(jminus1)] +
            ctx[par_in::At(jplus3)] + ctx[par_in::At(jplus2)] + ctx[par_in::At(jplus1)] +
            ctx[par_in::At(kminus3)] + ctx[par_in::At(kminus2)] + ctx[par_in::At(kminus1)] +
            ctx[par_in::At(kplus3)] + ctx[par_in::At(kplus2)] + ctx[par_in::At(kplus1)];
    }
};

// test the predefined offsets
TEST(StencilParameterUnittest, Offsets)
{
    const int isize = 17;
    const int jsize = 21;
    const int ksize = 31;

    IJKSize size;
    size.Init(isize, jsize, ksize);
    KBoundary boundary;
    boundary.Init(-3, 3);

    SwapDataField<IJKRealField> par;
    par.Init("par", size, boundary);

    // init the fields to 0
    for(int i = -cNumBoundaryLines; i < isize + cNumBoundaryLines; ++i)
    {
        for(int j = -cNumBoundaryLines; j < jsize + cNumBoundaryLines; ++j)
        {
            for(int k = -3; k < ksize + 3; ++k)
            {
                par.in()(i,j,k) = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                par.out()(i,j,k) = 0.0;
            }
        }
    }

    // apply offset update
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<par_in, cIn>(par.in()),
            Param<par_out, cInOut>(par.out())
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<OffsetUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >() 
                )
            )  
        )
    );
    testStencil.Apply();

    // verify the result
    for(int i = 0; i < isize; ++i)
    {
        for(int j = 0; j < jsize; ++j)
        {
            for(int k = 0; k < ksize; ++k)
            {
                Real expected = 
                    par.in()(i, j, k) +
                    par.in()(i-3, j, k) + par.in()(i-2, j, k) + par.in()(i-1, j, k) + 
                    par.in()(i+1, j, k) + par.in()(i+2, j, k) + par.in()(i+3, j, k) + 
                    par.in()(i, j-3, k) + par.in()(i, j-2, k) + par.in()(i, j-1, k) + 
                    par.in()(i, j+1, k) + par.in()(i, j+2, k) + par.in()(i, j+3, k) +
                    par.in()(i, j, k-3) + par.in()(i, j, k-2) + par.in()(i, j, k-1) + 
                    par.in()(i, j, k+1) + par.in()(i, j, k+2) + par.in()(i, j, k+3); 

                ASSERT_DOUBLE_EQ(expected, par.out()(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
            }
        }
    }
}
