#include "gtest/gtest.h"
#include <boost/mpl/at.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/equal.hpp>
#include "StencilFramework.h"
#include "StencilSweepFunctorOperationGroup.h"

enum { data, data_in, data_out, weight, intermediate };

// function used in order to compute the average field
Real compute_average(IJKRealField& data, IJKRealField weight, int i, int j, int k)
{
    return (1.0/3.0) * (
        data(i,j,k+1) * weight(i,j,k+1) +
        data(i,j,k) * weight(i,j,k) +
        data(i,j,k-1) * weight(i,j,k-1)
    ); 
}

// stage computing a weighted average
template<typename TEnv>
struct AveragingStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data)
    STAGE_PARAMETER(FullDomain, weight)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data::Center()] =  
            (1.0/3.0) * (
            ctx[data::At(kplus1)] * ctx[weight::At(kplus1)] +
            ctx[data::At(center)] * ctx[weight::At(center)] +
            ctx[data::At(kminus1)] * ctx[weight::At(kminus1)] 
        );
    }
};

// stage computing a weighted forward average using
template<typename TEnv>
struct ForwardAveragingStage
{
    STENCIL_STAGE(TEnv)
    
    STAGE_PARAMETER(FullDomain, intermediate)
    STAGE_PARAMETER(FullDomain, data)
    STAGE_PARAMETER(FullDomain, weight)
    
    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        computeAverage(ctx);

        // finalize the loop
        ctx[data::At(kplus2)] = 2.0;
        ctx[data::At(kplus1)] = 1.0;
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        computeAverage(ctx);
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        // prepare the loop
        ctx[intermediate::At(center)] = 
            ctx[data::At(center)] *
            ctx[weight::At(center)];
        ctx[intermediate::At(kminus1)] = 
            ctx[data::At(kminus1)] *
            ctx[weight::At(kminus1)];
        ctx[data::At(kminus1)] = -1.0;

        computeAverage(ctx);
    }

private:
    __ACC__
    static void computeAverage(Context ctx)
    {
        // multiply data and factor
        ctx[intermediate::At(kplus1)] = 
            ctx[data::At(kplus1)] *
            ctx[weight::At(kplus1)];

        ctx[data::Center()] =  
            (1.0/3.0) * (
            ctx[intermediate::At(kplus1)] +
            ctx[intermediate::At(center)] +
            ctx[intermediate::At(kminus1)]
        );
    }
};

// test computing a forward average using multiple caches with extended domain
TEST(StencilCacheUnittest, KCachingForwardLargeKRange)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);
    
    // define a result field
    IJKRealField data, weight, result;
    data.Init("data", size, boundary);
    weight.Init("weight", size, boundary);
    result.Init("result", size, boundary);

    // initialize data fields randomly
    for(int i=-cNumBoundaryLines; i<size.iSize()+cNumBoundaryLines; ++i)
    {
        for(int j=-cNumBoundaryLines; j<size.jSize()+cNumBoundaryLines; ++j)
        {
            for(int k=-1; k<size.kSize()+2; ++k)
            {
                data(i, j, k) = (Real)10.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                weight(i, j, k) = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
            }
        }
    }

    // compute the expected result
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); ++k)
            {
                result(i, j, k) = compute_average(data, weight, i, j, k);
            }
            result(i, j, -1) = -1.0;
            result(i, j, size.kSize()) = 1.0;
            result(i, j, size.kSize() + 1) = 2.0;
        }
    }

    // setup the stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::data, cInOut>(data),
            Param< ::weight, cIn>(weight)
        ),
        define_temporaries(
            StencilBuffer<intermediate, Real, KRange<FullDomain,-1,1> >()
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_caches(
                    KCache< ::data, cFillAndFlush, KWindow<-1,2>, KRange<FullDomain,-1,2> >(),
                    KCache< ::weight, cFill, KWindow<-1,1>, KRange<FullDomain,-1,1> >(),
                    KCache<intermediate, cLocal, KWindow<-1,1>, KRange<FullDomain,-1,1> >()
                ),
                define_stages(
                    StencilStage<ForwardAveragingStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // verify the result
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=-1; k<size.kSize() + 2; ++k)
            {
                ASSERT_DOUBLE_EQ(result(i, j, k), data(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 
            }
        }
    }
}

// test computing a forward average using multiple caches with extended domain
// (only loop over the flat coordinates and the k maximum domain)
TEST(StencilCacheUnittest, KCachingBackwardMultipleKLoopLegs)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);
    
    // define a result field
    IJKRealField data, weight, result;
    data.Init("data", size, boundary);
    weight.Init("weight", size, boundary);
    result.Init("result", size, boundary);

    // initialize data fields randomly
    for(int i=-cNumBoundaryLines; i<size.iSize()+cNumBoundaryLines; ++i)
    {
        for(int j=-cNumBoundaryLines; j<size.jSize()+cNumBoundaryLines; ++j)
        {
            for(int k=-1; k<size.kSize()+2; ++k)
            {
                data(i, j, k) = (Real)10.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                weight(i, j, k) = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
            }
        }
    }

    // compute the expected result do not work in the terrain coordinate domain
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); ++k)
            {
                if(k < cFlatLimit || k == size.kSize()-1)
                {
                    result(i, j, k) = compute_average(data, weight, i, j, k);
                }
                else
                {
                    result(i, j, k) = data(i, j, k); 
                }
            }
            result(i, j, -1) = -1.0;
            result(i, j, size.kSize()) = data(i, j, size.kSize());
            result(i, j, size.kSize() + 1) = data(i, j, size.kSize() + 1);
        }
    }

    // setup the stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::data, cInOut>(data),
            Param< ::weight, cIn>(weight)
        ),
        define_temporaries(
            StencilBuffer<intermediate, Real, KRange<FullDomain,-1,1> >()
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_caches(
                    KCache< ::data, cFillAndFlush, KWindow<-1,2>, KRange<FullDomain,-1,2> >(),
                    KCache< ::weight, cFill, KWindow<-1,1>, KRange<FullDomain,-1,1> >(),
                    KCache<intermediate, cLocal, KWindow<-1,1>, KRange<FullDomain,-1,1> >()
                ),
                define_stages(
                    StencilStage<ForwardAveragingStage, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,0,0> >(),
                    StencilStage<AveragingStage, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // verify the result
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=-1; k<size.kSize() + 2; ++k)
            {
                ASSERT_DOUBLE_EQ(result(i, j, k), data(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 
            }
        }
    }
}

// stage computing a weighted backward average using
template<typename TEnv>
struct BackwardAveragingStage
{
    STENCIL_STAGE(TEnv)
    
    STAGE_PARAMETER(FullDomain, intermediate)
    STAGE_PARAMETER(FullDomain, data)
    STAGE_PARAMETER(FullDomain, weight)
    
    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        // prepare the loop
        ctx[intermediate::At(kplus1)] = 
            ctx[data::At(kplus1)] *
            ctx[weight::At(kplus1)];
        ctx[intermediate::At(center)] = 
            ctx[data::At(center)] *
            ctx[weight::At(center)];
        ctx[data::At(kplus2)] = 2.0;
        ctx[data::At(kplus1)] = 1.0;

        computeAverage(ctx);
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        computeAverage(ctx);
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        computeAverage(ctx);
        
        // finalize the loop
        ctx[data::At(kminus1)] = -1.0;
    }

private:
    __ACC__
    static void computeAverage(Context ctx)
    {
        // multiply data and factor
        ctx[intermediate::At(kminus1)] = 
            ctx[data::At(kminus1)] *
            ctx[weight::At(kminus1)];

        ctx[data::Center()] =  
            (1.0/3.0) * (
            ctx[intermediate::At(kplus1)] +
            ctx[intermediate::At(center)] +
            ctx[intermediate::At(kminus1)]
        );
    }
};

// test computing a backward average using multiple caches with extended domain
TEST(StencilCacheUnittest, KCachingBackwardLargeKRange)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);
    
    // define a result field
    IJKRealField data, weight, result;
    data.Init("data", size, boundary);
    weight.Init("weight", size, boundary);
    result.Init("result", size, boundary);

    // initialize data fields randomly
    for(int i=-cNumBoundaryLines; i<size.iSize()+cNumBoundaryLines; ++i)
    {
        for(int j=-cNumBoundaryLines; j<size.jSize()+cNumBoundaryLines; ++j)
        {
            for(int k=-1; k<size.kSize()+2; ++k)
            {
                data(i, j, k) = (Real)10.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                weight(i, j, k) = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
            }
        }
    }

    // compute the expected result
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); ++k)
            {
                result(i, j, k) = compute_average(data, weight, i, j, k);
            }
            result(i, j, -1) = -1.0;
            result(i, j, size.kSize()) = 1.0;
            result(i, j, size.kSize() + 1) = 2.0;
        }
    }

    // setup the stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::data, cInOut>(data),
            Param< ::weight, cIn>(weight)
        ),
        define_temporaries(
            StencilBuffer<intermediate, Real, KRange<FullDomain,-1,1> >()
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_caches(
                    KCache< ::data, cFillAndFlush, KWindow<-1,2>, KRange<FullDomain,-1,2> >(),
                    KCache< ::weight, cFill, KWindow<-1,1>, KRange<FullDomain,-1,1> >(),
                    KCache<intermediate, cLocal, KWindow<-1,1>, KRange<FullDomain,-1,1> >()
                ),
                define_stages(
                    StencilStage<BackwardAveragingStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // verify the result
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=-1; k<size.kSize() + 2; ++k)
            {
                ASSERT_DOUBLE_EQ(result(i, j, k), data(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 
            }
        }
    }
}

// test computing a backward average using multiple caches with extended domain
// (only loop over the terrain coordinates and the k minimum domain)
TEST(StencilCacheUnittest, KCachingForwardMultipleKLoopLegs)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);
    
    // define a result field
    IJKRealField data, weight, result;
    data.Init("data", size, boundary);
    weight.Init("weight", size, boundary);
    result.Init("result", size, boundary);

    // initialize data fields randomly
    for(int i=-cNumBoundaryLines; i<size.iSize()+cNumBoundaryLines; ++i)
    {
        for(int j=-cNumBoundaryLines; j<size.jSize()+cNumBoundaryLines; ++j)
        {
            for(int k=-1; k<size.kSize()+2; ++k)
            {
                data(i, j, k) = (Real)10.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                weight(i, j, k) = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
            }
        }
    }

    // compute the expected result do not work in the terrain coordinate domain
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); ++k)
            {
                if(k >= cFlatLimit || k == 0)
                {
                    result(i, j, k) = compute_average(data, weight, i, j, k);
                }
                else
                {
                    result(i, j, k) = data(i, j, k); 
                }
            }
            result(i, j, -1) = data(i, j, -1);
            result(i, j, size.kSize()) = 1.0;
            result(i, j, size.kSize() + 1) = 2.0;
        }
    }

    // setup the stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::data, cInOut>(data),
            Param< ::weight, cIn>(weight)
        ),
        define_temporaries(
            StencilBuffer<intermediate, Real, KRange<FullDomain,-1,1> >()
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_caches(
                    KCache< ::data, cFillAndFlush, KWindow<-1,2>, KRange<FullDomain,-1,2> >(),
                    KCache< ::weight, cFill, KWindow<-1,1>, KRange<FullDomain,-1,1> >(),
                    KCache<intermediate, cLocal, KWindow<-1,1>, KRange<FullDomain,-1,1> >()
                ),
                define_stages(
                    StencilStage<BackwardAveragingStage, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,0> >(),
                    StencilStage<AveragingStage, IJRange<cComplete,0,0,0,0>, KRange<KMinimumCenter,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // verify the result
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=-1; k<size.kSize() + 2; ++k)
            {
                ASSERT_DOUBLE_EQ(result(i, j, k), data(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 
            }
        }
    }
}

// stage computing a weighted intermediate field
template<typename TEnv>
struct IntermediateStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data)
    STAGE_PARAMETER(FullDomain, weight)
    STAGE_PARAMETER(FullDomain, intermediate)

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        // do nothing
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[intermediate::At(kplus1)] = ctx[weight::At(kplus1)] * ctx[data::At(kplus1)];
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        ctx[intermediate::At(kplus1)] = ctx[weight::At(kplus1)] * ctx[data::At(kplus1)];
        ctx[intermediate::At(center)] = ctx[weight::At(center)] * ctx[data::At(center)];
    }
};

// average the intermediate field
template<typename TEnv>
struct AverageIntermediateStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data)
    STAGE_PARAMETER(FullDomain, intermediate)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data::Center()] =  
            (1.0/3.0) * (
            ctx[intermediate::At(kplus1)] +
            ctx[intermediate::At(center)] +
            ctx[intermediate::At(kminus1)] 
        );
    }
};

// compute the average using two stages with different k range
TEST(StencilCacheUnittest, KCachingDifferentStageKRange)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);
    
    // define a result field
    IJKRealField data, weight, result;
    data.Init("data", size, boundary);
    weight.Init("weight", size, boundary);
    result.Init("result", size, boundary);

    // initialize data fields randomly
    for(int i=-cNumBoundaryLines; i<size.iSize()+cNumBoundaryLines; ++i)
    {
        for(int j=-cNumBoundaryLines; j<size.jSize()+cNumBoundaryLines; ++j)
        {
            for(int k=-1; k<size.kSize()+2; ++k)
            {
                data(i, j, k) = (Real)10.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                weight(i, j, k) = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
            }
        }
    }

    // compute the expected result do not work in the terrain coordinate domain
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=1; k<size.kSize()-1; ++k)
            {
                result(i, j, k) = compute_average(data, weight, i, j, k);
            }
            result(i, j, -1) = data(i, j, -1);
            result(i, j, 0) = data(i, j, 0);
            result(i, j, size.kSize() - 1) = data(i, j, size.kSize() - 1);
            result(i, j, size.kSize()) = data(i, j, size.kSize());
            result(i, j, size.kSize() + 1) = data(i, j, size.kSize() + 1);
        }
    }

    // setup the stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::data, cInOut>(data),
            Param< ::weight, cIn>(weight)
        ),
        define_temporaries(
            StencilBuffer<intermediate, Real, KRange<FullDomain,0,0> >()
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_caches(
                    KCache< ::data, cFillAndFlush, KWindow<-1,1>, KRange<FullDomain,0,0> >(),
                    KCache< ::weight, cFill, KWindow<-1,1>, KRange<FullDomain,0,0> >(),
                    KCache<intermediate, cLocal, KWindow<-1,1>, KRange<FullDomain,0,0> >()
                ),
                define_stages(
                    StencilStage<IntermediateStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<AverageIntermediateStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-1> >()
                )
            )
        )
    );
    testStencil.Apply();

    // verify the result
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=-1; k<size.kSize() + 2; ++k)
            {
                ASSERT_DOUBLE_EQ(result(i, j, k), data(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 
            }
        }
    }
}

// function used in order to compute the 3d average field
Real compute_3d_average(IJKRealField& data, int i, int j, int k)
{
    return (1.0/8.0) * (
        data(i+1,j,k-1) + data(i+1,j,k+1) + 
        data(i-1,j,k-1) + data(i-1,j,k+1) + 
        data(i,j+1,k-1) + data(i,j+1,k+1) + 
        data(i,j-1,k-1) + data(i,j-1,k+1) 
    ); 
}

// average data in in k
template<typename TEnv>
struct AverageInKStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data_in)
    STAGE_PARAMETER(FullDomain, intermediate)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[intermediate::Center()] =  
            ctx[data_in::At(kplus1)] +
            ctx[data_in::At(kminus1)] ;
    }
};

// average intermediate in i and j
template<typename TEnv>
struct AverageInIJStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data_out)
    STAGE_PARAMETER(FullDomain, intermediate)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data_out::Center()] =  (1.0/8.0) * (
            ctx[intermediate::At(iplus1)] + ctx[intermediate::At(iminus1)] +
            ctx[intermediate::At(jplus1)] + ctx[intermediate::At(jminus1)]
        );
    }
};

// compute 3d average using an ij and a k cache
TEST(StencilCacheUnittest, IJCachesAndKCaches)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);
    
    // define a result field
    IJKRealField data_in, data_out, result;
    data_in.Init("data_in", size, boundary);
    data_out.Init("data_out", size, boundary);
    result.Init("result", size, boundary);

    // initialize data fields randomly
    for(int i=-cNumBoundaryLines; i<size.iSize()+cNumBoundaryLines; ++i)
    {
        for(int j=-cNumBoundaryLines; j<size.jSize()+cNumBoundaryLines; ++j)
        {
            for(int k=-1; k<size.kSize()+2; ++k)
            {
                data_in(i, j, k) = (Real)10.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                data_out(i, j, k) = data_in(i, j, k);
                result(i, j, k) = data_in(i, j, k);
            }
        }
    }

    // compute the expected result of the stencil
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=1; k<size.kSize()-1; ++k)
            {
                result(i, j, k) = compute_3d_average(data_in, i, j, k);
            }
        }
    }

    // setup the stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::data_out, cInOut>(data_out),
            Param< ::data_in, cIn>(data_in)
        ),
        define_temporaries(
            StencilBuffer<intermediate, Real, KRange<FullDomain,1,-1> >()
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_caches(
                    KCache< ::data_in, cFill, KWindow<-1,1>, KRange<FullDomain,0,0> >(),
                    IJCache< ::intermediate, KRange<FullDomain,0,0> >()
                ),
                define_stages(
                    StencilStage<AverageInKStage, IJRange<cComplete,-1,1,-1,1>, KRange<FullDomain,1,-1> >(),
                    StencilStage<AverageInIJStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-1> >()
                )
            )
        )
    );
    testStencil.Apply();

    // verify the result
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=-1; k<size.kSize() + 2; ++k)
            {
                ASSERT_DOUBLE_EQ(result(i, j, k), data_out(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 
            }
        }
    }
}

// compute the intermediate using data in
template<typename TEnv>
struct IntermediateUsingDataInStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data_in)
    STAGE_PARAMETER(FullDomain, weight)
    STAGE_PARAMETER(FullDomain, intermediate)

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        // do nothing
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[intermediate::At(kplus1)] = ctx[weight::At(kplus1)] * ctx[data_in::At(kplus1)];
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        ctx[intermediate::At(kplus1)] = ctx[weight::At(kplus1)] * ctx[data_in::At(kplus1)];
        ctx[intermediate::At(center)] = ctx[weight::At(center)] * ctx[data_in::At(center)];
    }
};

// average the intermediate field to the data out stage
template<typename TEnv>
struct AverageIntermediateToDataOutStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, data_out)
    STAGE_PARAMETER(FullDomain, intermediate)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[data_out::Center()] =  
            (1.0/3.0) * (
            ctx[intermediate::At(kplus1)] +
            ctx[intermediate::At(center)] +
            ctx[intermediate::At(kminus1)] 
        );
    }
};

// compute the average with some ij dependencies
TEST(StencilCacheUnittest, KCachesAndIJRanges)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(4, 4, 40);
    boundary.Init(-1, 2);
    
    // define a result field
    IJKRealField data_in, data_out, weight, intermediate, result;
    data_in.Init("data_in", size, boundary);
    data_out.Init("data_out", size, boundary);
    weight.Init("weight", size, boundary);
    intermediate.Init("intermediate", size, boundary);
    result.Init("result", size, boundary);

    // initialize data fields randomly
    for(int i=-cNumBoundaryLines; i<size.iSize()+cNumBoundaryLines; ++i)
    {
        for(int j=-cNumBoundaryLines; j<size.jSize()+cNumBoundaryLines; ++j)
        {
            for(int k=-1; k<size.kSize()+2; ++k)
            {
                data_in(i, j, k) = (Real)10.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                data_out(i, j, k) = data_in(i, j, k);
                result(i, j, k) = data_in(i, j, k);
                weight(i, j, k) = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                intermediate(i, j, k) = 7.0;
            }
        }
    }

    // compute the expected result do not work in the terrain coordinate domain
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=1; k<size.kSize()-1; ++k)
            {
                result(i, j, k) = compute_average(data_in, weight, i, j, k);
            }            
        }
    }

    // setup the stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<double, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::data_out, cInOut>(data_out),
            Param< ::data_in, cIn>(data_in),
            Param< ::weight, cIn>(weight),
            Param< ::intermediate, cInOut>(intermediate)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_caches(
                    KCache< ::data_in, cFill, KWindow<-1,1>, KRange<FullDomain,0,0> >(),
                    KCache< ::data_out, cFlush, KWindow<0,0>, KRange<FullDomain,1,-1> >(),
                    KCache< ::weight, cFill, KWindow<-1,1>, KRange<FullDomain,0,0> >(),
                    KCache< ::intermediate, cFlush, KWindow<-1,1>, KRange<FullDomain,0,0> >()
                ),
                define_stages(
                    StencilStage<IntermediateUsingDataInStage, IJRange<cComplete,-1,1,-1,1>, KRange<FullDomain,0,0> >(),
                    StencilStage<AverageIntermediateToDataOutStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-1> >()
                )
            )
        )
    );
    testStencil.Apply();

    // verify the result
    for(int i=-cNumBoundaryLines; i<size.iSize()+cNumBoundaryLines; ++i)
    {
        for(int j=-cNumBoundaryLines; j<size.jSize()+cNumBoundaryLines; ++j)
        {
            for(int k=-1; k<size.kSize() + 2; ++k)
            {
                ASSERT_DOUBLE_EQ(result(i, j, k), data_out(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 

                // additionally test that the intermediate value is set in the right domain!
                if( i >= -1 && i < size.iSize() + 1 && 
                    j >= -1 && j < size.jSize() + 1 && 
                    k >= 0 && k < size.kSize() )
                {
                    // compute the expected intermediate value
                    Real expected = data_in(i, j, k) * weight(i, j, k);

                    ASSERT_DOUBLE_EQ(expected, intermediate(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 
                }
                else
                {
                    ASSERT_DOUBLE_EQ(7.0, intermediate(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k; 
                }
            }
        }
    }
}

// test if we can check properly if a parameter is accessed
TEST(StencilCacheUnittest, IsParameterAccessed)
{
    // define stencil sweep group descriptors without cache
    typedef boost::mpl::vector1<
        StencilSweepGroupDescriptor<
            boost::mpl::vector1<
                StencilSweepDescriptor<
                    boost::mpl::void_,
                    boost::mpl::vector2<
                        StencilStage<AverageInKStage, IJRange<cComplete,-1,1,-1,1>, KRange<FullDomain,1,-1> >,
                        StencilStage<AverageInIJStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-1> >
                    >,
                    cKIncrement
                >
            >,
            boost::mpl::void_,
            -1
        >
    > SweepGroupDescriptorsNoCache;
    // define stencil sweep group descriptors with cache
    typedef boost::mpl::vector1<
        StencilSweepGroupDescriptor<
            boost::mpl::vector1<
                StencilSweepDescriptor<
                    boost::mpl::vector3<
                        KCache<data_in, cFill, KWindow<-1,1>, KRange<FullDomain,0,0> >,
                        KCache<data_out, cLocal, KWindow<-1,1>, KRange<FullDomain,0,0> >,
                        IJCache<intermediate, KRange<FullDomain,0,0> >
                    >,
                    boost::mpl::vector2<
                        StencilStage<AverageInKStage, IJRange<cComplete,-1,1,-1,1>, KRange<FullDomain,1,-1> >,
                        StencilStage<AverageInIJStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-1> >
                    >,
                    cKIncrement
                >
            >,
            boost::mpl::void_,
            -1
        >
    > SweepGroupDescriptorsWithCache;

    // make sure usage is detected at all (no caches)
    ASSERT_TRUE(
        (
            !stencil_sweep_group_descriptors_is_parameter_accessed<
                SweepGroupDescriptorsNoCache,
                FullDomain,
                boost::mpl::integral_c<int, weight>
            >::value
        )
    );
    ASSERT_TRUE(
        (
            stencil_sweep_group_descriptors_is_parameter_accessed<
                SweepGroupDescriptorsNoCache,
                FullDomain,
                boost::mpl::integral_c<int, intermediate>
            >::value
        )
    );

    // make sure caches are considered correctly
    ASSERT_TRUE(
        (
            !stencil_sweep_group_descriptors_is_parameter_accessed<
                SweepGroupDescriptorsWithCache,
                FullDomain,
                boost::mpl::integral_c<int, weight>
            >::value
        )
    );
    ASSERT_TRUE(
        (
            !stencil_sweep_group_descriptors_is_parameter_accessed<
                SweepGroupDescriptorsWithCache,
                FullDomain,
                boost::mpl::integral_c<int, intermediate>
            >::value
        )
    );
    // data in is the only one which is accessed (the others are cached locally)
    ASSERT_TRUE(
        (
            stencil_sweep_group_descriptors_is_parameter_accessed<
                SweepGroupDescriptorsWithCache,
                FullDomain,
                boost::mpl::integral_c<int, data_in>
            >::value
        )
    );
    ASSERT_TRUE(
        (
            !stencil_sweep_group_descriptors_is_parameter_accessed<
                SweepGroupDescriptorsWithCache,
                FullDomain,
                boost::mpl::integral_c<int, data_out>
            >::value
        )
    );
}

// meta function verifying two operation groups are identical
template<
    typename TExpectedOperationGroup,
    typename TActualOperationGroup>
struct check_stencil_sweep_functor_operation_group : boost::mpl::false_ {};

template<
    typename TExpectedOperations,
    typename TActualOperations,
    typename TIJRange>
struct check_stencil_sweep_functor_operation_group<
    StencilSweepFunctorOperationGroup<TExpectedOperations, TIJRange>,
    StencilSweepFunctorOperationGroup<TActualOperations, TIJRange> > : 
    boost::mpl::equal<TExpectedOperations, TActualOperations> 
{};

// test sweep function operation group computation
TEST(StencilCacheUnittest, ComputeStencilSweepFunctorOperationGroups)
{
    // define individual stages
    typedef StencilStage<IntermediateUsingDataInStage, IJRange<cComplete,-1,1,-1,1>, KRange<FullDomain,0,0> > Stage1;
    typedef StencilStage<AverageIntermediateToDataOutStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-1> > Stage2;
    
    // define individual caches
    typedef KCache<data_out, cFillAndFlush, KWindow<0,0>, KRange<FullDomain,1,-1> > DataOutCache;
    typedef KCache<intermediate, cFlush, KWindow<-1,1>, KRange<FullDomain,0,0> > IntermediateCache;
    typedef KCache<intermediate, cFlush, KWindow<-1,1>, KRange<FullDomain,0,0> > IntermediateCache;

    // define stage and cache vectors
    typedef boost::mpl::vector2<Stage1, Stage2> StencilStages;
    typedef boost::mpl::vector2<DataOutCache, IntermediateCache> Caches;

    // compute the full domain operation groups
    typedef compute_stencil_sweep_functor_operation_groups<
        StencilStages,
        Caches,
        FullDomain
    >::type FullDomainGroups;

    // define cache masks for intermediate and data_out
    // (using push back we get the right vector implementation in order to pass an is_same comparison)
    typedef boost::mpl::push_back<
        boost::mpl::vector0<>,
        boost::mpl::integral_c<int, intermediate>
    >::type IntermediateMask;
    typedef boost::mpl::push_back<
        boost::mpl::vector0<>,
        boost::mpl::integral_c<int, data_out>
    >::type DataOutMask;

    // assert three are 3 groups
    ASSERT_TRUE(boost::mpl::size<FullDomainGroups>::value == 3);

    // test the three groups individually
    ASSERT_TRUE(
        (
            check_stencil_sweep_functor_operation_group<
                StencilSweepFunctorOperationGroup<
                    boost::mpl::vector2<
                        StencilSweepFunctorOperationSlide<
                            IntermediateMask,
                            IJRange<cComplete,-1,1,-1,1>
                        >,
                        StencilSweepFunctorOperationStage<Stage1> 
                    >,
                    IJRange<cComplete,-1,1,-1,1>
                >,
                boost::mpl::at_c<FullDomainGroups, 0>::type
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_stencil_sweep_functor_operation_group<
                StencilSweepFunctorOperationGroup<
                    boost::mpl::vector3<
                        StencilSweepFunctorOperationFill<
                            DataOutMask,
                            IJRange<cComplete,0,0,0,0>
                        >,
                        StencilSweepFunctorOperationStage<Stage2>,
                        StencilSweepFunctorOperationFlush<
                            DataOutMask,
                            IJRange<cComplete,0,0,0,0>
                        >
                    >,
                    IJRange<cComplete,0,0,0,0>
                >,
                boost::mpl::at_c<FullDomainGroups, 1>::type
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_stencil_sweep_functor_operation_group<
                StencilSweepFunctorOperationGroup<
                    boost::mpl::vector1<
                        StencilSweepFunctorOperationFlush<
                            IntermediateMask,
                            IJRange<cComplete,-1,1,-1,1>
                        >
                    >,
                    IJRange<cComplete,-1,1,-1,1>
                >,
                boost::mpl::at_c<FullDomainGroups, 2>::type
            >::value
        )
    );

    // compute the k maximum operation groups
    typedef compute_stencil_sweep_functor_operation_groups<
        StencilStages,
        Caches,
        KMaximumCenter
    >::type KMaximumGroups;

    // assert three are 3 groups
    ASSERT_TRUE(boost::mpl::size<KMaximumGroups>::value == 3);

    // test the three groups individually
    // fill the data buffer in a separate group as it has a smaller ij range
    ASSERT_TRUE(
        (
            check_stencil_sweep_functor_operation_group<
                StencilSweepFunctorOperationGroup<
                    boost::mpl::vector1<
                        StencilSweepFunctorOperationFill<
                            DataOutMask,
                            IJRange<cComplete,0,0,0,0>
                        >
                    >,
                    IJRange<cComplete,0,0,0,0>
                >,
                boost::mpl::at_c<KMaximumGroups, 0>::type
            >::value
        )
    );
    ASSERT_TRUE(
        (
            check_stencil_sweep_functor_operation_group<
                StencilSweepFunctorOperationGroup<
                    boost::mpl::vector3<
                        StencilSweepFunctorOperationSlide<
                            IntermediateMask,
                            IJRange<cComplete,-1,1,-1,1>
                        >,
                        StencilSweepFunctorOperationStage<Stage1>,
                        StencilSweepFunctorOperationFlush<
                            IntermediateMask,
                            IJRange<cComplete,-1,1,-1,1>
                        >
                    >,
                    IJRange<cComplete,-1,1,-1,1>
                >,
                boost::mpl::at_c<KMaximumGroups, 1>::type
            >::value
        )
    );
    // flush data again in a separate group as it has a different ij range
    ASSERT_TRUE(
        (
            check_stencil_sweep_functor_operation_group<
                StencilSweepFunctorOperationGroup<
                    boost::mpl::vector1<
                        StencilSweepFunctorOperationFlush<
                            DataOutMask,
                            IJRange<cComplete,0,0,0,0>
                        >
                    >,
                    IJRange<cComplete,0,0,0,0>
                >,
                boost::mpl::at_c<KMaximumGroups, 2>::type
            >::value
        )
    );
}


