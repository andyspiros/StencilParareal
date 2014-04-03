#include <boost/mpl/at.hpp>
#include <boost/mpl/integral_c.hpp>
#include "gtest/gtest.h"
#include "StencilFramework.h"

// define parameters
enum { KMin, KMax, Full, Terrain, Flat, Result, UpSum, DownSum, Last, Buffer1, Buffer2 };

class StencilComplexLoopUnittest : public ::testing::Test
{
protected:  
  
    // input fields
    IJKRealField KMin_;
    IJKRealField KMax_;
    IJKRealField Full_;
    IJKRealField Terrain_;
    IJKRealField TerrainMin_;
    IJKRealField Flat_;
    IJKRealField Buffer1_;
    IJKRealField Buffer2_;
    // field setup by the test stencils
    IJKRealField Result_; 
    // 2d field for loop direction test
    IJRealField UpSum_;
    IJRealField DownSum_;
    IJRealField Last_;
    IJKSize size_;

    virtual void SetUp()
    {
        size_.Init(8, 12, 40);

        KBoundary boundary;
        boundary.Init(-1,1); 

        KMin_.Init("KMin", size_, boundary);
        KMax_.Init("KMax", size_, boundary);
        Full_.Init("Full", size_, boundary);
        Terrain_.Init("Terrain", size_, boundary);
        Flat_.Init("Flat", size_, boundary);
        Buffer1_.Init("Buffer1", size_, boundary);
        Buffer2_.Init("Buffer2", size_, boundary);
        Result_.Init("Result", size_, boundary);
        UpSum_.Init("UpSum", size_, boundary);
        DownSum_.Init("DownSum", size_, boundary);
        Last_.Init("Last", size_, boundary);
        
        for(int i=-cNumBoundaryLines; i<size_.iSize()+cNumBoundaryLines; ++i)
        {
            for(int j=-cNumBoundaryLines; j<size_.jSize()+cNumBoundaryLines; ++j)
            {
                for(int k=0; k<size_.kSize(); ++k)
                {
                    KMin_(i,j,k) = (Real)10.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                    Full_(i,j,k) = (Real)20.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                    KMax_(i,j,k) = (Real)30.0 +  (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                    Terrain_(i,j,k) = (Real)40.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                    Flat_(i,j,k) = (Real)50.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
                }
            }
        }

        for(int i=-cNumBoundaryLines; i<size_.iSize()+cNumBoundaryLines; ++i)
        {
            for(int j=-cNumBoundaryLines; j<size_.jSize()+cNumBoundaryLines; ++j)
            {
                for(int k=0; k<size_.kSize(); ++k)
                {
                    Result_(i,j,k) = (Real)-1.0;
                }
            }
        }
    }
     
    void VerifyRange(const IJKIndex& from, const IJKIndex& to, IJKRealField& expected, IJKRealField& actual)
    {
        for(int i=from.iIndex(); i<to.iIndex(); ++i)
        {
            for(int j=from.jIndex(); j<to.jIndex(); ++j)
            {
                for(int k=from.kIndex(); k<to.kIndex(); ++k)
                {
                    ASSERT_DOUBLE_EQ(expected(i,j,k), actual(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
                }
            }
        }
    }
    
    void VerifyRange(const IJKIndex& from, const IJKIndex& to, Real expected, IJKRealField& actual)
    {
        for(int i=from.iIndex(); i<to.iIndex(); ++i)
        {
            for(int j=from.jIndex(); j<to.jIndex(); ++j)
            {
                for(int k=from.kIndex(); k<to.kIndex(); ++k)
                {
                    ASSERT_DOUBLE_EQ(expected, actual(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k;
                }
            }
        }
    }
    
    void VerifyRangeNotModified(const IJKIndex& from, const IJKIndex& to, IJKRealField& actual)
    {
        for(int i=from.iIndex(); i<to.iIndex(); ++i)
        {
            for(int j=from.jIndex(); j<to.jIndex(); ++j)
            {
                for(int k=from.kIndex(); k<to.kIndex(); ++k)
                {
                    ASSERT_DOUBLE_EQ(-1.0, actual(i,j,k))  << "i=" << i << "; j=" << j << "; k=" << k;
                }
            }
        }
    }
};

// stencil stage updating the full domain
// including boundary specific updates
template<typename TEnv>
struct FullDomainUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)
    STAGE_PARAMETER(FullDomain, KMin)
    STAGE_PARAMETER(FullDomain, KMax)
    STAGE_PARAMETER(FullDomain, Full)

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        ctx[Result::Center()] = ctx[KMax::Center()];
    }
    
    __ACC__
    static void Do(Context ctx, KMaximumMinus2)
    {
        ctx[Result::Center()] = ctx[KMax::Center()];
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[Result::Center()] = ctx[Full::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMinimumPlus2)
    {
        ctx[Result::Center()] = ctx[KMin::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        ctx[Result::Center()] = ctx[KMin::Center()];
    }
};

// test a full domain update with custom boundaries in k increment direction
TEST_F(StencilComplexLoopUnittest, FullDomainKIncrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // test full domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,-1> >(), // update last level using a 2d stage
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >()
                )
            )           
        )
    );
    testStencil.Apply();
    
    // test k minimum levels
    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, KMin_, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), 2);
    VerifyRange(from, to, Full_, Result_);
    
    from.Init(0, 0, 2);
    to.Init(size_.iSize(), size_.jSize(), 3);
    VerifyRange(from, to, KMin_, Result_);

    // test full domain
    from.Init(0, 0, 3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-3);
    VerifyRange(from, to, Full_, Result_);
    
    // test k maximum levels
    from.Init(0, 0, size_.kSize()-3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, KMax_, Result_);

    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// test a full domain update with custom boundaries in k decrement direction
TEST_F(StencilComplexLoopUnittest, FullDomainKDecrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // test full domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,-1> >(), // update last level using a 2d stage
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >()
                )
            )           
        )
    );
    testStencil.Apply();
    
    // test k minimum levels
    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, KMin_, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), 2);
    VerifyRange(from, to, Full_, Result_);
    
    from.Init(0, 0, 2);
    to.Init(size_.iSize(), size_.jSize(), 3);
    VerifyRange(from, to, KMin_, Result_);

    // test full domain
    from.Init(0, 0, 3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-3);
    VerifyRange(from, to, Full_, Result_);
    
    // test k maximum levels
    from.Init(0, 0, size_.kSize()-3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, KMax_, Result_);

    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// stencil stage performing an atmosphere update 
// providing custom update functionality at the boundaries
template<typename TEnv>
struct AtmosphereUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)
    STAGE_PARAMETER(FullDomain, KMin)
    STAGE_PARAMETER(FullDomain, KMax)
    STAGE_PARAMETER(FullDomain, Terrain)
    STAGE_PARAMETER(FullDomain, Flat)

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        ctx[Result::Center()] = ctx[KMax::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMaximumMinus1)
    {
        ctx[Result::Center()] = ctx[KMax::Center()];
    }

    __ACC__
    static void Do(Context ctx, TerrainCoordinates)
    {
        ctx[Result::Center()] = ctx[Terrain::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMinimumTerrainCenter)
    {
        ctx[Result::Center()] = ctx[KMin::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMaximumFlatCenter)
    {
        ctx[Result::Center()] = ctx[KMax::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMaximumFlatMinus2)
    {
        ctx[Result::Center()] = ctx[KMax::Center()];
    }

    __ACC__
    static void Do(Context ctx, FlatCoordinates)
    {
        ctx[Result::Center()] = ctx[Flat::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMinimumPlus1)
    {
        ctx[Result::Center()] = ctx[KMin::Center()];
    }
};

// test an atmosphere update with custom boundaries in k increment direction
TEST_F(StencilComplexLoopUnittest, AtmosphereKIncrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // test full domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,-1> >()
                )
            )           
        )
    );
    testStencil.Apply();
    
    // test k minimum levels
    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), 2);
    VerifyRange(from, to, KMin_, Result_);
    
    // test flat
    from.Init(0, 0, 2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-3);
    VerifyRange(from, to, Flat_, Result_);

    // test flat maximum
    from.Init(0, 0, cFlatLimit-3);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-2);
    VerifyRange(from, to, KMax_, Result_);

    from.Init(0, 0, cFlatLimit-2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-0);
    VerifyRange(from, to, KMax_, Result_);

    // test terrain minimum
    from.Init(0, 0, cFlatLimit-0);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit+1);
    VerifyRange(from, to, KMin_, Result_);

    // test terrain
    from.Init(0, 0, cFlatLimit+1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, Terrain_, Result_);
    
    // test k maximum level
    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, KMax_, Result_);

    // make sure the top level was not touched
    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test an atmosphere update with custom boundaries in k decrement direction
TEST_F(StencilComplexLoopUnittest, AtmosphereKDecrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // test full domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,-1> >()
                )
            )           
        )
    );
    testStencil.Apply();
    
    // test k minimum levels
    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), 2);
    VerifyRange(from, to, KMin_, Result_);
    
    // test flat
    from.Init(0, 0, 2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-3);
    VerifyRange(from, to, Flat_, Result_);

    // test flat maximum
    from.Init(0, 0, cFlatLimit-3);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-2);
    VerifyRange(from, to, KMax_, Result_);

    from.Init(0, 0, cFlatLimit-2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-0);
    VerifyRange(from, to, KMax_, Result_);

    // test terrain minimum
    from.Init(0, 0, cFlatLimit-0);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit+1);
    VerifyRange(from, to, KMin_, Result_);

    // test terrain
    from.Init(0, 0, cFlatLimit+1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, Terrain_, Result_);
    
    // test k maximum level
    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, KMax_, Result_);

    // make sure the top level was not touched
    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test an atmosphere update in 2 steps with custom boundaries in k increment direction
TEST_F(StencilComplexLoopUnittest, Atmosphere2StepsKIncrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // test full domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,2,-1> >(),
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainCenter> >(),
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >()
                )
            )           
        )
    );
    testStencil.Apply();
    
    // make sure there are no minimum levels
    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 2);
    VerifyRangeNotModified(from, to, Result_);
    
    // test flat
    from.Init(0, 0, 2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-3);
    VerifyRange(from, to, Flat_, Result_);

    // test flat maximum
    from.Init(0, 0, cFlatLimit-3);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-2);
    VerifyRange(from, to, KMax_, Result_);

    from.Init(0, 0, cFlatLimit-2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-0);
    VerifyRangeNotModified(from, to, Result_);

    // test terrain minimum
    from.Init(0, 0, cFlatLimit-0);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit+1);
    VerifyRange(from, to, KMin_, Result_);

    // test terrain
    from.Init(0, 0, cFlatLimit+1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRangeNotModified(from, to, Result_);

    // test k maximum level
    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// test an atmosphere update in 2 steps with custom boundaries in k increment direction
TEST_F(StencilComplexLoopUnittest, Atmosphere2StepsKDecrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // test full domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,2,-1> >(),
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainCenter> >(),
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >()
                )
            )           
        )
    );
    testStencil.Apply();
    
    // make sure there are no minimum levels
    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 2);
    VerifyRangeNotModified(from, to, Result_);
    
    // test flat
    from.Init(0, 0, 2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-3);
    VerifyRange(from, to, Flat_, Result_);

    // test flat maximum
    from.Init(0, 0, cFlatLimit-3);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-2);
    VerifyRange(from, to, KMax_, Result_);

    from.Init(0, 0, cFlatLimit-2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-0);
    VerifyRangeNotModified(from, to, Result_);

    // test terrain minimum
    from.Init(0, 0, cFlatLimit-0);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit+1);
    VerifyRange(from, to, KMin_, Result_);

    // test terrain
    from.Init(0, 0, cFlatLimit+1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRangeNotModified(from, to, Result_);

    // test k maximum level
    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// stencil stage updating boundaries outside calculation domain
template<typename TEnv>
struct NegativeBoundaryUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)
    STAGE_PARAMETER(FullDomain, KMin)
    STAGE_PARAMETER(FullDomain, KMax)

    __ACC__
    static void Do(Context ctx, KMaximumPlus2)
    {
        ctx[Result::At(kminus2)] = ctx[KMax::At(kminus2)];
    }

    __ACC__
    static void Do(Context ctx, KMaximumPlus1)
    {
        ctx[Result::At(kminus2)] = ctx[KMax::At(kminus2)];
    }

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        ctx[Result::At(kminus2)] = ctx[KMax::At(kminus2)];
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        // do nothing
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        ctx[Result::At(kplus2)] = ctx[KMin::At(kplus2)];
    }

    __ACC__
    static void Do(Context ctx, KMinimumMinus1)
    {
        ctx[Result::At(kplus2)] = ctx[KMin::At(kplus2)];
    }

    __ACC__
    static void Do(Context ctx, KMinimumMinus2)
    {
        ctx[Result::At(kplus2)] = ctx[KMin::At(kplus2)];
    }
};

// test negative boundaries in k increment direction
TEST_F(StencilComplexLoopUnittest, NegativeBoundariesKIncrement)
{
    IJKIndex from, to;

    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<NegativeBoundaryUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,-2,0> >(),
                    StencilStage<NegativeBoundaryUpdate, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,0> >(),
                    StencilStage<NegativeBoundaryUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumPlus2> >()
                )
            ) 
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 3);
    VerifyRange(from, to, KMin_, Result_);

    from.Init(0, 0, 3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-3);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, size_.kSize()-3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, KMax_, Result_);
    
    // one k maximum level is not updated
    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// test negative boundaries in k decrement direction
TEST_F(StencilComplexLoopUnittest, NegativeBoundariesKDecrement)
{
    IJKIndex from, to;

    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<NegativeBoundaryUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,-2,0> >(),
                    StencilStage<NegativeBoundaryUpdate, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,0> >(),
                    StencilStage<NegativeBoundaryUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumPlus2> >()
                )
            ) 
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 3);
    VerifyRange(from, to, KMin_, Result_);

    from.Init(0, 0, 3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-3);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, size_.kSize()-3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, KMax_, Result_);
    
    // one k maximum level is not updated
    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// stencil stage summing up in loop direction
template<typename TEnv>
struct FullDomainSumUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)
    STAGE_PARAMETER(FullDomain, KMin)
    STAGE_PARAMETER(FullDomain, KMax)
    STAGE_PARAMETER(FullDomain, Full)
    STAGE_PARAMETER(FullDomain, UpSum)
    STAGE_PARAMETER(FullDomain, DownSum)
    STAGE_PARAMETER(FullDomain, Last)

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        ctx[Result::Center()] = ctx[KMax::Center()];

        if(ctx[Result::At(kminus1)] == ctx[Last::Center()])
        {
            ctx[UpSum::Center()] += (Real)1.0; 
        }
        ctx[Last::Center()] = ctx[Result::Center()];
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[Result::Center()] = ctx[Full::Center()];

        if(ctx[Result::At(kminus1)] == ctx[Last::Center()])
        {
            ctx[UpSum::Center()] += (Real)1.0; 
        }
        if(ctx[Result::At(kplus1)] == ctx[Last::Center()])
        {
            ctx[DownSum::Center()] += (Real)1.0;
        } 
        ctx[Last::Center()] = ctx[Result::Center()];
    }

    __ACC__
    static void Do(Context ctx, KMinimumCenter)
    {
        ctx[Result::Center()] = ctx[KMin::Center()];
        
        if(ctx[Result::At(kplus1)] == ctx[Last::Center()])
        {
            ctx[DownSum::Center()] += (Real)1.0;
        }        
        ctx[Last::Center()] = ctx[Result::Center()];
    }
};

// test verifying 2d reduction works in forward direction 
TEST_F(StencilComplexLoopUnittest, ReductionKIncrement)
{
    IJKIndex from, to;    
    for(int i=0; i<size_.iSize(); ++i)
    {
        for(int j=0; j<size_.jSize(); ++j)
        {
            UpSum_(i,j,0) = 0.0;
            DownSum_(i,j,0) = 0.0;
            Last_(i,j,0) = 0.0;
        }
    }

    // test full domain using 2 loops in order to make sure loops are ordered correctly
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Full, cIn>(Full_),
            Param<UpSum, cInOut>(UpSum_),
            Param<DownSum, cInOut>(DownSum_),  
            Param<Last, cInOut>(Last_)  
        ),
        define_loops(
        define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainSumUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,0,0> >(),
                    StencilStage<FullDomainSumUpdate, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,0> >()
                )
            ) 
        )
    );
    testStencil.Apply();

    // check the sum was done right
    for(int i=0; i<size_.iSize(); ++i)
    {
        for(int j=0; j<size_.jSize(); ++j)
        {
            ASSERT_DOUBLE_EQ((double)(size_.kSize()-1), UpSum_(i,j,0)) << "i=" << i << "; j=" << j;
        }
    }

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, KMin_, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// test verifying 2d reduction works in backward direction 
TEST_F(StencilComplexLoopUnittest, ReductionKDecrement)
{
    IJKIndex from, to;    
    for(int i=0; i<size_.iSize(); ++i)
    {
        for(int j=0; j<size_.jSize(); ++j)
        {
            UpSum_(i,j,0) = 0.0;
            DownSum_(i,j,0) = 0.0;
            Last_(i,j,0) = 0.0;
        }
    }

    // test full domain using 2 loops in order to make sure loops are ordered correctly
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Full, cIn>(Full_),
            Param<UpSum, cInOut>(UpSum_),
            Param<DownSum, cInOut>(DownSum_),  
            Param<Last, cInOut>(Last_)  
        ),
        define_loops(
        define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainSumUpdate, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,0> >(),    
                    StencilStage<FullDomainSumUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,0,0> >()
                )
            ) 
        )
    );
    testStencil.Apply();

    // check the sum was done right
    for(int i=0; i<size_.iSize(); ++i)
    {
        for(int j=0; j<size_.jSize(); ++j)
        {
            ASSERT_DOUBLE_EQ((double)(size_.kSize()-1), DownSum_(i,j,0)) << "i=" << i << "; j=" << j;
        }
    }

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, KMin_, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// stencil stage reseting last and result value
template<typename TEnv>
struct ResetResults
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)
    STAGE_PARAMETER(KMaximumCenter, Last)

    __ACC__
    static void Do(Context ctx, KMaximumCenter)
    {
        ctx[Last::Center()] = (Real)0.0;
        ctx[Result::Center()] = (Real)-1.0;
    }

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[Result::Center()] = (Real)-1.0;
    }
};

// test verifying 2d reduction works in both directions at once
TEST_F(StencilComplexLoopUnittest, Reduction2Way)
{
    IJKIndex from, to;    
    for(int i=0; i<size_.iSize(); ++i)
    {
        for(int j=0; j<size_.jSize(); ++j)
        {
            UpSum_(i,j,0) = 0.0;
            DownSum_(i,j,0) = 0.0;
            Last_(i,j,0) = 0.0;
        }
    }

    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<KMin, cIn>(KMin_),
            Param<KMax, cIn>(KMax_),
            Param<Full, cIn>(Full_),
            Param<UpSum, cInOut>(UpSum_),
            Param<DownSum, cInOut>(DownSum_),  
            Param<Last, cInOut>(Last_)  
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainSumUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            ),
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<ResetResults, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            ),
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainSumUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )            
        )
    );
    testStencil.Apply();

    // check the sum was done right
    for(int i=0; i<size_.iSize(); ++i)
    {
        for(int j=0; j<size_.jSize(); ++j)
        {
            ASSERT_DOUBLE_EQ((double)(size_.kSize()-1), UpSum_(i,j,0)) << "i=" << i << "; j=" << j;
            ASSERT_DOUBLE_EQ((double)(size_.kSize()-1), DownSum_(i,j,0)) << "i=" << i << "; j=" << j;
        }
    }

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, KMin_, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, KMax_, Result_);
}

// stencil stage pre-multiplying full by 2
template<typename TEnv>
struct PreMultiplyFull
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Buffer1)
    STAGE_PARAMETER(FullDomain, Full)
    STAGE_PARAMETER(FullDomain, KMin)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[Buffer1::Center()] = ctx[KMin::Center()] * ctx[Full::Center()];
    }
};

// stencil stage computing the average in i direction
template<typename TEnv>
struct AverageInI
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Buffer1)
    STAGE_PARAMETER(FullDomain, Buffer2)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[Buffer2::Center()] = 
            (T)1.0 / (T)3.0 * (
                ctx[Buffer1::At(iminus1)] +
                ctx[Buffer1::Center()] +
                ctx[Buffer1::At(iplus1)]
            );
    }
};

// stencil stage computing the average in j direction
template<typename TEnv>
struct AverageInJ
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)
    STAGE_PARAMETER(FullDomain, Buffer2)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[Result::Center()] = 
            (T)1.0 / (T)3.0 * (
                ctx[Buffer2::At(jminus1)] +
                ctx[Buffer2::Center()] +
                ctx[Buffer2::At(jplus1)]
            );
    }
};

// test 2 stage averaging in k increment direction
TEST_F(StencilComplexLoopUnittest, AverageKIncrement)
{
    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Buffer1, cInOut>(Buffer1_),
            Param<Buffer2, cInOut>(Buffer2_),
            Param<Full, cIn>(Full_),  
            Param<KMin, cIn>(KMin_)  
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<PreMultiplyFull, IJRange<cComplete,-1,1,-1,1>, KRange<FullDomain,0,0> >(),
                    StencilStage<AverageInI, IJRange<cComplete,0,0,-1,1>, KRange<FullDomain,0,0> >(),
                    StencilStage<AverageInJ, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )            
        )
    );
    testStencil.Apply();

    // check the sum was done right
    for(int i=0; i<size_.iSize(); ++i)
    {
        for(int j=0; j<size_.jSize(); ++j)
        {
            for(int k=0; k<size_.kSize(); ++k)
            {
                // compute average in i j
                Real average = (Real)1.0 / (Real)9.0 * 
                    (
                        KMin_(i-1,j-1,k) * Full_(i-1,j-1,k) + KMin_(i,j-1,k) * Full_(i,j-1,k) + KMin_(i+1,j-1,k) * Full_(i+1,j-1,k) +
                        KMin_(i-1,j,k) * Full_(i-1,j,k) + KMin_(i,j,k) * Full_(i,j,k) + KMin_(i+1,j,k) * Full_(i+1,j,k) +
                        KMin_(i-1,j+1,k) * Full_(i-1,j+1,k) + KMin_(i,j+1,k) * Full_(i,j+1,k) + KMin_(i+1,j+1,k) *Full_(i+1,j+1,k) 
                    );
                
                ASSERT_DOUBLE_EQ(average, Result_(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k; 
            }
            
        }
    }
}

// test 2 stage averaging in k decrement direction
TEST_F(StencilComplexLoopUnittest, AverageKDecrement)
{
    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Buffer1, cInOut>(Buffer1_),
            Param<Buffer2, cInOut>(Buffer2_),
            Param<Full, cIn>(Full_),  
            Param<KMin, cIn>(KMin_)  
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<PreMultiplyFull, IJRange<cComplete,-1,1,-1,1>, KRange<FullDomain,0,0> >(),
                    StencilStage<AverageInI, IJRange<cComplete,0,0,-1,1>, KRange<FullDomain,0,0> >(),
                    StencilStage<AverageInJ, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )            
        )
    );
    testStencil.Apply();

    // check the sum was done right
    for(int i=0; i<size_.iSize(); ++i)
    {
        for(int j=0; j<size_.jSize(); ++j)
        {
            for(int k=0; k<size_.kSize(); ++k)
            {
                // compute average in i j
                Real average = (Real)1.0 / (Real)9.0 * 
                    (
                        KMin_(i-1,j-1,k) * Full_(i-1,j-1,k) + KMin_(i,j-1,k) * Full_(i,j-1,k) + KMin_(i+1,j-1,k) * Full_(i+1,j-1,k) +
                        KMin_(i-1,j,k) * Full_(i-1,j,k) + KMin_(i,j,k) * Full_(i,j,k) + KMin_(i+1,j,k) * Full_(i+1,j,k) +
                        KMin_(i-1,j+1,k) * Full_(i-1,j+1,k) + KMin_(i,j+1,k) * Full_(i,j+1,k) + KMin_(i+1,j+1,k) *Full_(i+1,j+1,k) 
                    );
                
                ASSERT_DOUBLE_EQ(average, Result_(i,j,k)) << "i=" << i << "; j=" << j << "; k=" << k; 
            }
            
        }
    }
}

// stencil stage adding two to all elements of the full domain
template<typename TEnv>
struct AddTwo
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        if(ctx[Result::Center()] == (T)-1.0)
        {
            ctx[Result::Center()] = (Real)2.0;
        }
        else
        {
            ctx[Result::Center()] += (Real)2.0;
        }
    }
};

// stencil stage subtracting two if the value is two
template<typename TEnv>
struct MinusTwo
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        if(ctx[Result::Center()] == (T)2.0)
        {
            ctx[Result::Center()] -= (Real)2.0;
        }
    }
};

// test merging flat domain with various 2d domains in k increment direction
TEST_F(StencilComplexLoopUnittest, MergeFlatAndTerrainBoundaryKIncrement)
{
    IJKIndex from, to;  

    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,2,0> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus1> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainCenter> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumMinus1> >(),
                    StencilStage<MinusTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumMinus1> >()
                )
            )            
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit+1);
    VerifyRange(from, to, 2.0, Result_);
    
    from.Init(0, 0, cFlatLimit+1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRangeNotModified(from, to, Result_);
    
    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, 0.0, Result_); // note two should be subtracted again

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test merging flat domain with various 2d domains in k decrement direction
TEST_F(StencilComplexLoopUnittest, MergeFlatAndTerrainBoundaryKDecrement)
{
    IJKIndex from, to;  

    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,2,0> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus1> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainCenter> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumMinus1> >(),
                    StencilStage<MinusTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumMinus1> >()
                )
            )            
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit+1);
    VerifyRange(from, to, 2.0, Result_);
    
    from.Init(0, 0, cFlatLimit+1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRangeNotModified(from, to, Result_);
    
    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, 0.0, Result_); // note two should be subtracted again

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test splitting full domain with 2d domain in k increment direction
TEST_F(StencilComplexLoopUnittest, SplitFullDomainKIncrement)
{
    IJKIndex from, to;  

    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,-2> >(),
                    
                    StencilStage<MinusTwo, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus1> >(),

                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >()

                    

                )
            )            
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, 2.0, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), 2);
    VerifyRange(from, to, 0.0, Result_);
    
    from.Init(0, 0, 2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRange(from, to, 2.0, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, 4.0, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, 2.0, Result_);
    
    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, 2.0, Result_);
}

// test splitting full domain with 2d domain in k decrement direction
TEST_F(StencilComplexLoopUnittest, SplitFullDomainKDecrement)
{
    IJKIndex from, to;  

    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,-2> >(),
                    StencilStage<MinusTwo, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus1> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >()
                )
            )            
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRange(from, to, 2.0, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), 2);
    VerifyRange(from, to, 0.0, Result_);
    
    from.Init(0, 0, 2);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRange(from, to, 2.0, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, 4.0, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, 2.0, Result_);
    
    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, 2.0, Result_);
}

// merge different sized 3d domains in k increment direction
TEST_F(StencilComplexLoopUnittest, FullDomainAndTerrainKIncrement)
{
    IJKIndex from, to;  

    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-2> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,0> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >()
                )
            )            
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRange(from, to, 2.0, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, 4.0, Result_);

    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, 2.0, Result_);
}

// merge different sized 3d domains in k decrement direction
TEST_F(StencilComplexLoopUnittest, FullDomainAndTerrainKDecrement)
{
    IJKIndex from, to;  

    // test up and down sweep
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-2> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,0> >(),
                    StencilStage<AddTwo, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >()
                )
            )            
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRange(from, to, 2.0, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, 4.0, Result_);

    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, 2.0, Result_);
}

