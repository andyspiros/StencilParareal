#include <boost/mpl/at.hpp>
#include <boost/mpl/integral_c.hpp>
#include "gtest/gtest.h"
#include "StencilFramework.h"

// define parameters
enum { Full, Terrain, Flat, Result };

class StencilSimpleLoopUnittest : public ::testing::Test
{
protected:  
  
    // input fields
    IJKRealField Full_;
    IJKRealField Terrain_;
    IJKRealField Flat_;
    // field setup by the test stencils
    IJKRealField Result_; 
    IJKSize size_;

    virtual void SetUp()
    {
        size_.Init(8, 12, 40);

        KBoundary boundary;
        boundary.Init(/*-3,3,-2,2,*/0,0); 

        Full_.Init("Full", size_, boundary);
        Terrain_.Init("Terrain", size_, boundary);
        Flat_.Init("Flat", size_, boundary);
        Result_.Init("Result", size_, boundary);
        
        // init everything to 0
        for(int i=-3; i<size_.iSize()+3; ++i)
        {
            for(int j=-2; j<size_.jSize()+2; ++j)
            {
                for(int k=0; k<size_.kSize(); ++k)
                {
                    Full_(i,j,k) = (Real)0.0;
                    Terrain_(i,j,k) = (Real)0.0;
                    Flat_(i,j,k) = (Real)0.0;
                }
            }
        }

        for(int i=-cNumBoundaryLines; i<size_.iSize()+cNumBoundaryLines; ++i)
        {
            for(int j=-cNumBoundaryLines; j<size_.jSize()+cNumBoundaryLines; ++j)
            {
                for(int k=0; k<size_.kSize(); ++k)
                {
                    Full_(i,j,k) = (Real)20.0 + (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX));
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

template<typename TEnv>
struct FullDomainUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)
    STAGE_PARAMETER(FullDomain, Full)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[Result::Center()] = ctx[Full::Center()];
    }
};

template<typename TEnv>
struct AtmosphereUpdate
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, Result)
    STAGE_PARAMETER(TerrainCoordinates, Terrain)
    STAGE_PARAMETER(FlatCoordinates, Flat)

    __ACC__
    static void Do(Context ctx, FlatCoordinates)
    {
        ctx[Result::Center()] = ctx[Flat::Center()];
    }

    __ACC__
    static void Do(Context ctx, TerrainCoordinates)
    {
        ctx[Result::Center()] = ctx[Terrain::Center()];
    }
};

// test full domain stage in k inc direction
TEST_F(StencilSimpleLoopUnittest, FullDomainKIncrement)
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
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);
}

// test full domain stage in k dec direction
TEST_F(StencilSimpleLoopUnittest, FullDomainKDecrement)
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
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);
}

// test full domain stage with adapted k range in k inc direction
TEST_F(StencilSimpleLoopUnittest, FullDomainCustomRangeKIncrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // update flat coordinates
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,1,-1> >()
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
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    // update terrain coordinates
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,-1> >()
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
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test full domain stage with adapted k range in k dec direction
TEST_F(StencilSimpleLoopUnittest, FullDomainCustomRangeKDecrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // update flat coordinates
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,1,-1> >()
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
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    // update terrain coordinates
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,-1> >()
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
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, size_.kSize()-1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test atmosphere stage in k inc direction
TEST_F(StencilSimpleLoopUnittest, AtmosphereKIncrement)
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
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, Terrain_, Result_);
}

// test atmosphere stage in k dec direction
TEST_F(StencilSimpleLoopUnittest, AtmosphereKDecrement)
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
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, Terrain_, Result_);
}

// test atmosphere stage with custom k range in k inc direction
TEST_F(StencilSimpleLoopUnittest, AtmosphereCustomKIncrement)
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
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-2> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, Terrain_, Result_);

    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test atmosphere stage with custom k range in k dec direction
TEST_F(StencilSimpleLoopUnittest, AtmosphereCustomKDecrement)
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
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,1,-2> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, Terrain_, Result_);

    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test atmosphere stage with custom k range in k inc direction
// execute it in two steps for terrain and flat domain
TEST_F(StencilSimpleLoopUnittest, Atmosphere2StepsKIncrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // update flat domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,1,0> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    // update terrain domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,-2> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, Terrain_, Result_);

    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test atmosphere stage with custom k range in k dec direction
// execute it in two steps for terrain and flat domain
TEST_F(StencilSimpleLoopUnittest, Atmosphere2StepsKDecrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // update flat domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<FlatCoordinates,1,0> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    // update terrain domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Terrain, cIn>(Terrain_),
            Param<Flat, cIn>(Flat_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<AtmosphereUpdate, IJRange<cComplete,0,0,0,0>, KRange<TerrainCoordinates,0,-2> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, 1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit);
    VerifyRange(from, to, Flat_, Result_);

    from.Init(0, 0, cFlatLimit);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-2);
    VerifyRange(from, to, Terrain_, Result_);

    from.Init(0, 0, size_.kSize()-2);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
}

// test working on two dimensional ranges
TEST_F(StencilSimpleLoopUnittest, TwoDimensionalRangesKIncrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // update flat domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumCenter> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus1> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus2> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainCenter> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumMinus2> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumMinus1> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 3);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, 3);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit+1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, cFlatLimit+1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-3);
    VerifyRangeNotModified(from, to, Result_);
    
    from.Init(0, 0, size_.kSize()-3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);
}

// test working on two dimensional ranges
TEST_F(StencilSimpleLoopUnittest, TwoDimensionalRangesKDecrement)
{
    Stencil testStencil;
    IJKIndex from, to;

    // update flat domain
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),  
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumCenter> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus1> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumPlus2> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumFlatCenter> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMinimumTerrainCenter> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumMinus2> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumMinus1> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,0,0,0>, KRange<KMaximumCenter> >()
                )
            )           
        )
    );
    testStencil.Apply();

    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), 3);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, 3);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit-1);
    VerifyRangeNotModified(from, to, Result_);

    from.Init(0, 0, cFlatLimit-1);
    to.Init(size_.iSize(), size_.jSize(), cFlatLimit+1);
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, 0, cFlatLimit+1);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize()-3);
    VerifyRangeNotModified(from, to, Result_);
    
    from.Init(0, 0, size_.kSize()-3);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);
}

// test different ij ranges
TEST_F(StencilSimpleLoopUnittest, DifferentIJRangesKIncrement)
{
    IJKIndex from, to;

    // test no corners
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(        
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cIndented,-1,1,-1,1>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();
    
    // ceck the calculation domain
    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);
    
    // check the corners are not initialized
    from.Init(-1, -1, 0);
    to.Init(0, 0, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
    
    from.Init(size_.iSize(), -1, 0);
    to.Init(size_.iSize()+1, 0, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
    
    from.Init(-1, size_.jSize(), 0);
    to.Init(0, size_.jSize()+1, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    from.Init(size_.iSize(), size_.jSize(), 0);
    to.Init(size_.iSize()+1, size_.jSize()+1, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    // check the boundaries without corners are set to Full
    from.Init(-1, 0, 0);
    to.Init(0, size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, -1, 0);
    to.Init(size_.iSize(), 0, size_.kSize());
    VerifyRange(from, to, Full_, Result_);

    from.Init(size_.iSize(), 0, 0);
    to.Init(size_.iSize()+1, size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, size_.jSize(), 0);
    to.Init(size_.iSize(), size_.jSize()+1, size_.kSize());
    VerifyRange(from, to, Full_, Result_);

    // check plus boundary
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,2,0,2>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // check corners not initialized
    from.Init(-1, -1, 0);
    to.Init(0, 0, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    from.Init(size_.iSize(), -1, 0);
    to.Init(size_.iSize()+1, 0, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    from.Init(-1, size_.jSize(), 0);
    to.Init(0, size_.jSize()+1, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
    
    // except for this one
    from.Init(size_.iSize(), size_.jSize(), 1);
    to.Init(size_.iSize()+2, size_.jSize()+2, size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    // fill the remaining boundaries
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,-2,0,-2,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,-2,0,0,2>, KRange<FullDomain,0,0> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,2,-2,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // check the full region is initialized
    from.Init(-2, -2, 0);
    to.Init(size_.iSize()+2, size_.jSize()+2, size_.kSize());
    VerifyRange(from, to, Full_, Result_);
}

// test different ij ranges
TEST_F(StencilSimpleLoopUnittest, DifferentIJRangesKDecrement)
{
    IJKIndex from, to;

    // test no corners
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(        
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cIndented,-1,1,-1,1>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();
    
    // ceck the calculation domain
    from.Init(0, 0, 0);
    to.Init(size_.iSize(), size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);
    
    // check the corners are not initialized
    from.Init(-1, -1, 0);
    to.Init(0, 0, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
    
    from.Init(size_.iSize(), -1, 0);
    to.Init(size_.iSize()+1, 0, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
    
    from.Init(-1, size_.jSize(), 0);
    to.Init(0, size_.jSize()+1, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    from.Init(size_.iSize(), size_.jSize(), 0);
    to.Init(size_.iSize()+1, size_.jSize()+1, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    // check the boundaries without corners are set to Full
    from.Init(-1, 0, 0);
    to.Init(0, size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, -1, 0);
    to.Init(size_.iSize(), 0, size_.kSize());
    VerifyRange(from, to, Full_, Result_);

    from.Init(size_.iSize(), 0, 0);
    to.Init(size_.iSize()+1, size_.jSize(), size_.kSize());
    VerifyRange(from, to, Full_, Result_);

    from.Init(0, size_.jSize(), 0);
    to.Init(size_.iSize(), size_.jSize()+1, size_.kSize());
    VerifyRange(from, to, Full_, Result_);

    // check plus boundary
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,2,0,2>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // check corners not initialized
    from.Init(-1, -1, 0);
    to.Init(0, 0, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    from.Init(size_.iSize(), -1, 0);
    to.Init(size_.iSize()+1, 0, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);

    from.Init(-1, size_.jSize(), 0);
    to.Init(0, size_.jSize()+1, size_.kSize());
    VerifyRangeNotModified(from, to, Result_);
    
    // except for this one
    from.Init(size_.iSize(), size_.jSize(), 1);
    to.Init(size_.iSize()+2, size_.jSize()+2, size_.kSize()-1);
    VerifyRange(from, to, Full_, Result_);

    // fill the remaining boundaries
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size_,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param<Result, cInOut>(Result_),
            Param<Full, cIn>(Full_)
        ),
        define_loops(
            define_sweep<cKDecrement>(
                define_stages(
                    StencilStage<FullDomainUpdate, IJRange<cComplete,-2,0,-2,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,-2,0,0,2>, KRange<FullDomain,0,0> >(),
                    StencilStage<FullDomainUpdate, IJRange<cComplete,0,2,-2,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    testStencil.Apply();

    // check the full region is initialized
    from.Init(-2, -2, 0);
    to.Init(size_.iSize()+2, size_.jSize()+2, size_.kSize());
    VerifyRange(from, to, Full_, Result_);
}