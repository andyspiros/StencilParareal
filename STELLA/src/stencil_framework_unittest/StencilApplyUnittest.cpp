#include "gtest/gtest.h"
#include "StencilFramework.h"

enum { par, buffer };

// define specific data field boundaries for size tests!
const int cIMinusBoundaryLines = -3;
const int cIPlusBoundaryLines = 2;
const int cJMinusBoundaryLines = -5;
const int cJPlusBoundaryLines = 4;
typedef DataFieldIJBoundary<cIMinusBoundaryLines, cIPlusBoundaryLines, cJMinusBoundaryLines, cJPlusBoundaryLines> ApplyBoundary;

#ifdef __CUDA_BACKEND__
typedef DataFieldCUDA<Real, DataFieldStorageFormat<ApplyBoundary, StorageOrder::KJI, CUDARealAlignment> > ApplyField;
#else
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<ApplyBoundary, StorageOrder::JIK, OpenMPAlignment> > ApplyField;
#endif

// apply test stencils

template<typename TEnv>
struct ParameterInit
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, par)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[par::Center()] = (T)1.0;
    }
};

template<typename TEnv>
struct BufferInit
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, buffer)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[buffer::Center()] = (T)1.0;
    }
};

template<typename TEnv>
struct ParameterInitFromBuffer
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, buffer)
    STAGE_PARAMETER(FullDomain, par)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[par::Center()] = (ctx[buffer::At(iplus1)] + ctx[buffer::At(iminus1)] + ctx[buffer::At(jplus1)] + ctx[buffer::At(jminus1)]) / (T)4.0;
    }
};

// method testing apply for different calculation domain and k boundary sizes
void apply_size_round(const int isize, const int jsize, const int ksize, const int kMinusBoundary, const int kPlusBoundary)
{
    IJKSize size;
    size.Init(isize, jsize, ksize);
    KBoundary boundary;
    boundary.Init(kMinusBoundary, kPlusBoundary);

    ApplyField par;
    par.Init("par", size, boundary);

    // init the field with 0.0
    for(int i = cIMinusBoundaryLines; i < isize + cIPlusBoundaryLines; ++i)
    {
        for(int j = cJMinusBoundaryLines; j < jsize + cJPlusBoundaryLines; ++j)
        {
            for(int k = kMinusBoundary; k < ksize + kPlusBoundary; ++k)
            {
                par(i,j,k) = 0.0;
            }
        }
    }
    
    // setup a test stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::par, cInOut>(par)
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

    // check the calculation domain is set to 1.0
    for(int i = cIMinusBoundaryLines; i < isize + cIPlusBoundaryLines; ++i)
    {
        for(int j = cJMinusBoundaryLines; j < jsize + cJPlusBoundaryLines; ++j)
        {
            for(int k = kMinusBoundary; k < ksize + kPlusBoundary; ++k)
            {
                if(i >= 0 && i < isize && j >= 0 && j < jsize && k >= 0 && k < ksize)
                {
                    ASSERT_DOUBLE_EQ(1.0, par(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
                }
                else
                {
                    ASSERT_DOUBLE_EQ(0.0, par(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
                }
            }
        }
    }
}

// test if apply works for different calculation domain and k boundary sizes
TEST(StencilApplyUnittest, Size)
{
    // test for different sizes
    for(int iSize = 1; iSize < 51; iSize += 14)
    {
        for(int jSize = 1; jSize < 33; jSize += 6)
        {
            for(int kSize = 2; kSize < 13; kSize += 5) // the minimum k loop size is 2
            {
                for(int kMinusBoundary = -1; kMinusBoundary <= 0; kMinusBoundary += 1)
                {
                    for(int kPlusBoundary = 1; kPlusBoundary >= 0; kPlusBoundary -= 1)
                    {
                        apply_size_round(iSize, jSize, kSize, kMinusBoundary, kPlusBoundary);
                    }
                }
            }
        }
    }
}

// method testing apply for different boundary sizes
void apply_boundary_round(Stencil& testStencil, ApplyField& par, const int isize, const int jsize, const int ksize, const int iStart, const int iEnd, const int jStart, const int jEnd)
{
    IJBoundary applyBoundary;
    applyBoundary.Init(iStart, iEnd-isize, jStart, jEnd-jsize);

    // init the field with 0.0
    for(int i = cIMinusBoundaryLines; i < isize + cIPlusBoundaryLines; ++i)
    {
        for(int j = cJMinusBoundaryLines; j < jsize + cJPlusBoundaryLines; ++j)
        {
            for(int k = 0; k < ksize; ++k)
            {
                par(i,j,k) = 0.0;
            }
        }
    }

    testStencil.Apply(applyBoundary);

    // check the calculation domain + apply boundary is set to 1.0
    for(int i = cIMinusBoundaryLines; i < isize + cIPlusBoundaryLines; ++i)
    {
        for(int j = cJMinusBoundaryLines; j < jsize + cJPlusBoundaryLines; ++j)
        {
            for(int k = 0; k < ksize; ++k)
            {
                if(i >= iStart && i < iEnd && j >= jStart && j < jEnd && k >= 0 && k < ksize)
                {
                    ASSERT_DOUBLE_EQ(1.0, par(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
                }
                else
                {
                    ASSERT_DOUBLE_EQ(0.0, par(i, j, k)) << "i=" << i << "; j=" << j << "; k=" << k;
                }
            }
        }
    }
}

// test if apply works for different boundary sizes
TEST(StencilApplyUnittest, Boundary)
{
    IJKSize size;
    size.Init(10, 11, 4);
    KBoundary boundary;
    boundary.Init(0, 0);

    ApplyField par;
    par.Init("par", size, boundary);

    // setup a test stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::par, cInOut>(par)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<ParameterInit, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )  
        )
    );

    // test outer boundary
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -1, 11, -2, 15);

    // test inner boundary
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 2, 7, 1, 9);

    // top left corner
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -2, 3, -1, 2);
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -2, -1, -2, -1); // corner only

    // top right corner
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 9, 12, -1, 2);
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 10, 12, -2, -1); // corner only

    // bottom left corner
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -2, 3, 10, 13);
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -2, -1, 12, 13); // corner only

    // bottom right corner
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 9, 12, 10, 13);
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 10, 12, 12, 13); // corner only
}

// test if apply works for different boundary sizes
TEST(StencilApplyUnittest, BoundaryWithStages)
{
    IJKSize size;
    size.Init(10, 11, 4);
    KBoundary boundary;
    boundary.Init(0, 0);

    ApplyField par;
    par.Init("par", size, boundary);

    // setup a test stencil
    Stencil testStencil;
    StencilCompiler::Build(
        testStencil,
        "TestStencil",
        size,
        StencilConfiguration<Real, BlockSize<32,2> >(),
        pack_parameters(
            Param< ::par, cInOut>(par)
        ),
        define_temporaries(
            StencilBuffer<buffer, Real, KRange<FullDomain,0,0> >()
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<BufferInit, IJRange<cComplete,-1,1,-1,1>, KRange<FullDomain,0,0> >(),
                    StencilStage<ParameterInitFromBuffer, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )  
        )
    );

    // test outer boundary
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -1, 11, -2, 13);

    // test inner boundary
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 2, 7, 1, 9);

    // top left corner
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -2, 3, -1, 2);
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -2, -1, -2, -1); // corner only

    // top right corner
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 9, 12, -1, 2);
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 10, 12, -2, -1); // corner only

    // bottom left corner
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -2, 3, 10, 13);
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), -2, -1, 12, 13); // corner only

    // bottom right corner
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 9, 12, 10, 13);
    apply_boundary_round(testStencil, par, size.iSize(), size.jSize(), size.kSize(), 10, 12, 12, 13); // corner only
}


