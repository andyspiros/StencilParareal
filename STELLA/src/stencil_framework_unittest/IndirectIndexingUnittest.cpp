#include <boost/mpl/integral_c.hpp>
#include "gtest/gtest.h"
#include "StencilFramework.h"

#include <cstdlib>

// define the parameter enum
enum { in, outStatic, outDynamic, scalarOffset, i_idx, j_idx, k_idx, success};

namespace NinePoint{
    /******************************************************************************
     * performance tests for a simple 9 point stencil
     * implemented with both static and dynamic indexng
     ******************************************************************************/
    template<typename TEnv>
    struct StaticFilter9Stage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, in)
        STAGE_PARAMETER(FullDomain, outStatic)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[outStatic::Center()] = T(1./2.)*ctx[in::Center()] +
                                       T(1./16.) *
                                       (
                                            ctx[in::At(Offset<-1,-1,0>())] + ctx[in::At(Offset<-1,0,0>())] + ctx[in::At(Offset<-1,1,0>())] +
                                            ctx[in::At(Offset< 0,-1,0>())]                               + ctx[in::At(Offset< 0,1,0>())] +
                                            ctx[in::At(Offset< 1,-1,0>())] + ctx[in::At(Offset< 1,0,0>())] + ctx[in::At(Offset< 1,1,0>())]
                                       );
        }
    };

    template<typename TEnv>
    struct DynamicFilter9Stage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, in)
        STAGE_PARAMETER(FullDomain, outDynamic)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[outDynamic::Center()] = T(1./2.)*ctx[in::At(0,0,0)] +
                                       T(1./16.) *
                                       (
                                            ctx[in::At(-1,-1,0)] + ctx[in::At(-1,0,0)] + ctx[in::At(-1,1,0)] +
                                            ctx[in::At( 0,-1,0)]                       + ctx[in::At( 0,1,0)] +
                                            ctx[in::At( 1,-1,0)] + ctx[in::At( 1,0,0)] + ctx[in::At( 1,1,0)]
                                       );
        }
    };
}

// check that indirect indexing works for 9-point stencil
TEST(IndirectIndexingUnittest, Filter9Point)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(256, 256, 60);
    //size.Init(32, 4, 1);
    boundary.Init(0, 0);

    // define the fields
    IJKRealField inField;
    IJKRealField outFieldStatic;
    IJKRealField outFieldDynamic;
    inField.Init("input", size, boundary);
    outFieldStatic.Init("outStatic", size, boundary);
    outFieldDynamic.Init("outDynamic", size, boundary);

    // initialise the input field to random values
    for(int k=0; k<size.kSize(); k++)
        for(int i=-3; i<size.iSize()+3; ++i)
            for(int j=-3; j<size.jSize()+3; ++j)
                inField(i,j,k) = outFieldDynamic(i,j,k) = outFieldStatic(i,j,k) = 0.;

    for(int i=0; i<size.iSize(); ++i)
        for(int j=0; j<size.jSize(); ++j)
            for(int k=0; k<size.kSize(); k++)
                inField(i,j,k) = double(rand()%100)/100.;

    //inField.SynchronizeHostStorage();
    //outFieldStatic.SynchronizeHostStorage();
    //outFieldDynamic.SynchronizeHostStorage();

    using namespace NinePoint;
    Stencil staticStencil;
    StencilCompiler::Build(
        staticStencil,
        "StaticStencilFilter9",
        size,
#ifdef __CUDA_BACKEND__
        StencilConfiguration<double, BlockSize<32,4> >(),
#else
        StencilConfiguration<double, BlockSize<8,8> >(),
#endif
        pack_parameters(
            Param<in, cIn>(inField),
            Param<outStatic, cInOut>(outFieldStatic)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StaticFilter9Stage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    Stencil dynamicStencil;
    StencilCompiler::Build(
        dynamicStencil,
        "DynamicStencilFilter9",
        size,
#ifdef __CUDA_BACKEND__
        StencilConfiguration<double, BlockSize<32,4> >(),
#else
        StencilConfiguration<double, BlockSize<8,8> >(),
#endif
        pack_parameters(
            Param<in, cIn>(inField),
            Param<outDynamic, cInOut>(outFieldDynamic)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<DynamicFilter9Stage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );


    // apply once, then reset meters
    dynamicStencil.Apply();
    staticStencil.Apply();
    dynamicStencil.ResetMeters();
    staticStencil.ResetMeters();


    // apply the stencils multiple times
    int nSteps=25;
    for(int i=0; i<nSteps; i++)
    {
        staticStencil.Apply();
        dynamicStencil.Apply();
    }

    //outFieldStatic.SynchronizeDeviceStorage();
    //outFieldDynamic.SynchronizeDeviceStorage();

    std::cout << staticStencil.MetersToString() << std::endl;
    std::cout << dynamicStencil.MetersToString() << std::endl;

    // check the sum was done right
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(outFieldStatic(i,j,k), outFieldDynamic(i,j,k));
            }
        }
    }
}

namespace Laplace3D{
    /******************************************************************************
     * performance tests for a simple 3D Laplacian (7 points)
     * implemented with both static and dynamic indexng
     ******************************************************************************/
    template<typename TEnv>
    struct StaticLaplacianStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, in)
        STAGE_PARAMETER(FullDomain, outStatic)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[outStatic::Center()] =  T(4.)*ctx[in::At(Offset< 0, 0, 0>())]
                                             -ctx[in::At(Offset<-1, 0, 0>())]
                                             -ctx[in::At(Offset< 1, 0, 0>())]
                                             -ctx[in::At(Offset< 0,-1, 0>())]
                                             -ctx[in::At(Offset< 0, 1, 0>())]
                                             -ctx[in::At(Offset< 0, 0,-1>())]
                                             -ctx[in::At(Offset< 0, 0, 1>())];
        }
    };

    template<typename TEnv>
    struct DynamicLaplacianStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, in)
        STAGE_PARAMETER(FullDomain, outDynamic)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[outDynamic::Center()] = T(4.)*ctx[in::At( 0, 0, 0)]
                                             -ctx[in::At(-1, 0, 0)]
                                             -ctx[in::At( 1, 0, 0)]
                                             -ctx[in::At( 0,-1, 0)]
                                             -ctx[in::At( 0, 1, 0)]
                                             -ctx[in::At( 0, 0,-1)]
                                             -ctx[in::At( 0, 0, 1)];
        }
    };
}

// check that indirect indexing works for 9-point stencil
TEST(IndirectIndexingUnittest, Laplacian3D)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(256, 256, 60);
    boundary.Init(-1, 1);

    // define the fields
    IJKRealField inField;
    IJKRealField outFieldStatic;
    IJKRealField outFieldDynamic;
    inField.Init("input", size, boundary);
    outFieldStatic.Init("outStatic", size, boundary);
    outFieldDynamic.Init("outDynamic", size, boundary);

    // initialise the input field to random values
    for(int k=-1; k<size.kSize()+1; k++)
        for(int i=-3; i<size.iSize()+3; ++i)
            for(int j=-3; j<size.jSize()+3; ++j)
                inField(i,j,k) = outFieldDynamic(i,j,k) = outFieldStatic(i,j,k) = 0.;

    for(int i=0; i<size.iSize(); ++i)
        for(int j=0; j<size.jSize(); ++j)
            for(int k=0; k<size.kSize(); k++)
                inField(i,j,k) = double(rand()%100)/100.;

    using namespace Laplace3D;
    Stencil staticStencil;
    StencilCompiler::Build(
        staticStencil,
        "StaticStencilLaplacian",
        size,
#ifdef __CUDA_BACKEND__
        StencilConfiguration<double, BlockSize<32,4> >(),
#else
        StencilConfiguration<double, BlockSize<8,8> >(),
#endif
        pack_parameters(
            Param<in, cIn>(inField),
            Param<outStatic, cInOut>(outFieldStatic)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StaticLaplacianStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    Stencil dynamicStencil;
    StencilCompiler::Build(
        dynamicStencil,
        "DynamicStencilLaplacian",
        size,
#ifdef __CUDA_BACKEND__
        StencilConfiguration<double, BlockSize<32,4> >(),
#else
        StencilConfiguration<double, BlockSize<8,8> >(),
#endif
        pack_parameters(
            Param<in, cIn>(inField),
            Param<outDynamic, cInOut>(outFieldDynamic)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<DynamicLaplacianStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    // apply once, then reset meters
    dynamicStencil.Apply();
    staticStencil.Apply();
    dynamicStencil.ResetMeters();
    staticStencil.ResetMeters();


    // apply the stencils multiple times
    int nSteps=25;
    for(int i=0; i<nSteps; i++)
    {
        staticStencil.Apply();
        dynamicStencil.Apply();
    }
    /*
    outFieldStatic.SynchronizeDeviceStorage();
    outFieldDynamic.SynchronizeDeviceStorage();
    std::cout << "================" << "output dynamic" << "================" << std::endl;
    for(int k=0; k<size.kSize(); k++)
    {
        for(int i=-3; i<size.iSize()+3; ++i)
        {
            for(int j=-3; j<size.jSize()+3; ++j)
            {
                std::cout << outFieldDynamic(i,j,k) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "================" << "output static" << "================" << std::endl;
    for(int k=0; k<size.kSize(); k++)
    {
        for(int i=-3; i<size.iSize()+3; ++i)
        {
            for(int j=-3; j<size.jSize()+3; ++j)
            {
                std::cout << outFieldStatic(i,j,k) << " " ;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    */

    std::cout << staticStencil.MetersToString() << std::endl;
    std::cout << dynamicStencil.MetersToString() << std::endl;

    // check the sum was done right
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(outFieldStatic(i,j,k), outFieldDynamic(i,j,k));
            }
        }
    }
}

namespace RandAccess{
    /******************************************************************************
     * performance tests for a simple 9 point stencil
     * implemented with both static and dynamic indexng
     ******************************************************************************/
    template<typename TEnv>
    struct StaticRandStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, in)
        STAGE_PARAMETER(FullDomain, outStatic)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            if( ctx[in::Center()] > 0.5 )
            {
                ctx[outStatic::Center()] =
                    T(0.5) * (ctx[in::At(Offset<-1,0,0>())] + ctx[in::At(Offset<0,-1,0>())]);
            }else{
                ctx[outStatic::Center()] =
                    T(0.5) * (ctx[in::At(Offset<1,0,0>())] + ctx[in::At(Offset<0,1,0>())]);
            }
        }
    };

    template<typename TEnv>
    struct DynamicRandStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, in)
        STAGE_PARAMETER(FullDomain, outDynamic)
        STAGE_PARAMETER(FullDomain, scalarOffset)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            const int index = (ctx[in::Center()] > 0.5) ? -1 : 1;
            ctx[outDynamic::Center()] =
                T(0.5) * (
                    ctx[in::At(index, 0, 0)] +
                    ctx[in::At(0 ,index, 0)]
                );
        }
    };
}

// check that indirect indexing works for 9-point stencil
TEST(IndirectIndexingUnittest, RandAccess)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(256, 256, 60);
    boundary.Init(-1, 1);

    int myOffset=0;

    // define the fields
    IJKRealField inField;
    IJKRealField outFieldStatic;
    IJKRealField outFieldDynamic;
    inField.Init("input", size, boundary);
    outFieldStatic.Init("outStatic", size, boundary);
    outFieldDynamic.Init("outDynamic", size, boundary);

    // initialise the input field to random values
    for(int k=-1; k<size.kSize()+1; k++)
        for(int i=-3; i<size.iSize()+3; ++i)
            for(int j=-3; j<size.jSize()+3; ++j)
                inField(i,j,k) = outFieldDynamic(i,j,k) = outFieldStatic(i,j,k) = 0.;

    for(int i=0; i<size.iSize(); ++i)
        for(int j=0; j<size.jSize(); ++j)
            for(int k=0; k<size.kSize(); k++)
                inField(i,j,k) = double(rand()%100)/100.;

    using namespace RandAccess;
    Stencil staticStencil;
    StencilCompiler::Build(
        staticStencil,
        "StaticStencilRandAccess",
        size,
#ifdef __CUDA_BACKEND__
        StencilConfiguration<double, BlockSize<32,4> >(),
#else
        StencilConfiguration<double, BlockSize<8,8> >(),
#endif
        pack_parameters(
            Param<in, cIn>(inField),
            Param<outStatic, cInOut>(outFieldStatic),
            Param<scalarOffset, cScalar>(myOffset)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<StaticRandStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    Stencil dynamicStencil;
    StencilCompiler::Build(
        dynamicStencil,
        "DynamicStencilRandAccess",
        size,
#ifdef __CUDA_BACKEND__
        StencilConfiguration<double, BlockSize<32,4> >(),
#else
        StencilConfiguration<double, BlockSize<8,8> >(),
#endif
        pack_parameters(
            Param<in, cIn>(inField),
            Param<outDynamic, cInOut>(outFieldDynamic),
            Param<scalarOffset, cScalar>(myOffset)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<DynamicRandStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    // apply once, then reset meters
    dynamicStencil.Apply();
    staticStencil.Apply();
    dynamicStencil.ResetMeters();
    staticStencil.ResetMeters();


    // apply the stencils multiple times
    int nSteps=25;
    for(int i=0; i<nSteps; i++)
    {
        staticStencil.Apply();
        dynamicStencil.Apply();
    }

    ////outFieldStatic.SynchronizeDeviceStorage();
    //outFieldDynamic.SynchronizeDeviceStorage();

    std::cout << staticStencil.MetersToString() << std::endl;
    std::cout << dynamicStencil.MetersToString() << std::endl;

    // check the sum was done right
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(outFieldStatic(i,j,k), outFieldDynamic(i,j,k));
            }
        }
    }
}


namespace AbsoluteIndex {


    template<typename TEnv>
    struct AbsoluteIndexStage 
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, i_idx)
        STAGE_PARAMETER(FullDomain, j_idx)
        STAGE_PARAMETER(FullDomain, k_idx)
        STAGE_PARAMETER(FullDomain, success)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            int myi = ctx[i_idx::Center()];
            int myj = ctx[j_idx::Center()];
            int myk = ctx[k_idx::Center()];

            bool success = true;
            
            for (int i = -3; i <= 3; ++i)
                for (int j = -3; j <= 3; ++j)
                    for (int k = -1; k <= 1; ++k)
                    {
                        success = success && 
                            (ctx[i_idx::At(i, j, k)] == myi + i) &&
                            (ctx[j_idx::At(i, j, k)] == myj + j) &&
                            (ctx[k_idx::At(i, j, k)] == myk + k);
                    }

            ctx[success::Center()] = success ? 1 : 0;
        }
    };
}

TEST(IndirectIndexingUnittest, AbsoluteIndex)
{
    IJKSize size;
    KBoundary boundary;
    size.Init(256, 256, 60);
    boundary.Init(-1, 1);

    // define the fields
    IRealField iField;
    JRealField jField;
    KRealField kField;
    IJKRealField successField;
    iField.Init("i", size, boundary);
    jField.Init("j", size, boundary);
    kField.Init("k", size, boundary);
    successField.Init("success", size, boundary);

    // initialise the index fields
    for(int k=-1; k<size.kSize()+1; k++)
        for(int i=-3; i<size.iSize()+3; ++i)
            for(int j=-3; j<size.jSize()+3; ++j)
            {
                iField(i,j,k) = i;
                jField(i,j,k) = j;
                kField(i,j,k) = k;
            }

    using namespace AbsoluteIndex;
    Stencil absStencil;
    StencilCompiler::Build(
        absStencil,
        "AbsoluteIndexStencil",
        size,
#ifdef __CUDA_BACKEND__
        StencilConfiguration<double, BlockSize<32,4> >(),
#else
        StencilConfiguration<double, BlockSize<8,8> >(),
#endif
        pack_parameters(
            Param<i_idx, cIn>(iField),
            Param<j_idx, cIn>(jField),
            Param<k_idx, cIn>(kField),
            Param<success, cInOut>(successField)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<AbsoluteIndexStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    // apply once, then reset meters
    absStencil.Apply();

    // check the success field
    for(int i=0; i<size.iSize(); ++i)
    {
        for(int j=0; j<size.jSize(); ++j)
        {
            for(int k=0; k<size.kSize(); k++)
            {
                ASSERT_DOUBLE_EQ(successField(i,j,k), 1.);
            }
        }
    }
}

