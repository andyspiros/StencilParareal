#include <iostream>
#include <cmath>
#include "MatFile.h"
#include "SharedInfrastructure.h"
#include "StencilFramework.h"
#include "mpi.h"

#include "Convection.h"
#include "ConvectionFunctions.h"

enum
{
    nu, cx, cy, cz,
    dx, dx2, dt, dteuler,
    k, k1, k2, k3, k4,
    q, qeuler, qrhs
};


namespace ConvectionStages
{

    template<typename TEnv>
    struct RHSStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, qrhs)
        STAGE_PARAMETER(FullDomain, k)
        STAGE_PARAMETER(FullDomain, nu)
        STAGE_PARAMETER(FullDomain, cx)
        STAGE_PARAMETER(FullDomain, cy)
        STAGE_PARAMETER(FullDomain, cz)
        STAGE_PARAMETER(FullDomain, dx)
        STAGE_PARAMETER(FullDomain, dx2)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[k::Center()] = 

                // Compute advection
                ctx[Call<Advection>::With(
                        qrhs::Center(),
                        dx::Center(),
                        cx::Center(),
                        cy::Center(),
                        cz::Center()
                    )]

                +
                
                // Compute laplacian
                ctx[nu::Center()] * ctx[Call<Laplace>::With(
                        qrhs::Center(),
                        dx2::Center()
                )];
        }
    };

    template<typename TEnv>
    struct EulerStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, q)
        STAGE_PARAMETER(FullDomain, qeuler)
        STAGE_PARAMETER(FullDomain, k)
        STAGE_PARAMETER(FullDomain, dteuler)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[qeuler::Center()] =
                  ctx[q::Center()]
                + ctx[dteuler::Center()]*ctx[k::Center()];
        }
    };


    template<typename TEnv>
    struct RKStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, q)
        STAGE_PARAMETER(FullDomain, k1)
        STAGE_PARAMETER(FullDomain, k2)
        STAGE_PARAMETER(FullDomain, k3)
        STAGE_PARAMETER(FullDomain, k4)
        STAGE_PARAMETER(FullDomain, dt)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[q::Center()] = ctx[q::Center()] + ctx[dt::Center()]/6. *
                (
                         ctx[k1::Center()]
                   + 2.* ctx[k2::Center()]
                   + 2.* ctx[k3::Center()]
                   +     ctx[k4::Center()]
                );
        }
    };

}

Convection::Convection(ConvectionRealField& data,
                       Real nucoeff, Real cxcoeff, Real cycoeff, Real czcoeff,
                       Real stepsize, Real timestepsize, MPI_Comm comm)
    : q_(data)
    , nu_(nucoeff), cx_(cxcoeff), cy_(cycoeff), cz_(czcoeff)
    , dx_(stepsize), dx2_(stepsize*stepsize)
    , dt_(timestepsize), dthalf_(.5 * timestepsize)
    , heQ_(true, true, true, comm)
    , heQinternal_(true, true, true, comm)
    , matq("q.mat"), matk1("k1.mat"), matk2("k2.mat"), matk3("k3.mat"), matk4("k4.mat")
{
#ifdef __CUDA_BACKEND__
    typedef BlockSize<32, 4> BSize;
#else
    typedef BlockSize<8, 8> BSize;
#endif

    using namespace ConvectionStages;

    stencil1_ = new Stencil;
    stencil2_ = new Stencil;
    stencil3_ = new Stencil;
    stencil4_ = new Stencil;
    
    IJKBoundary inputBoundary = data.boundary();
    assert(inputBoundary.iMinusOffset() == -2);
    assert(inputBoundary.iPlusOffset() == 2);
    assert(inputBoundary.jMinusOffset() == -2);
    assert(inputBoundary.jPlusOffset() == 2);
    assert(inputBoundary.kMinusOffset() == -2);
    assert(inputBoundary.kPlusOffset() == 2);

    IJKSize calculationDomain = data.calculationDomain();
    KBoundary kboundary;
    kboundary.Init(-convectionBoundaryLines, convectionBoundaryLines);

    // Initialize temporary fields
    qinternal_.Init("qinternal", calculationDomain, kboundary);
    k1_.Init("k1", calculationDomain, kboundary);
    k2_.Init("k2", calculationDomain, kboundary);
    k3_.Init("k3", calculationDomain, kboundary);
    k4_.Init("k4", calculationDomain, kboundary);

    // Initialize joker fields
    //qtmp_.Init("qtmp", q_);
    //k_.Init("k", k1_);

    StencilCompiler::Build(
        *stencil1_,
        "Convection1",
        calculationDomain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Main field
            Param<q, cIn>(q_),
            Param<qrhs, cIn>(q_),
            Param<qeuler, cInOut>(qinternal_),

            // k fields
            Param<k, cInOut>(k1_),

            // Scalars
            Param<nu, cScalar>(nu_),
            Param<cx, cScalar>(cx_),
            Param<cy, cScalar>(cy_),
            Param<cz, cScalar>(cz_),
            Param<dx, cScalar>(dx_),
            Param<dx2, cScalar>(dx2_),
            Param<dteuler, cScalar>(dthalf_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<RHSStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<EulerStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        *stencil2_,
        "Convection2",
        calculationDomain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Main field
            Param<q, cIn>(q_),
            Param<qrhs, cIn>(qinternal_),
            Param<qeuler, cInOut>(qinternal_),

            // k fields
            Param<k, cInOut>(k2_),

            // Scalars
            Param<nu, cScalar>(nu_),
            Param<cx, cScalar>(cx_),
            Param<cy, cScalar>(cy_),
            Param<cz, cScalar>(cz_),
            Param<dx, cScalar>(dx_),
            Param<dx2, cScalar>(dx2_),
            Param<dteuler, cScalar>(dthalf_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<RHSStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<EulerStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        *stencil3_,
        "Convection3",
        calculationDomain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Main field
            Param<q, cIn>(q_),
            Param<qrhs, cIn>(qinternal_),
            Param<qeuler, cInOut>(qinternal_),

            // k fields
            Param<k, cInOut>(k3_),

            // Scalars
            Param<nu, cScalar>(nu_),
            Param<cx, cScalar>(cx_),
            Param<cy, cScalar>(cy_),
            Param<cz, cScalar>(cz_),
            Param<dx, cScalar>(dx_),
            Param<dx2, cScalar>(dx2_),
            Param<dteuler, cScalar>(dt_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<RHSStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<EulerStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        *stencil4_,
        "Convection4",
        calculationDomain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Main field
            Param<q, cInOut>(q_),
            Param<qrhs, cInOut>(qinternal_),

            // k fields
            Param<k, cInOut>(k4_),
            Param<k1, cIn>(k1_),
            Param<k2, cIn>(k2_),
            Param<k3, cIn>(k3_),
            Param<k4, cIn>(k4_),

            // Scalars
            Param<nu, cScalar>(nu_),
            Param<cx, cScalar>(cx_),
            Param<cy, cScalar>(cy_),
            Param<cz, cScalar>(cz_),
            Param<dx, cScalar>(dx_),
            Param<dx2, cScalar>(dx2_),
            Param<dt, cScalar>(dt_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<RHSStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >(),
                    StencilStage<RKStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );


    heQ_.registerField(q_);
    heQinternal_.registerField(qinternal_);

    // MAT files
    //matq.startCell("q", 32);
    //matk1.startCell("k1", 32);
    //matk2.startCell("k2", 32);
    //matk3.startCell("k3", 32);
    //matk4.startCell("k4", 32);
}

Convection::~Convection()
{
    delete stencil1_;
    delete stencil2_;
    delete stencil3_;
    delete stencil4_;
}

void Convection::DoTimeStep()
{
    heQ_.exchange();
    stencil1_->Apply();

    heQinternal_.exchange();
    stencil2_->Apply();

    heQinternal_.exchange();
    stencil3_->Apply();

    heQinternal_.exchange();
    stencil4_->Apply();
}

