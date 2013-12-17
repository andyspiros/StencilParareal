#include <iostream>
#include "Heat.h"

enum
{
    q, qstart, qend, qtens, coeff, di2, dj2, dk2, dt,
    k1, k2, k3, k4
};


namespace HeatStages
{
    template<typename TEnv>
    struct DiffStageAndUpdate
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, q)
        STAGE_PARAMETER(FullDomain, qstart)
        STAGE_PARAMETER(FullDomain, qend)
        STAGE_PARAMETER(FullDomain, qtens)
        STAGE_PARAMETER(FullDomain, coeff)
        STAGE_PARAMETER(FullDomain, di2)
        STAGE_PARAMETER(FullDomain, dj2)
        STAGE_PARAMETER(FullDomain, dk2)
        STAGE_PARAMETER(FullDomain, dt)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            T d2q_di2 =
                - 1. * ctx[q::At(Offset<-2, 0, 0>())]
                +16. * ctx[q::At(Offset<-1, 0, 0>())]
                -30. * ctx[q::At(Offset< 0, 0, 0>())]
                +16. * ctx[q::At(Offset< 1, 0, 0>())]
                - 1. * ctx[q::At(Offset< 2, 0, 0>())];
            T d2q_dj2 =
                - 1. * ctx[q::At(Offset<0, -2, 0>())]
                +16. * ctx[q::At(Offset<0, -1, 0>())]
                -30. * ctx[q::At(Offset<0,  0, 0>())]
                +16. * ctx[q::At(Offset<0,  1, 0>())]
                - 1. * ctx[q::At(Offset<0,  2, 0>())];
            T d2q_dk2 =
                - 1. * ctx[q::At(Offset<0, 0, -2>())]
                +16. * ctx[q::At(Offset<0, 0, -1>())]
                -30. * ctx[q::At(Offset<0, 0,  0>())]
                +16. * ctx[q::At(Offset<0, 0,  1>())]
                - 1. * ctx[q::At(Offset<0, 0,  2>())];
            T tens = ctx[coeff::Center()] * (
                   d2q_di2 / ctx[di2::Center()]
                 + d2q_dj2 / ctx[dj2::Center()]
                 + d2q_dk2 / ctx[dk2::Center()]
            ) / 12.;
            ctx[qtens::Center()] = tens;
            ctx[qend::Center()] = ctx[qstart::Center()] + ctx[dt::Center()] * tens;
        }
    };

    template<typename TEnv>
    struct DiffStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, q)
        STAGE_PARAMETER(FullDomain, qtens)
        STAGE_PARAMETER(FullDomain, coeff)
        STAGE_PARAMETER(FullDomain, di2)
        STAGE_PARAMETER(FullDomain, dj2)
        STAGE_PARAMETER(FullDomain, dk2)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            T d2q_di2 =
                - 1. * ctx[q::At(Offset<-2, 0, 0>())]
                +16. * ctx[q::At(Offset<-1, 0, 0>())]
                -30. * ctx[q::At(Offset< 0, 0, 0>())]
                +16. * ctx[q::At(Offset< 1, 0, 0>())]
                - 1. * ctx[q::At(Offset< 2, 0, 0>())];
            T d2q_dj2 =
                - 1. * ctx[q::At(Offset<0, -2, 0>())]
                +16. * ctx[q::At(Offset<0, -1, 0>())]
                -30. * ctx[q::At(Offset<0,  0, 0>())]
                +16. * ctx[q::At(Offset<0,  1, 0>())]
                - 1. * ctx[q::At(Offset<0,  2, 0>())];
            T d2q_dk2 =
                - 1. * ctx[q::At(Offset<0, 0, -2>())]
                +16. * ctx[q::At(Offset<0, 0, -1>())]
                -30. * ctx[q::At(Offset<0, 0,  0>())]
                +16. * ctx[q::At(Offset<0, 0,  1>())]
                - 1. * ctx[q::At(Offset<0, 0,  2>())];
            T tens = ctx[coeff::Center()] * (
                   d2q_di2 / ctx[di2::Center()]
                 + d2q_dj2 / ctx[dj2::Center()]
                 + d2q_dk2 / ctx[dk2::Center()]
            ) / 12.;
            ctx[qtens::Center()] = tens;
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

Heat::Heat(IJKRealField& data, Real heatcoeff, Real stepi, Real stepj, Real stepk, Real timestepsize)
    : q_(data), coeff_(heatcoeff), di2_(stepi*stepi),  dj2_(stepj*stepj), dk2_(stepk*stepk),
      dt_(timestepsize), dthalf_(0.5*dt_)
{
    IJKSize calculationDomain = q_.calculationDomain();
    KBoundary kboundary;
    kboundary.Init(-3, 3);

    qtempin_.Init("qtemp", calculationDomain, kboundary);
    qtempout_.Init("qtemp", calculationDomain, kboundary);
    k1_.Init("k1", calculationDomain, kboundary);
    k2_.Init("k2", calculationDomain, kboundary);
    k3_.Init("k3", calculationDomain, kboundary);
    k4_.Init("k4", calculationDomain, kboundary);

    StencilCompiler::Build(
        stencil1_,
        "Heat",
        calculationDomain,
        StencilConfiguration<Real, BlockSize<8, 8> >(),
        pack_parameters(
            // Input fields
            Param<q, cIn>(q_),
            Param<qstart, cIn>(qtempin_),
            // Output fields
            Param<qend, cInOut>(qtempout_),
            Param<qtens, cInOut>(k1_),
            // Scalars
            Param<coeff, cScalar>(coeff_),
            Param<di2, cScalar>(di2_),
            Param<dj2, cScalar>(dj2_),
            Param<dk2, cScalar>(dk2_),
            Param<dt, cScalar>(dthalf_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<HeatStages::DiffStageAndUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        stencil2_,
        "Heat",
        calculationDomain,
        StencilConfiguration<Real, BlockSize<8, 8> >(),
        pack_parameters(
            // Input fields
            Param<q, cIn>(q_),
            Param<qstart, cIn>(qtempin_),
            // Output fields
            Param<qend, cInOut>(qtempout_),
            Param<qtens, cInOut>(k2_),
            // Scalars
            Param<coeff, cScalar>(coeff_),
            Param<di2, cScalar>(di2_),
            Param<dj2, cScalar>(dj2_),
            Param<dk2, cScalar>(dk2_),
            Param<dt, cScalar>(dthalf_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<HeatStages::DiffStageAndUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        stencil3_,
        "Heat",
        calculationDomain,
        StencilConfiguration<Real, BlockSize<8, 8> >(),
        pack_parameters(
            // Input fields
            Param<q, cIn>(q_),
            Param<qstart, cIn>(qtempin_),
            // Output fields
            Param<qend, cInOut>(qtempout_),
            Param<qtens, cInOut>(k3_),
            // Scalars
            Param<coeff, cScalar>(coeff_),
            Param<di2, cScalar>(di2_),
            Param<dj2, cScalar>(dj2_),
            Param<dk2, cScalar>(dk2_),
            Param<dt, cScalar>(dt_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<HeatStages::DiffStageAndUpdate, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        stencil4_,
        "Heat",
        calculationDomain,
        StencilConfiguration<Real, BlockSize<8, 8> >(),
        pack_parameters(
            // Input fields
            Param<q, cIn>(q_),
            // Output fields
            Param<qtens, cInOut>(k3_),
            // Scalars
            Param<coeff, cScalar>(coeff_),
            Param<di2, cScalar>(di2_),
            Param<dj2, cScalar>(dj2_),
            Param<dk2, cScalar>(dk2_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<HeatStages::DiffStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        stencilRK_,
        "HeatRK",
        calculationDomain,
        StencilConfiguration<Real, BlockSize<8, 8> >(),
        pack_parameters(
            // Input fields
            Param<k1, cIn>(k1_),
            Param<k2, cIn>(k2_),
            Param<k3, cIn>(k3_),
            Param<k4, cIn>(k4_),
            // Output fields
            Param<q, cInOut>(q_),
            // Scalars
            Param<dt, cScalar>(dt_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<HeatStages::RKStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
}

void Heat::DoTimeStep()
{
    stencil1_.Apply();
    qtempin_.SwapWith(qtempout_);
    stencil2_.Apply();
    qtempin_.SwapWith(qtempout_);
    stencil3_.Apply();
    qtempin_.SwapWith(qtempout_);
    stencil4_.Apply();

    stencilRK_.Apply();
}
