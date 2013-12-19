#include <iostream>
#include <cmath>
#include "MatFile.h"
#include "Heat.h"

const double pi = 3.14159265358979;

Real printNorm(const IJKRealField& field, const std::string& name)
{
    const IJKSize size =  field.calculationDomain();
    Real norm = 0.;
    for (int i = 0; i < size.iSize(); ++i)
        for (int j = 0; j < size.jSize(); ++j)
            for (int k = 0; k < size.kSize(); ++k)
            {
                norm = std::max(std::abs(field(i,j,k)), norm);
            }
    std::cout << "Norm of " << name << " is " << norm << "\n";
}

enum
{
    coeff, dx2, dt,
    k, k1, k2, k3, k4,
    qmain, qlapl, qtemp
};


namespace HeatStages
{

    template<typename TEnv>
    struct LaplaceStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, qlapl)
        STAGE_PARAMETER(FullDomain, k)
        STAGE_PARAMETER(FullDomain, coeff)
        STAGE_PARAMETER(FullDomain, dx2)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
#ifdef NINETEEN_POINT_STENCIL
            T tens = -24. * ctx[qlapl::Center()]
                + 2. * (
                        ctx[qlapl::At(Offset<-1,  0,  0>())]
                     +  ctx[qlapl::At(Offset<+1,  0,  0>())]
                     +  ctx[qlapl::At(Offset< 0, -1,  0>())]
                     +  ctx[qlapl::At(Offset< 0, +1,  0>())]
                     +  ctx[qlapl::At(Offset< 0,  0, -1>())]
                     +  ctx[qlapl::At(Offset< 0,  0, +1>())]
               )
                     + ctx[qlapl::At(Offset< 0, -1, -1>())]
                     + ctx[qlapl::At(Offset< 0, -1, +1>())]
                     + ctx[qlapl::At(Offset< 0, +1, -1>())]
                     + ctx[qlapl::At(Offset< 0, +1, +1>())]
                     + ctx[qlapl::At(Offset<-1,  0, -1>())]
                     + ctx[qlapl::At(Offset<-1,  0, +1>())]
                     + ctx[qlapl::At(Offset<+1,  0, -1>())]
                     + ctx[qlapl::At(Offset<+1,  0, +1>())]
                     + ctx[qlapl::At(Offset<-1, -1,  0>())]
                     + ctx[qlapl::At(Offset<-1, +1,  0>())]
                     + ctx[qlapl::At(Offset<+1, -1,  0>())]
                     + ctx[qlapl::At(Offset<+1, +1,  0>())];
            tens /= 6. * ctx[dx2::Center()];
#else
# ifdef FIVE_POINT_STENCIL
            T d2q_di2 =
                - 1. * ctx[qlapl::At(Offset<-2, 0, 0>())]
                +16. * ctx[qlapl::At(Offset<-1, 0, 0>())]
                -30. * ctx[qlapl::At(Offset< 0, 0, 0>())]
                +16. * ctx[qlapl::At(Offset< 1, 0, 0>())]
                - 1. * ctx[qlapl::At(Offset< 2, 0, 0>())];
            T d2q_dj2 =
                - 1. * ctx[qlapl::At(Offset<0, -2, 0>())]
                +16. * ctx[qlapl::At(Offset<0, -1, 0>())]
                -30. * ctx[qlapl::At(Offset<0,  0, 0>())]
                +16. * ctx[qlapl::At(Offset<0,  1, 0>())]
                - 1. * ctx[qlapl::At(Offset<0,  2, 0>())];
            T d2q_dk2 =
                - 1. * ctx[qlapl::At(Offset<0, 0, -2>())]
                +16. * ctx[qlapl::At(Offset<0, 0, -1>())]
                -30. * ctx[qlapl::At(Offset<0, 0,  0>())]
                +16. * ctx[qlapl::At(Offset<0, 0,  1>())]
                - 1. * ctx[qlapl::At(Offset<0, 0,  2>())];
            T tens = ctx[coeff::Center()] * (d2q_di2 + d2q_dj2 + d2q_dk2);
            tens /= ctx[dx2::Center()] * 12.;
# else
            T d2q_di2 =
                +1. * ctx[qlapl::At(Offset<-1, 0, 0>())]
                -2. * ctx[qlapl::At(Offset< 0, 0, 0>())]
                +1. * ctx[qlapl::At(Offset< 1, 0, 0>())];
            T d2q_dj2 =
                +1. * ctx[qlapl::At(Offset<0, -1, 0>())]
                -2. * ctx[qlapl::At(Offset<0,  0, 0>())]
                +1. * ctx[qlapl::At(Offset<0,  1, 0>())];
            T d2q_dk2 =
                +1. * ctx[qlapl::At(Offset<0, 0, -1>())]
                -2. * ctx[qlapl::At(Offset<0, 0,  0>())]
                +1. * ctx[qlapl::At(Offset<0, 0,  1>())];
            T tens = ctx[coeff::Center()] * (d2q_di2 + d2q_dj2 + d2q_dk2);
            tens /= ctx[dx2::Center()];
# endif
#endif
            ctx[k::Center()] = tens;
        }
    };

    template<typename TEnv>
    struct EulerStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, qmain)
        STAGE_PARAMETER(FullDomain, qtemp)
        STAGE_PARAMETER(FullDomain, k)
        STAGE_PARAMETER(FullDomain, dt)

        static void Do(Context ctx, FullDomain)
        {
            const T dt_ = ctx[dt::Center()];
            ctx[qtemp::Center()] = ctx[qmain::Center()]
                                 + dt_*ctx[k::Center()];
        }
    };

    template<typename TEnv>
    struct RKStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, qmain)
        STAGE_PARAMETER(FullDomain, k1)
        STAGE_PARAMETER(FullDomain, k2)
        STAGE_PARAMETER(FullDomain, k3)
        STAGE_PARAMETER(FullDomain, k4)
        STAGE_PARAMETER(FullDomain, dt)

        static void Do(Context ctx, FullDomain)
        {
            ctx[qmain::Center()] = ctx[qmain::Center()] + ctx[dt::Center()]/6. *
                (
                         ctx[k1::Center()]
                   + 2.* ctx[k2::Center()]
                   + 2.* ctx[k3::Center()]
                   +     ctx[k4::Center()]
                );
        }
    };

}

Heat::Heat(IJKRealField& data, Real heatcoeff, Real stepsize, Real timestepsize, MPI_Comm comm)
    : qs_(data), ks_(data)
    , coeff_(heatcoeff), dx2_(stepsize*stepsize)
    , dt_(timestepsize), dthalf_(.5 * timestepsize)
    , he_(true, true, true, comm)
{
    IJKSize calculationDomain = data.calculationDomain();
    StencilCompiler::Build(
        laplaceStencil_,
        "Laplace",
        calculationDomain,
        StencilConfiguration<Real, BlockSize<8, 8> >(),
        pack_parameters(
            // Input fields
            Param<qlapl, cIn>(qs_.qlaplace()),
            // Output fields
            Param<k, cInOut>(ks_.k()),
            // Scalars
            Param<coeff, cScalar>(coeff_),
            Param<dx2, cScalar>(dx2_),
            Param<dt, cScalar>(dtparam_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<HeatStages::LaplaceStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
    std::cout << "Registerig qlapl as " << &qs_.qlaplace() << "\n";

    StencilCompiler::Build(
        eulerStencil_,
        "Euler",
        calculationDomain,
        StencilConfiguration<Real, BlockSize<8, 8> >(),
        pack_parameters(
            // Input fields
            Param<qmain, cIn>(qs_.qmain()),
            Param<k, cIn>(ks_.k()),
            // Output fields
            Param<qtemp, cInOut>(qs_.qtemp()),
            // Scalars
            Param<dt, cScalar>(dtparam_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<HeatStages::EulerStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        rkStencil_,
        "RungeKutta",
        calculationDomain,
        StencilConfiguration<Real, BlockSize<8, 8> >(),
        pack_parameters(
            // Input fields
            Param<k1, cIn>(ks_.k1()),
            Param<k2, cIn>(ks_.k2()),
            Param<k3, cIn>(ks_.k3()),
            Param<k4, cIn>(ks_.k4()),
            // Output fields
            Param<qmain, cInOut>(qs_.qmain()),
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

    he_.registerField(qs_.qlaplace());
}

void Heat::DoTimeStep()
{
    /* First RK timestep */
    ks_.set(1);
    qs_.setLaplaceMain();
    he_.exchange();
    laplaceStencil_.Apply();
    // Now: k1 = laplace(qmain)

    /* Second RK timestep */
    qs_.restore();
    dtparam_ = dthalf_;
    eulerStencil_.Apply();
    // Now: qtemp = qmain + dthalf_*k1

    ks_.set(2);
    qs_.setLaplaceTemp();
    he_.exchange();
    laplaceStencil_.Apply();
    // Now: k2 = laplace(qmain + dthalf_*k1)
    
    /* Third RK timestep */
    qs_.restore();
    eulerStencil_.Apply();
    // Now: qtemp = qmain + dthalf_*k2

    ks_.set(3);
    qs_.setLaplaceTemp();
    he_.exchange();
    laplaceStencil_.Apply();
    // Now: k3 = laplace(qmain + dthalf_*k2)
    
    /* Fourth RK timestep */
    qs_.restore();
    dtparam_ = dt_;
    eulerStencil_.Apply();
    // Now: qtemp = qmain + dt_*k3

    ks_.set(4);
    qs_.setLaplaceTemp();
    he_.exchange();
    laplaceStencil_.Apply();
    // Now: k4 = laplace(qmain + dt_*k3)

    /* Final RK stage: put things together */
    ks_.restore();
    qs_.restore();
    rkStencil_.Apply();
}

void Heat::DoTest()
{
    // Stencil teststencil;
    // IJKSize domain = q_.calculationDomain();
    // KBoundary kboundary;
    // kboundary.Init(-3,3);
    // IJKRealField test;
    // test.Init("test", domain, kboundary);

    // double one = 1.;

    // StencilCompiler::Build(
    //     teststencil,
    //     "Test",
    //     domain,
    //     StencilConfiguration<Real, BlockSize<8,8> >(),
    //     pack_parameters(
    //         Param<q, cIn>(q_),
    //         Param<qtens, cInOut>(test),
    //         Param<coeff, cScalar>(one),
    //         Param<dx2, cScalar>(dx2_)
    //     ),
    //     define_loops(
    //         define_sweep<cKIncrement>(
    //             define_stages(
    //                 StencilStage<HeatStages::DiffStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
    //             )
    //         )
    //     )
    // );

    // teststencil.Apply();

    // MatFile errfile("error.mat");
    // IJKRealField errfield;
    // errfield.Init("error", domain, kboundary);

    // // Test
    // for (int i = 0; i < domain.iSize(); ++i)
    //     for (int j = 0; j < domain.jSize(); ++j)
    //         for (int k = 0; k < domain.kSize(); ++k)
    //         {
    //             double exact = -3*pi*pi*q_(i, j,k);

    //             double actual = test(i, j, k);

    //             double error = std::abs(exact - actual);

    //             errfield(i, j, k) = error;
    //             if (error > 1e-10)
    //             {
    //                 std::cerr << "Error at ("
    //                     << i+1 << "," << j+1 << "," << k+1 << ") : "
    //                     << actual << " instead of " << exact << "\n";
    //             }
    //         }

    // errfile.addField(errfield, 0);
}
