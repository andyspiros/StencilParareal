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
    nu, cx, cy, cz, dx, dx2, dt,
    k, k1, k2, k3, k4,
    qmain, qrhs, qtemp
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
            // Compute advection
            T tensadv = ctx[Call<Advection>::With(qrhs::Center(), dx::Center(), cx::Center(), cy::Center(), cz::Center())];
            
            // Compute laplacian
            T tenslapl = ctx[nu::Center()] * ctx[Call<Laplace>::With(qrhs::Center(), dx2::Center())];

            ctx[k::Center()] = tensadv + tenslapl;
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

        __ACC__
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

        __ACC__
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

Convection::Convection(IJKRealField& data,
                       Real nucoeff, Real cxcoeff, Real cycoeff, Real czcoeff,
                       Real stepsize, Real timestepsize, MPI_Comm comm)
    : nu_(nucoeff), cx_(cxcoeff), cy_(cycoeff), cz_(czcoeff)
    , dx_(stepsize), dx2_(stepsize*stepsize)
    , dt_(timestepsize), dthalf_(.5 * timestepsize)
    , qMain_(data)
    , he1_(true, true, true, comm)
    , he2_(true, true, true, comm)
    , he3_(true, true, true, comm)
    , he4_(true, true, true, comm)
{
#ifdef __CUDA_BACKEND__
    typedef BlockSize<32, 4> BSize;
#else
    typedef BlockSize<8, 8> BSize;
#endif

    using namespace ConvectionStages;

    rhsStencil_ = new Stencil;
    eulerStencil_ = new Stencil;
    rkStencil_ = new Stencil;
    IJKSize calculationDomain = data.calculationDomain();
    KBoundary kboundary;
    kboundary.Init(-3, 3);

    // Initialize internal fields
    qInternal_.Init("qInternal", calculationDomain, kboundary);
    k1_.Init("k1", calculationDomain, kboundary);
    k2_.Init("k2", calculationDomain, kboundary);
    k3_.Init("k3", calculationDomain, kboundary);
    k4_.Init("k4", calculationDomain, kboundary);

    // Initialize joker fields
    qrhs_.Init("qrhs", qInternal_);
    k_.Init("k", k1_);

    StencilCompiler::Build(
        *rhsStencil_,
        "RightHandSide",
        calculationDomain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Input fields
            Param<qrhs, cIn>(qrhs_),
            // Output fields
            Param<k, cInOut>(k_),
            // Scalars
            Param<nu, cScalar>(nu_),
            Param<cx, cScalar>(cx_),
            Param<cy, cScalar>(cy_),
            Param<cz, cScalar>(cz_),
            Param<dx, cScalar>(dx_),
            Param<dx2, cScalar>(dx2_),
            Param<dt, cScalar>(dtparam_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<RHSStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        *eulerStencil_,
        "Euler",
        calculationDomain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Input fields
            Param<qmain, cIn>(qMain_),
            Param<k, cIn>(k_),
            // Output fields
            Param<qtemp, cInOut>(qInternal_),
            // Scalars
            Param<dt, cScalar>(dtparam_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<EulerStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    StencilCompiler::Build(
        *rkStencil_,
        "RungeKutta",
        calculationDomain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Input fields
            Param<k1, cIn>(k1_),
            Param<k2, cIn>(k2_),
            Param<k3, cIn>(k3_),
            Param<k4, cIn>(k4_),
            // Output fields
            Param<qmain, cInOut>(qMain_),
            // Scalars
            Param<dt, cScalar>(dt_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<RKStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    he1_.registerField(qMain_);
    he2_.registerField(qInternal_);
    he3_.registerField(k3_);
    he4_.registerField(k4_);
}

Convection::~Convection()
{
    delete rhsStencil_;
    delete eulerStencil_;
    delete rkStencil_;
}

void Convection::DoTimeStep()
{
    /* First RK timestep */
    qrhs_.set_dataField(qMain_);
    k_.set_dataField(k1_);
    he1_.exchange();
    rhsStencil_->Apply();

    /* Second RK timestep */
    dtparam_ = dthalf_;
    eulerStencil_->Apply();

    k_.set_dataField(k2_);
    qrhs_.set_dataField(qInternal_);
    he2_.exchange();
    rhsStencil_->Apply();
    
    /* Third RK timestep */
    dtparam_ = dthalf_;
    eulerStencil_->Apply();

    k_.set_dataField(k3_);
    qrhs_.set_dataField(qInternal_);
    he2_.exchange();
    rhsStencil_->Apply();
    
    /* Fourth RK timestep */
    dtparam_ = dt_;
    eulerStencil_->Apply();

    k_.set_dataField(k4_);
    qrhs_.set_dataField(qInternal_);
    he2_.exchange();
    rhsStencil_->Apply();

    /* Final RK stage: put things together */
    rkStencil_->Apply();
}

