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
    : qs_(data), ks_(data)
    , nu_(nucoeff), cx_(cxcoeff), cy_(cycoeff), cz_(czcoeff)
    , dx_(stepsize), dx2_(stepsize*stepsize)
    , dt_(timestepsize), dthalf_(.5 * timestepsize)
    , he_(true, true, true, comm)
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
    StencilCompiler::Build(
        *rhsStencil_,
        "RightHandSide",
        calculationDomain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Input fields
            Param<qrhs, cIn>(qs_.qrhs()),
            // Output fields
            Param<k, cInOut>(ks_.k()),
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
                    StencilStage<RKStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    he_.registerField(qs_.qrhs());
}

Convection::~Convection()
{
    delete rhsStencil_;
    delete eulerStencil_;
    delete rkStencil_;
}

void Convection::DoTimeStep()
{
    double erhs, eeuler, erk, ecomm;
    DoTimeStep(erhs, eeuler, erk, ecomm);
}

void Convection::DoTimeStep(double& erhs, double& eeuler, double& erk, double& ecomm)
{
    double e;
    erhs = 0., eeuler = 0., erk = 0., ecomm = 0.;
    
    /* First RK timestep */
    ks_.set(1);
    qs_.setRHSMain();
    e = MPI_Wtime(); he_.exchange(); ecomm += MPI_Wtime() - e;
    e = MPI_Wtime(); rhsStencil_->Apply(); erhs += MPI_Wtime()-e;
    // Now: k1 = laplace(qmain)

    /* Second RK timestep */
    qs_.restore();
    dtparam_ = dthalf_;
    e = MPI_Wtime(); eulerStencil_->Apply(); eeuler += MPI_Wtime()-e;
    // Now: qtemp = qmain + dthalf_*k1

    ks_.set(2);
    qs_.setRHSTemp();
    e = MPI_Wtime(); he_.exchange(); ecomm += MPI_Wtime()-e;
    e = MPI_Wtime(); rhsStencil_->Apply(); erhs += MPI_Wtime()-e;
    // Now: k2 = laplace(qmain + dthalf_*k1)
    
    /* Third RK timestep */
    qs_.restore();
    e = MPI_Wtime(); eulerStencil_->Apply(); eeuler += MPI_Wtime()-e;
    // Now: qtemp = qmain + dthalf_*k2

    ks_.set(3);
    qs_.setRHSTemp();
    e = MPI_Wtime(); he_.exchange(); ecomm += MPI_Wtime()-e;
    e = MPI_Wtime(); rhsStencil_->Apply(); erhs += MPI_Wtime()-e;
    // Now: k3 = laplace(qmain + dthalf_*k2)
    
    /* Fourth RK timestep */
    qs_.restore();
    dtparam_ = dt_;
    e = MPI_Wtime(); eulerStencil_->Apply(); eeuler += MPI_Wtime()-e;
    // Now: qtemp = qmain + dt_*k3

    ks_.set(4);
    qs_.setRHSTemp();
    e = MPI_Wtime(); he_.exchange(); ecomm += MPI_Wtime()-e;
    e = MPI_Wtime(); rhsStencil_->Apply(); erhs += MPI_Wtime()-e;
    // Now: k4 = laplace(qmain + dt_*k3)

    /* Final RK stage: put things together */
    ks_.restore();
    qs_.restore();
    e = MPI_Wtime(); rkStencil_->Apply(); erk += MPI_Wtime()-e;
}

