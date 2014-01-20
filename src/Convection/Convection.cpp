#include <iostream>
#include <cmath>
#include "MatFile.h"
#include "SharedInfrastructure.h"
#include "StencilFramework.h"

#include "Convection.h"
#include "ConvectionFunctions.h"

#include "Periodicity.h"

#ifdef __CUDA_BACKEND__
    typedef BlockSize<32, 4> BSize;
#else
    typedef BlockSize<8, 8> BSize;
#endif


enum
{
    nu, cx, cy, cz, dx, dx2, dt,
    k, k1, k2, k3, k4,
    qmainin, qmainout, qrhsin, qeulerout
};


namespace ConvectionStages
{

    template<typename TEnv>
    struct RHS2Stage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, qrhsin)
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
            T tensadv = ctx[Call<Advection2>::With(qrhsin::Center(), dx::Center(), cx::Center(), cy::Center(), cz::Center())];

            // Compute laplacian
            T tenslapl = ctx[nu::Center()] * ctx[Call<Laplace2>::With(qrhsin::Center(), dx2::Center())];

            ctx[k::Center()] = tensadv + tenslapl;
        }
    };

    template<typename TEnv>
    struct RHS4Stage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, qrhsin)
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
            T tensadv = ctx[Call<Advection4>::With(qrhsin::Center(), dx::Center(), cx::Center(), cy::Center(), cz::Center())];

            // Compute laplacian
            T tenslapl = ctx[nu::Center()] * ctx[Call<Laplace4>::With(qrhsin::Center(), dx2::Center())];

            ctx[k::Center()] = tensadv + tenslapl;
        }
    };

    template<typename TEnv>
    struct EulerStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, qmainin)
        STAGE_PARAMETER(FullDomain, qeulerout)
        STAGE_PARAMETER(FullDomain, k)
        STAGE_PARAMETER(FullDomain, dt)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            const T dt_ = ctx[dt::Center()];
            ctx[qeulerout::Center()] = ctx[qmainin::Center()]
                                 + dt_*ctx[k::Center()];
        }
    };

    template<typename TEnv>
    struct RKStage
    {
        STENCIL_STAGE(TEnv)

        STAGE_PARAMETER(FullDomain, qmainin)
        STAGE_PARAMETER(FullDomain, qmainout)
        STAGE_PARAMETER(FullDomain, k1)
        STAGE_PARAMETER(FullDomain, k2)
        STAGE_PARAMETER(FullDomain, k3)
        STAGE_PARAMETER(FullDomain, k4)
        STAGE_PARAMETER(FullDomain, dt)

        __ACC__
        static void Do(Context ctx, FullDomain)
        {
            ctx[qmainout::Center()] = ctx[qmainin::Center()] + ctx[dt::Center()]/6. *
                (
                         ctx[k1::Center()]
                   + 2.* ctx[k2::Center()]
                   + 2.* ctx[k3::Center()]
                   +     ctx[k4::Center()]
                );
        }
    };

}

Convection::Convection(int isize, int jsize, int ksize, Real dx,
                       Real nucoeff, Real cxcoeff, Real cycoeff, Real czcoeff)
    : nu_(nucoeff), cx_(cxcoeff), cy_(cycoeff), cz_(czcoeff)
    , dx_(dx), dx2_(dx*dx)
{
    IJKSize domain;
    domain.Init(isize, jsize, ksize);
    InitializeStencils(domain);
}

Convection::~Convection()
{
    delete rhs2Stencil_;
    delete rhs4Stencil_;
    delete eulerStencil_;
    delete rkStencil_;
}

void Convection::InitializeStencils(const IJKSize& domain)
{
    // Initialize temporary fields
    KBoundary kboundary;
    kboundary.Init(-convectionBoundaryLines, convectionBoundaryLines);

    qInternal_.Init("InternalField", domain, kboundary);
    k1_.Init("k1", domain, kboundary);
    k2_.Init("k2", domain, kboundary);
    k3_.Init("k3", domain, kboundary);
    k4_.Init("k4", domain, kboundary);

    // Initialize joker data fields
    //qInput_.Init("InputField", qInternal_);
    //qOutput_.Init("OutputField", qInternal_);
    qMainIn_.Init("InputMainField", qInternal_);
    qMainOut_.Init("OutputMainField", qInternal_);
    qRHSIn_.Init("RHSField", qInternal_);
    qEulerOut_.Init("EulerField", qInternal_);
    k_.Init("KField", k1_);

    // Allocate stencils
    rhs2Stencil_ = new Stencil;
    rhs4Stencil_ = new Stencil;
    eulerStencil_ = new Stencil;
    rkStencil_ = new Stencil;

    using namespace ConvectionStages;

    // Initialize RHS2 stencil
    StencilCompiler::Build(
        *rhs2Stencil_,
        "RightHandSide",
        domain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Input fields
            Param<qrhsin, cIn>(qRHSIn_),
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
                    StencilStage<RHS2Stage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    // Initialize RHS4 stencil
    StencilCompiler::Build(
        *rhs4Stencil_,
        "RightHandSide",
        domain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Input fields
            Param<qrhsin, cIn>(qRHSIn_),
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
                    StencilStage<RHS4Stage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    // Initialize Euler stencil
    StencilCompiler::Build(
        *eulerStencil_,
        "Euler",
        domain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // Input fields
            Param<qmainin, cIn>(qMainIn_),
            Param<k, cIn>(k_),
            // Output fields
            Param<qeulerout, cInOut>(qEulerOut_),
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

    // Initialize RK stencil
    StencilCompiler::Build(
        *rkStencil_,
        "RungeKutta",
        domain,
        StencilConfiguration<Real, BSize>(),
        pack_parameters(
            // K fields
            Param<k1, cIn>(k1_),
            Param<k2, cIn>(k2_),
            Param<k3, cIn>(k3_),
            Param<k4, cIn>(k4_),
            // Solution fields
            Param<qmainin, cIn>(qMainIn_),
            Param<qmainout, cInOut>(qMainOut_),
            // Scalars
            Param<dt, cScalar>(dtparam_)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<RKStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );
}

double Convection::DoRK4(ConvectionField& inputField, ConvectionField& outputField,
                       double dt, unsigned timesteps)
{
    // First timestep, read from input, write into output
    if (timesteps > 0)
        DoRK4Timestep(inputField, outputField, dt);

    // All other timesteps
    for (--timesteps; timesteps > 0; --timesteps)
        DoRK4Timestep(outputField, outputField, dt);

    return 0.;
}

double Convection::DoEuler(ConvectionField& inputField, ConvectionField& outputField,
                       double dt, unsigned timesteps)
{
    // First timestep, read from input, write into output
    if (timesteps > 0)
        DoEulerTimestep(inputField, outputField, dt);

    // All other timesteps
    for (--timesteps; timesteps > 0; --timesteps)
        DoEulerTimestep(outputField, outputField, dt);

    return 0.;
}

void Convection::DoRK4Timestep(ConvectionField& inputField, ConvectionField& outputField, double dt)
{
    double dthalf = dt * .5;
    qEulerOut_.set_dataField(qInternal_);
    qMainIn_.set_dataField(inputField);
    qMainOut_.set_dataField(outputField);

    qRHSIn_.set_dataField(inputField);
    Periodicity<ConvectionField> periodicity(qRHSIn_);

#ifndef NDEBUG
    static int timestep = 0;
    std::ostringstream fname;
    fname << "trallallero_" << (++timestep) << ".mat";
    MatFile mat(fname.str().c_str());
    mat.addField("Input", qRHSIn_.dataField());
#endif

    /* First RK stage */
    // k1 = RHS(inputField)
    k_.set_dataField(k1_);
    periodicity.Apply();
#ifndef NDEBUG
    mat.addField("InputPeriodic", qRHSIn_.dataField());
#endif
    rhs4Stencil_->Apply();
#ifndef NDEBUG
    mat.addField("k1", k_.dataField());
#endif

    /* Second RK timestep */
    // Euler: qInternal = qInput + dthalf * k1
    dtparam_ = dthalf;
    eulerStencil_->Apply();
#ifndef NDEBUG
    mat.addField("qTmp1", qInternal_);
#endif

    // k2 = RHS(qInternal)
    k_.set_dataField(k2_);
    qRHSIn_.set_dataField(qInternal_);
    periodicity.Apply();
#ifndef NDEBUG
    mat.addField("qTmp1Periodic", qRHSIn_.dataField());
#endif
    rhs4Stencil_->Apply();
#ifndef NDEBUG
    mat.addField("k2", k_.dataField());
#endif

    /* Third RK timestep */
    // Euler: qInternal = qInput + dthalf * k2
    dtparam_ = dthalf;
    eulerStencil_->Apply();
#ifndef NDEBUG
    mat.addField("qTmp2", qInternal_);
#endif

    // k3 = RHS(qInternal)
    k_.set_dataField(k3_);
    periodicity.Apply();
#ifndef NDEBUG
    mat.addField("qTmp2Periodic", qRHSIn_.dataField());
#endif
    rhs4Stencil_->Apply();
#ifndef NDEBUG
    mat.addField("k3", k_.dataField());
#endif

    /* Fourth RK timestep */
    // Euler: qInternal = qInput + dt * k3
    dtparam_ = dt;
    eulerStencil_->Apply();
#ifndef NDEBUG
    mat.addField("qTmp3", qInternal_);
#endif

    // k4 = RHS(qInternal)
    k_.set_dataField(k4_);
    periodicity.Apply();
#ifndef NDEBUG
    mat.addField("qTmp3Periodic", qRHSIn_.dataField());
#endif
    rhs4Stencil_->Apply();
#ifndef NDEBUG
    mat.addField("k4", k_.dataField());
#endif

    /* Final RK stage: put things together */
    // outputField = inputField + dt/6 * (k1 + 2 k2 + 2 k3 + k4)
    rkStencil_->Apply();
}

void Convection::DoEulerTimestep(ConvectionField& inputField, ConvectionField& outputField, double dt)
{
    // Compute RHS into k1
    qEulerOut_.set_dataField(outputField);
    qMainIn_.set_dataField(inputField);
    qMainOut_.set_dataField(outputField);

    Periodicity<ConvectionField> periodicity(qRHSIn_);

    // k1 = RHS(inputField)
    qRHSIn_.set_dataField(inputField);
    k_.set_dataField(k1_);
    periodicity.Apply();
    cudaDeviceSynchronize();
    rhs2Stencil_->Apply();

    // Euler: outputField = inputField + dt * k1
    dtparam_ = dt;
    eulerStencil_->Apply();
}

