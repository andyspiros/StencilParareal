#ifndef HEAT_H_
#define HEAT_H_
#include "SharedInfrastructure.h"
#include "StencilFramework.h"
#include "KField.h"
#include "QField.h"
#include "HaloExchange3D.h"

class Heat
{
public:
    Heat(IJKRealField& q, Real coeff, Real dx, Real dt, MPI_Comm comm);

    void DoTimeStep();

    void DoTest();

private:
    // Scalars
    Real coeff_, dx2_;
    Real dt_, dthalf_, dtparam_;

    // Data fields
    QField<IJKRealField> qs_;
    KField<IJKRealField> ks_;

    // Stencils
    Stencil laplaceStencil_, eulerStencil_, rkStencil_;

    // HaloExchange
    HaloExchange3D<IJKRealField> he_;
};

#endif // HEAT_H_
