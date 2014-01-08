#ifndef CONVECTION_H_
#define CONVECTION_H_

#include "SharedInfrastructure.h"
#include "KField.h"
#include "QField.h"
#include "HaloExchange3D.h"

class Stencil;

class Convection
{
public:
    Convection(IJKRealField& q, Real nu, Real cx, Real cy, Real cz, Real dx, Real dt, MPI_Comm comm);
    ~Convection();

    void DoTimeStep();
    void DoTimeStep(double&, double&, double&, double&);

private:
    // Scalars
    Real nu_, cx_, cy_, cz_, dx_, dx2_;
    Real dt_, dthalf_, dtparam_;

    // Data fields
    QField<IJKRealField> qs_;
    KField<IJKRealField> ks_;

    // Stencils
    Stencil *rhsStencil_, *eulerStencil_, *rkStencil_;

    // HaloExchange
    HaloExchange3D<IJKRealField> he_;
};

#endif // CONVECTION_H_
