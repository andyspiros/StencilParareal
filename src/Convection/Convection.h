#ifndef CONVECTION_H_
#define CONVECTION_H_

#include "SharedInfrastructure.h"
#include "KField.h"
#include "QField.h"
#include "HaloExchange3D.h"

class Stencil;

typedef IJKRealField ConvectionField;

class Convection
{
public:
    Convection(ConvectionField& q, Real nu, Real cx, Real cy, Real cz, Real dx, Real dt, MPI_Comm comm);
    ~Convection();

    void DoTimeStep();
    void DoTimeStep(double&, double&, double&, double&);

private:
    // Scalars
    Real nu_, cx_, cy_, cz_, dx_, dx2_;
    Real dt_, dthalf_, dtparam_;

    // Data fields
    ConvectionField& qMain_;
    ConvectionField qInternal_;
    ConvectionField k1_, k2_, k3_, k4_;

    JokerDataField<ConvectionField> qrhs_, k_;

    // Stencils
    Stencil *rhsStencil_, *eulerStencil_, *rkStencil_;

    // HaloExchange
    HaloExchange3D<ConvectionField> he1_, he2_, he3_, he4_;
};

#endif // CONVECTION_H_
