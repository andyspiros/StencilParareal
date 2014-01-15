#ifndef CONVECTION_H_
#define CONVECTION_H_

#include "SharedInfrastructure.h"
#include "KField.h"
#include "QField.h"
#include "HaloExchange3D.h"

class Stencil;

// Define data field
static const int convectionBoundaryLines = 2;
typedef DataFieldIJBoundary<-convectionBoundaryLines, convectionBoundaryLines, -convectionBoundaryLines, convectionBoundaryLines> ConvectionIJBoundary;

#ifdef CUDA_BACKEND
typedef DataFieldCUDA<Real, DataFieldStorageFormat<ConvectionIJBoundary, StorageOrder::KJI, CUDARealAlignment> > ConvectionField;
#else
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<ConvectionIJBoundary, StorageOrder::KJI, OpenMPAlignment> > ConvectionField;
#endif

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
    HaloExchange3D<ConvectionField> he1_, he2_;
};

#endif // CONVECTION_H_

