#ifndef CONVECTION_H_
#define CONVECTION_H_

#include "SharedInfrastructure.h"
#include "KField.h"
#include "QField.h"
#include "HaloExchange3D.h"
#include "MatFile.h"

class Stencil;

// Define data field
static const int convectionBoundaryLines = 2;
typedef DataFieldIJBoundary<-convectionBoundaryLines, convectionBoundaryLines, -convectionBoundaryLines, convectionBoundaryLines> ConvectionIJBoundary;

#ifdef CUDA_BACKEND
typedef DataFieldCUDA<Real, DataFieldStorageFormat<ConvectionIJBoundary, StorageOrder::KJI, CUDARealAlignment> > ConvectionRealField;
#else
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<ConvectionIJBoundary, StorageOrder::KJI, OpenMPAlignment> > ConvectionRealField;
#endif

class Convection
{
public:
    Convection(ConvectionRealField& q, Real nu, Real cx, Real cy, Real cz, Real dx, Real dt, MPI_Comm comm);
    ~Convection();

    void DoTimeStep();
    void DoTimeStep(double&, double&, double&, double&);

private:
    // Scalars
    Real nu_, cx_, cy_, cz_, dx_, dx2_;
    Real dt_, dthalf_, dtparam_;

    ConvectionRealField& q_;
    ConvectionRealField qinternal_;

    ConvectionRealField k1_, k2_, k3_, k4_;

    // Stencils
    Stencil* stencil1_;
    Stencil* stencil2_;
    Stencil* stencil3_;
    Stencil* stencil4_;

    // HaloExchange
    HaloExchange3D<ConvectionRealField> heQ_, heQinternal_;

    MatFile matq, matk1, matk2, matk3, matk4;
};

#endif // CONVECTION_H_

