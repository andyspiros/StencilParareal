#include "SharedInfrastructure.h"
#include "StencilFramework.h"

class Heat
{
public:
    Heat(IJKRealField& q, Real coeff, Real di, Real dj, Real dk, Real dt);

    void DoTimeStep();

private:
    // Scalars
    Real coeff_, di2_, dj2_, dk2_;
    Real dt_, dthalf_;

    // Data fields
    IJKRealField& q_;
    IJKRealField qtempin_, qtempout_, k1_, k2_, k3_, k4_;

    // Stencils
    Stencil stencil1_;
    Stencil stencil2_;
    Stencil stencil3_;
    Stencil stencil4_;
    Stencil stencilRK_;
};
