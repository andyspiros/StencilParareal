#include <iostream>
#include <cmath>
#include "Heat.h"
#include "MatFile.h"

const double pi = 3.14159265358979;

void initializeQ(IJKRealField& q, Real dx, Real dy, Real dz)
{
    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();

    double x, y, z;

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                x = (i+1)*dx;
                y = (j+1)*dy;
                z = (k+1)*dz;
                q(i, j, k) = sin(pi*x)*sin(pi*y)*sin(pi*z);
            }
}

double computeError(const IJKRealField& q, double t, Real dx, Real dy, Real dz)
{
    double error = 0.;

    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();

    double x, y, z;
    double errinf = 0.;
    double exact, e;

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                x = (i+1)*dx;
                y = (j+1)*dy;
                z = (k+1)*dz;

                exact = sin(pi*x)*sin(pi*y)*sin(pi*z);
                exact *= std::exp(-3.*pi*pi*t);

                e = q(i,j,k) - exact;
                errinf = std::max(errinf, std::abs(exact));
                error = std::max(error, std::abs(e));
            }
    return error / errinf;
}

int main(int argc, char **argv)
{
    // Read command-line arguemnts
    int size = 64;
    double dt, dx, dy, dz;
    int timesteps = 100;
    if (argc > 1)
    {
        if (argc > 1)
        {
            size = std::atoi(argv[1]);
        }
        if (argc > 2)
        {
            timesteps = std::atoi(argv[2]);
        }
        if (argc > 3)
        {
            std::cout << "Usage: heat size timesteps\n";
            std::cout << "Aborting\n";
            return -1;
        }
    }
    timesteps = (timesteps + 19)/20*20;
    dt = 0.05 / timesteps;
    dx = dy = dz = 1. / (size + 1);

    double cfl = 1. * dt / (dx*dx);

    // Output
    std::cout << "Running with:\n"
        << " - domain size: " << size << "x" << size << "x" << size << "\n"
        << " - spatial discretization step: " << dx << "\n"
        << " - timestep size: " << dt << "\n"
        << " - timesteps: " << timesteps << "\n"
        << " - final time: " << timesteps*dt << "\n"
        << " - CFL: " << cfl << "\n\n";

    // Create q
    IJKSize calculationDomain;
    calculationDomain.Init(size, size, size);
    KBoundary kboundary;
    kboundary.Init(-3, 3);
    IJKRealField q, errfield;
    q.Init("q", calculationDomain, kboundary);
    errfield.Init("error", calculationDomain, kboundary);

    // Initialize content of q
    initializeQ(q, dx, dy ,dz);
 
    // Initialize stencil
    Heat heat(q, 1., dx, dy, dz, dt);

    // Initialize MAT file
    MatFile mat("heat.mat");
    mat.startCell("q", timesteps/20+1);
    mat.addField(q, 0);

    double error;

    for (int t = 1; t <= timesteps; ++t) {
        heat.DoTimeStep();
        error = computeError(q, dt * t, dx, dy ,dz);
        if (t % 20 == 0) {
            mat.addField(q, t/20);
        }
    }
    
    std::cout << "Error at end: " << computeError(q, timesteps*dt, dx, dy, dz) << "\n";

    mat.endCell();

}
