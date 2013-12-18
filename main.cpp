#include <iostream>
#include <cmath>
#include "Heat.h"
#include "MatFile.h"

const double pi = 3.14159265358979;

void initializeQ(IJKRealField& q, Real dx, Real dy, Real dz, Real t = 0.)
{
    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();

    double x, y, z;

    for (int i = -1; i <= iSize; ++i)
        for (int j = -1; j <= jSize; ++j)
            for (int k = -1; k <= kSize; ++k)
            {
                x = (i+1)*dx;
                y = (j+1)*dy;
                z = (k+1)*dz;
                q(i, j, k) = sin(pi*x)*sin(pi*y)*sin(pi*z)
                    * exp(-3.*pi*pi*t);
                //q(i, j, k) = sin(x);
            }
}

double computeError(const IJKRealField& q, Real t, Real dx, Real dy, Real dz, IJKRealField* errfield=0)
{
    double error = 0.;

    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();

    Real x, y, z;
    Real errinf = 0.;
    Real exact, e;
    Real exactinf;

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                // Coordinates
                x = (i+1)*dx;
                y = (j+1)*dy;
                z = (k+1)*dz;

                // Exact solution
                exact = sin(pi*x)*sin(pi*y)*sin(pi*z);
                exact *= std::exp(-3.*pi*pi*t);
                //exact = x*x;
                exactinf = std::max(std::abs(exact), exactinf);

                // Error
                e = q(i,j,k) - exact;
                errinf = std::max(std::abs(e), errinf);

                // Error field
                if (errfield)
                    (*errfield)(i, j, k) = e;
            }
    std::cout << "Infinity norm of exact solution: " << exactinf << "\n";
    std::cout << "Infinity norm of error: " << errinf << "\n";
    return errinf / exactinf;
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
        //if (argc > 2)
        //{
        //    timesteps = std::atoi(argv[2]);
        //}
        //if (argc > 3)
        if (argc > 2)
        {
            std::cout << "Usage: heat size timesteps\n";
            std::cout << "Aborting\n";
            return -1;
        }
    }

    // Keep CFL constant at 0.1
    const double endtime = 0.05;
    double cfl = 0.1;
    timesteps = endtime * (size+1)*(size+1) / cfl + 0.5;
    timesteps = 000;

    if (size == 0) {
        std::cin >> size;
        timesteps = 0;
    }

    //timesteps = (timesteps + 19)/20*20;
    dt = 0.05 / timesteps;
    dx = dy = dz = 1. / (size + 1);

    cfl = 1. * dt / (dx*dx);

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
    IJKRealField q, errfield, exactfield;
    q.Init("q", calculationDomain, kboundary);
    errfield.Init("error", calculationDomain, kboundary);
    exactfield.Init("exact", calculationDomain, kboundary);

    // Initialize content of q
    initializeQ(q, dx, dy ,dz);
 
    // Initialize stencil
    Heat heat(q, 1., dx, dt);

    if (timesteps == 0)
    {
        // Do test, exit
        heat.DoTest();
        std::cout << "Exiting\n";
        return 0;
    }

    // Initialize MAT file
    // MatFile mat("heat.mat");
    // MatFile errfile("error.mat");
    // mat.startCell("q", timesteps/20+1);
    // errfile.startCell("q", timesteps/20+1);
    // mat.addField(q, 0);

    double error;
    // error = computeError(q, 0., dx, dy ,dz, &errfield);
    // errfile.addField(errfield, 0);

    for (int t = 1; t <= timesteps; ++t) {
        heat.DoTimeStep();
        //error = computeError(q, dt * t, dx, dy ,dz, &errfield);

         if (t % 20 == 0) {
             std::cout << "Timestep " << t << " done" << std::endl;
        //     mat.addField(q, t/20);
        //     errfile.addField(errfield, timesteps/20);
         }
    }
    
    std::cout << "Error at end: " << computeError(q, endtime, dx, dy, dz, &errfield) << "\n";
    initializeQ(exactfield, dx, dy, dz, timesteps*dt);

    MatFile matfile("result.mat");
    matfile.addField(q, timesteps);
    matfile.addField(errfield, timesteps);
    matfile.addField(exactfield, timesteps);

    // mat.endCell();
    // errfile.endCell();
}
