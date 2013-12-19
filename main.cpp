#include <iostream>
#include <cmath>
#include "Heat.h"
#include "MatFile.h"
#include "HaloExchange3D.h"

#include "GCL.h"
#include "utils/layout_map.h"
#include "utils/boollist.h"
#include "L2/include/halo_exchange.h"
#include "L2/include/descriptors.h"
#include "L3/include/proc_grids_3D.h"
#include "L3/include/Halo_Exchange_3D.h"

const double pi = 3.14159265358979;

void initializeQ(IJKRealField& q, Real dx, Real dy, Real dz, Real nu = 1., Real t = 0.)
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
                x = (2*i+1)*dx / 2.;
                y = (2*j+1)*dy / 2.;
                z = (2*k+1)*dz / 2.;
                q(i, j, k) = sin(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*z)
                    * exp(-12.*pi*pi*nu*t);
            }
}

double computeError(const IJKRealField& q, Real t, Real dx, Real dy, Real dz, Real nu = 1., IJKRealField* errfield=0)
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
                x = (2*i+1)*dx / 2.;
                y = (2*j+1)*dy / 2.;
                z = (2*k+1)*dz / 2.;

                // Exact solution
                exact = sin(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*z);
                exact *= std::exp(-12.*pi*pi*nu*t);
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
    typedef GCL::MPI_3D_process_grid_t<GCL::gcl_utils::boollist<3> > GridType;

    // Initialize GCL
    GCL::GCL_Init();

    // Initialize grid
    int commsize;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int dims[3] = {0, 0, 1};
    int periods[3] = {1, 1, 1};
    MPI_Dims_create(GCL::PROCS, 3, dims);
    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm);

    // Read command-line arguemnts
    int size = 64;
    double dt, dx, dy, dz;
    int timesteps = 100;
    bool forceTimesteps = false;
    if (argc > 1)
    {
        if (argc > 1)
        {
            size = std::atoi(argv[1]);
        }
        if (argc > 2)
        {
            timesteps = std::atoi(argv[2]);
            forceTimesteps = true;
        }
        if (argc > 3)
        {
            std::cout << "Usage: heat [size [timesteps]]\n";
            std::cout << "Aborting\n";
            return -1;
        }
    }

    // Keep CFL constant at 0.1
    const double endtime = 0.05;
    double cfl = 0.2;
    if (!forceTimesteps)
        timesteps = endtime * (size+1)*(size+1) / cfl + 0.5;

    if (size == 0) {
        std::cin >> size;
        timesteps = 0;
    }

    //timesteps = (timesteps + 19)/20*20;
    dt = 0.05 / timesteps;
    dx = dy = dz = 1. / (size);
    double nu = dx;

    cfl = 1. * dt / (dx*dx);

    // Output configuration
    std::cout << "Running with:\n"
        << " - domain size: " << size << "x" << size << "x" << size << "\n"
        << " - spatial discretization step: " << dx << "\n"
        << " - timestep size: " << dt << "\n"
        << " - timesteps: " << timesteps << "\n"
        << " - final time: " << timesteps*dt << "\n"
        << " - heat coefficient: " << nu << "\n"
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
    initializeQ(q, dx, dy, dz, nu);
 
    // Initialize stencil
    Heat heat(q, nu, dx, dt, comm);

    if (timesteps == 0)
    {
        // Do test, exit
        heat.DoTest();
        std::cout << "Exiting\n";
        return 0;
    }

    // Initialize MAT file
    MatFile mat("heat.mat");
    // MatFile errfile("error.mat");
    mat.startCell("q", timesteps/20+1);
    // errfile.startCell("q", timesteps/20+1);
    mat.addField(q, 0);

    double error;
    // error = computeError(q, 0., dx, dy ,dz, &errfield);
    // errfile.addField(errfield, 0);

    for (int t = 1; t <= timesteps; ++t) {
        heat.DoTimeStep();
        //error = computeError(q, dt * t, dx, dy ,dz, &errfield);

         //if (t % 20 == 0) {
         //    std::cout << "Timestep " << t << " done" << std::endl;
         //    mat.addField(q, t/20);
         //    errfile.addField(errfield, timesteps/20);
         //}
    }
    
    std::cout << "Error at end: " << computeError(q, endtime, dx, dy, dz, nu, &errfield) << "\n";
    initializeQ(exactfield, dx, dy, dz, nu, timesteps*dt);

    MatFile matfile("result.mat");
    matfile.addField(q, timesteps);
    matfile.addField(errfield, timesteps);
    matfile.addField(exactfield, timesteps);

    mat.endCell();
    // errfile.endCell();

    // Finalize GCL
    GCL::GCL_Finalize();
}
