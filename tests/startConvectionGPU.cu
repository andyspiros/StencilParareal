#define __CUDA_BACKEND__
#include "Convection.h"
#include "mpi.h"
#include <iostream>
#include <iomanip>

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    const int s = 64;
    Convection convection(s, s, s, 1./s, 0.1, 1., 1., 1.);
    ConvectionField qin, qout;
    IJKSize domain; domain.Init(s, s, s);
    KBoundary kboundary; kboundary.Init(-2, 2);
    qin.Init("qin", domain, kboundary);
    qout.Init("qout", domain, kboundary);

    double e;

    // RK
    for (int i = 1; i <= 1024; i*=2)
    {
        cudaDeviceSynchronize();
        e = MPI_Wtime();
        convection.DoRK4(qin, qout, 0.001, i);
        e = MPI_Wtime() - e;
        std::cout << std::setw(4) << i;
        std::cout << "      " << std::fixed << std::setprecision(7) << e << "\n";
    }

    // Euler
    for (int i = 1; i <= 1024; i*=2)
    {
        cudaDeviceSynchronize();
        double e = MPI_Wtime();
        convection.DoEuler(qin, qout, 0.001, i);
        e = MPI_Wtime() - e;
        std::cout << std::setw(4) << i;
        std::cout << "      " << std::fixed << std::setprecision(7) << e << "\n";
    }

    cudaDeviceReset();
    MPI_Finalize();
}

