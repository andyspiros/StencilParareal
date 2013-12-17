#include <iostream>
#include <cmath>
#include "Heat.h"

void initializeQ(IJKRealField& q)
{
    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();
    const double pi = 3.1415926536;

    const double dx = 1./(iSize-1);
    const double dy = 1./(jSize-1);
    const double dz = 1./(kSize-1);

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                q(i, j, k) = sin(pi*i*dx)*sin(pi*j*dy)*sin(pi*k*dz);
            }
}

int main(int argc, char **argv)
{
    // Read command-line arguemnts
    int size = 64;
    int timesteps = 100;
    if (argc > 1)
    {
        if (argc > 2)
        {
            size = std::atoi(argv[1]);
        }
        if (argc > 3)
        {
            timesteps = std::atoi(argv[2]);
        }
        if (argc > 4)
        {
            std::cout << "Usage: heat isize jsize ksize\n";
            std::cout << "Aborting\n";
            return -1;
        }
    }

    // Create q
    IJKSize calculationDomain;
    calculationDomain.Init(size, size, size);
    KBoundary kboundary;
    kboundary.Init(-3, 3);
    IJKRealField q;
    q.Init("q", calculationDomain, kboundary);

    // Initialize content of q
    initializeQ(q);
    
    // Initialize stencil
    Heat heat(q, 1., 1., 1., 1., 0.1);

    for (int t = 0; t < 10; ++t) {
        std::cout << "Doing timestep " << t << "\n";
        heat.DoTimeStep();
    }

}
