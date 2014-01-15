#include <iostream>
#include "RuntimeConfiguration.h"
#include "ConvectionSolution.h"

int main(int argc, char **argv)
{
    RuntimeConfiguration conf(argc, argv);

    std::cout << "Running with:\n"
        << " - spatial discretization step: " << conf.dx() << "\n"
        << " - timestep size: " << conf.dt() << "\n"
        << " - endtime: " << conf.endtime() << "\n"
        << " - heat coefficient: " << conf.nu() << "\n"
        << " - advection velocity in x: " << conf.cx() << "\n"
        << " - advection velocity in y: " << conf.cy() << "\n"
        << " - advection velocity in z: " << conf.cz() << "\n"
        << " - CFL fine: " << conf.cflFine() << "\n"
        << " - CFL coarse: " << conf.cflCoarse() << "\n"
        << " - time slices: " << conf.timeSlices() << "\n"
        << "\n";

    // Calculation domain and boundaries
    IJKSize domain; domain.Init(conf.gridSize(), conf.gridSize(), conf.gridSize()); 
    KBoundary kboundary; kboundary.Init(-convectionBoundaryLines, convectionBoundaryLines);
}

