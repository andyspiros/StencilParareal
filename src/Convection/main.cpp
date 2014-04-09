#include <iostream>
#include <cmath>
#include <cstdlib>

#include "boost/program_options.hpp"

#include "Convection.h"
#include "ConvectionSolution.h"
#include "MatFile.h"
//#include "HaloExchange3D.h"

//#include "GCL.h"
//#include "utils/layout_map.h"
//#include "utils/boollist.h"
//#include "L2/include/halo_exchange.h"
//#include "L2/include/descriptors.h"
//#include "L3/include/proc_grids_3D.h"
//#include "L3/include/Halo_Exchange_3D.h"

#include "mpi.h"

#ifndef __CUDA_BACKEND__
# include <omp.h>
#endif


namespace po = boost::program_options;

void computeLocalSizes(int globalSize, int nprocs,
        std::vector<int>& localSizes,
        std::vector<double>& xstarts, std::vector<double>& xends)
{
    localSizes.resize(nprocs);
    xstarts.resize(nprocs);
    xends.resize(nprocs);

    const float dx = 1. / globalSize;

    double l = static_cast<double>(globalSize) / nprocs;
    int lp = static_cast<int>(std::ceil(l));
    int lm = static_cast<int>(std::floor(l));

    int pp = globalSize - lm*nprocs;
    int p;

    int start = 0;
    for (p = 0; p < pp; ++p)
    {
        localSizes[p] = lp;
        xstarts[p] = start*dx;
        start += lp;
        xends[p] = start*dx;
    }
    for (; p < nprocs; ++p)
    {
        localSizes[p] = lm;
        xstarts[p] = start*dx;
        start += lm;
        xends[p] = start*dx;
    }
}


struct HeatConfiguration
{
    double nu0, nufreq;
    double cx, cy, cz;
    int gridsize;
    int timesteps;
    double endtime;
    double cfl;
    double dx;
    double dt;
    bool mat;
    bool rk;
};

HeatConfiguration parseCommandLine(int argc, char **argv)
{
    HeatConfiguration conf;
    conf.nu0 = 1.;
    conf.nufreq = 0.;
    conf.cx = 1.;
    conf.cy = 1.;
    conf.cz = 1.;
    conf.gridsize = 32;
    conf.timesteps = 512;
    conf.endtime = 0.05;
    conf.dx = 1. / conf.gridsize;
    conf.dt = conf.endtime / conf.timesteps;
    conf.cfl = conf.dt / (conf.dx*conf.dx);
    conf.mat = false;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produce this help message")
        ("nu0", po::value<double>(), "Initial diffusion coefficient")
        ("nufreq", po::value<double>(), "Frequency of diffusion coefficient")
        ("cx", po::value<double>(), "Advection velocity in x direction")
        ("cy", po::value<double>(), "Advection velocity in y direction")
        ("cz", po::value<double>(), "Advection velocity in z direction")
        ("gridSize", po::value<int>(), "Number of gridpoints along each direction")
        ("timeSteps", po::value<int>(), "Number of timesteps")
        ("endTime", po::value<double>(), "Time to simulate")
        ("cfl", po::value<double>(), "Fixed CFL condition")
        ("order", po::value<int>(), "Order of timestepping (valid values are 1 and 4)")
        ("mat", "Output intermediate states")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        std::exit(0);
    }

    conf.mat = vm.count("mat");

    if (vm.count("order"))
    {
        if (vm["order"].as<int>() == 1)
            conf.rk = false;
        else if (vm["order"].as<int>() == 4)
            conf.rk = true;
        else
        {
            std::cerr << "Order " << vm["order"].as<int>()<< " is not supported\n";
            std::cerr << "Aborting.\n";
            std::exit(-1);
        }
    }
    else
        conf.rk = true;

    if (vm.count("nu0"))
    {
        conf.nu0 = vm["nu0"].as<double>();
    }

    if (vm.count("nufreq"))
    {
        conf.nufreq = vm["nufreq"].as<double>();
    }

    if (vm.count("cx"))
    {
        conf.cx = vm["cx"].as<double>();
    }
    if (vm.count("cy"))
    {
        conf.cy = vm["cy"].as<double>();
    }
    if (vm.count("cz"))
    {
        conf.cz = vm["cz"].as<double>();
    }

    if (vm.count("endTime"))
    {
        conf.endtime = vm["endTime"].as<double>();
        conf.dt = conf.endtime / conf.timesteps;
    }

    if (vm.count("cfl") && vm.count("gridSize") && vm.count("timeSteps"))
    {
        std::cerr << "Error: setting cfl, gridSize and timeSteps not allowed\n";
        std::cerr << "Aborting\n";
        std::exit(1);
    }

    bool gridsizeSet = false;
    bool timestepsSet = false;
    if (vm.count("gridSize"))
    {
        conf.gridsize = vm["gridSize"].as<int>();
        conf.dx = 1. / conf.gridsize;
        gridsizeSet = true;
    }

    if (vm.count("timeSteps"))
    {
        conf.timesteps = vm["timeSteps"].as<int>();
        conf.dt = conf.endtime / conf.timesteps;
        timestepsSet = true;
    }

    if (vm.count("cfl"))
    {
        double cfl = vm["cfl"].as<double>();
        if (gridsizeSet)
        {
            double dt = conf.dx*conf.dx * cfl / (3./2.*conf.nu0);
            conf.timesteps = conf.endtime / dt + .5;
            conf.dt = conf.endtime / conf.timesteps;
            conf.cfl = conf.dt * 3./2.*conf.nu0 / (conf.dx*conf.dx);
        }
        else // Also apply if nothing else is specified
        {
            double dx = std::sqrt(conf.dt * conf.nu0 / cfl);
            conf.gridsize = 1. / dx + .5;
            conf.dx = 1. / conf.gridsize;
            conf.cfl = conf.dt * 3./2.*conf.nu0 / (conf.dx*conf.dx);
        }
    }

    // Compute the CFL in any case
    conf.cfl = conf.dt * 3./2.*conf.nu0 / (conf.dx * conf.dx);

    return conf;
}

int main(int argc, char **argv)
{
    // Parse command-line arguments
    HeatConfiguration conf = parseCommandLine(argc, argv);

    // Initialize GCL
    MPI_Init(&argc, &argv);
    //GCL::GCL_Init();
    //typedef GCL::gcl_utils::boollist<3> CyclicType;
    //typedef GCL::MPI_3D_process_grid_t<CyclicType> GridType;
    const bool isRoot = 1;

    // Initialize grid
    int commsize = 1;
    //MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int dims[3] = {1, 1, 1};
    int periods[3] = {1, 1, 1};
    //MPI_Dims_create(GCL::PROCS, 3, dims);
    //MPI_Comm comm;
    //MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm);

    std::vector<int> localSizesI, localSizesJ, localSizesK;
    std::vector<double> xstarts, xends, ystarts, yends, zstarts, zends;
    computeLocalSizes(conf.gridsize, dims[0], localSizesI, xstarts, xends);
    computeLocalSizes(conf.gridsize, dims[1], localSizesJ, ystarts, yends);
    computeLocalSizes(conf.gridsize, dims[2], localSizesK, zstarts, zends);


    // Get position on grid
    //GridType grid(CyclicType(true, true, true), comm);
    int myPI = 0, myPJ = 0, myPK = 0;
    //grid.coords(myPI, myPJ, myPK);

    // Owned physical space
    const double Dx = 1. / dims[0];
    const double Dy = 1. / dims[1];
    const double Dz = 1. / dims[2];
    const double xstart = xstarts[myPI];
    const double ystart = ystarts[myPJ];
    const double zstart = zstarts[myPK];
    const double xend = xends[myPI];
    const double yend = yends[myPJ];
    const double zend = zends[myPK];

    // Output configuration
    if (isRoot)
    {
        std::cout << "Running with:\n"
            << " - grid size: " << conf.gridsize << "x" << conf.gridsize << "x" << conf.gridsize << "\n"
            << " - spatial discretization step: " << conf.dx << "\n"
            << " - timestep size: " << conf.dt << "\n"
            << " - timesteps: " << conf.timesteps << "\n"
            << " - endtime: " << conf.endtime << "\n"
            << " - initial diffusion coefficient: " << conf.nu0 << "\n"
            << " - frequency of diffusion coefficient: " << conf.nufreq << "\n"
            << " - advection velocity in x: " << conf.cx << "\n"
            << " - advection velocity in y: " << conf.cy << "\n"
            << " - advection velocity in z: " << conf.cz << "\n"
            << " - CFL: " << conf.cfl << "\n"
            << " - solver: " << (conf.rk ? "Runge-Kutta 4" : "explicit Euler") << "\n"
            << "\n";

        std::cout << "Spatial configuration:\n";
        for (int p = 0; p < 1; ++p)
        {
            int coords[3] = {0, 0, 0};
            //MPI_Cart_coords(comm, p, 3, coords);
            std::cout << " - Process " << p << " is "
                << "(" << coords[0] << ", " << coords[1] << ", " << coords[2] << ")\n"
                << "   -- Local grid: " << localSizesI[coords[0]] << "x"
                << localSizesJ[coords[1]] << "x" << localSizesK[coords[2]] << "\n"
                << "   -- Starting coordinates: (" << xstarts[coords[0]] << ", "
                << ystarts[coords[1]] << ", " << zstarts[coords[2]] << ")\n"
                << "   -- Ending coordinates: (" << xends[coords[0]] << ", "
                << yends[coords[1]] << ", " << zends[coords[2]] << ")\n"
                ;
        }
    }

    // Create q
    IJKSize calculationDomain;
    calculationDomain.Init(localSizesI[myPI], localSizesJ[myPJ], localSizesK[myPK]);
    KBoundary kboundary;
    kboundary.Init(-convectionBoundaryLines, convectionBoundaryLines);
    ConvectionField qIn, qOut, errfield, exactfield;
    qIn.Init("q", calculationDomain, kboundary);
    qOut.Init("q", calculationDomain, kboundary);
    errfield.Init("error", calculationDomain, kboundary);
    exactfield.Init("exact", calculationDomain, kboundary);

    // Initialize content of q
    fillQ(qIn, conf.nu0, conf.nufreq, conf.cx, conf.cy, conf.cz, 0., xstart, xend, ystart, yend, zstart, zend);
#ifdef __CUDA_BACKEND__
    qIn.SynchronizeDeviceStorage();
    qOut.SynchronizeDeviceStorage();
#endif
    MatFile initmat("init.mat");
    initmat.addField("qIn", qIn);
    initmat.close();

    // Initialize stencil
    Convection convection(localSizesI[myPI], localSizesJ[myPJ], localSizesK[myPK], conf.dx, conf.nu0, conf.nufreq, conf.cx, conf.cy, conf.cz);

    // Initialize MAT file
    MatFile *mat;
    if (conf.mat)
    {
        std::ostringstream fnameMat;
        fnameMat << "convection_" << myPI << "_" << myPJ << "_" << myPK << ".mat";
        mat = new MatFile(fnameMat.str());
        mat->startCell("q", conf.timesteps);
        mat->addField(qIn, 0);
    }

    // Solve equation
    //double e = MPI_Wtime();
    //double erhsTot = 0., eeulerTot = 0., erkTot = 0., ecommTot = 0.;
    //double erhs, eeuler, erk, ecomm;
    //for (int t = 1; t <= conf.timesteps; ++t) {
    //    if (t == 1)
    //        convection.DoRK4Timestep(qIn, qOut, conf.dt, qInHE);
    //    else
    //        convection.DoRK4Timestep(qOut, qOut, conf.dt, qOutHE);

    //    if (conf.mat)
    //        mat->addField(qIn, t);
    //}
    //e = MPI_Wtime() - e;
    double e = MPI_Wtime();
    if (conf.rk)
        convection.DoRK4(qIn, qOut, 0., conf.dt, conf.timesteps);
    else
        convection.DoEuler(qIn, qOut, 0., conf.dt, conf.timesteps);
    e = MPI_Wtime() - e;

    // Close MAT file
    if (conf.mat)
    {
        mat->endCell();
        delete mat;
    }

    // Compute error
    double infAtEnd[2];
    double errorAtEnd = computeError(qOut, conf.nu0, conf.nufreq, conf.cx, conf.cy, conf.cz, conf.endtime,
                                     xstart, xend, ystart, yend,
                                     zstart, zend, infAtEnd[0], infAtEnd[1], &errfield);
    std::vector<double> infsAtEnd(2);
    MPI_Gather(infAtEnd, 2, MPI_DOUBLE, &infsAtEnd[0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    infsAtEnd[0] = infAtEnd[0];
    infsAtEnd[1] = infAtEnd[1];

    if (isRoot)
    {
        std::cout << "Timing results:\n"
            //<< " - right hand side stencil: " << erhsTot * 1000. << " msec\n"
            //<< " - euler stencil: " << eeulerTot * 1000. << " msec\n"
            //<< " - rk stencil: " << erkTot * 1000. << " msec\n"
            //<< " - halo exchange: " << ecommTot * 1000. << " msec\n"
            << " - total time: " << e * 1000. << " msec\n"
            << "\n";

        std::cout << "Errors on all processes:\n";
        double totExactInf = 0.;
        double totErrorInf = 0.;
        for (int p = 0; p < 1; ++p)
        {
            totExactInf = std::max(totExactInf, infsAtEnd[2*p]);
            totErrorInf = std::max(totErrorInf, infsAtEnd[2*p+1]);
            std::cout << " - Process " << p << " has exact=" << infsAtEnd[2*p]
                << ", error=" << infsAtEnd[2*p+1] << "\n";
        }
        std::cout << "\n Total relative error: "
            << totErrorInf/totExactInf << "\n";
    }


    if (conf.mat)
    {
        fillQ(exactfield, conf.nu0, conf.nufreq, conf.cx, conf.cy, conf.cz, conf.endtime, xstart, xend, ystart, yend, zstart, zend);

        std::ostringstream fnameResult;
        fnameResult << "result_" << myPI << "_" << myPJ << "_" << myPK << ".mat";

        MatFile matfile(fnameResult.str());
        matfile.addField(qIn, -1);
        matfile.addField(qOut, -1);
        matfile.addField(errfield, -1);
        matfile.addField(exactfield, -1);
    }

    // Finalize GCL
    MPI_Finalize();
}

