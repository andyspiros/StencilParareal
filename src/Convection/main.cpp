#include <iostream>
#include <cmath>
#include <cstdlib>

#include "boost/program_options.hpp"

#include "Convection.h"
#include "MatFile.h"
#include "HaloExchange3D.h"

#include "GCL.h"
#include "utils/layout_map.h"
#include "utils/boollist.h"
#include "L2/include/halo_exchange.h"
#include "L2/include/descriptors.h"
#include "L3/include/proc_grids_3D.h"
#include "L3/include/Halo_Exchange_3D.h"

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

inline double exactQ(double nu, double cx, double cy, double cz,
                     double x, double y, double z, double t)
{
    const double pi = 3.14159265358979;
    return sin(2.*pi*(x-cx*t))*sin(2.*pi*(y-cy*t))*sin(2.*pi*(z-cz*t)) * exp(-12.*pi*pi*nu*t);
}

void fillQ(ConvectionField& q,
                double nu, double cx, double cy, double cz,
                double t,
                double xstart, double xend,
                double ystart, double yend,
                double zstart, double zend
            )
{
    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();

    const double dxhalf = (xend-xstart)/iSize / 2.;
    const double dyhalf = (yend-ystart)/jSize / 2.;
    const double dzhalf = (zend-zstart)/kSize / 2.;

    double x, y, z;

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                x = xstart + (2*i+1)*dxhalf;
                y = ystart + (2*j+1)*dyhalf;
                z = zstart + (2*k+1)*dzhalf;
                q(i, j, k) = exactQ(nu, cx, cy, cz, x, y, z, t);
            }
}

template<typename TDataField>
double computeError(const TDataField& q,
                double nu, double cx, double cy, double cz,
                double t,
                double xstart, double xend,
                double ystart, double yend,
                double zstart, double zend
,
                double& exact_inf, double &error_inf,
                TDataField* errfield=0
            )
{
    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();

    const double dxhalf = (xend-xstart)/iSize / 2.;
    const double dyhalf = (yend-ystart)/jSize / 2.;
    const double dzhalf = (zend-zstart)/kSize / 2.;

    double x, y, z;

    error_inf = 0.;
    exact_inf = 0.;
    double exact, e;

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                // Coordinates
                x = xstart + (2*i+1)*dxhalf;
                y = ystart + (2*j+1)*dyhalf;
                z = zstart + (2*k+1)*dzhalf;

                // Exact solution
                exact = exactQ(nu, cx, cy, cz, x, y, z, t);
                exact_inf = std::max(std::abs(exact), exact_inf);

                // Error
                e = q(i,j,k) - exact;
                error_inf = std::max(std::abs(e), error_inf);

                // Error field
                if (errfield)
                    (*errfield)(i, j, k) = e;
            }
    return error_inf / exact_inf;
}

struct HeatConfiguration
{
    double nu;
    double cx, cy, cz;
    int gridsize;
    int timesteps;
    double endtime;
    double cfl;
    double dx;
    double dt;
    bool mat;
};

HeatConfiguration parseCommandLine(int argc, char **argv)
{
    HeatConfiguration conf;
    conf.nu = 1.;
    conf.cx = 1.;
    conf.cy = 1.;
    conf.cz = 1.;
    conf.gridsize = 32;
    conf.timesteps = 500;
    conf.endtime = 0.05;
    conf.dx = 1. / conf.gridsize;
    conf.dt = conf.endtime / conf.timesteps;
    conf.cfl = conf.dt / (conf.dx*conf.dx);
    conf.mat = false;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produce this help message")
        ("nu", po::value<double>(), "Heat coefficient")
        ("cx", po::value<double>(), "Advection velocity in x direction")
        ("cy", po::value<double>(), "Advection velocity in y direction")
        ("cz", po::value<double>(), "Advection velocity in z direction")
        ("gridsize", po::value<int>(), "Number of gridpoints along each direction")
        ("timesteps", po::value<int>(), "Number of timesteps")
        ("endtime", po::value<double>(), "Time to simulate")
        ("cfl", po::value<double>(), "Fixed CFL condition")
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

    if (vm.count("nu"))
    {
        conf.nu = vm["nu"].as<double>();
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

    if (vm.count("endtime"))
    {
        conf.endtime = vm["endtime"].as<double>();
        conf.dt = conf.endtime / conf.timesteps;
    }

    if (vm.count("cfl") && vm.count("gridsize") && vm.count("timesteps"))
    {
        std::cerr << "Error: setting cfl, gridsize and timesteps not allowed\n";
        std::cerr << "Aborting\n";
        std::exit(1);
    }

    bool gridsizeSet = false;
    bool timestepsSet = false;
    if (vm.count("gridsize"))
    {
        conf.gridsize = vm["gridsize"].as<int>();
        conf.dx = 1. / conf.gridsize;
        gridsizeSet = true;
    }

    if (vm.count("timesteps"))
    {
        conf.timesteps = vm["timesteps"].as<int>();
        conf.dt = conf.endtime / conf.timesteps;
        timestepsSet = true;
    }

    if (vm.count("cfl"))
    {
        double cfl = vm["cfl"].as<double>();
        if (gridsizeSet)
        {
            double dt = conf.dx*conf.dx * cfl / conf.nu;
            conf.timesteps = conf.endtime / dt + .5;
            conf.dt = conf.endtime / conf.timesteps;
            conf.cfl = conf.dt * conf.nu / (conf.dx*conf.dx);
        }
        else // Also apply if nothing else is specified
        {
            double dx = std::sqrt(conf.dt * conf.nu / cfl);
            conf.gridsize = 1. / dx + .5;
            conf.dx = 1. / conf.gridsize;
            conf.cfl = conf.dt * conf.nu / (conf.dx*conf.dx);
        }
    }

    // Compute the CFL in any case
    conf.cfl = conf.dt * conf.nu / (conf.dx * conf.dx);

    return conf;
}

int main(int argc, char **argv)
{
    // Parse command-line arguments
    HeatConfiguration conf = parseCommandLine(argc, argv);

    // Initialize GCL
    GCL::GCL_Init();
    typedef GCL::gcl_utils::boollist<3> CyclicType;
    typedef GCL::MPI_3D_process_grid_t<CyclicType> GridType;
    const bool isRoot = GCL::PID == 0;

    // Initialize grid
    int commsize;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int dims[3] = {0, 0, 0};
    int periods[3] = {1, 1, 1};
    MPI_Dims_create(GCL::PROCS, 3, dims);
    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm);

    std::vector<int> localSizesI, localSizesJ, localSizesK;
    std::vector<double> xstarts, xends, ystarts, yends, zstarts, zends;
    computeLocalSizes(conf.gridsize, dims[0], localSizesI, xstarts, xends);
    computeLocalSizes(conf.gridsize, dims[1], localSizesJ, ystarts, yends);
    computeLocalSizes(conf.gridsize, dims[2], localSizesK, zstarts, zends);


    // Get position on grid
    GridType grid(CyclicType(true, true, true), comm);
    int myPI, myPJ, myPK;
    grid.coords(myPI, myPJ, myPK);

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
            << " - heat coefficient: " << conf.nu << "\n"
            << " - advection velocity in x: " << conf.cx << "\n"
            << " - advection velocity in y: " << conf.cy << "\n"
            << " - advection velocity in z: " << conf.cz << "\n"
            << " - CFL: " << conf.cfl << "\n"
            << "\n";

        std::cout << "Spatial configuration:\n";
        for (int p = 0; p < GCL::PROCS; ++p)
        {
            int coords[3];
            MPI_Cart_coords(comm, p, 3, coords);
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
    ConvectionField q, errfield, exactfield;
    q.Init("q", calculationDomain, kboundary);
    errfield.Init("error", calculationDomain, kboundary);
    exactfield.Init("exact", calculationDomain, kboundary);

    // Initialize content of q
    fillQ(q, conf.nu, conf.cx, conf.cy, conf.cz, 0., xstart, xend, ystart, yend, zstart, zend);

    // Initialize stencil
    Convection convection(q, conf.nu, conf.cx, conf.cy, conf.cz, conf.dx, conf.dt, comm);

    // Initialize MAT file
    MatFile *mat;
    if (conf.mat)
    {
        std::ostringstream fnameMat;
        fnameMat << "convection_" << myPI << "_" << myPJ << "_" << myPK << ".mat";
        mat = new MatFile(fnameMat.str());
        mat->startCell("q", conf.timesteps);
        mat->addField(q, 0);
    }

    // Solve equation
    double e = MPI_Wtime();
    double erhsTot = 0., eeulerTot = 0., erkTot = 0., ecommTot = 0.;
    double erhs, eeuler, erk, ecomm;
    for (int t = 1; t <= conf.timesteps; ++t) {
        convection.DoTimeStep();

        if (conf.mat)
            mat->addField(q, t);
    }
    e = MPI_Wtime() - e;

    // Close MAT file
    if (conf.mat)
    {
        mat->endCell();
        delete mat;
    }

    // Compute error
    double infAtEnd[2];
    double errorAtEnd = computeError(q, conf.nu, conf.cx, conf.cy, conf.cz, conf.endtime,
                                     xstart, xend, ystart, yend,
                                     zstart, zend, infAtEnd[0], infAtEnd[1], &errfield);
    std::vector<double> infsAtEnd(2*GCL::PROCS);
    MPI_Gather(infAtEnd, 2, MPI_DOUBLE, &infsAtEnd[0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (isRoot)
    {
        std::cout << "Timing results:\n"
            << " - right hand side stencil: " << erhsTot * 1000. << " msec\n"
            << " - euler stencil: " << eeulerTot * 1000. << " msec\n"
            << " - rk stencil: " << erkTot * 1000. << " msec\n"
            << " - halo exchange: " << ecommTot * 1000. << " msec\n"
            << " - total time: " << e * 1000. << " msec\n"
            << "\n";

        std::cout << "Errors on all processes:\n";
        double totExactInf = 0.;
        double totErrorInf = 0.;
        for (int p = 0; p < GCL::PROCS; ++p)
        {
            totExactInf = std::max(totExactInf, infsAtEnd[2*p]);
            totErrorInf = std::max(totErrorInf, infsAtEnd[2*p+1]);
            std::cout << " - Process " << p << " has exact=" << infsAtEnd[2*p]
                << ", error=" << infsAtEnd[2*p+1] << "\n";
        }
        std::cout << "\n Total relative error: "
            << totErrorInf/totExactInf << "\n";
    }

    fillQ(exactfield, conf.nu, conf.cx, conf.cy, conf.cz, conf.endtime, xstart, xend, ystart, yend, zstart, zend);

    std::ostringstream fnameResult;
    fnameResult << "result_" << myPI << "_" << myPJ << "_" << myPK << ".mat";

    MatFile matfile(fnameResult.str());
    matfile.addField(q, -1);
    matfile.addField(errfield, -1);
    matfile.addField(exactfield, -1);

    // Finalize GCL
    GCL::GCL_Finalize();
}

