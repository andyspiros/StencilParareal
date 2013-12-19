#include <iostream>
#include <cmath>
#include <cstdlib>

#include "boost/program_options.hpp"

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

namespace po = boost::program_options;


inline double exactQ(double nu, double x, double y, double z, double t)
{
    const double pi = 3.14159265358979;
    return sin(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*z) * exp(-12.*pi*pi*nu*t);
}

void fillQ(IJKRealField& q, double nu, double t,
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
                q(i, j, k) = exactQ(nu, x, y, z, t);
            }
}

double computeError(const IJKRealField& q, double nu, double t,
                double xstart, double xend,
                double ystart, double yend,
                double zstart, double zend,
                IJKRealField* errfield=0
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

    double errinf = 0.;
    double exact, e;
    double exactinf;

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                // Coordinates
                x = xstart + (2*i+1)*dxhalf;
                y = ystart + (2*j+1)*dyhalf;
                z = zstart + (2*k+1)*dzhalf;

                // Exact solution
                exact = exactQ(nu, x, y, z, t);
                exactinf = std::max(std::abs(exact), exactinf);

                // Error
                e = q(i,j,k) - exact;
                errinf = std::max(std::abs(e), errinf);

                // Error field
                if (errfield)
                    (*errfield)(i, j, k) = e;
            }
    return errinf / exactinf;
}

struct HeatConfiguration
{
    double nu;
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
        ("coeff", po::value<double>(), "Heat coefficient")
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

    if (vm.count("coeff"))
    {
        conf.nu = vm["coeff"].as<double>();
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

    if (dims[0] != dims[1] || dims[0] != dims[2])
    {
        if (isRoot)
        {
            std::cerr << "Non-cubic distribution is not allowed\n";
            std::cerr << "Aborting\n";
        }
        return 2;
    }

    // Get position on grid
    GridType grid(CyclicType(true, true, true), comm);
    int myPI, myPJ, myPK;
    grid.coords(myPI, myPJ, myPK);

    // Owned physical space
    const double Dx = 1. / dims[0];
    const double Dy = 1. / dims[1];
    const double Dz = 1. / dims[2];
    const double xstart = Dx * myPI;
    const double ystart = Dy * myPJ;
    const double zstart = Dz * myPK;
    const double xend = Dx * (myPI + 1);
    const double yend = Dy * (myPJ + 1);
    const double zend = Dz * (myPK + 1);

    // Output configuration
    if (isRoot)
    std::cout << "Running with:\n"
        << " - grid size: " << conf.gridsize << "x" << conf.gridsize << "x" << conf.gridsize << "\n"
        << " - spatial discretization step: " << conf.dx << "\n"
        << " - timestep size: " << conf.dt << "\n"
        << " - timesteps: " << conf.timesteps << "\n"
        << " - endtime: " << conf.endtime << "\n"
        << " - heat coefficient: " << conf.nu << "\n"
        << " - CFL: " << conf.cfl << "\n\n";

    // Create q
    IJKSize calculationDomain;
    calculationDomain.Init(conf.gridsize/dims[0], conf.gridsize/dims[1], conf.gridsize/dims[2]);
    KBoundary kboundary;
    kboundary.Init(-3, 3);
    IJKRealField q, errfield, exactfield;
    q.Init("q", calculationDomain, kboundary);
    errfield.Init("error", calculationDomain, kboundary);
    exactfield.Init("exact", calculationDomain, kboundary);

    // Initialize content of q
    fillQ(q, conf.nu, 0., xstart, xend, ystart, yend, zstart, zend);

    // Initialize stencil
    Heat heat(q, conf.nu, conf.dx, conf.dt, comm);

    // Initialize MAT file
    MatFile *mat;
    if (conf.mat)
    {
        std::ostringstream fnameMat;
        fnameMat << "heat_" << myPI << "_" << myPJ << "_" << myPK << ".mat";
        mat = new MatFile(fnameMat.str());
        mat->startCell("q", conf.timesteps);
        mat->addField(q, 0);
    }

    double error;

    for (int t = 1; t <= conf.timesteps; ++t) {
        heat.DoTimeStep();

        if (conf.mat)
            mat->addField(q, t);
    }

    // Close MAT file
    if (conf.mat)
    {
        mat->endCell();
        delete mat;
    }

    double errorAtEnd = computeError(q, conf.nu, conf.endtime,
                                     xstart, xend, ystart, yend,
                                     zstart, zend, &errfield);
    std::vector<double> errorsAtEnd(GCL::PROCS);
    MPI_Gather(&errorAtEnd, 1, MPI_DOUBLE, &errorsAtEnd[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (isRoot)
    for (int p = 0; p < GCL::PROCS; ++p)
    {
        std::cout << "Error in proc " << p << ": " << errorsAtEnd[p] << "\n";
    }

    fillQ(exactfield, conf.nu, conf.endtime, xstart, xend, ystart, yend, zstart, zend);

    std::ostringstream fnameResult;
    fnameResult << "result_" << myPI << "_" << myPJ << "_" << myPK << ".mat";

    MatFile matfile(fnameResult.str());
    matfile.addField(q, conf.timesteps);
    matfile.addField(errfield, conf.timesteps);
    matfile.addField(exactfield, conf.timesteps);

    // Finalize GCL
    GCL::GCL_Finalize();
}
