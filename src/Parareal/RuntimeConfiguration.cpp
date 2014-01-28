#include "RuntimeConfiguration.h"
#include "boost/program_options.hpp"
#include "mpi.h"

RuntimeConfiguration::RuntimeConfiguration(int argc, char **argv)
    : nu_(1.), cx_(1.), cy_(1.), cz_(1.)
    , gridSize_(32), endTime_(0.05)
    , timeStepsFine_(128), timeStepsCoarse_(32)
    , kmax_(5), mat_(false), async_(true)
{
    // Retrieve MPI information
    MPI_Comm_size(MPI_COMM_WORLD, &timeSlices_);

    // Create description of allowed options
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produce this help message")
        ("nu", po::value<double>(), "Heat coefficient")
        ("cx", po::value<double>(), "Advection velocity in x direction")
        ("cy", po::value<double>(), "Advection velocity in y direction")
        ("cz", po::value<double>(), "Advection velocity in z direction")
        ("gridSize", po::value<int>(), "Number of grid points along each direction")
        ("endTime", po::value<double>(), "Time to simulate")
        ("timeStepsFine", po::value<int>(), "Number of fine timesteps per timeslice")
        ("timeStepsCoarse", po::value<int>(), "Number of coarse timesteps per timeslice")
        ("kmax", po::value<int>(), "Number of parareal iterations")
        ("mat", "Serialize intermediate fields")
        ("noasync", "Avoid asynchronous communication")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Parse options

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        std::exit(0);
    }

    if (vm.count("nu"))
    {
        this->set_nu(vm["nu"].as<double>());
    }

    if (vm.count("cx"))
    {
        this->set_cx(vm["cx"].as<double>());
    }

    if (vm.count("cy"))
    {
        this->set_cy(vm["cy"].as<double>());
    }

    if (vm.count("cz"))
    {
        this->set_cz(vm["cz"].as<double>());
    }

    if (vm.count("gridSize"))
    {
        this->set_gridSize(vm["gridSize"].as<int>());
    }

    if (vm.count("endTime"))
    {
        this->set_endTime(vm["endTime"].as<double>());
    }

    if (vm.count("timeStepsCoarse"))
    {
        this->set_timeStepsCoarse(vm["timeStepsCoarse"].as<int>());
    }

    if (vm.count("timeStepsFine"))
    {
        this->set_timeStepsFine(vm["timeStepsFine"].as<int>());
    }

    if (vm.count("kmax"))
    {
        this->set_kmax(vm["kmax"].as<int>());
    }

    if (vm.count("mat"))
    {
        this->set_mat(true);
    }

    if (vm.count("noasync"))
    {
        this->set_async(false);
    }
}

